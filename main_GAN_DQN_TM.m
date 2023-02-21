% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This version builds a Q-network so that it can predict all traffic elements, 
% i.e. train a Q-network with all OD flows.
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % clear variables in base space
close all % close all Figure windows
clc % clear the contents of the command window
%% Parameter configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dataset parameterss
length_training = 500; % Training data set length and number of samples
length_data = 2016;
precision = 5; % The extent to which the traffic size narrows the data range

% Parameters for GAN
K = 100; % Number of updates of discriminator network D

num_G = 1000; % Number of samples generated each time

maxEpoch_G = 1;
miniBatchSize_G = 15000;
InitialLearnRate_G = 0.001;

maxEpoch_D = 1;
miniBatchSize_D = 15000;
InitialLearnRate_D = 0.001;

% Parameters for DQN
epsilon_threshold = 0.5;
C = 5; % Counting, for Target network updates
Iteration = 100; % Number of iterations

%% Loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./result/TM_Abilene.dat');
load('./result/A_Abilene.dat');
load('./result/link_Abilene.dat');

R = A_Abilene;                         % Routing Matrix
TM_Abilene = TM_Abilene.';
% Column coordinates indicate the length of time
TM = TM_Abilene(:, 1:length_data);      % Traffic Matrix
Link = link_Abilene(:, 1:length_data);  % Link loads

% Narrowing the data range of the traffic size, for better reinforcement learning calculations
TM = floor(TM ./ (10 ^ precision));

% Get the maximum value of the traffic size
max_TM = max(max(TM));

% TM_pre indicates the estimation matrix
TM_pre = zeros(size(TM, 1), (length_data - length_training));
% r_TM indicates the number of od streams in the network
r_TM = size(TM, 1);

% TM_prior indicates the traffic matrix used for training
TM_prior = TM(:, 1:length_training);
TM_prior_vector = reshape(TM_prior.', 1, r_TM * length_training);

%% Q matrix of the sample, in order to calculate the number of frequency of traffic size jumps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising the QS matrix, Rows and columns indicate the extent of network traffic
Q_S = zeros(max_TM + 1, max_TM + 1);

% Enumerate each od stream
for kkk = 1:r_TM
    % Enumerating the range of network traffic
    for i = 1:(max_TM + 1)
        % Get the moment when the current od traffic flow size is equal to this traffic value
        loc_c = find(TM(kkk, :) == i - 1);
        if isempty(loc_c) == 0
            % Since there is no information about the traffic size at the next moment at the last moment,
            % the traffic size at the current moment is removed
            loc_max = find(loc_c == size(TM, 2));
            if isempty(loc_max) == 0
                loc_c(loc_max) = [];
            end
            
            % Get the location of the next moment
            loc_c_next = loc_c + 1;
            
            % Iterate through the moments when the current traffic size is i
            for j = 1:size(loc_c, 2)
                % Due to the size of the flow from i to the next moment, the reward value plus one
                Q_S(i, TM(kkk, loc_c_next(j)) + 1) = Q_S(i, TM(kkk, loc_c_next(j)) + 1) + 1;
            end
        end
    end
end

%% Constructing discriminative network datasets=======================================================
% Input to the D network: a traffic size, batch size Total time for training minus one
Input_training_D = zeros(1, (length_training - 1) * r_TM);
% Output of the D network: rewards for each traffic size, batch size Total time for training minus one
Output_training_D = zeros(size(Q_S, 2), (length_training - 1) * r_TM);

% Based on TM matrix and Q_S matrix assign values to the training data, j Indexes assigned to the training data
j = 1;
% Enumerates the current TM traffic size matrix, kkk is the number of od streams, and i is the moment the network is running
for kkk = 1:r_TM
    for i = 2:length_training
        % Input is the current traffic size of the visit
        Input_training_D(j) = TM(kkk, (i - 1));  % Input: traffic size * 1
        Output_training_D(:, j) = Q_S(Input_training_D(1, i - 1) + 1, :); %  Output: Reward value for each range of traffic * 1
        j = j + 1;
    end
end

% Data standardisation ===============================================================
% Mean and variance of Input_training_D data
Mu_Input_D = mean(Input_training_D);
sig_Input_D = std(Input_training_D);

% Mean and variance of Output_training_D data
Mu_Output_D = mean(mean(Output_training_D));
% Stretching the Output_training_D data into a vector form
reshape_Output_training_D = reshape(Output_training_D, 1, size(Output_training_D, 1) * size(Output_training_D, 2));
% Get the standard deviation of the custom size vector
sig_Output_D = std(reshape_Output_training_D);

% Constructing datasets =============================================================== 
% Initialisation
% Using the mean and standard deviation data obtained above to obtain the actual data set for training
num_obv = size(Input_training_D, 2); % num_obv is the total number of current training sessions, i.e. the total number of moments
XTrain_D = cell(num_obv, 1);
YTrain_D = cell(num_obv, 1);

% Building data
for j = 1:num_obv
    % Normalise Input_training_D data and assign values
    XTrain_D{j} = (Input_training_D(j) - Mu_Output_D) ./ sig_Output_D;      % Input: traffic 2 * 1
    % Normalise Output_training_D data and assign values
    YTrain_D{j} = (Output_training_D(:,j) - Mu_Input_D) ./ sig_Input_D;     % Output: reward value 1 * 1
    % Since the output of the D network not only inputs the reward value, but also outputs the probability of whether the current is real data, the label 1 is added to the YTrain_D output at this point
    YTrain_D{j} = [YTrain_D{j}; 1];
end

% save('./result/XTrain_D.mat','XTrain_D')
% save('./result/YTrain_D.mat','YTrain_D')

%% Constructing generative network datasets =======================================================
% Construct the input and output of the G-network, the input is a random
% number and the number is the length of the training time (num_obv)
Input_training_G = rand(1, num_obv); % Input: random number
% The output of the G network is the input to the D network, so splice Input_training_D and Output_training_D into each other
Output_training_G = [Input_training_D; Output_training_D]; % Output: traffic value + reward value

% Standardisation?=============================================================== 
% Similar to above, the mean and standard deviation of the input and output of the G network are normalized
Mu_Input_G = mean(Input_training_G);
sig_Input_G = std(Input_training_G);

Mu_Output_G = mean(mean(Output_training_G));
reshape_Output_training_G = reshape(Output_training_G, 1, size(Output_training_G, 1) * size(Output_training_G, 2));
sig_Output_G = std(reshape_Output_training_G);

% Constructing datasets =============================================================== 

% Initialisation
% num_obv = length_training - 1;
% Calculate the true input and output of the G network from the mean and standard deviation
XTrain_G = cell(num_obv, 1);
YTrain_G = cell(num_obv, 1);

% Building data

for j = 1:num_obv
    XTrain_G{j} = (Input_training_G(j) - Mu_Input_G) ./ sig_Input_G;
    YTrain_G{j} = (Output_training_G(:, j) - Mu_Output_G) ./ sig_Output_G;
end

% Get the length of the output vectors of the D and G networks, plus one because the value of the current label is to be recorded
length_output_D = size(Output_training_D, 1) + 1;
length_output_G = length_output_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The training dataset was constructed above
%  XTrain_D traffic value 1 * 1
%  YTrain_D reward value + data truth (max_TM + 0or1) * 1
%  XTrain_G Random number 1 * 1
%  YTrain_G traffic value + reward value (max_TM + 1) * 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising the network
net_D = function_Train_D(XTrain_D, YTrain_D, miniBatchSize_D, maxEpoch_D, InitialLearnRate_D, length_output_D);
net_G = function_Train_G(XTrain_G, YTrain_G, miniBatchSize_G, maxEpoch_G, InitialLearnRate_G, length_output_G);

% Since the DQN algorithm requires two networks, a deep copy of the current discriminator network
net_D_Target = net_D;

state_space = 0:max_TM; % Status space, The range of the state space in reinforcement learning, which we define as the range of traffic size data
c = 1; % Counting, for Target network updates

%% Training GAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial state? ================================================================================
% XTrain_D_pool, YTrain_D_pool stores the data in the D network training experience pool, 
% first adding the original training data to it, and then adding the new data to the experience pool after the generation of new data by the network
XTrain_D_pool = cell([], 1);
YTrain_D_pool = cell([], 1);

XTrain_D_pool = XTrain_D;
YTrain_D_pool = YTrain_D;

%   =====================================================================================
% Total loop iteration of the reinforcement learning algorithm
for ite = 1:Iteration
    % Random enumeration of the current traffic size value X_D
    X_D = randi([0 max_TM], 1, 1);
    % The greed rate is obtained from the rand function
    epsilon = rand(1);
    % space derived for the state X_D_next at the next moment and Y_D (reward value + probability of whether it is true data) output by the D network
    X_D_next = cell(miniBatchSize_G, 1);
    Y_D = cell(miniBatchSize_G, 1);
    
    % When the greed rate reaches a threshold, the state X_D_next is obtained randomly at the next moment
    if epsilon > epsilon_threshold
        for kkk = 1:miniBatchSize_G
            X_D_next{kkk} = randi([0 max_TM], 1, 1); % Select a random action
        end
    else
        % Output the maximum reward value according to net_D and get the state at the next moment
        for kk = 1:miniBatchSize_G
            % Use the net_D network to give a prediction of the next reward value for the X_D data (traffic size)
            Y_D{kk} = predict(net_D, X_D, 'MiniBatchSize', miniBatchSize_G);
        end
        for kkk = 1:miniBatchSize_G
            % Enumerate the current reward value prediction results
            Y_D_max = Y_D{kkk};
            % Get the reward and location information for the one with the largest reward value
            [R_max, max_location] = max(Y_D_max(1:(size(Y_D_max, 1) - 1)));
            % The position with the largest reward value is the size of the traffic at the next moment
            X_D_next{kkk} = state_space(max_location);
        end
    end
    
    % Training Discriminative network D ========================================================================
    % Predicted reward for the next moment traffic size based on net_D_Target network and X_D_next traffic size
    Y_D_next = predict(net_D_Target, X_D_next, 'MiniBatchSize', miniBatchSize_D);  
    % num_obv: number of real samples
    % num_G: number of samples generated at a time
    % Put the X_D_next and Y_D_next information into the experience pool of the training D network at this point
    XTrain_D_pool = [XTrain_D_pool; X_D_next];
    YTrain_D_pool = [YTrain_D_pool; Y_D_next];
    
    sample_num = size(XTrain_D_pool, 1); % Get the entire size of the current experience pool
    sampling = randperm(sample_num); % Use the randperm function to generate a random permutation of the length of sample_num
    sampling = sampling(1:miniBatchSize_D);
    
    % Select the X_D_next_training, Y_D_next_training data currently used to train the D network
    % Allocation of space followed by assignment based on random sampling
    X_D_next_training = cell(miniBatchSize_D, 1); 
    Y_D_next_training = cell(miniBatchSize_D, 1);
    for kkkk = 1:miniBatchSize_D
        X_D_next_training{kkkk} = XTrain_D_pool{sampling(kkkk)};
        Y_D_next_training{kkkk} = YTrain_D_pool{sampling(kkkk)};
    end
    
    % Train the net_D network according to the parameters with a training count of K
    for k = 1:K
        % The loss function in the training process consists of two parts, one is the reward value after the jump and the other is the probability of whether the data is currently true or not
        net_D = function_Train_D(X_D_next_training, Y_D_next_training, miniBatchSize_D, maxEpoch_D, InitialLearnRate_D, length_output_D); %%%%%%%%%%%%%%%%%%%%%%%%%%% X_D_next
    end
    
    % Training Generator Network G ===========================================================================
    % Train the G-network based on the collated XTrain_G and YTrain_G information
    net_G = function_Train_G(XTrain_G, YTrain_G, miniBatchSize_G, maxEpoch_G, InitialLearnRate_G, length_output_G);
    
    % Generating samples =======================================================================
    % Create a random number based on the number of num_G to be generated
    Input_Testing_G = rand(1, num_G);
    % Allocating space for input data to the G network
    XTest_G = cell(num_G,1);
    % Since Mu_Input_G and sig_Input_G are used in the standardization above, they are directly used here
    for j = 1:num_G
        XTest_G{j} = (Input_Testing_G(j) - Mu_Input_G) ./ sig_Input_G;
    end
    % Get the input to the generated network, i.e. the fake data we created
    YTest_G = predict(net_G, XTest_G, 'MiniBatchSize', miniBatchSize_G);
    
    % Split the generated data into the input and output data of the D network 
    % and add it to the current experience pool
    XTrain_D_G = cell(num_G, 1);
    YTrain_D_G = cell(num_G, 1);
    for j = 1:num_G
        s = YTest_G{j};
        % The first figure is the traffic size
        s_input = s(1);
        % Other data are rewards after jumping to other traffic values
        s_output = s(2:end);
        XTrain_D_G{j} = s_input;
        YTrain_D_G{j} = [s_output; 0]; % Generated data, so the label = 0; 
    end
    
    % Add to XTrain_D_pool, YTrain_D_pool experience pool
    XTrain_D_pool = [XTrain_D_pool; XTrain_D_G]; % Update the sample pool
    YTrain_D_pool = [YTrain_D_pool; YTrain_D_G]; % Update the sample pool
    
    % Update the target discriminator network when the target discriminator update time is reached, 
    % this network will be more stable in terms of reward value
    if c == C
        net_D_Target = net_D;
        c = 1;
    else
        c = c + 1;
    end
    disp(ite);
end
%% Traffic estimates®¡ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use TM_pre to record all traffic sizes predicted, first adding the training set (the first half of the dataset) to the final result set
TM_pre = zeros(r_TM, length_data);
TM_pre(:, 1:length_training) = TM(:, 1:length_training);


input_pre = cell(1, 1);
% Enumerate all od streams
for num_OD = 1:r_TM
    % Enumerate all time periods of the current run
    for num_time = length_training: (length_data - 1)
        % Get the size of the flow at the moment and normalise it
        input_pre = (TM_pre(num_OD, num_time) - Mu_Input_D) ./ sig_Input_D;
        % Use discriminator net_D to predict the traffic size at the next moment
        TM_pre_next_R = predict(net_D, input_pre, 'MiniBatchSize', miniBatchSize_D);
        Q_TM_pre_next_R = TM_pre_next_R(1:(size(TM_pre_next_R, 1) - 1));  % Eliminate label
        % Get the size of the traffic with the largest reward value at the next moment, 
        % i.e. the most likely traffic jump at the next moment
        [max_report_R, location_max_report_R] = max(Q_TM_pre_next_R);
        % Since the position index starts at 1 and our traffic size starts at 0, 
        % it is subtracted by one. Assign the traffic size at the next moment.
        TM_pre(num_OD, num_time + 1) = location_max_report_R - 1;
    end
end

% Return to original flow size range
TM_pre = TM_pre * (10 ^ precision);
TM_Pre_GAN_DQN = TM_pre;

% close(h)
%% Save results ==========================================================

save('./result/TM_Pre_GAN_DQN.mat', 'TM_Pre_GAN_DQN')

%% The following IPFP algorithm is used to correct the x estimates?================================================
% The IPFP algorithm picks up where our previous work on traffic prediction left off
TM_Pre_GAN_DQN_IPFP = function_IPFP(Link, R, TM_Pre_GAN_DQN);

%% Save results ===================================================================

save('./result/TM_Pre_GAN_DQN_IPFP.mat', 'TM_Pre_GAN_DQN_IPFP');
disp('Done')
