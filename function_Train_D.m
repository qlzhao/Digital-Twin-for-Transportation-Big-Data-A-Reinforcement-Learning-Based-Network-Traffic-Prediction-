function net = function_Train_D(XTrain, YTrain, miniBatchSize, maxEpoch, InitialLearnRate, length_output_D)

% Construction of discriminative networks combining GAN and DQN, number of network layers, initialization of various parameters during training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers = [
    sequenceInputLayer(1, "Name", "sequence")
    fullyConnectedLayer(1000, "Name", "fc_1")
    fullyConnectedLayer(1000, "Name", "fc_2")
    % length_output_D is the size of the vector to be output from network D 
    fullyConnectedLayer(length_output_D, "Name", "fc_3")
    regressionLayer("Name", "regressionoutput")];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adam is a neural network gradient descent method
% maxEpoch is the maximum number of rounds to be used for training
% miniBatchSize is the size of the smallest batch used for each training iteration
% GradientThreshold is the positive threshold for gradient gradients.
% InitialLearnRate is the initial learning rate used for training. If the learning rate is too low, the training will take a long time, but if the learning rate is too high, the training may fall into sub-optimal results.
% LearnRateSchedule is used to specify a method for reducing the overall learning rate during training. piecewise method for multiplying the learning rate by a factor every time a certain number of periods have elapsed.
% The effect of LearnRateDropPeriod is that each time this period is elapsed, a learning rate drop factor is applied to the global learning rate.
% The LearnRateDropFactor is a multiplicative factor that is applied to the learning rate each time a certain number of epochs pass
% A Verbose of 0 means that training progress information is not printed.
options = trainingOptions('adam', ...
    'MaxEpochs', maxEpoch, ...
    'MiniBatchSize', miniBatchSize, ...
    'GradientThreshold', 1, ...
    'InitialLearnRate', InitialLearnRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 125, ...
    'LearnRateDropFactor', 0.2, ...
    'Verbose', 0); 
%      'Plots','training-progress',...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
net = trainNetwork(XTrain, YTrain, layers, options);

end
