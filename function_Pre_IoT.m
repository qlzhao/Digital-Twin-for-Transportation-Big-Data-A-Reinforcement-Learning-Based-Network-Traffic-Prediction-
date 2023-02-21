function A = function_Pre_IoT

%==========================================================================
% Read real traffic matrix data
load('./result/TM_IIoT.mat'); % real_od
real_od  = TM_IIoT;

% Reading estimated data
% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_IPFP_IIoT.mat'));
RL_TM = TM_Pre_GAN_DQN_IPFP_IIoT;

% PCA method
load(sprintf('./result/TM_PCA_Prediction_IPFP_IIoT.mat'));
PCA_TM = TM_PCA_Prediction_IPFP_IIoT;

% SRMF method
load(sprintf('./result/TM_SRMF_Prediction_IPFP_IIoT.mat'));
SRMF_TM = TM_SRMF_Prediction_IPFP_IIoT;

% ==========================================================================

% Uniform data
tm_start = 1700;
tm_end = 2015;

x_axes = [1 : tm_end - tm_start + 1];
% c = 800;
% Intercepted partial estimates =============================================================
real_od = real_od(:, tm_start:tm_end);
RL_TM = RL_TM(:, tm_start:tm_end);
PCA_TM = PCA_TM(:, tm_start:tm_end);
SRMF_TM = SRMF_TM(:, tm_start:tm_end);

% ==========================================================================
% Traversing the matrix
display(sum(abs(PCA_TM - real_od)));
display(sum(abs(RL_TM - real_od)));

result1 = [];
for i = 1:size(real_od, 1)
    all = 0.0;
    for j = 1:size(real_od, 2)
        all = all + abs(RL_TM(i, j) - PCA_TM(i, j));
    end
    result1 = [result1 all];
end
display(result1);

% Calculation of the data by means of the improvement rate formula
RL_PCA_A = (sum(sum(abs(PCA_TM - real_od))) - sum(sum(abs(RL_TM - real_od)))) / sum(sum(abs(PCA_TM - real_od)));
RL_SRMF_A = (sum(sum(abs(SRMF_TM - real_od))) - sum(sum(abs(RL_TM - real_od)))) / sum(sum(abs(SRMF_TM - real_od)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [RL_PCA_A, RL_SRMF_A];
