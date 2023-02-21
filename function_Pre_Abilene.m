function A = function_Pre_Abilene

% ==========================================================================
% Read real traffic matrix data
load('./result/TM_Abilene.dat'); % real_od
real_od  = TM_Abilene.';

% Reading estimated data
% RL method
load('./result/TM_Pre_GAN_DQN_IPFP.mat');
RL_TM = TM_Pre_GAN_DQN_IPFP;

% PCA method
load('./result/PCA_Abilene.mat');
PCA_TM = x_avg;

% SRMF method
load('./result/TM_SRMF_Prediction_Abilene.mat');
SRMF_TM = TM_SRMF_Prediction_Abilene;

%==========================================================================

% Uniform data
tm_start = 501;
tm_end = 2015;

x_axes = [1:tm_end - tm_start + 1];
% c = 800;
% Intercepted partial estimates =============================================================
real_od = real_od(:, tm_start : tm_end);
RL_TM = RL_TM(:, tm_start : tm_end);
PCA_TM = PCA_TM(:, tm_start : tm_end);
SRMF_TM = SRMF_TM(:, tm_start : tm_end);

% ==========================================================================
% Calculation of the data by means of the improvement rate formula
RL_PCA_A = (sum(sum(abs(PCA_TM - real_od))) - sum(sum(abs(RL_TM - real_od)))) / sum(sum(abs(PCA_TM - real_od)));
RL_SRMF_A = (sum(sum(abs(SRMF_TM - real_od))) - sum(sum(abs(RL_TM - real_od)))) / sum(sum(abs(SRMF_TM - real_od)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [RL_PCA_A, RL_SRMF_A];
