function A = function_Pre_GEANT

%==========================================================================
% Read real traffic matrix data
load('./result/tm_GEANT.mat');% real_od
real_od_GEANT  = tm_GEANT;

% Reading estimated data
% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_GEANT_IPFP.mat'));
RL_TM_GEANT = TM_Pre_GAN_DQN_GEANT_IPFP;

% PCA method
load('./result/TM_PCA_Prediction_GEANT.mat');
PCA_TM_GEANT = TM_PCA_GEANT;

% SRMF method
load('./result/TM_SRMF_Prediction_GEANT.mat');
SRMF_TM_GEANT = x_srmf;


% ==========================================================================

% Uniform data
tm_start_GEANT = 2001;
tm_end_GEANT = 3360;

x_axes = [1:tm_end_GEANT - tm_start_GEANT + 1];
% c = 800;
% Intercepted partial estimates =============================================================
real_od_GEANT = real_od_GEANT(:, tm_start_GEANT : tm_end_GEANT);
RL_TM_GEANT = RL_TM_GEANT(:, tm_start_GEANT : tm_end_GEANT);
PCA_TM_GEANT = PCA_TM_GEANT(:, tm_start_GEANT : tm_end_GEANT);
SRMF_TM_GEANT = SRMF_TM_GEANT(:, tm_start_GEANT : tm_end_GEANT);

%==========================================================================
% Calculation of the data by means of the improvement rate formula
RL_PCA_A = (sum(sum(abs(PCA_TM_GEANT - real_od_GEANT))) - sum(sum(abs(RL_TM_GEANT - real_od_GEANT)))) / sum(sum(abs(PCA_TM_GEANT - real_od_GEANT)));
RL_SRMF_A = (sum(sum(abs(SRMF_TM_GEANT - real_od_GEANT))) - sum(sum(abs(RL_TM_GEANT - real_od_GEANT)))) / sum(sum(abs(SRMF_TM_GEANT - real_od_GEANT)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [RL_PCA_A, RL_SRMF_A];
