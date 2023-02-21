%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias
%
% In the Abilene network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
%==========================================================================
% Read real traffic matrix data
load('./result/TM_Abilene.dat');% real_od
real_od = TM_Abilene.';

% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_IPFP.mat'));
RL_TM = TM_Pre_GAN_DQN_IPFP;

% PCA method
load('./result/PCA_Abilene.mat');
PCA_TM = x_avg;

% SRMF method
load('./result/TM_SRMF_Prediction_Abilene.mat');
SRSVD_TM = TM_SRMF_Prediction_Abilene;


% ==========================================================================

% Uniform data
tm_start = 501;
tm_end = 2015;

% Customize the current x-axis
x_axes = [1:tm_end - tm_start + 1];
% c = 800;
% Intercepted partial estimates =============================================================

% Reads actual data and predicted data from three methods
real_od = real_od(:, tm_start:tm_end);
RL_TM = RL_TM(:, tm_start:tm_end);
SRSVD_TM = SRSVD_TM(:, tm_start:tm_end);
PCA_TM = PCA_TM(:, tm_start:tm_end);


% OD Stream Sorting ==================================================================
% Calculate the mean of each matrix and reorder the positions by the mean
aver_RL_TM = mean(RL_TM');
aver_SRSVD_TM = mean(SRSVD_TM');
aver_PCA_TM = mean(PCA_TM');

[val_RL_TM, loc_RL_TM] = sort(aver_RL_TM, 'descend');
[val_SRSVD_TM, loc_SRSVD_TM] = sort(aver_SRSVD_TM, 'descend');
[val_PCA_TM, loc_PCA_TM] = sort(aver_PCA_TM, 'descend');

% Calculation bias ===================================================================
RL_TM_Bias = mean(RL_TM' - real_od');
PCA_TM_Bias = mean(PCA_TM' - real_od');
SRSVD_TM_Bias = mean(SRSVD_TM' - real_od');

% calculate standard deviaton 
SD_RL_TM = var(RL_TM');
SD_PCA_TM = var(PCA_TM');
SD_SRSVD_TM = var(SRSVD_TM');

% Deviation after sorting ===============================================================
Des_RL_TM_Bias = RL_TM_Bias(loc_RL_TM);
Des_SRSVD_TM_Bias = SRSVD_TM_Bias(loc_SRSVD_TM);
Des_PCA_TM_Bias = PCA_TM_Bias(loc_PCA_TM);

% ==========================================================================
% Graphing of Des_RL_TM_Bias, Des_PCA_TM_Bias, Des_SRSVD_TM_Bias data
figure(1);
subplot(2, 1, 1);
plot(Des_RL_TM_Bias, 'or');
hold on;
plot(Des_PCA_TM_Bias, '*b');
plot(Des_SRSVD_TM_Bias, 'Pg');
legend('Our method',  'PCA' , 'SRMF')
% Set the font for the x-axis display to Times New Roman and the font size to 14
xlabel('Flow ID', 'FontSize', 14, 'FontName', 'Times New Roman');
% Set the font for the y-axis display to Times New Roman and the font size to 14
ylabel('Bias', 'FontSize', 14, 'FontName', 'Times New Roman');
hold off;
grid off;
box on;
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
% Set the title of the current figure
title('The Bias of OD Flows in Abilene')

% ==========================================================================
% Draw a second graph with the horizontal axis set to the standard deviation of the matrix
figure(1);
subplot(2, 1, 2);
plot(SD_RL_TM, RL_TM_Bias, 'or');
hold on;
plot(SD_PCA_TM, PCA_TM_Bias, '*b');
plot(SD_SRSVD_TM, SRSVD_TM_Bias, 'Pg');
legend('Our method', 'PCA', 'SRMF')

% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
xlabel('SD in Error', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Bias', 'FontSize', 14, 'FontName', 'Times New Roman');
hold off;
grid off;
box on;
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
title('The Bias and SD in Abilene')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias
%
% In the GEANT network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear 
% close all
% clc
%==========================================================================
% Read real traffic matrix data
load('./result/tm_GEANT.mat');
real_od = tm_GEANT;

% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_GEANT_IPFP.mat'));
RL_TM  = TM_Pre_GAN_DQN_GEANT_IPFP;

% PCA method
load('./result/TM_PCA_Prediction_GEANT.mat');
PCA_TM = TM_PCA_GEANT;

% SRMF method
load('./result/TM_SRMF_Prediction_GEANT.mat')
SRSVD_TM = x_srmf;

%==========================================================================

% Uniform data, Current time range of the GENAT network is [tm_start, tm_end]
tm_start = 2001;
tm_end = 3360;

% Customize the current x-axis
x_axes = [1:tm_end - tm_start + 1];
% c = 800;

% Intercepted partial estimates =============================================================
% Intercept data based on tm_start and tm_end, including real traffic data and predicted estimates
real_od = real_od(:, tm_start:tm_end);
RL_TM = RL_TM(:, tm_start:tm_end);
PCA_TM = PCA_TM(:, tm_start:tm_end);
SRSVD_TM = SRSVD_TM(:, tm_start:tm_end);

% OD Stream Sorting==================================================================
aver_RL_TM = mean(RL_TM');
aver_PCA_TM = mean(PCA_TM');
aver_SRSVD_TM = mean(SRSVD_TM');

% Sorted by mean value from largest to smallest
[val_RL_TM, loc_RL_TM] = sort(aver_RL_TM, 'descend');
[val_PCA_TM, loc_PCA_TM] = sort(aver_PCA_TM, 'descend');
[val_SRSVD_TM, loc_SRSVD_TM] = sort(aver_SRSVD_TM, 'descend');


% Calculation bias ===================================================================
RL_TM_Bias = mean(RL_TM' - real_od');
PCA_TM_Bias = mean(PCA_TM' - real_od');
SRSVD_TM_Bias = mean(SRSVD_TM'- real_od')';


% calculate standard deviaton 
SD_RL_TM = var(RL_TM');
SD_PCA_TM = var(PCA_TM');
SD_SRSVD_TM = var(SRSVD_TM');


% Deviation after sorting ===============================================================
Des_RL_TM_Bias = RL_TM_Bias(loc_RL_TM);
Des_SRSVD_TM_Bias = SRSVD_TM_Bias(loc_SRSVD_TM);
Des_PCA_TM_Bias = PCA_TM_Bias(loc_PCA_TM);

% ==========================================================================
% Graphing of Des_RL_TM_Bias, Des_PCA_TM_Bias, Des_SRSVD_TM_Bias data
figure(2)
subplot(2, 1, 1)
plot(Des_RL_TM_Bias, 'or');
hold on;
plot(Des_PCA_TM_Bias, '*b');
plot(Des_SRSVD_TM_Bias, 'Pg');
legend('Our method', 'PCA', 'SRMF')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
xlabel('Flow ID', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Bias','FontSize', 14, 'FontName', 'Times New Roman');
hold off;
grid off;
box on;
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
title('The Bias of OD Flows in GEANT')

%==========================================================================
% Draw a second graph with the horizontal axis set to the standard deviation of the matrix
figure(2)
subplot(2, 1, 2)
plot(SD_RL_TM, RL_TM_Bias, 'or');
hold on
plot(SD_PCA_TM, PCA_TM_Bias, '*b');
plot(SD_SRSVD_TM, SRSVD_TM_Bias, 'Pg');

legend('Our method', 'PCA', 'SRMF')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
xlabel('SD in Error', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Bias', 'FontSize', 14, 'FontName', 'Times New Roman');
hold off;
grid off;
box on;
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
title('The Bias and SD in GEANT')


%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias
% 
% In the IoT network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear 
% close all
% clc
%==========================================================================
% Read real traffic matrix data
load('./result/TM_IIoT.mat'); % real_od
real_od = TM_IIoT;

% Reading estimated data
% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_IPFP_IIoT.mat'));
RL_TM  = TM_Pre_GAN_DQN_IPFP_IIoT;

% PCA method
load(sprintf('./result/TM_PCA_Prediction_IPFP_IIoT.mat'));
PCA_TM = TM_PCA_Prediction_IPFP_IIoT;

% SRMF method
load(sprintf('./result/TM_SRMF_Prediction_IPFP_IIoT.mat'));
SRSVD_TM = TM_SRMF_Prediction_IPFP_IIoT;

%==========================================================================

% Uniform data
tm_start = 1700;
tm_end = 2015;

x_axes = [1:tm_end - tm_start + 1];
% c = 800;
% Intercepted partial estimate =============================================================
% Reads actual data and predicted data from three methods in IoT
real_od = real_od(:, tm_start:tm_end);
RL_TM = RL_TM(:, tm_start:tm_end);
SRSVD_TM = SRSVD_TM(:, tm_start:tm_end);
PCA_TM = PCA_TM(:, tm_start:tm_end);


% OD Stream Sorting ==================================================================
aver_RL_TM = mean(RL_TM');
aver_SRSVD_TM = mean(SRSVD_TM');
aver_PCA_TM = mean(PCA_TM');

[val_RL_TM, loc_RL_TM] = sort(aver_RL_TM, 'descend');
[val_SRSVD_TM, loc_SRSVD_TM] = sort(aver_SRSVD_TM, 'descend');
[val_PCA_TM, loc_PCA_TM] = sort(aver_PCA_TM, 'descend');

% Calculation bias ===================================================================
RL_TM_Bias = mean(RL_TM' - real_od');
PCA_TM_Bias = mean(PCA_TM' - real_od');
SRSVD_TM_Bias = mean(SRSVD_TM' - real_od');

% Calulate standard deviaton 
SD_RL_TM = var(RL_TM');
SD_PCA_TM = var(PCA_TM');
SD_SRSVD_TM = var(SRSVD_TM');


% Deviation after sorting===============================================================
Des_RL_TM_Bias = RL_TM_Bias(loc_RL_TM);
Des_SRSVD_TM_Bias = SRSVD_TM_Bias(loc_SRSVD_TM);
Des_PCA_TM_Bias = PCA_TM_Bias(loc_PCA_TM);

% ==========================================================================
figure(3)
subplot(2,1,1)
plot(Des_RL_TM_Bias, 'or');
hold on;
plot(Des_PCA_TM_Bias, '*b');
plot(Des_SRSVD_TM_Bias, 'Pg');

legend('Our method', 'PCA', 'SRMF')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
xlabel('Flow ID', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Bias', 'FontSize', 14, 'FontName', 'Times New Roman');
hold off;
grid off;
box on;
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
title('The Bias of OD Flows in IoT')
ax = gca;
ax.YAxis.Exponent = 2;
% ==========================================================================
% Draw a second graph with the horizontal axis set to the standard deviation of the matrix
figure(3)
subplot(2, 1, 2)
plot(SD_RL_TM, RL_TM_Bias, 'or');
hold on;
plot(SD_PCA_TM, PCA_TM_Bias, '*b');
plot(SD_SRSVD_TM, SRSVD_TM_Bias, 'Pg');

legend('Our method',  'PCA' , 'SRMF')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
xlabel('SD in Error', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Bias', 'FontSize', 14, 'FontName', 'Times New Roman');
hold off;
grid off;
box on;
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')
title('The Bias and SD in IoT')
ax = gca;
ax.YAxis.Exponent = 2;
