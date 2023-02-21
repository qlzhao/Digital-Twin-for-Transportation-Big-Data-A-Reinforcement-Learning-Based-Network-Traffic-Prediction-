clear all
close all
clc
% ==========================================================================
% Set the numbering of od flows to visualise the predicted effect of the
% traffic size in GEANT.
ID1 = 33;
ID2 = 44;

tm_start = 2001;
end_time = 3360;
%==========================================================================

% Reading estimated data
% RL method
load('./result/tm_GEANT.mat'); % real_od
real_od_GEANT = tm_GEANT;

% Reading estimated data

load(sprintf('./result/TM_Pre_GAN_DQN_GEANT_IPFP.mat'));
RL_TM_GEANT  = TM_Pre_GAN_DQN_GEANT_IPFP;

load('./result/TM_PCA_Prediction_GEANT.mat');
PCA_TM_GEANT = TM_PCA_GEANT;

load('./result/TM_SRMF_Prediction_GEANT.mat')
SRSVD_TM_GEANT = x_srmf;
%==========================================================================

% Uniform data
% c = 800;
% Intercepted partial estimates =============================================================

real_od = real_od_GEANT(:, tm_start : end_time);
RL_TM_GEANT = RL_TM_GEANT(:, tm_start : end_time);
PCA_TM_GEANT = PCA_TM_GEANT(:, tm_start : end_time);
% ==========================================================================
% Set the [xmin, xmax], [ymin, ymax] ranges and graph using the true and predicted values
figure(1)

xmin = 0;
% xmax = 100;
% xmax = 1000;
xmax = size(real_od, 2) + 1;
% xmax = max(size(real_od), len(FSR_TM_GEANT));

ymin = min(real_od(ID1,:)) * 0.9;
ymax = max(real_od(ID1,:)) * 1.4;

subplot(2,1,1)
plot(real_od(ID1,:), 'b')
hold on
plot(RL_TM_GEANT(ID1, :), 'r-.')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
ylabel('Volume of Traffic Flow', 'FontName', 'Times New Roman', 'FontSize',14);
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 14)
title(sprintf('OD %d', ID1), 'FontName', 'Times New Roman', 'FontSize', 14)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
legend('Real traffic data', 'Predicted traffic data')
axis([xmin xmax ymin ymax]);

%=========================================
% Set the [xmin, xmax], [ymin, ymax] ranges and graph using the true and predicted values
xmin = 0;
% xmax = 100;
xmax = size(real_od, 2) + 1;
ymin = min(real_od(ID2, :)) * 0.9;
ymax = max(real_od(ID2, :)) * 1.7;

subplot(2, 1, 2)
plot(real_od(ID2,:), 'b')
hold on
plot(RL_TM_GEANT(ID2,:), 'r-.')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
ylabel('Volume of Traffic Flow', 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 14)
title(sprintf('OD %d', ID2), 'FontName', 'Times New Roman', 'FontSize', 14)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
legend('Real traffic data', 'Predicted traffic data')
axis([xmin xmax ymin ymax]);
