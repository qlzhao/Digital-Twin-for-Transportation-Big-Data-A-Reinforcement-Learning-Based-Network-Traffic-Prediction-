clear all
close all
clc
%==========================================================================
% Set the numbering of od flows to visualise the predicted effect of the
% traffic size in IoT.
ID1 = 30;
ID2 = 78;
tm_start = 1700;
end_time = 2015;
%==========================================================================

% Reading estimated data
% RL method
load('./result/TM_IIoT.mat'); % real_od
real_od = TM_IIoT;

% Reading estimated data
load(sprintf('./result/TM_Pre_GAN_DQN_IPFP_IIoT.mat'));
OM_TM  = TM_Pre_GAN_DQN_IPFP_IIoT;

% PCA method
load(sprintf('./result/TM_PCA_Prediction_IPFP_IIoT.mat'));
PCA_TM = TM_PCA_Prediction_IPFP_IIoT;

% SRSVD method
load(sprintf('./result/TM_SRMF_Prediction_IPFP_IIoT.mat'));
SRSVD_TM = TM_SRMF_Prediction_IPFP_IIoT;
%==========================================================================

% Uniform data
% c = 800;
% Intercepted partial estimates =============================================================

real_od = real_od(:, tm_start : end_time);
OM_TM = OM_TM(:, tm_start : end_time);
PCA_TM = PCA_TM(:, tm_start : end_time);
% ==========================================================================
% Set the [xmin, xmax], [ymin, ymax] ranges and graph using the true and predicted values
figure(1)
xmin = 0;
% xmax = 100;
xmax = size(real_od, 2) + 1;
ymin = min(real_od(ID1,:)) * 0.9;
ymax = max(real_od(ID1,:)) * 1.4;

subplot(2, 1, 1)
plot(real_od(ID1, :), 'b')
hold on
plot(OM_TM(ID1, :), 'r-.')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
ylabel('Volume of Traffic Flow', 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 14)
title(sprintf('OD %d', ID1), 'FontName', 'Times New Roman', 'FontSize', 14)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% ytickformat('%.1e');
legend('Real traffic data', 'Predicted traffic data')
axis([xmin xmax ymin ymax]);
ax = gca;
ax.YAxis.Exponent = 2;

%=========================================
% Set the [xmin, xmax], [ymin, ymax] ranges and graph using the true and predicted values
xmin = 0;
xmax = size(real_od, 2) + 1;
ymin = min(real_od(ID2, :)) * 0.9;
ymax = max(real_od(ID2, :)) * 1.7;

subplot(2, 1, 2)
plot(real_od(ID2, :), 'b')
hold on
plot(OM_TM(ID2, :), 'r-.')
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
ylabel('Volume of Traffic Flow', 'FontName', 'Times New Roman', 'FontSize', 14);
% ytickformat('%.1e');
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 14)
title(sprintf('OD %d', ID2), 'FontName', 'Times New Roman', 'FontSize', 14)
set(gca,'FontName', 'Times New Roman', 'FontSize', 14)
legend('Real traffic data', 'Predicted traffic data')
axis([xmin xmax ymin ymax]);
ax = gca;
% Set the index scale for the vertical axis
ax.YAxis.Exponent = 2;
