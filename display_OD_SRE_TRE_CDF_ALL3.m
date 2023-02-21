% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is mainly used to calculate the relative error between the predicted and true values of the traffic matrix
%
% Creator:  Laisen Nie
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear variables in base space, close all Figure windows, clear the
% contents of the command window.
clear
close all
clc

% Read the estimated data from the Abilene network
load('./result/TM_Abilene.dat'); % real_od
real_od = TM_Abilene.'; % Invert the current estimation matrix, the rows indicate the number of the od stream and the columns indicate the time period

% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_IPFP.mat'));
RL_TM  = TM_Pre_GAN_DQN_IPFP;

% PCA method This is the earliest version of the PCA algorithm, with a typical training data set length of 500.
load('./result/PCA_Abilene.mat');
PCA_TM = x_avg;

% SRMF method
load('./result/TM_SRMF_Prediction_Abilene.mat');
SRSVD_TM = TM_SRMF_Prediction_Abilene;

% ==========================================================================

% Uniform data, test set start position and end position
tm_start = 501;
end_time = 2015;

% Customise the current axis position
x_axes = [1 : end_time - tm_start + 1];
% c = 800;

%==========================================================================
% Uniform data, Take the real predicted traffic data according to tm_start and end_time
real_od = real_od(:, tm_start : end_time);
RL_TM = RL_TM(:, tm_start : end_time);
PCA_TM = PCA_TM(:, tm_start : end_time);
SRSVD_TM = SRSVD_TM(:, tm_start : end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial relative error in Abilene %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rows of RL_TM represent the number of the od stream and the columns represent the length of time
[r, c] = size(RL_TM);

% The spatial relative error between the predicted flow size and the true flow size for the three methods is obtained
SRE_RL_TM = function_SRE(RL_TM, real_od);
SRE_PCA_TM = function_SRE(PCA_TM, real_od);
SRE_SRSVD_TM = function_SRE(SRSVD_TM, real_od);

% Set the range of the x and y axes of the spatial relative error figure for the three methods in the Abilene network
len = r;
xmin = 0; xmax = len;
ymin = 0;
ymax = max(SRE_SRSVD_TM);
% disp(SRE_SRSVD_TM);

% Set the various properties of the chart's display
figure(1);
subplot(2, 1, 1)
% -r means the color of the display line is red
plot(SRE_RL_TM, '-r', 'linewidth', 2);
hold on;
% -b means the color of the display line is blue
plot(SRE_PCA_TM, '-b');
% -g means the color of the display line is green
plot(SRE_SRSVD_TM, '-g');
legend('Our method', 'PCA', 'SRMF')
% Set the display font of the y-axis labels to Times New Roman, with a font size of 14
ylabel('SRE', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin, xmax, ymin, ymax]);
% Set the display font of the x-axis labels to Times New Roman, with a font size of 14
xlabel( 'Flow ID', 'FontName', 'Times New Roman', 'FontSize', 14);
title('SREs in Abilene')
ax = gca;
ax.YAxis.Exponent = 1;

% The average spatial relative error of the three methods in the Abilene network was obtained
mean_SRE_RL = mean(SRE_RL_TM);
% disp(mean_SRE_RL);
mean_SRE_SRSVD = mean(SRE_SRSVD_TM);
mean_SRE_PCA = mean(SRE_PCA_TM);

%%%%%%%%%%%%%%%%%%%%%%%%% Temporal Relative error in Abilene %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rows of real_od represent the number of the od stream and the columns represent the length of time
[r, c] = size(real_od);

% Customise the current axis position
x_axes = 1:(end_time - tm_start + 1);

% The Temporal relative error between the predicted traffic size and the real traffic size for the three methods is obtained
TRE_RL_TM = function_TRE(RL_TM, real_od);
TRE_PCA_TM = function_TRE(PCA_TM, real_od);
TRE_SRSVD_TM = function_TRE(SRSVD_TM, real_od);

% Obtain the average Temporal relative error of the three methods
mean_TRE_RL = mean(TRE_RL_TM);
% disp(mean_TRE_RL);
mean_TRE_SRSVD = mean(TRE_SRSVD_TM);
mean_TRE_PCA = mean(TRE_PCA_TM);

% Set the current position range of the x and y axes
len = c;
xmin = 1; xmax = 1 + len;
ymin = 0; ymax = max(TRE_PCA_TM) * 2;

% Superimpose the current figure under the previous one
figure(1);
subplot(2, 1, 2);
% Set the various label attribute information for the figure TREs in Abilene
plot(x_axes, TRE_RL_TM, '-r', 'linewidth', 2);
hold on;
plot(x_axes, TRE_PCA_TM, '-b');
plot(x_axes, TRE_SRSVD_TM, '-g');
legend('Our method', 'PCA', 'SRMF')
% Set the display font of the x-axis labels to Times New Roman, with a font size of 14
xlabel('Time Slot Order', 'FontName', 'Times New Roman', 'FontSize', 14);
% Set the display font of the y-axis labels to Times New Roman, with a font size of 14
ylabel('TRE', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin, xmax, ymin, ymax]);
title('TREs in Abilene')

%%%%%%%%%%%%%%%%%%%%%%    Spatial Relative Error CDF in Abilene   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculation of the cumulative distribution function of the spatial relative errors of the three methods
[errsp_vector_RL_TM, errsp_vector_RL_TM_cdf] = function_SRE_CDF(SRE_RL_TM);
[errsp_vector_PCA_TM, errsp_vector_PCA_TM_cdf] = function_SRE_CDF(SRE_PCA_TM);
[errsp_vector_SRSVD_TM, errsp_vector_SRSVD_TM_cdf] = function_SRE_CDF(SRE_SRSVD_TM);

% Set graph properties for displaying relative spatial error probability density functions
% Find [xmin, xmax] [ymin, ymax] range
r1 = max(errsp_vector_RL_TM);
r2 = max(errsp_vector_PCA_TM);
r4 = max(errsp_vector_SRSVD_TM);

len = max(r1, r2);
% len = max(len, r3);
% len = max(len, r4);
xmin = 0; xmax = len;
c1 = max(errsp_vector_RL_TM_cdf);
c3 = max(errsp_vector_PCA_TM_cdf);
len = max(c1, c3);
ymin = 0; ymax = len;

% As above, set the colour of the line graph, font size, etc.
figure(2);
subplot(2, 1, 1)
plot(errsp_vector_RL_TM, errsp_vector_RL_TM_cdf, '-r', 'linewidth', 2);
% disp(errsp_vector_RL_TM);
% disp(errsp_vector_RL_TM_cdf);
hold on;
plot(errsp_vector_PCA_TM, errsp_vector_PCA_TM_cdf, '-b');
% disp(errsp_vector_PCA_TM);
% disp(errsp_vector_PCA_TM_cdf);

plot(errsp_vector_SRSVD_TM, errsp_vector_SRSVD_TM_cdf, '-g');
% disp(errsp_vector_SRSVD_TM);
% disp(errsp_vector_SRSVD_TM_cdf);

legend('Our method', 'PCA', 'SRMF')
xlabel( 'SRE', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel( 'CDF', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
axis([xmin, xmax, ymin, ymax]);
title('CDF of SRE in Abilene.')
%%%%%%%%%%%%%%%%%%%%%   Temporal Relative Error CDF in Abilene   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the cumulative distribution function of the Temporal relative
% errors of the three methods in Abilene
[errsp_vector_val_RL_TM, errsp_vector_val_RL_TM_cdf] = function_TRE_CDF(TRE_RL_TM);
[errsp_vector_val_PCA_TM, errsp_vector_val_PCA_TM_cdf] = function_TRE_CDF(TRE_PCA_TM);
[errsp_vector_val_SRSVD_TM, errsp_vector_val_SRSVD_TM_cdf] = function_TRE_CDF(TRE_SRSVD_TM);

% Find [xmin, xmax] [ymin, ymax] range
r1 = max(errsp_vector_val_RL_TM);
r3 = max(errsp_vector_val_PCA_TM);
len = max(r1, r3);
xmin = 0; xmax = len * 1.5;
c1 = max(errsp_vector_val_RL_TM_cdf);
c3 = max(errsp_vector_val_PCA_TM_cdf);
len = max(c1, c3);
ymin = 0; ymax = len;

figure(2)
subplot(2, 1, 2)
% As above, set the colour of the line graph, font size, etc.
plot(errsp_vector_val_RL_TM, errsp_vector_val_RL_TM_cdf, '-r', 'linewidth', 2);
% disp(errsp_vector_val_RL_TM);
% disp(errsp_vector_val_RL_TM_cdf);

axis([xmin, xmax, ymin, ymax]);
hold on;
plot(errsp_vector_val_PCA_TM, errsp_vector_val_PCA_TM_cdf, '-b');
% disp(errsp_vector_val_PCA_TM);
% disp(errsp_vector_val_PCA_TM_cdf);

plot(errsp_vector_val_SRSVD_TM, errsp_vector_val_SRSVD_TM_cdf, '-g');
% disp(errsp_vector_val_SRSVD_TM);
% disp(errsp_vector_val_SRSVD_TM_cdf);

% Set the display font to Times New Roman and the font size to 14
legend('Our method', 'PCA', 'SRMF')
xlabel('TRE', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('CDF', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
set(0, 'defaultfigurecolor', 'w') 
title('CDF of TRE in Abilene.')
% set(gcf, 'position', [500, 500, 500, 500])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is mainly used to calculate the relative error between the predicted and true values of the traffic matrix
% The comments in this section are similar to the above in that they both logically read the data and later draw the figure
%
% Creator:  Laisen Nie
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all
% clc

% Reading estimated data
load('./result/tm_GEANT.mat');
real_od_GEANT = tm_GEANT;

% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_GEANT_IPFP.mat'));
RL_TM_GEANT = TM_Pre_GAN_DQN_GEANT_IPFP;

% PCA method
load('./result/TM_PCA_Prediction_GEANT.mat');
PCA_TM_GEANT = TM_PCA_GEANT;

% SRMF method
load('./result/TM_SRMF_Prediction_GEANT.mat')
SRSVD_TM_GEANT = x_srmf;

%==========================================================================
% Uniform data, test set start position and end position
tm_start_GEANT = 2001;
end_time_GEANT = 3360;

% Customise the current axis position
x_axes_GEANT = [1:end_time_GEANT - tm_start_GEANT + 1];
% c = 800;

%==========================================================================
% Uniform data, Take the real predicted traffic data according to tm_start_GEANT and end_time_GEANT
real_od_GEANT = real_od_GEANT(:, tm_start_GEANT:end_time_GEANT);
RL_TM_GEANT = RL_TM_GEANT(:, tm_start_GEANT:end_time_GEANT);
PCA_TM_GEANT = PCA_TM_GEANT(:, tm_start_GEANT:end_time_GEANT);
SRSVD_TM_GEANT = SRSVD_TM_GEANT(:, tm_start_GEANT:end_time_GEANT);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Relative spatial error in GEANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rows of RL_TM_GEANT represent the number of the od stream and the columns represent the length of time
[r_GEANT, c_GEANT] = size(RL_TM_GEANT);

% The spatial relative error between the predicted flow size and the true flow size for the three methods is obtained
SRE_RL_TM_GEANT = function_SRE(RL_TM_GEANT, real_od_GEANT);
SRE_PCA_TM_GEANT = function_SRE(PCA_TM_GEANT, real_od_GEANT);
SRE_SRSVD_TM_GEANT = function_SRE(SRSVD_TM_GEANT, real_od_GEANT);

% Find [xmin, xmax] [ymin, ymax] range
len_GEANT = r_GEANT;
xmin_GEANT = 0; xmax_GEANT_GEANT = len_GEANT;
ymin_GEANT = 0; ymax_GEANT = norm(SRE_SRSVD_TM_GEANT, inf);
ymax_GEANT = ymax_GEANT + (ymax_GEANT / 10);

% Set the display font to Times New Roman and the font size to 14
figure(3);
subplot(2, 1, 1)
plot(SRE_RL_TM_GEANT, '-r', 'linewidth', 2);
hold on;
plot(SRE_PCA_TM_GEANT, '-b');
plot(SRE_SRSVD_TM_GEANT, '-g');

legend('Our method', 'PCA', 'SRMF')
xlabel( 'Flow ID', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel( 'SRE', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin_GEANT, xmax_GEANT_GEANT, ymin_GEANT, ymax_GEANT]);
title('SREs in GEANT')

% Obtain the average spatial relative error of the three methods
mean_SRE_RL = mean(SRE_RL_TM_GEANT);
% disp(SRE_RL_TM_GEANT);
mean_SRE_SRSVD = mean(SRE_SRSVD_TM_GEANT);
mean_SRE_PCA = mean(SRE_PCA_TM_GEANT);

%%%%%%%%%%%%%%%%%%%%%%%%% Temporal Relative Error in GEANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal Relative Error
[r_GEANT, c_GEANT] = size(real_od_GEANT);
x_axes_GEANT = 1:(end_time_GEANT - tm_start_GEANT + 1);

% The Temporal relative error between the predicted traffic size and the real traffic size for the three methods is obtained
TRE_RL_TM = function_TRE(RL_TM_GEANT, real_od_GEANT);
TRE_PCA_TM = function_TRE(PCA_TM_GEANT, real_od_GEANT);
TRE_SRSVD_TM = function_TRE(SRSVD_TM_GEANT, real_od_GEANT);

% Obtain the average Temporal relative error of the three methods
mean_TRE_RL = mean(TRE_RL_TM);
mean_TRE_SRSVD = mean(TRE_SRSVD_TM);
mean_TRE_PCA = mean(TRE_PCA_TM);

% Find [xmin_GEANT, xmax_GEANT_GEANT] [ymin_GEANT, ymax_GEANT] range
len_GEANT = c_GEANT;
xmin_GEANT = 1; xmax_GEANT_GEANT = 1 + len_GEANT;
ymin_GEANT = 0; ymax_GEANT = max(TRE_PCA_TM) * 2;
% ymax = max(ymax, 1.2);

% Graphs based on the various parameters calculated above
figure(3);
subplot(2, 1, 2)
plot(x_axes_GEANT, TRE_RL_TM, '-r', 'linewidth', 2);
hold on;
plot(x_axes_GEANT, TRE_PCA_TM, '-b');
plot(x_axes_GEANT, TRE_SRSVD_TM, '-g');
legend('Our method', 'PCA', 'SRMF')
xlabel('Time Slot Order', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('TRE','FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin_GEANT, xmax_GEANT_GEANT, ymin_GEANT, ymax_GEANT]);
% title('(d)', 'position', [23,-4], 'FontSize', 12)
title('TREs in GEANT')
%%%%%%%%%%%%%%%%%%%%%%    Spatial Relative Error CDF in GEANT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the cumulative distribution function of the spatial relative errors of the three methods
[errsp_vector_RL_TM, errsp_vector_RL_TM_cdf] = function_SRE_CDF(SRE_RL_TM_GEANT);
[errsp_vector_PCA_TM, errsp_vector_PCA_TM_cdf] = function_SRE_CDF(SRE_PCA_TM_GEANT);
[errsp_vector_SRSVD_TM, errsp_vector_SRSVD_TM_cdf] = function_SRE_CDF(SRE_SRSVD_TM_GEANT);
% Find [xmin_GEANT, xmax_GEANT_GEANT] [ymin_GEANT, ymax_GEANT] range
r1 = max(errsp_vector_RL_TM);
r2 = max(errsp_vector_PCA_TM);
r4 = max(errsp_vector_SRSVD_TM);

len_GEANT = max(r1, r2);
% len = max(len, r3);
% len = max(len, r4);
xmin_GEANT = 0; xmax_GEANT_GEANT = len_GEANT;
c1 = max(errsp_vector_RL_TM_cdf);
c3 = max(errsp_vector_PCA_TM_cdf);
len_GEANT = max(c1, c3);
ymin_GEANT = 0; ymax_GEANT = len_GEANT;

% As above, set the colour of the line graph, font size, etc.
figure(4);
subplot(2, 1, 1)
plot(errsp_vector_RL_TM, errsp_vector_RL_TM_cdf, '-r' , 'linewidth', 2);
hold on;
plot(errsp_vector_PCA_TM, errsp_vector_PCA_TM_cdf, '-b');
plot(errsp_vector_SRSVD_TM, errsp_vector_SRSVD_TM_cdf, '-g');

legend('Our method', 'PCA', 'SRMF')
xlabel('SRE', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('CDF', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin_GEANT, xmax_GEANT_GEANT, ymin_GEANT, ymax_GEANT]);
title('CDF of SRE in GEANT.')
%%%%%%%%%%%%%%%%%%%%%    Temporal Relative Error CDF in GEANT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the cumulative distribution function of the Temporal relative
% errors of the three methods in GEANT
[errsp_vector_val_RL_TM, errsp_vector_val_RL_TM_cdf] = function_TRE_CDF(TRE_RL_TM);
[errsp_vector_val_PCA_TM, errsp_vector_val_PCA_TM_cdf] = function_TRE_CDF(TRE_PCA_TM);
[errsp_vector_val_SRSVD_TM, errsp_vector_val_SRSVD_TM_cdf] = function_TRE_CDF(TRE_SRSVD_TM);

r1 = max(errsp_vector_val_RL_TM);
r3 = max(errsp_vector_val_PCA_TM);
len_GEANT = max(r1, r3);
xmin_GEANT = 0; xmax_GEANT_GEANT = len_GEANT * 1.5;
c1 = max(errsp_vector_val_RL_TM_cdf);
c3 = max(errsp_vector_val_PCA_TM_cdf);
len_GEANT = max(c1, c3);
ymin_GEANT = 0; ymax_GEANT = len_GEANT;
% Find [xmin_GEANT, xmax_GEANT_GEANT] [ymin_GEANT, ymax_GEANT] range

figure(4)
subplot(2, 1, 2)
plot(errsp_vector_val_RL_TM, errsp_vector_val_RL_TM_cdf, '-r', 'linewidth', 2);
axis([xmin_GEANT, xmax_GEANT_GEANT, ymin_GEANT, ymax_GEANT]);
hold on;
plot(errsp_vector_val_PCA_TM, errsp_vector_val_PCA_TM_cdf, '-b');
plot(errsp_vector_val_SRSVD_TM, errsp_vector_val_SRSVD_TM_cdf, '-g');

% Set the display font to Times New Roman and the font size to 14
legend('Our method', 'PCA', 'SRMF')
xlabel('TRE', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('CDF', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize',14)
set(0, 'defaultfigurecolor', 'w') 
title('CDF of TRE in GEANT.')

%==========================================================================


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is mainly used to calculate the relative error between the predicted and true values of the traffic matrix
%
% Creator:  Laisen Nie
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all
% clc

% Read the estimated data from the IoT network
load('./result/TM_IIoT.mat'); % real_od
real_od = TM_IIoT;

% RL method
load(sprintf('./result/TM_Pre_GAN_DQN_IPFP_IIoT.mat'));
RL_TM  = TM_Pre_GAN_DQN_IPFP_IIoT;

% PCA method
load(sprintf('./result/TM_PCA_Prediction_IPFP_IIoT.mat'));
PCA_TM = TM_PCA_Prediction_IPFP_IIoT;

% SRSVD method
load(sprintf('./result/TM_SRMF_Prediction_IPFP_IIoT.mat'));
SRSVD_TM = TM_SRMF_Prediction_IPFP_IIoT;

%==========================================================================

% Uniform data, test set start position and end position
tm_start = 1700;
end_time = 2015;

% Customise the current axis position
x_axes = [1:end_time - tm_start + 1];
% c = 800;

%==========================================================================
% Uniform data, Take the real predicted traffic data according to tm_start and end_time

real_od = real_od(:, tm_start:end_time);
RL_TM = RL_TM(:, tm_start:end_time);
PCA_TM = PCA_TM(:, tm_start:end_time);
SRSVD_TM = SRSVD_TM(:, tm_start:end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial relative error in IoT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rows of RL_TM represent the number of the od stream and the columns represent the length of time
[r, c] = size(RL_TM);

% The spatial relative error between the predicted flow size and the true flow size for the three methods is obtained
SRE_RL_TM = function_SRE(RL_TM, real_od);
SRE_PCA_TM = function_SRE(PCA_TM, real_od);
SRE_SRSVD_TM = function_SRE(SRSVD_TM, real_od);

% Set the range of the x and y axes of the spatial relative error figure for the three methods in the IOT network
len = r;
xmin = 0; xmax = len;
ymin = 0; ymax = norm(SRE_PCA_TM, inf);
ymax = max(SRE_RL_TM) + max(SRE_RL_TM) / 10;

figure(5);
subplot(2, 1, 1)
% Set the various label attribute information for the figure SREs in IoT
plot(SRE_RL_TM, '-r', 'linewidth', 2);
hold on;
plot(SRE_PCA_TM, '-b');
plot(SRE_SRSVD_TM, '-g');

legend('Our method', 'PCA', 'SRMF')

ylabel('SRE', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin, xmax, ymin, ymax]);
xlabel('Flow ID', 'FontName', 'Times New Roman', 'FontSize', 14);
title('SREs in IoT')
mean_SRE_RL = mean(SRE_RL_TM)
mean_SRE_SRSVD = mean(SRE_SRSVD_TM)
mean_SRE_PCA = mean(SRE_PCA_TM)

%%%%%%%%%%%%%%%%%%%%%%%%% Temporal Relative error in IoT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal Relative error
% The rows of real_od represent the number of the od stream and the columns represent the length of time
[r, c] = size(real_od);

% Customise the current axis position
x_axes = 1:(end_time - tm_start + 1);

TRE_RL_TM = function_TRE(RL_TM, real_od);
TRE_PCA_TM = function_TRE(PCA_TM, real_od);
TRE_SRSVD_TM = function_TRE(SRSVD_TM, real_od);

% Find [xmin, xmax] [ymin, ymax] range
len = c;
xmin = 1; xmax = 1 + len;
ymin = 0; ymax = max(TRE_PCA_TM) * 1.5;
% ymax = max(ymax, 1.2);

% Obtain the average Temporal relative error of the three methods
mean_TRE_RL = mean(TRE_RL_TM)
mean_TRE_SRSVD = mean(TRE_SRSVD_TM)
mean_TRE_PCA = mean(TRE_PCA_TM)

figure(5);
subplot(2, 1, 2)
% Set the display font to Times New Roman and the font size to 14
plot(x_axes, TRE_RL_TM, '-r', 'linewidth', 2);
hold on;
plot(x_axes, TRE_PCA_TM, '-b');
plot(x_axes, TRE_SRSVD_TM, '-g');
legend('Our method', 'PCA', 'SRMF')
xlabel( 'Time Slot Order', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel( 'TRE', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
axis([xmin, xmax, ymin, ymax]);
title('TREs in IoT')

%%%%%%%%%%%%%%%%%%%%%%    Spatial Relative Error CDF in IoT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculation of the cumulative distribution function of the spatial relative errors of the three methods
[errsp_vector_RL_TM, errsp_vector_RL_TM_cdf] = function_SRE_CDF(SRE_RL_TM);
[errsp_vector_PCA_TM, errsp_vector_PCA_TM_cdf] = function_SRE_CDF(SRE_PCA_TM);
[errsp_vector_SRSVD_TM, errsp_vector_SRSVD_TM_cdf] = function_SRE_CDF(SRE_SRSVD_TM);

% Find [xmin, xmax] [ymin, ymax] range
r1 = max(errsp_vector_RL_TM);
r2 = max(errsp_vector_PCA_TM);
r4 = max(errsp_vector_SRSVD_TM);

len = max(r1, r2);
% len = max(len, r3);
% len = max(len, r4);
xmin = 0; xmax = len;
c1 = max(errsp_vector_RL_TM_cdf);
c3 = max(errsp_vector_PCA_TM_cdf);
len = max(c1, c3);
ymin = 0; ymax = len;

figure(6);
subplot(2, 1, 1);
% Set the display font to Times New Roman and the font size to 14
plot(errsp_vector_RL_TM, errsp_vector_RL_TM_cdf, '-r' ,'linewidth', 2);
hold on;
plot(errsp_vector_PCA_TM, errsp_vector_PCA_TM_cdf, '-b');
plot(errsp_vector_SRSVD_TM, errsp_vector_SRSVD_TM_cdf, '-g');
legend('Our method', 'PCA', 'SRMF')
xlabel('SRE', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('CDF', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca,'FontName','Times New Roman','FontSize',14)
axis([xmin, xmax, ymin, ymax]);
title('CDF of SRE in IoT')
%%%%%%%%%%%%%%%%%%%%%     Temporal Relative Error CDF in IoT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the cumulative distribution function of the Temporal relative
% errors of the three methods in IoT
[errsp_vector_val_RL_TM, errsp_vector_val_RL_TM_cdf] = function_TRE_CDF(TRE_RL_TM);
[errsp_vector_val_PCA_TM, errsp_vector_val_PCA_TM_cdf] = function_TRE_CDF(TRE_PCA_TM);
[errsp_vector_val_SRSVD_TM, errsp_vector_val_SRSVD_TM_cdf] = function_TRE_CDF(TRE_SRSVD_TM);

% Find [xmin, xmax] [ymin, ymax] range
r1 = max(errsp_vector_val_RL_TM);
r3 = max(errsp_vector_val_PCA_TM);
len = max( r1, r3);
xmin = 0; xmax = len * 1.5;
c1 = max(errsp_vector_val_RL_TM_cdf);
c3 = max(errsp_vector_val_PCA_TM_cdf);
len = max(c1, c3);
ymin = 0; ymax = len;

figure(6);
subplot(2, 1, 2);
% As above, set the color of the line graph, font size, etc.
plot(errsp_vector_val_RL_TM, errsp_vector_val_RL_TM_cdf, '-r', 'linewidth',2);
axis([xmin, xmax, ymin, ymax]);
hold on;
plot(errsp_vector_val_PCA_TM, errsp_vector_val_PCA_TM_cdf, '-b');
plot(errsp_vector_val_SRSVD_TM, errsp_vector_val_SRSVD_TM_cdf, '-g');
disp(errsp_vector_val_RL_TM);
disp(errsp_vector_val_RL_TM_cdf);
disp(errsp_vector_val_PCA_TM);
disp(errsp_vector_val_PCA_TM_cdf);
disp(errsp_vector_val_SRSVD_TM);
disp(errsp_vector_val_SRSVD_TM_cdf);

legend('Our method', 'PCA', 'SRMF')
xlabel('TRE', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('CDF', 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
set(0, 'defaultfigurecolor', 'w') 
title('CDF of TRE in IoT')

% set(gcf, 'position', [500, 500, 500, 500])
