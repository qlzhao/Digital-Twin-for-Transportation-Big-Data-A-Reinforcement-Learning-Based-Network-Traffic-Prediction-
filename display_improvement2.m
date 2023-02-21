% clear variables in base space, close all Figure windows, clear the
% contents of the command window.
clear
close all
clc
% ==========================================================================
% The ratio of the improvement rate between our method and the other two
% prediction methods is obtained for the three networks.
A = function_Pre_Abilene;
G = function_Pre_GEANT;
I = function_Pre_IoT;
% The first scenario is RL:PCA and the second is RL:SRMF
a = [A(1) A(2); G(1) G(2); I(1) I(2)];
disp(a);
% Drawings figure ==========================================================================
figure1 = figure();
c = categorical({'Abilene', 'GEANT', 'IoT'});
bar(c, a, 1, 'grouped')
legend('RL:PCA', 'RL:SRMF');
% Set the [xmin, xmax], [ymin, ymax] ranges and graph using a structure
xmin = 0.5;
xmax = 3.5;
ymin = -1;
ymax = 1;
% axis( [ xmin, xmax, ymin, ymax ] );
ylim([ymin, ymax]);
% set(gca, 'xTicklabel', {})
% xlabel( 'In Abilene, GEANT and IoT Network', 'FontSize', 14, 'FontName', 'Times New Roman' );
% Set the font for the x-axis and y-axis display to Times New Roman and the font size to 14
ylabel( 'Performance Improvement Ratio', 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
% set(gca, 'yTicklabel', {'0', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})

% ==========================================================================
