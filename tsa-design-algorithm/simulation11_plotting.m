%% Setup
clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);
load('MAT Files/sim11738576599664.mat');

fig=figure;
subplot(2,1,1);
histogram(T,'Normalization','probability');
title('Random-Varying Force Profiles, $n=1,000$');
xlabel('Min. Required Stall Torque');
subplot(2,1,2);
histogram(W1,'Normalization','probability');
xlabel('Min. Required Free-Run Speed');
han=axes(fig,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Relative Frequency');

A1 = categorical(top_item);
B1 = categories(A1);
B = categorical(top_vendor);

% figure
% histogram(top_min_dist);

figure
subplot(2,1,1);
h1 = boxplot(T,'Orientation','horizontal');
set(findobj(gca,'type','line'),'linew',2)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
xlabel('Stall Torque, $\tau_{s,d}$')
sgtitle('Random-Varying Force Profiles, $n=1,000$','interpreter','latex');

subplot(2,1,2)
h2 = boxplot(W1,'Orientation','horizontal');
set(findobj(gca,'type','line'),'linew',2)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
xlabel('Free-Run Speed, $\omega_{NL,d}$')


