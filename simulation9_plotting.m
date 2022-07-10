clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);
load('MAT Files/sim9--738577555451.mat');
T = real(T);
W1 = real(T);

figure
histogram(tEnd(tEnd<=1),40);
ylabel('Frequency');
xlabel('Compuation Time (s)');

figure
histogram(tEnd(tEnd>1),40);
xlabel('Compuation Time (s)');
ylabel('Frequency');

figure
subplot(2,1,1);
histogram(W1(W1<=50));
% subplot(2,1,2);
% histogram(W1(W1>50));
xlabel('Req. Free-Run Speed [RPM]');
ylabel('Frequency');

subplot(2,1,2);
histogram(top_act_NLRPM,40);
ylabel('Frequency');
xlabel('Selected Free-Run Speed [RPM]');

figure
subplot(2,1,1);
histogram(T(T<=50),40);
xlabel('Req. Stall Torque (kg$\cdot$cm)');
ylabel('Frequency');

subplot(2,1,2);
histogram(top_act_torque,40);
xlabel('Selected Stall Torque (kg$\cdot$cm)');
ylabel('Frequency');