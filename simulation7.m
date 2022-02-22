clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);
dt = 2;
dX = 0.1;
Fz0 = 10;
Fz_MAX = 25;
Si = 1000;
N = 2;
strain_ideal = 0.3;
R0 = 0.001;
FP = 1;
[L0,X0, K] = dX2X0(dX, Fz0, Si, N,strain_ideal);
b = [1,100,200,300,400,500,600,700];
m = length(b);
l = zeros(m,1);
tEnd = zeros(m,1);
tau = 0.1;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_ = [bX, bO, Kr, J, Fc, Tc, M];
% Compute dOdt
% fig = figure;
% for i = 1:length(b)
% [dOdt, Tm, tf, P, s, Xdiff{i}, tEnd(i)] = ...
%     coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt, dX, b(i),tau,dynamic_constants_);
%     l(i) = length(Xdiff{i});
%     if b(i) < 400
%         s1 = figure(1);
%     else
%         s2 = figure(2);
%     end
%     g{i} = stairs((1:l(i)), Xdiff{i}); 
%     g{i}.LineWidth = 2;
%     hold on;
% end
% save('sim7.mat');'=
load('sim7.mat');
close all;
clearvars -except g Xdiff l
%% Figure 1
figure(1)
for i = 1:4
    g{i} = stairs(1:l(i), Xdiff{i}); 
    hold on;
    g{i}.LineWidth = 2;
end
leg1 = legend('$\rho = 1$', '$\rho = 100$', '$\rho = 200$', '$\rho = 300$');
leg1.Interpreter = 'latex';
leg1.NumColumns = 1;
set(gca,'XScale','log');
xlabel('Current Iteration, $l$');
ylabel('$X_0 - X(\mathrm{end}) - \Delta X$ [m]');
g{1}.Marker = '+';
g{2}.Marker = '*';
g{3}.Marker = '^';
g{4}.Marker = 'd'; 
g{1}.LineStyle = '-';
g{2}.LineStyle = ':';
g{3}.LineStyle = '-.';
g{4}.LineStyle = '-.';
%% Figure 2
figure(2)
for i = 8:-1:5
    g{i} = stairs((1:l(i)), Xdiff{i});
    g{i}.LineWidth = 10-i;
    hold on;
end
set(gca,'XScale','log');
leg2 = legend('$\rho = 700$', '$\rho = 600$', '$\rho = 500$', '$\rho = 400$');
xlabel('Current Iteration, $l$');
ylabel('$X_0 - X(\mathrm{end}) - \Delta X$ [m]');
xlim([0,100]);
leg2.Interpreter = 'latex';
leg2.NumColumns = 4;
g{5}.LineStyle = '-';
g{6}.LineStyle = '-.';
g{7}.LineStyle = '--';
g{8}.LineStyle = ':';
g{5}.Marker = '+';
g{6}.Marker = '*';
g{7}.Marker = '^';
g{8}.Marker = 'd';