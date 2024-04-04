clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);
load('sim7.mat');
close all;
clearvars -except g Xdiff l b tEnd
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
for i = 9:-1:6
    switch i
        case 9
            figure(3)
            
        otherwise 
            figure(2)
    end
    g{i} = stairs((1:l(i)), Xdiff{i});
    if i == 9
        xlabel('Current Iteration, $l$');
        ylabel('$X_0 - X(\mathrm{end}) - \Delta X$ [m]');
        g{i}.Color = [0.4940, 0.1840, 0.5560];
        legend('$\rho = 10,000$');
    end
    g{i}.LineWidth = 2;
    hold on;
end
set(gca,'XScale','log');
leg2 = legend('$\rho = 500$', '$\rho = 100$', '$\rho = 10$');
xlabel('Current Iteration, $l$');
ylabel('$X_0 - X(\mathrm{end}) - \Delta X$ [m]');
xlim([0,100]);
leg2.Interpreter = 'latex';
leg2.NumColumns = 4;
g{9}.LineStyle = '-';
g{6}.LineStyle = '-.';
g{7}.LineStyle = '--';
g{8}.LineStyle = ':';
g{9}.Marker = '+';
g{6}.Marker = '*';
g{7}.Marker = '^';
g{8}.Marker = 'd';