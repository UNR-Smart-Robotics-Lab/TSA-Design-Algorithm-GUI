% David Bombara
% July 4, 2022
% Simulation 14

clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);
%%
m_vec = [100,200,300,400,500,1000,1500,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000]';

ddotX_uf_ = cell(length(m_vec),1);
ddotX_ = cell(length(m_vec),1);
tf_vec_ = cell(length(m_vec),1);
for iii = 1:length(m_vec)
    str = ['designAlgoResults/design_algorithm_results_',num2str(m_vec(iii)),'_.mat'];
    load(str);
    tf_vec_{iii,1} = tf_vec;
    ddotX_uf_{iii,1} = ddotX_uf;
    ddotX_{iii,1} = ddotX;
end

%% Figure 1: 3D Plot
figure
for iii = 1:length(m_vec)
    plot3(m_vec(iii)*ones(length(tf_vec_{iii,1}),1),tf_vec_{iii,1}, ddotX_{iii,1})
    hold on
end
xlabel('Discretization Level, $m$')
ylabel('Time, $t$ [s]')
zlabel('Linear Acceleration, $\ddot{X}$ [m/s$^2$]')

%% Figure 2: 2D Plot
figure
for iii = [1,6,17]
    switch iii
        case 1
            lnstl = '-';
        case 6
            lnstl = '--';
        case 17
            lnstl = ':';
    end
    plot(tf_vec_{iii,1}, ddotX_{iii,1},'LineStyle',lnstl,"LineWidth",3)
    hold on;
end
xlabel('Time, $t$ [s]')
ylabel('Linear Acceleration, $\ddot{X}$ [m/s$^2$]');
legend('$m=100$','$m=1000$','$m=10000$')