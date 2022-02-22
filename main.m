%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the ICRA code, this file was named "Pt1_main2.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all; format compact; format long;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   378   560   380]);
set(0,'defaultAxesFontSize',14);
set(0, 'DefaultLineLineWidth',3);
[dt, dX, Fz0, Fz_MAX, Si, N, R0, FP, strain_ideal] = sim_inputs();
% Same as P1_main except with NO OPTIMIZATION
[L0,X0,K] = dX2X0(dX, Fz0, Si, N,strain_ideal);
[dOdt, Tm, tf, P, s, Xdiff,O] = ...
    coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt, dX, 50);
T = max(Tm)*5; % the stall torque should not exceed 20% of the torque the
% TSA will require (conservative estimate since it's using the max torque)
% compute free-run speed based on the actual speed and stall torque
W1 = (max(dOdt)/(max(Tm) - T))*(0 - T) + 0;
% I think there's a mistake here
% W2 = dOdt/(1-T);
% assume a linear relationship between speed and torque.

W1 = W1*60/2/pi; % convert rad/s to RPM
% W2 = W2*60/2/pi;

T = T*10.197162129779; % convert (N*m) to (kg*cm);
% Source: https://www.convertunits.com/from/N-m/to/kg-cm

ideal_volume = 5000;
%max_volume = 100000;
max_volume = Inf;
min_volume = 5000;
% [item_final, price_final, vol_final] = euclid_dist_fnctn(W1, T, Input_Volume);
[item_final, vendor,act_torque, torque, RPM] = ...
    euclid_dist_fnctn2(W1, T, ideal_volume, max_volume,12);
item_final % Goal: 4698
% 4841: T = 1.7; W = 2200
% 4757:  T = 3, W = 1600
figure
plot(torque, RPM,'o','LineWidth',3); hold on;
% plot(3.5, 1600, 'o');
plot(T, W1, '*','MarkerSize',25);
plot(1.7, 2200, 'd','MarkerSize',20);
plot(3, 1600, 's','MarkerSize',20)
taumax = 15;
omeganlmax = 10000;
xlim([0, taumax])
ylim([0, omeganlmax])
xlabel('Stall Torque, $\tau_s$ [kg$\cdot$cm]')
ylabel('No-load Speed, $\omega_{NL}$ [RPM]')

% C = [0,0,0];
% ar2 = fill([T, T, 100, 100], [W1, 10000, 10000, W1], C);
% ar2.FaceAlpha = 0;
% ar2.FaceColor = [0.3010, 0.7450, 0.9330];
% ar2.FaceColor = [0,0,0];
% ar2.LineStyle = ':';
% ar2.LineWidth = 3; hold on;

p{1} = plot([T, taumax+10], [W1, W1],'Color',[0 0.5 0]); 
p{2} = plot([T, taumax+10], [omeganlmax+10, omeganlmax+10],'Color',[0 0.5 0]);
p{3} = plot([T, T],[W1, omeganlmax+10],'Color',[0 0.5 0]);
p{4} = plot([taumax+10,taumax+10],[W1, omeganlmax+10],'Color',[0 0.5 0]);

for i = 1:4
    p{i}.Color = [0.4660, 0.6740, 0.1880];
    p{i}.LineStyle = '-.';
end
leg = legend('Database Entries','($\tau_{s, des}, \omega_{NL, des}$)',...
    '$1^{\mathrm{st}}$ Motor Choice', '$2^{\mathrm{nd}}$ Motor Choice', 'Acceptable Region');
leg.NumColumns = 1;
leg.Interpreter = 'latex';

%%%%%%%%%%%%% 
% Figure 2
figure
% Curve 1 - our ideal values
plot(max(Tm)*10.197162129779, max(dOdt)*60/2/pi,'s','MarkerSize',15);
hold on;
plot([0, T], [W1, 0],'o-'); hold on;
%
% Curve 2 - the first choice motor
plot([0, 1.7], [2200, 0],'o--')
%
% Curve 3 - the second choice motor
plot([0, 3], [1600,0],'o-.');
%
leg2 = legend('$(\mathrm{max}(\tau_m), \dot{\theta}_p)$','Ideal Curve',...
    '1$^{\mathrm{st}}$ Motor Choice', '2$^{\mathrm{nd}}$ Motor Choice');
leg2.Interpreter = 'latex';
xlabel('Torque, $\tau$ [kg$\cdot$cm]');
ylabel('Angular Speed, $\omega$ [RPM]');