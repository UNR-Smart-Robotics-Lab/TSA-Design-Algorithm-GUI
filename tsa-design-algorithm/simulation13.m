%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the ICRA code, this file was named "Pt1_main2.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all; format compact; format long;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   378   560   380]);
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth',3);

%% Input Parameters
dt = 2;
dX = 0.1;
Fz0 = 10;
Fz_MAX = 25;
N = 2;
strain_ideal = 0.3;
FP = 1;
K__ = 500;
b=60;
tau=0.1;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_ = [bX, bO, Kr, J, Fc, Tc, M];
[R0_, string_mat, K_] = getR(K__);
R0 = R0_(1)/1000;
K = K_(1);

ideal_vol=5000;
max_vol=Inf;
Input_Voltage=0;
Price_Lim=Inf;
Mass_Lim=Inf;
algo_weights = [0.5,0.5,0.5];
%% Run the algorithm
[L0,X0,~]=dX2X0(dX,Fz0,K,N,strain_ideal);
[dOdt,Tm,tf,P,s,Xdiff,tEnd] = ...
coupledProps3(R0,Fz0,Fz_MAX,FP,N,X0,L0,K,dt,dX,b,tau,dynamic_constants_);
conv_len = length(Xdiff);
T = max(Tm)*5;
W1=(max(dOdt)/(max(Tm) - T))*(0 - T)+0;
W1=W1*60/2/pi;
T = T*10.197162129779; % convert (N*m) to (kg*cm)
[item_final, vendor, act_torque,torque,RPM,min_dist,act_vol,act_mass,act_NLRPM,act_price,act_voltage,act_type] = ...
    euclid_dist_fnctn2(W1,T,ideal_vol,max_vol,Input_Voltage,Price_Lim,Mass_Lim,algo_weights);

% 4841: T = 1.7; W = 2200
% 4757:  T = 3, W = 1600

%% Plot the results
figure
plot(torque, RPM,'o','LineWidth',3); hold on;
plot(T, W1, '*','MarkerSize',25/2);
plot(1.7, 2200, 'd','MarkerSize',20/2);
plot(3, 1600, 's','MarkerSize',20/2)
taumax = 15;
omeganlmax = 10000;
xlim([0, taumax])
ylim([0, omeganlmax])
xlabel('Stall Torque, $\tau_s$ [kg$\cdot$cm]')
ylabel('No-load Speed, $\omega_{NL}$ [RPM]')

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
plot(max(Tm)*10.197162129779, max(dOdt)*60/2/pi,'s','MarkerSize',15/2);
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

