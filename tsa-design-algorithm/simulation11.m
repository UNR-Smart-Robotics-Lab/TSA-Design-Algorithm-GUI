%% Simulation 11: Testing the Random Fz Values to see if they affect the motor selection
%% Formatting and Setup
clear; clc; close all; format compact; format long;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   380]);
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
%% Parameters
dt = 2;
dX = 0.1;
Fz0 = 10;
Fz_MAX = 25;
K = 500;
N = 2;
strain_ideal = 0.3;
FP = 3; % random
[R0, ~, K_] = getR(K);
R0 = R0(1)/1000;
b=200; % step size gain
tau=0.1;
n = 1000; % number of trials in the simulation
ideal_vol=0;
max_vol = Inf;
Input_Voltage = 0;
Price_Lim = Inf;
Mass_Lim = Inf;
min_volume = 0;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_ = [bX, bO, Kr, J, Fc, Tc, M];
[L0,X0, ~] = dX2X0(dX,Fz0,K_(1),N,strain_ideal); % does not depend on the force profile
%% Initialize Arrays
W1=zeros(n,1);
T=zeros(n,1);
dOdt = cell(n,1);
Tm=cell(n,1);
P=cell(n,1);
tf = zeros(n,1);
s=cell(n,1);
top_item = cell(n,1);
top_vendor = cell(n,1);
top_act_torque = zeros(n,1);
top_min_dist=zeros(n,1);
top_act_vol=zeros(n,1);
top_act_mass=zeros(n,1);
top_act_NLRPM=zeros(n,1);
top_act_price=zeros(n,1);
top_act_voltage=zeros(n,1);
top_act_type=cell(n,1);
conv_len = zeros(n,1);
Xdiff = cell(n,1);
tEnd = zeros(n,1);
formatSpec = 'dt: %4.2f, dX: %4.3f, Fz0: %4.3f, Fz_MAX: %4.3f, Si: %4.3f, FP: %4.3f, i: %4.3f, rho: %4.3f, tEnd: %4.3f \n';
%% Simulation
for i = 1:n
        [dOdt{i}, Tm{i}, tf(i), P{i}, s{i}, Xdiff{i}, tEnd(i)] = ...
        coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt, dX, b,tau,dynamic_constants_);
    conv_len(i) = length(Xdiff{i});
    T(i) = max(Tm{i})*5;
    W1(i) = (max(dOdt{i})/(max(Tm{i}) - T(i)))*(0 - T(i)) + 0;
    W1(i) = W1(i)*60/2/pi;
    T(i) = T(i)*10.197162129779; % convert (N*m) to (kg*cm)
    [item_final, vendor, act_torque,~,~, min_dist,act_vol,act_mass,act_NLRPM,act_price,act_voltage,act_type] = ...
    euclid_dist_fnctn2(W1(i), T(i), ideal_vol, max_vol, Input_Voltage,Price_Lim,Mass_Lim);
    top_item{i}=item_final{1};
    top_vendor{i}=vendor{1};
    top_act_torque(i)=act_torque{1};
    top_min_dist(i)=min_dist(1);
    top_act_vol(i)=act_vol{1};
    top_act_mass(i)=act_mass{1};
    top_act_NLRPM(i)=act_NLRPM{1};
    top_act_price(i)=act_price{1};
    top_act_voltage(i)=act_voltage{1};
    top_act_type{i}=act_type{1};
    fprintf(formatSpec,dt,dX,Fz0,Fz_MAX,K,FP,i,b,tEnd(i));
end
save(['MAT Files/sim11',num2str(round(now*1000000)),'.mat']);