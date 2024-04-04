%% Simulation 12: Checking how the Euclidean distance changes based on different perturbations
%% Formatting and Setup
clear; clc; close all; format compact; format long;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   380]);
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
%% Define our initial inputs
%% Parameters
dt = 2;dX = 0.1;Fz0 = 10;Fz_MAX = 25;
K = 500;N = 2;strain_ideal = 0.3;
FP = 1; % random
[R0, ~, K_] = getR(K);
R0 = R0(1)/1000;
b=60; % step size gain
tau=0.1;
ideal_vol=0;
max_vol = Inf;
Input_Voltage = 0;
Price_Lim = Inf;
Mass_Lim = Inf;
min_volume = 0;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_ = [bX, bO, Kr, J, Fc, Tc, M];
[L0,X0, ~] = dX2X0(dX,Fz0,K_(1),N,strain_ideal); % does not depend on the force profile
%var_vec = (0.95:0.001:1.05);
f = 0.5; % 50 percent
var_vec = linspace(1-f,1+f,100);
n=length(var_vec);

%% Initializing Arrays
X = cell(n,5);
Y = zeros(n,5);
dOdt = X;P = X;Tm = X;s = X;top_item=X;top_vendor=X;top_act_type=X;Xdiff = X;
top_act_torque=Y;top_min_dist=Y;top_act_vol=Y;top_act_mass=Y;top_act_NLRPM=Y;top_act_price=Y;
top_act_voltage=Y;tf = Y;tEnd = Y;conv_len = Y;T=Y;W1=Y;
%% Simulation 1: Vary dt
dt_var = dt*var_vec;
%dt_var = linspace(0.2,10,n);
j=1;
for i = 1:n
    [dOdt{i,j}, Tm{i,j}, tf(i,j), P{i,j}, s{i,j}, Xdiff{i,j}, tEnd(i,j)] = ...
        coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt_var(i), dX, b,tau,dynamic_constants_);
    disp(i);
end
disp('dt done');
%% Simulation 2: Vary dX
dX_var = dX*var_vec;
%dX_var = linspace(0.1,0.3,n);
j=2;
for i = 1:n
    [dOdt{i,j}, Tm{i,j}, tf(i,j), P{i,j}, s{i,j}, Xdiff{i,j}, tEnd(i,j)] = ...
        coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt, dX_var(i), b,tau,dynamic_constants_);
    disp(i);
end
disp('dX done')
%% Simulation 3: Vary Fz0
Fz0_var = Fz0*var_vec;
%Fz0_var = linspace(0,20,n);
j=3;
for i = 1:n
    [dOdt{i,j}, Tm{i,j}, tf(i,j), P{i,j}, s{i,j}, Xdiff{i,j}, tEnd(i,j)] = ...
        coupledProps3(R0, Fz0_var(i), Fz_MAX, FP, N, X0, L0, K, dt, dX, b,tau,dynamic_constants_);
    disp(i);
end
disp('Fz0 done')
%% Simulation 4: Vary Fz_{max}
Fz_MAX_var = Fz_MAX*var_vec;
%Fz_MAX_var = Fz0_var + linspace(0,10,n);
j=4;
for i = 1:n
    [dOdt{i,j}, Tm{i,j}, tf(i,j), P{i,j}, s{i,j}, Xdiff{i,j}, tEnd(i,j)] = ...
        coupledProps3(R0,Fz0,Fz_MAX_var(i),FP,N,X0,L0,K,dt,dX,b,tau,dynamic_constants_);
    disp(i);
end
disp('Fz_max done')
%% Simulation 5: Vary K
K_var = K*var_vec;
%K_var = linspace(100,500,n);
j=5;
for i = 1:n
    [dOdt{i,j}, Tm{i,j}, tf(i,j), P{i,j}, s{i,j}, Xdiff{i,j}, tEnd(i,j)] = ...
        coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K_var(i), dt, dX, b,tau,dynamic_constants_);
    disp(i);
end
disp('K done')
%% Do some computations on all the data
for j = 1:5
    for i = 1:n
        conv_len(i,j) = length(Xdiff{i,j});
        T(i,j) = max(Tm{i,j})*5;
        W1(i,j) = (max(dOdt{i,j})/(max(Tm{i,j}) - T(i,j)))*(0 - T(i,j)) + 0;
        W1(i,j) = W1(i,j)*60/2/pi;
        T(i,j) = T(i,j)*10.197162129779; % convert (N*m) to (kg*cm)
        [item_final, vendor, act_torque,~,~, min_dist,act_vol,act_mass,act_NLRPM,act_price,act_voltage,act_type] = ...
            euclid_dist_fnctn2(W1(i,j), T(i,j),ideal_vol,max_vol,Input_Voltage,Price_Lim,Mass_Lim,[0.5,0.5,0.5]);
        top_item{i,j}=item_final{1};
        top_vendor{i,j}=vendor{1};
        top_act_torque(i,j)=act_torque{1};
        top_min_dist(i,j)=min_dist(1);
        top_act_vol(i,j)=act_vol{1};
        top_act_mass(i,j)=act_mass{1};
        top_act_NLRPM(i,j)=act_NLRPM{1};
        top_act_price(i,j)=act_price{1};
        top_act_voltage(i,j)=act_voltage{1};
        top_act_type{i,j}=act_type{1};
        disp(i);
    end
end
save(['MAT Files/sim12--',num2str(round(now*1000000)),'.mat']);