%% Script Setup
clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);
%% Define our inputs and simulation parameters
n = 1000;
m = 100;
b = linspace(1,1000,m);
Fz0 = rand(n,1)*20;
Fz_MAX = Fz0 + rand(n,1)*10;
K__ = rand(n,1)*(500 - 100)+100;
N = 2;
strain_ideal = 0.3;
FP = randi(3,n,1);
dt = (rand(n,1)*(10-0.2) + 0.2);
dX = (rand(n,1)*(0.3-0.1) + 0.1);
%% Initialize Arrays
conv_len = zeros(n,1);
tEnd = zeros(n,1);
Xdiff = cell(n);
mean_conv_len = zeros(n,1);
formatSpec = 'dt: %4.2f, dX: %4.3f, Fz0: %4.3f, Fz_MAX: %4.3f, K: %4.3f, FP: %4.3f, i: %4.3f rho: %4.3f \n';
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_=[bX, bO, Kr, J, Fc, Tc, M];
tau=0.1;
L0=zeros(n,1);
X0=zeros(n,1);
K=zeros(n,1);
W1=zeros(n,1);
T=zeros(n,1);
dOdt = cell(n,1);
Tm=cell(n,1);
P=cell(n,1);
tf=zeros(n,1);
s=cell(n,1);
R0 = zeros(n,1);
%% Run the simulation
for i = 1:n
    [R0_, ~, K_] = getR(K__(i));
    R0(i) = R0_(1)/1000;
    K(i) = K_(1);
    [L0(i),X0(i), ~] = dX2X0(dX(i), Fz0(i), K(i), N,strain_ideal);
    for j = 1:m
        fprintf(formatSpec,dt(i),dX(i),Fz0(i),Fz_MAX(i),K(i),FP(i),i,b(j));
        [~, ~, tf(i,j), ~, s{i,j}, Xdiff{i,j}, tEnd(i,j)] = ...
            coupledProps3(R0(i), Fz0(i), Fz_MAX(i), FP(i), N, X0(i), L0(i), K(i), dt(i), dX(i), b(j),tau,dynamic_constants_);
        conv_len(i,j) = length(Xdiff{i,j});
        %T(i,j) = max(Tm{i,j})*5;
        %W1(i,j) = (max(dOdt{i,j})/(max(Tm{i,j}) - T(i,j)))*(0 - T(i,j)) + 0;
        %W1(i,j) = W1(i,j)*60/2/pi;
        %T(i,j) = T(i,j)*10.197162129779; % convert (N*m) to (kg*cm), source: https://www.convertunits.com/from/N-m/to/kg-cm
    end
end
%% Save the data
save(['MAT Files/sim10--',num2str(round(now*100000)),'.mat'])