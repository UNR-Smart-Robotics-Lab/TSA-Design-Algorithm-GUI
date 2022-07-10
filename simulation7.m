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
K = 500;
N = 2;
strain_ideal = 0.3;
FP = 1;
[R0, ~, K_] = getR(K);
R0 = R0(1)/1000;
[L0,X0, ~] = dX2X0(dX, Fz0, K_(1), N,strain_ideal);
b = [10^2,10^2,10^2,10^2,1,10,100,500,10000];
m = length(b);
l = zeros(m,1);
tEnd = zeros(m,1);
tau = 0.1;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_ = [bX, bO, Kr, J, Fc, Tc, M];
for i = 1:length(b)
[dOdt, Tm, tf, P, s, Xdiff{i}, tEnd(i)] = ...
    coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K_(1), dt, dX, b(i),tau,dynamic_constants_);
    l(i) = length(Xdiff{i});
    disp(tEnd(i));
    disp(i);
end
save('sim7.mat');
