clear; clc; close all; format compact; format long;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   378   560   380]);
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
% Inputs
dt = 2;
dX = 0.1;
Fz0 = 10;
Fz_MAX = 25;
N = 2;
strain_ideal = 0.3;
R0 = 0.001;
FP = 1;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
dynamic_constants_ = [bX, bO, Kr, J, Fc, Tc, M];
rho = [0.01:0.01:0.09,0.1:0.1:0.9,1:700];
Si = 1000;
[L0,X0,K] = dX2X0(dX, Fz0, Si, N,strain_ideal);
tau = 0.1;
% for i = 1:length(rho)
% [dOdt, Tm, tf, P, s, Xdiff{i}, tEnd(i)] = ...
%     coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt, dX, rho(i),tau,dynamic_constants_);
%     l{1}(i) = length(Xdiff{i});
%     disp(rho(i));
%     disp(tEnd(i));
% end
% save('sim5_data.mat');

load('sim5_data.mat');

figure
yyaxis left
semilogy(rho,l{1},'-','LineWidth',3);
ylabel('Iterations, $l$');
yyaxis right
semilogy(rho,tEnd,'-');
ylabel('Computation time, $t_c$ (s)');
xlabel('Step Size Gain, $\rho$ (rad/m/s)');
xlim([-100,850])