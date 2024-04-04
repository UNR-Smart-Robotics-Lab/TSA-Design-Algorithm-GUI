clear; clc; close all;
load('sim5_data.mat');
figure
yyaxis left
semilogy(rho,l{1},'-','LineWidth',3);
ylabel('Iterations, $l$');
yyaxis right
semilogy(rho,tEnd,'-.');
ylabel('Computation time, $t_c$ (s)');
xlabel('Step Size Gain, $\rho$ (rad/m/s)');
xlim([-100,600]);