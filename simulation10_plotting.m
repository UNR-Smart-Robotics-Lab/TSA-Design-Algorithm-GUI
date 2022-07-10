clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   420]);

load('MAT Files/sim10--73857666235.mat');
I1 = zeros(n,1);
I2 = zeros(n,1);
B1 = zeros(n,1);
B2=B1;
mintEnd = zeros(n,1);
minConvLen=zeros(n,1);
meantEnd = zeros(n,1);
meantEnd_rho = zeros(n,1);
maxtEnd_rho = zeros(n,1);
Bpass = zeros(n,m);
% Rows: differnet inputs
% Colums: different rho values
for i = 1:n 
    [mintEnd(i),I1(i)] = mink(tEnd(i,:),1);
    [minConvLen(i),I2(i)] = mink(conv_len(i,:),1);
    B1(i) = b(I1(i));
    B2(i) = b(I2(i));
    meantEnd(i) =  mean(tEnd(i,:));
end

% Psedocode
% for every b inputs
% count the number of inputs for each the algorithm converges
% whichever b value converges the most, you should use during the next
% simulations
ct = zeros(length(b),1);
for i = 1:length(b)
    tmp = tEnd(:,i);
    ct(i) = length(tmp(tmp<0.05));
    
    
end

figure
semilogx(b,ct./n);
xlabel('Step Size Gain, $\rho$');
ylabel('Probability of $t_c \leq 0.05$\,s');

figure
histogram(B1,20,'Normalization','probability');
xlabel('Optimal $\rho$');
ylabel('Relative Freqency');

% figure
% histogram(B2,20,'Normalization','countdensity');
% xlabel('$\rho$');
% ylabel('Relative Freqency');

figure
plot(dt,B2,'.');

figure
subplot(3,1,1);
plot(dX,B2,'.');
xlabel('Contraction Range, $\Delta X$ (m)');
ylabel('Optimal $\rho$');

subplot(3,1,2);
plot(Fz0, B2,'.');
ylabel('Optimal $\rho$');
xlabel('Min. Force, $F_{z0}$ (N)');

subplot(3,1,3);
plot(Fz_MAX,B2,'.');
ylabel('Optimal $\rho$');
xlabel('Max. Force, $F_{z,max}$ (N)');

figure
subplot(2,1,1);
plot(FP,B2,'.');
ylabel('Optimal $\rho$');
xlabel('Force Profile');

subplot(2,1,2);
plot(K,B2,'.');
ylabel('Optimal $\rho$');
xlabel('Norm. Stiffness, $K$ (N)');

figure
plot(dX, mintEnd,'.');
xlabel('Contraction Range, $\Delta X$ (m)');
ylabel('Computation time at the optimal $\rho$ value')

figure
plot(dX, meantEnd,'.');
xlabel('Contraction Range, $\Delta X$ (m)');
ylabel('Average computation time for all $\rho$ values')

% figure
% plot(b,meantEnd_rho,'.');
% xlabel('$\rho$');
% ylabel('Mean. Computation Time');
% 
% figure
% plot(b,maxtEnd_rho,'.');
% xlabel('$\rho$');
% ylabel('Max. Computation Time');