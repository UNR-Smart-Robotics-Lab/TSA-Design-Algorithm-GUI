%% Simulation 12: Checking how the Euclidean distance changes based on different perturbations
%% Formatting and Setup
clear; clc; close all; format compact; format long;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   278   560   380]);
set(0,'defaultAxesFontSize',16);
set(0, 'DefaultLineLineWidth', 2);
load('MAT Files/sim12--738581659454.mat');

%% dt
figure;
yyaxis left
plot(dt_var,T(:,1));
ylabel('Req. Stall Torque (kg$\cdot$cm)');

yyaxis right
plot(dt_var,W1(:,1),'--');
xlabel('Contraction Time, $\Delta t$ (m)')
ylabel('Req. Free-Run Speed (RPM)');
legend('Stall Torque', 'Free-Run Speed')

%% dX
figure
yyaxis left
plot(dX_var,T(:,2));
ylabel('Req. Stall Torque (kg$\cdot$cm)');

yyaxis right
plot(dX_var,W1(:,2),'--');
xlabel('Contraction Range, $\Delta X$ (m)')
ylabel('Req. Free-Run Speed (RPM)');
legend('Stall Torque', 'Free-Run Speed')

%% F_{z0}
fig = figure;
yyaxis left
plot(Fz0_var,T(:,3));
ylabel('Req. Stall Torque (kg$\cdot$cm)');

yyaxis right
plot(Fz0_var,W1(:,3),'--');
xlabel('Min. Axial Force, $F_{z,0}$');
ylabel('Req. Free-Run Speed (RPM)');
legend('Stall Torque', 'Free-Run Speed','location','best')

figure
yyaxis left
plot(Fz_MAX_var,T(:,4));
ylabel('Req. Stall Torque (kg$\cdot$cm)');

yyaxis right
plot(Fz_MAX_var,W1(:,4),'--');
xlabel('Max. Axial Force, $F_{z,max}$');
ylabel('Req. Free-Run Speed (RPM)');
legend('Stall Torque', 'Free-Run Speed','location','best')

figure
yyaxis left
plot(K_var,T(:,5));
ylabel('Req. Stall Torque (kg$\cdot$cm)');

yyaxis right
plot(K_var,W1(:,5),'--');
xlabel('Normalized Stiffness, $K$ (N)');
ylabel('Req. Free-Run Speed (RPM)');
legend('Stall Torque', 'Free-Run Speed','location','best')

%% Bar Graph
A = {'$\Delta t$', '$\Delta X$', '$F_{z,0}$', '$F_{z,max}$', '$K$'};
B = categorical(A);
B = reordercats(B,A);
x = zeros(1,5); a=x;b=x;c=x;d=x;e=x;f=x;
for i = 1:5
    a(i) = median(T(:,i));
    b(i) = min(T(:,i));
    c(i) = max(T(:,i));
    d(i) = median(W1(:,i));
    e(i) = min(W1(:,i));
    f(i) = max(W1(:,i));
end
% Data for the Upper
g(1,:)=(c-a)./a; % torque
g(2,:)=(f-d)./d; % speed
g(:,5)=-g(:,5);
g(:,1)=-g(:,1);
% Data for the lower
h(1,:)=((b-a)./a); % torque
h(2,:)=((e-d)./d); % speed
h(:,5)=-h(:,5);
h(:,1)=-h(:,1);
g=g';h=h';
figure
subplot(2,1,1)
bar(B,g*100)
title('$+50$\% Nominal Value');
ylabel('Percent Change')
legend('Stall Torque', 'No-Load Speed')
subplot(2,1,2);
bar(B,h*100)
title('$-50$\% Nominal Value');
ylabel('Percent Change')

legend('Stall Torque', 'No-Load Speed')




