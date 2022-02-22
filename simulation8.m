% David Bombara
% Simulation to see how theta, theta_dot, X, and dXdt change with time.
clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   378   560   380]);
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 3);
[dt, dX, Fz0, Fz_MAX, Si, N, R0, FP, strain_ideal] = sim_inputs();

% Same as P1_main except with NO OPTIMIZATION
[L0,X0, K] = dX2X0(dX, Fz0, Si, N,strain_ideal);

m = 1000;
[bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();


Xmin = X0 - dX; % unloaded contraction range

tf_vec = (linspace(0, dt, m))';

b =531;

j = 1; % initialized loop counter

% Initial guess for Omax using variable radius, finite stiffness
[O, ~] = dX2O(X0, dX, Fz_MAX, R0, N, L0, K);
Omax = max(O); 
clear O;

% Initial guess for Omax using constant radius, infinite stiffness.
% Omax = sqrt(dX./R0^2);
% initial guess for dOdt
dOdt_p = Omax./dt;
tStart = tic;
s = zeros(10000, 1);
Xdiff = zeros(10000,1);
while 1
    % For constant dOdt:
    % dOdt = (Omax / dt)*ones(m,1); % for now, assume dOdt == constant
    % O = (linspace(0, Omax, m))';
    % O = dOdt*tf_vec;
    
    % Initializing Arrays
    X        = ones(m, 1) * X0;
    R        = ones(m, 1) * R0;
    Fi       = ones(m,1) * Fz0./N;
    dXdt     = zeros(m,1);
    F        = ones(m,1) * Fz0;
    Ffr      = zeros(m,1);
    dRdt_num = zeros(m,1);
    dRdt_den = zeros(m,1);
    dRdt     = zeros(m,1);
    ddotX    = zeros(m,1);
    ddotX_uf = zeros(m,1);
    Fz       = ones(m,1) * Fz0./N;
    H        = zeros(m,1);
    
    % For non-constant dOdt:
    % dOdt_p = 0.95*max(O)/(0.9*max(tf_vec)); % peak velocity
%     dOdt1 = (linspace(0, dOdt_p, 0.05*m))';
%     dOdt2 = dOdt_p*ones(0.9*m, 1);
%     dOdt3 = (linspace(dOdt_p, 0, 0.05*m))';
%     dOdt = [dOdt1; dOdt2; dOdt3];
    tau = 0.1;
    dOdt = dOdt_p*(1- exp(-tf_vec/tau));
    O = cumtrapz(tf_vec, dOdt);
    
    % For i = 1
    H(1) = (O(1).*R(1).^2./sqrt((L0.^2).*(1 + Fi(1)./K).^2 - O(1).^2.*R(1).^2));
    dXdt(1) = ((dOdt(1).*O(1).*R(1).^2 + O(1).^2.*dRdt(1).*R(1))./...
        sqrt(L0.^2.*((1 + Fi(1)./K).^2) - O(1).^2.*R(1).^2));
    Ffr(1) = (Fc*sign(dXdt(1)) + bX*dXdt(1));
    
    dRdt_num(1) = real(R0*(X0^.5).*O(1).*R(1)^2.*dOdt(1));
    dRdt_den(1) = real(2*X(1).^(5/2) - R0.*X0.^(1/2).*R0.*O(1).^2);
    dRdt(1) = (dRdt_num(1)./dRdt_den(1));
    
    for i = 2:m
        X(i)  = (sqrt(L0^2.*(1 + Fi(i-1)/K).^2 - O(i).^2*R(i-1).^2));
        Fz(i) = (Fz_fnctn(Fz0, Fz_MAX, FP, X0, Xmin, X(i)));
        R(i) = (R0.*sqrt(X0./X(i))); % (m) variable radius
        dRdt_num(i) = real(R0*(X0^.5).*O(i).*R(i)^2.*dOdt(i));
        dRdt_den(i) = real(2*X(i).^(5/2) - R0.*X0.^(1/2).*R0.*O(i).^2);
        dRdt(i) = (dRdt_num(i)./dRdt_den(i));
        dXdt(i) = ((dOdt(i).*O(i).*R(i).^2 + O(i).^2.*dRdt(i).*R(i))./...
            sqrt(L0.^2.*((1 + Fi(i-1)./K).^2) - O(i).^2.*R(i).^2));
        Ffr(i) = (Fc*sign(dXdt(i)) + bX*dXdt(i));
        F(i) = (M*ddotX(i) + Fz(i) + Ffr(i));
        Fi(i) = (F(i)./N./X(i).*X0); % fiber tension
        H(i) = (O(i).*R(i).^2./sqrt((L0.^2).*(1 + Fi(i)./K).^2 - O(i).^2.*R(i).^2));
        ddotX_uf(i) = ((dXdt(i) - dXdt(i-1))./(tf_vec(i) - tf_vec(i-1)));
        if i <= 5
            ddotX(i) = ddotX_uf(i);
        else
            ddotX(i) = sum(ddotX_uf(i-4:i))./5;
        end
    end
    Xdiff(j) = real(X0 - min(X) - dX);
    if abs(Xdiff(j)) > 1e-4
        %  Omax = Omax - s * sign((X0 - min(X)) - dX);
        % Adjust dOdt_p instead
        s(j) = b*abs(Xdiff(j));
        if (X0 - min(X)) < dX
            dOdt_p = dOdt_p + s(j);
        else
            dOdt_p = dOdt_p - s(j);
        end
        
    else
        break
    end
    j = j + 1;
    % disp('Contraction range');
    % disp(real(X0 - min(X)));

end
tEnd = toc(tStart);
disp(tEnd);
s = nonzeros(s);
Xdiff = nonzeros(Xdiff);
% stiffness
dSdO = 2.*(R0./(L0.^2 - O.^2.*R0.^2)).^2.*...
    (2.*L0.^2 - R0.^2.*O.^2./Kr.*R0.*O.^3 + L0.^3 ./ K .* O);

tf = (max(X) - min(X))./mean(dXdt); % Time to complete actuation (actual)

Tfr = bO*dOdt + Tc*sign(dOdt);
alph = [0; diff(dOdt)./diff(tf_vec)];
Tm = J.*alph + (H .* 0.5 .* dSdO.*F).*F + Tfr;

% Power = torque * angular speed
P = Tm.*dOdt; % required power OUTPUT from motor (not including efficiency)

T = max(Tm)*5; % the stall torque should not exceed 20% of the torque the
% TSA will require (conservative estimate since it's using the max torque)

% compute free-run speed based on the actual speed and stall torque
W1 = (max(dOdt)/(max(Tm) - T))*(0 - T) + 0;

W1 = W1*60/2/pi; % convert rad/s to RPM

T = T*10.197162129779; % convert (N*m) to (kg*cm);
% Source: https://www.convertunits.com/from/N-m/to/kg-cm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
plot(tf_vec, O);
xlabel('Time (s)')
ylabel('$\theta$ [rad]')

subplot(2,1,2);
plot(tf_vec, dOdt);
xlabel('Time (s)')
ylabel('$\dot{\theta}$ [rad/s]')

figure
subplot(2,1,1)
plot(tf_vec, X0 - X);
xlabel('Time (s)')
ylabel('$X_0 - X$ [m]')

subplot(2,1,2)
plot(tf_vec, dXdt);
xlabel('Time (s)')
ylabel('$\dot{X}$ [m/s]')