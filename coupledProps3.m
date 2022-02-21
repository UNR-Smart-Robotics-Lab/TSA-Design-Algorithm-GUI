function [dOdt, Tm, tfinal, P, s, Xdiff, tEnd,O] = ...
    coupledProps3(R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dt, dX, b,tau,dynamic_constants_)
if nargin < 11
    b = 324;
    tau = 0.1;
end
% if nargin < 13
%     [bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants();
% end
bX = dynamic_constants_(1);
bO = dynamic_constants_(2);
Kr = dynamic_constants_(3);
J = dynamic_constants_(4);
Fc = dynamic_constants_(5);
Tc = dynamic_constants_(6);
M = dynamic_constants_(7);
m = 1000; % level of discretization

Xmin = X0 - dX; % unloaded contraction range
tf_vec = (linspace(0, dt, m))';
j = 1; % initialized loop counter
% Initial guess for Omax using variable radius, finite stiffness
[O_pre, ~] = dX2O(X0, dX, Fz_MAX, R0, N, L0, K);
Omax = max(O_pre); 
dOdt_p = Omax./dt;
tStart = tic;
s = zeros(10000, 1);
Xdiff = zeros(10000,1);
while 1
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
    ddotX_uf = zeros(m,1);
    ddotX    = zeros(m,1);
    Fz       = ones(m,1) * Fz0./N;
    H        = zeros(m,1);
    dFidt    = zeros(m,1);

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
        dRdt(i) = 0;
        dXdt(i) = ( ...
            dOdt(i).*O(i).*R(i).^2 + ...
            O(i).^2.*dRdt(i).*R(i)... % - (L0^2).*(1 + Fi(i-1)./K).*dFidt(i-1)./K...
            )./...
            sqrt(L0.^2.*((1 + Fi(i-1)./K).^2) - O(i).^2.*R(i).^2);
        Ffr(i) = (Fc*sign(dXdt(i)) + bX*dXdt(i));
        F(i) = (M*ddotX(i) + Fz(i) + Ffr(i));
        Fi(i) = (F(i)./N./X(i).*X0); % fiber tension
        % Fi(i) = K/L0*(X(i) - L0);
        % dFidt(i) = K./L0*dXdt(i);
        %dFidt(i) = (Fi(i) - Fi(i-1))/(tf_vec(i) - tf_vec(i-1));
        H(i) = (O(i).*R(i).^2./sqrt((L0.^2).*(1 + Fi(i)./K).^2 - O(i).^2.*R(i).^2));
        %H(i) = (O(i)*R(i)^2)*dOdt(i)./(X(i));
        % Unfiltered Acceleration
        ddotX_uf(i) = ((dXdt(i) - dXdt(i-1))./(tf_vec(i) - tf_vec(i-1)));
        if i <= 5
            ddotX(i) = ddotX_uf(i);
        else
            ddotX(i) = sum(ddotX_uf(i-4:i))./5;
        end
        % Filtered Acceleartin
        %ddotX_vec = sgolayfilt(ddotX_uf,3,5);
        %ddotX(i) = ddotX_vec(i);
    end
    Xdiff(j) = real(X(1) - X(end) - dX);
    strain = (X(1) - X(end))./X(1);
    %if abs(Xdiff(j)) > 1e-3*dX
    if abs(Xdiff(j)) > 1e-3*dX
        %  Omax = Omax - s * sign((X0 - min(X)) - dX);
        % Adjust dOdt_p instead
        s(j) = b*abs(Xdiff(j));
        if (X0 - min(real(X))) < dX
            dOdt_p = dOdt_p + s(j);
        else
            dOdt_p = dOdt_p - s(j);
        end
    else
        break
    end
    j = j + 1;
     %disp('Contraction range');
     %disp(abs(Xdiff(j-1)));
    if toc(tStart) >= 10
        break
    end
end
tEnd = toc(tStart);
% figure(1)
% plot(j,tEnd,'.k'); hold on;
s = nonzeros(s);
Xdiff = nonzeros(Xdiff);
% stiffness
dSdO = 2.*(R0./(L0.^2 - O.^2.*R0.^2)).^2.*...
    (2.*L0.^2 - R0.^2.*O.^2./Kr.*R0.*O.^3 + L0.^3 ./ K .* O);
%tfinal = (X(1) - X(end))./mean((dXdt)); % Time to complete actuation (actual)
tfinal = (max(O) - min(O))./mean(dOdt);
% figure(1)
% subplot(4,3,1); plot(O); title('$O$');
% subplot(4,3,2); plot(X); title('$X$');
% subplot(4,3,3); plot(R); title('$R$');
% subplot(4,3,4); plot(Fz); title('$F_z$');
% subplot(4,3,5); plot(F); title('$F$');
% subplot(4,3,6); plot(Ffr); title('$F_{fr}$');
% subplot(4,3,7); plot(Fz); title('$F_z$');
% subplot(4,3,8); plot(dRdt); title('$dR/dt$');
% subplot(4,3,9); plot(ddotX); title('$\ddot{X}$');
% subplot(4,3,10); plot(H); title('$H$');
% subplot(4,3,11); plot(dOdt); title('$\dot{O}$')
% subplot(4,3,12); plot(dXdt); title('$\dot{X}$')

Tfr = bO*dOdt + Tc*sign(dOdt);
alph = [0; diff(dOdt)./diff(tf_vec)];
Tm = J.*alph + (H .* 0.5 .* dSdO.*F).*F + Tfr;
% Power = torque * angular speed
P = Tm.*dOdt; % required power OUTPUT from motor (not including efficiency)
