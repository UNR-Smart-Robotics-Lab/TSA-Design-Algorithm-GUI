function [dt,s2] = forward_model(s1)
% T: stall torque of the motor from the database
% W: free-run speed of the motor from the database
% Omax: maximum rotations
% dt: actual attainable contraction time

%R0, Fz0, Fz_MAX, FP, N, X0, L0, K, dX, b,tau,dynamic_constants_
%max_Tm = s1.T/5;
%s1.dOdt_p = s1.W/s1.T*(s1.T-max_Tm);
dt = s1.Omax/s1.dOdt_p;
bX = s1.dynamic_constants_(1);
bO = s1.dynamic_constants_(2);
Kr = s1.dynamic_constants_(3);
J  = s1.dynamic_constants_(4);
Fc = s1.dynamic_constants_(5);
Tc = s1.dynamic_constants_(6);
M  = s1.dynamic_constants_(7);
m  = 1000; % level of discretization
Xmin = s1.X0 - s1.dX; % unloaded contraction range
s2.tf_vec = (linspace(0, dt, m))';
s2.X        = ones(m, 1) * s1.X0;
s2.R        = ones(m, 1) * s1.R0;
s2.Fi       = ones(m,1) * s1.Fz0./s1.N;
s2.dXdt     = zeros(m,1);
s2.F        = ones(m,1) * s1.Fz0;
s2.Ffr      = zeros(m,1);
s2.dRdt     = zeros(m,1);
s2.ddotX_uf = zeros(m,1);
s2.ddotX    = zeros(m,1);
s2.Fz       = ones(m,1) * s1.Fz0./s1.N;
s2.H        = zeros(m,1);
s2.dOdt = s1.dOdt_p*(1- exp(-s2.tf_vec/s1.tau));
s2.O = cumtrapz(s2.tf_vec,s2.dOdt);
s2.H(1) = ...
    (s2.O(1).*s2.R(1).^2./sqrt((s1.L0.^2).*(1+s2.Fi(1)./s1.K).^2-s2.O(1).^2.*s2.R(1).^2));
s2.dXdt(1) = ((s2.dOdt(1).*s2.O(1).*s2.R(1).^2 + s2.O(1).^2.*s2.dRdt(1).*s2.R(1))./...
    sqrt(s1.L0.^2.*((1+s2.Fi(1)./s1.K).^2)-s2.O(1).^2.*s2.R(1).^2));
s2.Ffr(1) = (Fc*sign(s2.dXdt(1)) + bX*s2.dXdt(1));
for i = 2:m
    s2.X(i)  = (sqrt(s1.L0^2.*(1 + s2.Fi(i-1)/s1.K).^2 - s2.O(i).^2*s2.R(i-1).^2));
    s2.Fz(i) = (Fz_fnctn(s1.Fz0,s1.Fz_MAX,s1.FP,s1.X0,Xmin,s2.X(i)));
    s2.R(i) = (s1.R0.*sqrt(s1.X0./s2.X(i))); % (m) variable radius
    s2.dRdt(i) = 0;
    s2.dXdt(i) = ( ...
        s2.dOdt(i).*s2.O(i).*s2.R(i).^2 + ...
        s2.O(i).^2.*s2.dRdt(i).*s2.R(i)... % - (L0^2).*(1 + Fi(i-1)./K).*dFidt(i-1)./K...
        )./...
        sqrt(s1.L0.^2.*((1+s2.Fi(i-1)./s1.K).^2)-s2.O(i).^2.*s2.R(i).^2);
    s2.Ffr(i) = (Fc*sign(s2.dXdt(i))+bX*s2.dXdt(i));
    s2.F(i) = (M*s2.ddotX(i) + s2.Fz(i) + s2.Ffr(i));
    s2.Fi(i) = (s2.F(i)./s1.N./s2.X(i).*s1.X0); % fiber tension
    s2.H(i) = (s2.O(i).*s2.R(i).^2./...
        sqrt((s1.L0.^2).*(1+s2.Fi(i)./s1.K).^2-s2.O(i).^2.*s2.R(i).^2));
    s2.ddotX_uf(i) =((s2.dXdt(i)-s2.dXdt(i-1))./(s2.tf_vec(i)-s2.tf_vec(i-1)));
    if i <= 5
        s2.ddotX(i) = s2.ddotX_uf(i);
    else
        s2.ddotX(i) = sum(s2.ddotX_uf(i-4:i))./5;
    end
end
s2.dSdO = 2.*(s1.R0./(s1.L0.^2 - s2.O.^2.*s1.R0.^2)).^2.*...
    (2.*s1.L0.^2-s1.R0.^2.*s2.O.^2./Kr.*s1.R0.*s2.O.^3+s1.L0.^3./s1.K .*s2.O);
s2.Tfr = bO*s2.dOdt + Tc*sign(s2.dOdt);
s2.alph = [0;diff(s2.dOdt)./diff(s2.tf_vec)];
s2.Tm = J.*s2.alph + (s2.H.*0.5.*s2.dSdO.*s2.F).*s2.F+s2.Tfr;
s2.P = s2.Tm.*s2.dOdt;
end