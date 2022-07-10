function [L0,X0,K] = dX2X0(dX, Fz, K, N, strain_ideal,Xp)
% dX: contraction range (m)
% Fz: axial force (N)
% N = number of strings
% strain_ideal = ideal strain
% K: normalized stiffness (N)
X0 = (dX + Xp)/strain_ideal*N/N;
L0 = X0/(1+Fz/(2*K));
end