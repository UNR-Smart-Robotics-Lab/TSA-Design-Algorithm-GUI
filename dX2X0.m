function [L0,X0,K] = dX2X0(dX, Fz, K, N, strain_ideal)
% dX: contraction range (m)
% Fz: axial force (N)
% N = number of strings
% strain_ideal = ideal strain
% K: normalized stiffness (N)
X0 = dX/strain_ideal;
L0 = X0/(1+Fz/(2*K));
end