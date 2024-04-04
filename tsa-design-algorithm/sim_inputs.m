function [dt, dX, Fz0, Fz_MAX, Si, N, R0, FP, strain_ideal] = sim_inputs()
% Inputs
dt = 2;
dX = 0.1;
Fz0 = 10;
Fz_MAX = 25;
Si = 500;
N = 2;
strain_ideal = 0.3;
R0 = 0.001;
FP = 1; % linear
% P = 200; % really high
% E = 200; % really high
end