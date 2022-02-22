function [bX, bO, Kr, J, Fc, Tc, M] = dynamic_constants()
    % bX: Viscous friction coefficient of the load on the slider (N*s/m)
    bX = 10;
    % bO: Viscous friction coefficient of the motor (N*m*s)
    bO = 2e-4;
    % Radial stiffness: Kr (N/m)
    Kr = 1.2e4;
    % J: otor moment of inertia
    J = 1e-6; % kg*m^2
    % Fc: Sliding Coulombic frictional force (N)
    Fc = 5;
    % tauc: Rotational Coulombic frictional force (N*m)
    Tc = 2e-3;
    % m: mass of payload
    M = 1; % kg
end