function [O, Xmin2] = dX2O(X0, dX, Fz_MAX, R0, N, L0, K)
Xmin = X0 - dX;
Rmax = R0*sqrt(X0/Xmin);
Fimax = Fz_MAX*X0/N/Xmin;
Omax = sqrt((L0^2*(1+Fimax/K)^2 - Xmin^2)/(Rmax^2));
O = transpose(linspace(0, Omax, 1000));
Xmin2 = Xmin;
end