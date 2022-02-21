function Fz = Fz_fnctn(Fz_ZERO, Fz_MAX, FP, X0, Xmin, X)
switch FP
    case 1 % Linear Profile
        dXmax = X0 - Xmin;
        dX = X0 - X;
        m = (Fz_MAX - Fz_ZERO)./(dXmax);
        Fz = m*dX + Fz_ZERO;
    case 2 % Parabolic profile
        % Assume the vertex of the parabola is at Fz_ZERO
        a = (Fz_MAX - Fz_ZERO)./((X0 - Xmin).^2);
        Fz = transpose(a*(X0 - X).^2 + Fz_ZERO);
    case 3
        Fz = (Fz_MAX - Fz_ZERO)*rand(length(X),1) + Fz_ZERO;
    otherwise % Assume linear force profile
        m = (Fz_MAX - Fz_ZERO)./(X0 - Xmin);
        Fz = m*X + Fz_ZERO;
end
end