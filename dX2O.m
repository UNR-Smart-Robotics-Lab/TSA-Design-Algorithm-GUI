function [O, Xmin2] = dX2O(X0, dX, Fz_MAX, R0, N, L0, K)

Xmin = X0 - dX; % unloaded contraction range
% Omax = sqrt((X0^2 - Xmin^2)/R0^2); % initial guess
Rmax = R0*sqrt(X0/Xmin);
Fimax = Fz_MAX*X0/N/Xmin;
Omax = sqrt((L0^2*(1+Fimax/K)^2 - Xmin^2)/(Rmax^2));
%k = 1; % initialize loop counter
% tic % start timer
% O_step = 0.1; 
% Xmin2 = Xmin;
% while 1
%     % Variable radius model
%     Rmax = R0*sqrt(X0/Xmin2);
%     
%     % Fiber Tension
%     Fimax = Fz_MAX*X0/N/Xmin2;
%     
%     % TSA Model accounting for finite stiffness
%     Xmin2 = sqrt(L0^2 * (1 + Fimax / K)^2 - Omax^2 * Rmax^2);
%     
%     if ~isreal(Xmin2) % If Xmin2 is complex
%         Omax = Omax - O_step; % decrease the turns
%     end
%     
%     % Checking if Xmin is within the tolerance specified by the user.
%     if abs(Xmin - Xmin2) < 0.00005
%         break
%     elseif Xmin2 > Xmin
%         Omax = Omax + O_step; % Change Omax based on the computed Xmin.
%     elseif Xmin2 < Xmin
%         Omax = Omax - O_step;
%     end
%     
%     % If the loop runs for more than 20 seconds, exit.
%     etime = toc;
%     if etime > 2
%         break
%     end
%     
%     % Loop counter
%     k = k + 1;
% end
O = transpose(linspace(0, Omax, 1000));
Xmin2 = Xmin;
end