function [L0,X0, K] = dX2X0(dX, Fz, Si, N, strain_ideal)
% Untwisted TSA Parameters
%i = 1; % initialize loop counter
%tic;
%strain_UL = strain_ideal;  % UL = unloaded

K = Si/N;
X0 = dX/strain_ideal;
L0 = X0 - Fz/2/K;
K = K*L0;

% Iterative Solving
% while 1
%     L0 = dX/strain_UL; % UNLOADED TSA length
%     
%     K = Si*L0/N; % NORMALIZED Stiffness of INDIVIDUAL string based on stiffness of TSA at
%     % ZERO turns (untwisted)
%     
%     X0 = Fz*L0/2/K + L0; % LOADED TSA length
%     
%     strain_loaded = dX/X0;
%     if (abs(strain_loaded - strain_ideal)) < 0.00000001
%         break
%     elseif strain_loaded < strain_ideal
%         strain_UL = strain_UL + 0.00000001; % change the unloaded string length
%     elseif strain_loaded > strain_ideal
%         strain_UL = strain_UL - 0.00000001; % change the unloaed string length
%     end
%     
%     etime = toc;
%     if etime > 1 % if the loop takes longer than 1s to finish, exit because something is wrong.
%         break
%     end
%     
%     i = i + 1; % loop counter
% end
% clear i; % free the variable to be used for something else
end