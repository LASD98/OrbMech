function [c,ceq] = nonlcon(X,h_atm_venus,T)
% non linear constraint function for fmincon
%   INPUTs: 
%       x[3]: mjd2000 departure, fly-by & and arrival transfer time mjd2000
%       h_atm_venus[1] km
%       T[3]: days  period of revolution of Mercury, Venus, NEO81
%   OUTPUTs:
%       c <= 0   constraint function
%       ceq == 0
%
% FUNCTIONS CALLED:
% cost_tot.m

[~,rp] = cost_tot(X,h_atm_venus,T);
if anynan(rp)
    c = [];
    ceq = [];
    return
end
R_V = astroConstants(22);

c = [R_V + h_atm_venus - rp; % non collision constraint 
     X(1) - X(2);            % dates order constraint
     X(2) - X(3);
     date2mjd2000([2028,1,1,0,0,0]) - X(1); % dates constraint
     X(3) - date2mjd2000([2058,1,1,0,0,0]);
     X(2) - X(1) - T(2);  % max ToF constraint
     X(3) - X(2) - T(3)]; 
                   % nonlinear inequality c <= 0
ceq = [];          % nonlinear equality c = 0
end
