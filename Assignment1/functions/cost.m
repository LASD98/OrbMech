function [Dv1,Dv2,TOF,cartesian,DV1,DV2] = cost(x,ID,mod,TOFMAX)
% Cost (objective) function for each arc, part of "time grid minimization", 
% and part of total cost function for fmincon (cost_tot)
% Can be upgraded so to receive gregorian date as input, which can be 
% converted here to mjd2000, rather than outside
%   INPUTs: 
%       x[2]: mjd2000 "initial tr. time" & "final transfer time or ToF"
%       ID[2]
%       mod[1]: 1 for tw2 input, 2 for ToF input
%       TOFMAX[1]: maximum ToF for the Lambert arc [days]
%   OUTPUTs:
%       Dv1[1]: cost of first maneuver [km/s]
%	    Dv2[1]: cost of second maneuver [km/s]
%       TOF[1]: ToF of arc [seconds]
%       cartesian[3x6]: [RI,RF,vi,vf,VI',VF']
%	    DV1[3]: vectorial cost of first maneuver [km/s]
%	    DV2[3]: vectorial cost of second maneuver [km/s]

% FUNCTIONS CALLED:
% kep2car.m

tw1 = x(1); 
if mod == 1
    tw2 = x(2);
    if tw2 <= tw1
        Dv1 = NaN;
        Dv2 = NaN;
        TOF = NaN;
        cartesian = NaN;
        DV1 = NaN;
        DV2 = NaN;
    return
    end
elseif mod == 2
    ToF = x(2);
    tw2 = tw1 + ToF*365.2417; 
end
MU = astroConstants(4);
TOF = (tw2-tw1)*3600*24; % [s]
if TOF/(3600*24) > TOFMAX
    Dv1 = NaN;
    Dv2 = NaN;
    TOF = NaN;
    cartesian = NaN;
    DV1 = NaN;
    DV2 = NaN;
    return 
end

mjd2000 = [tw1,tw2];
rvect=zeros(3,2); vvect = rvect;
for s = 1:2
    if ID(s) < 12
        kep = uplanet(mjd2000(s),ID(s));
    else
        kep = ephNEO(mjd2000(s),ID(s));
    end
    a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
    [rvect(:,s),vvect(:,s)] = kep2car(a,e,i,Omega,omega,theta,MU);
end

RI = rvect(:,1); % [km]     position 1
RF = rvect(:,2); % [km]     position 2
vi = vvect(:,1); % [km/s]   velocity before transfer (initial orbit)
vf = vvect(:,2); % [km/s]   velocity after transfer (final orbit)

orbitType = 0; % (0: prograde, 1: retrograde)
[~,~,~,~,VI,VF,~,~] = lambertMR(RI,RF,TOF,MU,orbitType);
cartesian = [RI,RF,vi,vf,VI',VF'];

% cost
DV1 = VI'-vi;
DV2 = VF'-vf;
Dv1 = norm(DV1);
Dv2 = norm(DV2);
% Dvtot = Dv1+Dv2;

end