clearvars ; clc; close all; format long g

MercuryID  = 1;                     % ID for uplanet.m
VenusID    = 2;                     % ID for uplanet.m
NeoID      = 81;                    % ID for ephNEO.m
mu_Sun     = astroConstants(4);     % Sun gravitational parameter   [km^3/s2]      
mu_Venus   = astroConstants(12);    % Venus gravitational parameter [km^3/s2]
R_Venus    = astroConstants(22);    % Venus radius                  [km]
h_atm      = 500;                   % Venus atmosphere's height     [km]

% Period of the planets 
kep_Mercury = uplanet(10837.4268,MercuryID);
kep_Venus   = uplanet(11716.0648,VenusID);
kep_81      = ephNEO(12981.0216,NeoID);

a_M  = kep_Mercury(1);
a_V  = kep_Venus(1);
a_81 = kep_81(1);
T_M  = 2*pi*sqrt(a_M^3/mu_Sun);
T_V  = 2*pi*sqrt(a_V^3/mu_Sun);
T_81 = 2*pi*sqrt(a_81^3/mu_Sun);

T_Mdays = T_M / (24*3600);
T_Vdays = T_V / (24*3600);
T_81days = T_81 / (24*3600);

%%
guess = [10827.4268 11706.0648 12991.0216];
lb = [];
ub = [];
costfun = @(t) deltaVtot(t,MercuryID,VenusID,NeoID,mu_Sun,mu_Venus,R_Venus);
constrfun = @(t) nonlincon(t,MercuryID,VenusID,NeoID,mu_Sun,mu_Venus,R_Venus,h_atm);
options = optimoptions("fmincon", 'Algorithm','sqp', 'Display','iter-detailed', ...
                       'MaxFunctionEvaluations',3000,'StepTolerance',1e-12, ...
                       'ConstraintTolerance',1e-10);
[bestdate, deltaVopt] = fmincon(costfun,guess,[],[],[],[],lb, ub,constrfun,options);

[~,rp_norm] = deltaVtot(bestdate,MercuryID,VenusID,NeoID,mu_Sun,mu_Venus,R_Venus);
%% FUNCTIONS

function [deltaV_dep,deltaV_arr,Vi_arc,Vf_arc,t_par] = deltaV_arc(Body1_ID, Body2_ID, t1, t2, mu_Sun)

% BODY 1: check if it's a planet or an asteroid, retrieve keplerian
% parameters (from uplanet.m|ephNEO.m|ephMoon.m) and retrieve its position r1
% and velocity v1 at time t1 [MJD2000, days]

if Body1_ID < 11 % It's a planet (except Moon)
    kep1    = uplanet(t1, Body1_ID);
    [r1,v1] = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),mu_Sun);
elseif Body1_ID > 11 % it's an asteroid
    kep1    = ephNEO(t1, Body1_ID);
    [r1,v1] = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),mu_Sun);
else % It's the Moon
    fprintf('strunz hai sbagliato');
end

% BODY 2: check if it's a planet or an asteroid, retrieve keplerian
% parameters (from uplanet.m|ephNEO.m|ephMoon.m) and retrieve its position r2
% and velocity v2 at time t2 [MJD2000, days]

if Body2_ID < 11 % It's a planet
    kep2    = uplanet(t2, Body2_ID);
    [r2,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),mu_Sun);
elseif Body2_ID > 11 % it's an asteroid
    kep2    = ephNEO(t2, Body2_ID);
    [r2,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),mu_Sun);
else % It's the Moon
    fprintf('strunz hai sbagliato');
end

% Get time of flight in seconds
days2sec = 24*3600;         % [s]
ToF = (t2-t1)*days2sec;     % [s]

% Lambert's arc given r1,r2 and ToF
if ToF >= 0
    [~,~,~,~,Vi_arc,Vf_arc,t_par] = lambertMR(r1,r2,ToF,mu_Sun,0,0,0,2);
    deltaV_dep = norm(Vi_arc - v1');
    deltaV_arr = norm(v2' - Vf_arc);
else
    deltaV_dep = NaN;
    deltaV_arr = NaN;
    Vf_arc = NaN;
    t_par  = NaN;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [deltaV,rp_norm] = deltaVtot(t,Body1_ID,Body2_ID,Body3_ID,mu_Sun,mu2,R2)

t1 = t(1); % departure date
t2 = t(2); % flyby date
t3 = t(3); % arrival date

% Second Body data
kep2   = uplanet(t2, Body2_ID);
[~,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),mu_Sun);


[deltaV_dep_arc1,~,~,Vf_arc1]  = deltaV_arc(Body1_ID, Body2_ID, t1, t2, mu_Sun);
[~,deltaV_arr_arc2,Vi_arc2]  = deltaV_arc(Body2_ID, Body3_ID, t2, t3, mu_Sun);


% Powered gravity assist
Vinf_minus = Vf_arc1' - v2;
Vinf_plus  = Vi_arc2' - v2;

% Turning angle
delta = acos(dot(Vinf_minus,Vinf_plus)/(norm(Vinf_minus)*norm(Vinf_plus)));

% Radius of pericentre
fun = @(x)  asin(1./(1+((x.*dot(Vinf_minus,Vinf_minus))./mu2))) + ...
            asin(1./(1+((x.*dot(Vinf_plus,Vinf_plus))  ./mu2))) - ...
            delta;

opt = optimset('Display','off', 'TolFun',1e-13);
rp_norm = fsolve(fun, R2, opt);


Vp_minus_norm = sqrt( (norm(Vinf_minus))^2 + (2*mu2)/rp_norm);
Vp_plus_norm  = sqrt( (norm(Vinf_plus))^2 + (2*mu2)/rp_norm);

% Magnitude of the manoeuvre at pericentre
deltaVp = norm(Vp_minus_norm - Vp_plus_norm);

% Total Cost function
deltaV  = deltaV_dep_arc1 + deltaV_arr_arc2 + deltaVp;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,Ceq] = nonlincon(t,Body1_ID,Body2_ID,Body3_ID,mu_Sun,mu2,R2,h_atm)

t1 = t(1); % departure date
t2 = t(2); % flyby date
t3 = t(3); % arrival date

% BODY 1: check if it's a planet or an asteroid, retrieve keplerian
% parameters (from uplanet.m|ephNEO.m|ephMoon.m) and retrieve its position r1 
% and velocity v1 at time t1 [MJD2000, days]

if Body1_ID < 11 % It's a planet (except Moon)
    kep1 = uplanet(t1, Body1_ID);
    r1   = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),mu_Sun);
elseif Body1_ID > 11 % it's an asteroid
    kep1 = ephNEO(t1, Body1_ID);
    r1   = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),mu_Sun);
else % It's the Moon
    fprintf('strunz hai sbagliato');
end

% BODY 2: check if it's a planet or an asteroid, retrieve keplerian
% parameters (from uplanet.m|ephNEO.m|ephMoon.m) and retrieve its position r2
% and velocity v2 at time t2 [MJD2000, days]

if Body2_ID < 11 % It's a planet
    kep2    = uplanet(t2, Body2_ID);
    [r2,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),mu_Sun);
elseif Body2_ID > 11 % it's an asteroid
    kep2    = ephNEO(t2, Body2_ID);
    [r2,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),mu_Sun);
else 
    fprintf('strunz hai sbagliato');
end

% BODY 3: check if it's a planet or an asteroid, retrieve keplerian
% parameters (from uplanet.m|ephNEO.m|ephMoon.m) and retrieve its position r3
% and velocity v3 at time t3 [MJD2000, days]

if Body3_ID < 11 % It's a planet
    kep3 = uplanet(t3, Body3_ID);
    r3   = kep2car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),mu_Sun);
elseif Body3_ID > 11 % it's an asteroid
    kep3 = ephNEO(t3, Body3_ID);
    r3   = kep2car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),mu_Sun);
else
    fprintf('strunz hai sbagliato');
end

% Time of flight for the two arcs
days2sec = 24*3600;

tof1 = (t2-t1)*days2sec;
tof2 = (t3-t2)*days2sec;

%First tranfer arc
[~,~,~,~,~,Vf_arc1] = lambertMR(r1,r2,tof1,mu_Sun,0,0,0,2);

% Second tranfer arc
[~,~,~,~,Vi_arc2]   = lambertMR(r2,r3,tof2,mu_Sun,0,0,0,2);

Vinf_minus = Vf_arc1' - v2;
Vinf_plus  = Vi_arc2' - v2;

% Turning angle
delta = acos(dot(Vinf_minus,Vinf_plus)/(norm(Vinf_minus)*norm(Vinf_plus)));

% Radius of pericentre
fun = @(x)  asin(1./(1+((x.*dot(Vinf_minus,Vinf_minus))./mu2))) + ...
            asin(1./(1+((x.*dot(Vinf_plus,Vinf_plus))  ./mu2))) - ...
            delta;

opt = optimset('Display','off', 'TolFun',1e-13);
rp_norm = fsolve(fun, R2, opt);

% Compute Body 2 sphere of influence radius
rSOI = kep2(1) * (mu2/mu_Sun)^(2/5);

Ceq = [];

C(1) = R2 + h_atm - rp_norm;
C(2) = rp_norm - rSOI;
C(3) = t1 - t2;
C(4) = t2 - t3;
end

