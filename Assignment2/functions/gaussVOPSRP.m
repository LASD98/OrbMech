function dkep = gaussVOPSRP(t,kep,Re,omE,body,mu,J2)

% Defines the two-body differential problem with SRP and J2 perturbations,
% to be solved with an ODE solver, like ode113.
% 
% PROTOTYPE
%   dkep = gaussVOPSRP(~,kep,Re,omE,body,mu,J2)
% 
% INPUT:
%   kep[6] = vector of Cartesian position and velocities, positions first.
%   Re[1] = radius of the parent body [km]
%   omE[1] = parent body rotation speed [rad/s]
%   body[2] = vector containing area to mass ratio [m^2/kg] and cR [-]
%   mu[1] = gravitational parameter [km^3/s^2]
%   J2[1] = second zonal harmonic of the zonal variations of the gravitational
%     field, due to the oblateness of the attractor [-]
% 
% OUTPUT:
%   dkep[6] = vector of derivatives of Keplerian parameters.
% 
% FUNCTIONS CALLED:
%   astroConstants
%   uplanet
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

if nargin < 7
    J2 = 0.00108263;
end

% SRP Constants
AU = astroConstants(2);
pSR = 4.5*1e-6; % SRP @ 1AU
mu_s = astroConstants(4);
A_m = body(1);
cr = body(2);
pl_ID = 3; % Earth
kepE = uplanet(t/86400, pl_ID);

a = kep(1);
e = kep(2);
i = kep(3);
Om = kep(4);
om = kep(5);
f = kep(6);

[rr,vv] = par2car(kep,mu);

h = norm(cross(rr,vv));
p = h^2/mu;
th = om + f;

xx = rr(1);
yy = rr(2);
zz = rr(3);
r = norm(rr);



accJ2 = [3/2*J2*mu*Re^2/r^4*(xx/r*(5*zz^2/r^2-1))
         3/2*J2*mu*Re^2/r^4*(yy/r*(5*zz^2/r^2-1))
         3/2*J2*mu*Re^2/r^4*(zz/r*(5*zz^2/r^2-3))];


[RR1, ~] = kep2car(kepE(1), kepE(2), kepE(3), kepE(4), kepE(5), kepE(6), mu_s);
rrSSC = RR1(:);

rSSC = norm(rrSSC);
aSRP = pSR*(AU^2/norm(rSSC)^2)*cr*(A_m)*1e-3;
a_SRP_x = -aSRP*rrSSC(1)/norm(rSSC);
a_SRP_y = -aSRP*rrSSC(2)/norm(rSSC);
a_SRP_z = -aSRP*rrSSC(3)/norm(rSSC);
accSRP = [a_SRP_x a_SRP_y a_SRP_z];


ROm = [cos(Om) -sin(Om) 0
       sin(Om)  cos(Om) 0
       0 0 1];
   
Rom = [cos(om) -sin(om) 0
       sin(om)  cos(om) 0
       0 0 1];

Ri = [1 0 0
      0 cos(i) -sin(i)
      0 sin(i)  cos(i)];

% T = rotz(rad2deg(Om))*rotx(rad2deg(i))*rotz(rad2deg(om));
T = ROm*Ri*Rom;

TT =[cos(f) -sin(f) 0
     sin(f)  cos(f) 0
     0 0 1];

R = T*TT;

acc = R\(accJ2 + accSRP);

da = 2*a^2/h*(e*sin(f)*acc(1)+p/r*acc(2));
de = 1/h*(p*sin(f)*acc(1)+((p+r)*cos(f)+r*e)*acc(2));
di = r*cos(th)/h*acc(3);
dOm = r*sin(th)/(h*sin(i))*acc(3);
dom = 1/(h*e)*(-p*cos(f)*acc(1)+(p+r)*sin(f)*acc(2))-r*sin(th)*cos(i)/h/sin(i)*acc(3);
df = h/r^2+1/(e*h)*(p*cos(f)*acc(1)-(p+r)*sin(f)*acc(2));

dkep = [da
        de
        di
        dOm
        dom
        df];
