function dy = ode_2bodyPerturb( t, y, mu, perturbations, method )

% ode_2body ODE system for the perturbed two-body problem
% PROTOTYPE
% dy = ode_2bodyPerturb( t, y, mu, perturbations, method )
%
% INPUT:
% t[1] = Time 
% y[6x1] = State vector:
%                   elems. 1:3 -> position vector [km]
%                   elems. 4:6 -> velocity vector [km/s]
%   
% mu[1] = Gravitational parameter [km^3/s^2]
% perturbations[3] = vector containing a term (1 or 0) to take into account 
%           J2 effects, and a 2x1 vector containing area to mass ratio [m^2/kg]
%           and reflectivity coefficient [-] for SRP
% method = 'cart' or 'gauss'
%
% OUTPUT:
%   dy[6x1] = derivatives of the state (Keplerian elements)
%
% FUNCTIONS CALLED:
%   astroConstants
%   uplanet
%   kep2car
%
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

% J2 Constants
    pertJ2  = perturbations{1};
    pl_ID = 3; % Earth
    R_e = 6378.137;
    J2 = 0.00108263;

% SRP Constants
    pertSRP = perturbations{2};
    AU = astroConstants(2);
    pSR = 4.5*1e-6; % SRP @ 1AU
    mu_s = astroConstants(4);

switch method
    case 'cart'
        rr = y(1:3);
        vv = y(4:6);
        
    case 'gauss'
        a  = y(1);   e  = y(2);   
%         i  = y(3);
%         OM = y(4);   om = y(5);   th = y(6);

        i  = mod(y(3),pi);
        OM = mod(y(4),2*pi);   om = mod(y(5),2*pi);   th = mod(y(6),2*pi);


        [rr, vv] = kep2car(a,e,i,OM,om,th, mu);
        
    otherwise
        error(sprintf('Method: %s NOT FOUND!', method))
end

r = norm(rr);
v = norm(vv);
x = rr(1);
y = rr(2);
z = rr(3);

% Initialize perturbing accelerations
a_J2_x = 0;   a_SRP_x = 0;
a_J2_y = 0;   a_SRP_y = 0;
a_J2_z = 0;   a_SRP_z = 0;

% J2 perturbing acceleration in cartesian coordinates:
if pertJ2
    kJ2 = 1.5*J2*mu*R_e^2/r^4;
    a_J2_x = kJ2 * x/r*(5*z^2/r^2-1);
    a_J2_y = kJ2 * y/r*(5*z^2/r^2-1);
    a_J2_z = kJ2 * z/r*(5*z^2/r^2-3);
end

% SRP perturbing acceleration in cartesian coordinates:
if pertSRP
    cr = pertSRP(1);
    A_m = pertSRP(2);
    kep = uplanet(t/86400, pl_ID);
%     [RR1, ~] = kep2car(kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu_s);
    [RR1, ~] = kep2car(kep(1), kep(2), 23.4*pi/180, kep(4), kep(5), kep(6), mu_s);

%     rrSSC = RR1(:);
    rrSSC = RR1(:) + [x;y;z];

    rSSC = norm(rrSSC);
    aSRP = pSR*(AU^2/norm(rSSC)^2)*cr*(A_m)*1e-3;
    a_SRP_x = -aSRP*rrSSC(1)/norm(rSSC);
    a_SRP_y = -aSRP*rrSSC(2)/norm(rSSC);
    a_SRP_z = -aSRP*rrSSC(3)/norm(rSSC);
    
end

% Overall perturbing acceleration in cartesian coordinates:
    a_p_x = a_J2_x + a_SRP_x;
    a_p_y = a_J2_y + a_SRP_y;
    a_p_z = a_J2_z + a_SRP_z;




if strcmp(method, 'cart')
    % Derivatives of the state vector
    dy = [  vv(1)                   ;
            vv(2)                   ;
            vv(3)                   ;
            -mu/r^3 * x  +  a_p_x  ;
            -mu/r^3 * y  +  a_p_y  ;
            -mu/r^3 * z  +  a_p_z  ];

elseif strcmp(method, 'gauss')
    % J2 perturbing acceleration in {t,n,h} coordinates:
    tt = vv/norm(vv);                 % Tangent  unit vector
    hh = cross(rr,vv); hh = hh/norm(hh); % Normal   unit vector
    nn = cross(hh,tt);                % Binormal unit vector
    ROT_tnh2xyz = [tt(:) nn(:) hh(:)];
    a_p_xyz = [a_p_x a_p_y a_p_z]';
    
    a_p_tnh = ROT_tnh2xyz' * a_p_xyz;
    a_t = a_p_tnh(1);  a_n = a_p_tnh(2);  a_h = a_p_tnh(3);
    
    
    % Variation of orbital elements:
    b = a * sqrt(1-e^2);
    p = b^2/a;
    n = sqrt(mu/a^3);
    h = n*a*b;
    th_star = th + om;

    a_dot  = 2*a^2*v/mu * a_t;
    e_dot  = 1/v* ( 2*(e+cos(th))*a_t - r/a*sin(th)*a_n );
    i_dot  = r*cos(th_star)/h * a_h;
    OM_dot = r*sin(th_star)/(h*sin(i)) * a_h;
    om_dot = 1/(e*v) * ( 2*sin(th)*a_t + (2*e + r/a*cos(th))*a_n ) - r*sin(th_star)*cos(i)/(h*sin(i))*a_h;
    th_dot = h/(r^2) - 1/(e*v) * ( 2*sin(th)*a_t + (2*e + r/a*cos(th))*a_n );


    % Set the derivatives of the state
    dy = [  a_dot ;
            e_dot ;
            i_dot ;
            OM_dot;
            om_dot;
            th_dot];
end

    
  
 
end

