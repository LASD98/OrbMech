%%   FUNCTIONS CALLED:
%   Earth3d
%   plotOrbit
%   groundTrackODE
%   groundTrackJ2
%   groundTrackJ2SRP
%   planisphere
%   plotGroundTrack
%   repeatingGTJ2
%   orbit
%   orbitJ2
%   orbitSRP
%   orbitSRPJ2
%   ode_2bodyPerturb
%   par2car
%   car2par
%   anglePostProcessingShift
%   anglePostProcessingMod
%   spectralAnalysis
%   MAF
%   colorPlot3 
%   OrbitEvolutionVideo

clear; close all; clc

addpath('functions')
addpath('horizons')

%% Orbit definition and representation

totalT = tic;

mu = 398600;

a = 2.7251*1e4;   
e = 0.4485;
i = deg2rad(11.3122);
Om = 0;
om = 0;
th = 0;

figure
Earth3d
orb = plotOrbit(a,e,i,Om,om,0,2*pi,2*pi/100);
view(40,20)
grid on
pbaspect([1 1 1])
title('Initial orbit representation')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

%% Ground tracks (no perturbation and J2 + SRP)

tic

C_R = 1;
A2m = 5;
perturbations = {1 [C_R A2m]};
body = [A2m,C_R];

kep0 = [a,e,i,Om,om,th];
lonG0 = 0;

t2 = 10;
t3 = 365; % time used for main analyses
npts1 = 1000;
deltath2 = 2; % 2 minutes to calculate every Δθ = 1°
npts2 = 10*24*60/deltath2;
deltath3 = 25; % 25 minutes to calculate every Δθ = 12° (or 5 minutes to calculate every Δθ = 2.4°)
npts3 = 365*24*60/deltath3;

time1 = linspace(0,2*pi*sqrt(a^3/mu),npts1); % one orbit
time2 = linspace(0,seconds(days(t2)),npts2); % t2 days
time3 = linspace(0,seconds(days(t3)),npts3); % t3 days

[GT1.alpha, GT1.delta, GT1.lon, GT1.lat] = groundTrackODE(kep0, lonG0, time1);
[GT2.alpha, GT2.delta, GT2.lon, GT2.lat] = groundTrackODE(kep0, lonG0, time2);
[GT3.alpha, GT3.delta, GT3.lon, GT3.lat] = groundTrackODE(kep0, lonG0, time3);

% J2
% [GT1J2.alpha, GT1J2.delta, GT1J2.lon, GT1J2.lat] = groundTrackJ2(kep0, lonG0, time1);
% [GT2J2.alpha, GT2J2.delta, GT2J2.lon, GT2J2.lat] = groundTrackJ2(kep0, lonG0, time2);
% [GT3J2.alpha, GT3J2.delta, GT3J2.lon, GT3J2.lat] = groundTrackJ2(kep0, lonG0, time3);

% J2 + SRP
[GT1J2.alpha, GT1J2.delta, GT1J2.lon, GT1J2.lat] = groundTrackJ2SRP(kep0, lonG0, time1, body);
[GT2J2.alpha, GT2J2.delta, GT2J2.lon, GT2J2.lat] = groundTrackJ2SRP(kep0, lonG0, time2, body);
[GT3J2.alpha, GT3J2.delta, GT3J2.lon, GT3J2.lat] = groundTrackJ2SRP(kep0, lonG0, time3, body);

planisphere
GT1.p = plotGroundTrack(rad2deg(GT1.lon),rad2deg(GT1.lat),'g');
GT1J2.p = plotGroundTrack(rad2deg(GT1J2.lon),rad2deg(GT1J2.lat),'r');
title('Ground track over one orbit')
legend('start','end','unperturbed','start','end','perturbed')

planisphere
% yticks(-90:1:90)
% xticks(-180:1:180)
GT2.p = plotGroundTrack(rad2deg(GT2.lon),rad2deg(GT2.lat),'g');
GT2J2.p = plotGroundTrack(rad2deg(GT2J2.lon),rad2deg(GT2J2.lat),'r');
title("Ground track over " + t2 + " days")
legend('start','end','unperturbed','start','end','perturbed')

planisphere
GT3.p = plotGroundTrack(rad2deg(GT3.lon),rad2deg(GT3.lat),'g');
GT3J2.p = plotGroundTrack(rad2deg(GT3J2.lon),rad2deg(GT3J2.lat),'r');
title("Ground track over " + t3 + " days")
legend('start','end','unperturbed','start','end','perturbed')

fprintf('\nInitial ground tracks evaluation: ')
toc

%% Repeating ground tracks

m = 1;
k = 2;
aRepJ2 = repeatingGTJ2(m,k,a,e,i);

kep1 = kep0;
kep1(1) = aRepJ2;

[GT3rep.alpha, GT3rep.delta, GT3rep.lon, GT3rep.lat] = groundTrackODE(kep1, lonG0, time3);
% [GT3J2rep.alpha, GT3J2rep.delta, GT3J2rep.lon, GT3J2rep.lat] = groundTrackJ2(kep1, lonG0, time3);
[GT3J2rep.alpha, GT3J2rep.delta, GT3J2rep.lon, GT3J2rep.lat] = groundTrackJ2SRP(kep1, lonG0, time3, body);

planisphere
GT3J2rep.p = plotGroundTrack(rad2deg(GT3J2rep.lon),rad2deg(GT3J2rep.lat),'r');
GT3rep.p = plotGroundTrack(rad2deg(GT3rep.lon),rad2deg(GT3rep.lat),'g');

title("Repeating ground track over " + t3 + " days")
legend('start','end','perturbed','start','end','unperturbed')

%%

kep1 = kep0;

[rr,vv] = par2car(kep1);

opts = odeset('Reltol',1e-12,'AbsTol',1e-15);

Re = 6378.137;
J2 = 0.00108263;
omE = deg2rad(15.04/3600);

time = linspace(0,seconds(days(t3)),npts3);
t1day = linspace(0,seconds(days(1)),npts1);

%% Orbital parameters propagation (every combination of perturbations)

% initializing structures to speed up the for cycle
case0.a = zeros(length(time),1);
case0.e = zeros(length(time),1);
case0.i = zeros(length(time),1);
case0.Omega = zeros(length(time),1);
case0.omega = zeros(length(time),1);
case0.theta = zeros(length(time),1);
case1 = case0;
case2 = case0;
case3 = case0;

[case0.time,case0.cart] = ode113(@(t,y) orbit(t,y,mu), time, [rr,vv], opts); % no perturbation

[case1.time,case1.cart] = ode113(@(t,y) orbitJ2(t,y,mu), time, [rr,vv], opts); % J2 only

[case2.time,case2.cart] = ode113(@(t,y) orbitSRP(t,y,mu,Re,omE,body), time, [rr,vv], opts); % SRP only

tic

[case31day.time,case31day.cart] = ode113(@(t,y) ode_2bodyPerturb(t, y, mu, perturbations, 'cart'), t1day, [rr,vv], opts); % SRP + J2, cartesian

[case3.time,case3.cart] = ode113(@(t,y) ode_2bodyPerturb(t, y, mu, perturbations, 'cart'), time, [rr,vv], opts); % SRP + J2, cartesian
fprintf('\nCartesian J2 + SRP: ')
toc

tic

[case41day.time,case41day.kep] = ode113(@(t,y) ode_2bodyPerturb(t, y, mu, perturbations, 'gauss'), t1day, kep1, opts); % SRP + J2, gauss

[case4.time,case4.kep] = ode113(@(t,y) ode_2bodyPerturb(t, y, mu, perturbations, 'gauss'), time, kep1, opts); % SRP + J2, gauss
fprintf('\nGauss VOP J2 + SRP: ')
toc

%% Converting everything into Keplerian parameters

tic

for jj = 1:length(case3.time)
    
    [case0.a(jj),case0.e(jj),case0.i(jj),case0.Omega(jj),case0.omega(jj),case0.theta(jj)] = car2par(case0.cart(jj,1:3),case0.cart(jj,4:6));    
    [case1.a(jj),case1.e(jj),case1.i(jj),case1.Omega(jj),case1.omega(jj),case1.theta(jj)] = car2par(case1.cart(jj,1:3),case1.cart(jj,4:6));    
    [case2.a(jj),case2.e(jj),case2.i(jj),case2.Omega(jj),case2.omega(jj),case2.theta(jj)] = car2par(case2.cart(jj,1:3),case2.cart(jj,4:6));    
    [case3.a(jj),case3.e(jj),case3.i(jj),case3.Omega(jj),case3.omega(jj),case3.theta(jj)] = car2par(case3.cart(jj,1:3),case3.cart(jj,4:6)); 
       
end

for jj = 1:length(case31day.time)
    [case31day.a(jj),case31day.e(jj),case31day.i(jj),case31day.Omega(jj),case31day.omega(jj),case31day.theta(jj)] = car2par(case31day.cart(jj,1:3),case31day.cart(jj,4:6));
end

case4.a = case4.kep(:,1);
case4.e = case4.kep(:,2);
case4.i = case4.kep(:,3);
case4.Omega = case4.kep(:,4);
case4.omega = case4.kep(:,5);
case4.theta = case4.kep(:,6);

case41day.a = case41day.kep(:,1);
case41day.e = case41day.kep(:,2);
case41day.i = case41day.kep(:,3);
case41day.Omega = case41day.kep(:,4);
case41day.omega = case41day.kep(:,5);
case41day.theta = case41day.kep(:,6);

case1.Omega = anglePostProcessingShift(case1.Omega);
case3.Omega = anglePostProcessingShift(case3.Omega);
case4.theta = anglePostProcessingMod(case4.theta);
case31day.Omega = anglePostProcessingShift(case31day.Omega);
case41day.theta = anglePostProcessingMod(case41day.theta);

fprintf('\nKeplerian parameters conversion: ')
toc

%% Comparing results for J2 + SRP, Cartesian and Gauss VOP (nominal values, one day)

tic

figure
sgtitle('Unfiltered data in time domain (nominal values, one day)')
subplot(3,2,1),
hold on
plot(case31day.time/86400,case31day.a,"LineWidth",2.5)
plot(case41day.time/86400,case41day.a,"LineWidth",1.5)
xlabel('Time [days]')
ylabel('a [km]')
axis tight, grid on, title('Semi-major axis')
legend('Cartesian','Gaussian')

subplot(3,2,3), hold on
plot(case31day.time/86400,case31day.e,"LineWidth",2.5)
plot(case41day.time/86400,case41day.e,"LineWidth",1.5)
xlabel('Time [days]')
ylabel('e [-]')
axis tight, grid on, title('Eccentricity')
legend('Cartesian','Gaussian')

subplot(3,2,5), hold on
plot(case31day.time/86400,rad2deg(case31day.i),"LineWidth",2.5)
plot(case41day.time/86400,rad2deg(case41day.i),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('i [deg]')
axis tight, grid on, title('Inclination')
legend('Cartesian','Gaussian')

subplot(3,2,2), hold on
plot(case31day.time/86400,rad2deg(case31day.Omega),"LineWidth",2.5)
plot(case41day.time/86400,rad2deg(case41day.Omega),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('\Omega [deg]')
axis tight, grid on, title('RAAN')
legend('Cartesian','Gaussian')

subplot(3,2,4), hold on
plot(case31day.time/86400,rad2deg(case31day.omega),"LineWidth",2.5)
plot(case41day.time/86400,rad2deg(case41day.omega),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('\omega [deg]')
axis tight, grid on, title('Argument of perigee')
legend('Cartesian','Gaussian')

subplot(3,2,6), hold on
plot(case31day.time/86400,rad2deg(case31day.theta),"LineWidth",2.5)
plot(case41day.time/86400,rad2deg(case41day.theta),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('\theta [deg]')
yticks([0 90 180 270 360])
axis tight, grid on, title('True anomaly')
legend('Cartesian','Gaussian')

%% Comparing results for J2 + SRP, Cartesian and Gauss VOP (nominal values)

figure
sgtitle('Unfiltered data in time domain (nominal values)')
subplot(3,2,1),
hold on
plot(case3.time/86400,case3.a)
plot(case4.time/86400,case4.a)
xlabel('Time [days]')
ylabel('a [km]')
axis tight, grid on, title('Semi-major axis')
legend('Cartesian','Gaussian')
xticks(0:30:360)

subplot(3,2,3), hold on
plot(case3.time/86400,case3.e)
plot(case4.time/86400,case4.e)
xlabel('Time [days]')
ylabel('e [-]')
axis tight, grid on, title('Eccentricity')
legend('Cartesian','Gaussian')
xticks(0:30:360)

subplot(3,2,5), hold on
plot(case3.time/86400,rad2deg(case3.i))
plot(case4.time/86400,rad2deg(case4.i))
xlabel('Time [days]')
ylabel('i [deg]')
axis tight, grid on, title('Inclination')
legend('Cartesian','Gaussian')
xticks(0:30:360)

subplot(3,2,2), hold on
plot(case3.time/86400,rad2deg(case3.Omega))
plot(case4.time/86400,rad2deg(case4.Omega))
xlabel('Time [days]')
ylabel('\Omega [deg]')
axis tight, grid on, title('RAAN')
legend('Cartesian','Gaussian')
xticks(0:30:360)

subplot(3,2,4), hold on
plot(case3.time/86400,rad2deg(case3.omega))
plot(case4.time/86400,rad2deg(case4.omega))
xlabel('Time [days]')
ylabel('\omega [deg]')
axis tight, grid on, title('Argument of perigee')
legend('Cartesian','Gaussian')
xticks(0:30:360)

subplot(3,2,6), hold on
plot(case3.time/86400,rad2deg(case3.theta))
plot(case4.time/86400,rad2deg(case4.theta))
xlabel('Time [days]')
ylabel('\theta [deg]')
axis tight, grid on, title('True anomaly')
yticks([0 90 180 270 360])
legend('Cartesian','Gaussian')
xticks(0:30:360)

%% Comparing results for J2 + SRP, Cartesian and Gauss VOP (relative values)

figure
sgtitle('Unfiltered data in time domain (relative values)')
subplot(3,2,1),
hold on
plot(case3.time/86400,case3.a-case4.a)
xlabel('Time [days]')
ylabel('a [km]')
axis tight, grid on, title('Semi-major axis error')
xticks(0:30:360)

subplot(3,2,3), hold on
plot(case3.time/86400,case3.e-case4.e)
xlabel('Time [days]')
ylabel('e [-]')
axis tight, grid on, title('Eccentricity error')
xticks(0:30:360)

subplot(3,2,5), hold on
plot(case3.time/86400,rad2deg(case3.i-case4.i))
xlabel('Time [days]')
ylabel('i [deg]')
axis tight, grid on, title('Inclination error')
xticks(0:30:360)

subplot(3,2,2), hold on
plot(case3.time/86400,rad2deg(case3.Omega-case4.Omega))
xlabel('Time [days]')
ylabel('\Omega [deg]')
axis tight, grid on, title('RAAN error')
xticks(0:30:360)

subplot(3,2,4), hold on
plot(case3.time/86400,rad2deg(case3.omega-case4.omega))
xlabel('Time [days]')
ylabel('\omega [deg]')
axis tight, grid on, title('Argument of perigee error')
xticks(0:30:360)

subplot(3,2,6), hold on
plot(case3.time(2:end)/86400,rad2deg(case3.theta(2:end)-case4.theta(2:end)))
xlabel('Time [days]')
ylabel('\theta [deg]')
axis tight, grid on, title('True anomaly error')
xticks(0:30:360)

%% Spectral analysis

figure
sgtitle('Unfiltered data in frequency domain')
subplot(3,2,1)
[case3.frequencyRange,case3.aAmplitude] = spectralAnalysis(case3.a,case3.time); hold on
[case4.frequencyRange,case4.aAmplitude] = spectralAnalysis(case4.a,case4.time);
title('Semi-major axis')
legend('Cartesian','Gaussian')
ylabel('Amplitude [km]')
%xlim([-6e-3,6e-3])

subplot(3,2,3)
[case3.frequencyRange,case3.eAmplitude] = spectralAnalysis(case3.e,case3.time); hold on
[case4.frequencyRange,case4.eAmplitude] = spectralAnalysis(case4.e,case4.time);
title('Eccentricity')
legend('Cartesian','Gaussian')
ylabel('Amplitude [-]')
%xlim([-6e-3,6e-3])

subplot(3,2,5)
[case3.frequencyRange,case3.iAmplitude] = spectralAnalysis(rad2deg(case3.i),case3.time); hold on
[case4.frequencyRange,case4.iAmplitude] = spectralAnalysis(rad2deg(case4.i),case4.time);
title('Inclination')
legend('Cartesian','Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

subplot(3,2,2)
[case3.frequencyRange,case3.OmegaAmplitude] = spectralAnalysis(rad2deg(case3.Omega),case3.time); hold on
[case4.frequencyRange,case4.OmegaAmplitude] = spectralAnalysis(rad2deg(case4.Omega),case4.time);
title('RAAN')
legend('Cartesian','Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

subplot(3,2,4)
[case3.frequencyRange,case3.omegaAmplitude] = spectralAnalysis(rad2deg(case3.omega),case3.time); hold on
[case4.frequencyRange,case4.omegaAmplitude] = spectralAnalysis(rad2deg(case4.omega),case4.time);
title('Argument of perigee')
legend('Cartesian','Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

subplot(3,2,6)
[case3.frequencyRange,case3.thetaAmplitude] = spectralAnalysis(rad2deg(case3.theta),case3.time); hold on
[case4.frequencyRange,case4.thetaAmplitude] = spectralAnalysis(rad2deg(case4.theta),case4.time);
title('True anomaly')
legend('Cartesian','Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

% % This shows that peaks correspond to multiples of the period.
% period = 2*pi*sqrt(case3.a(1)^3/mu);
% 
% for jj = 1:6
%     
%     subplot(3,2,jj)
%     for ii = 1:20
% 
%         plot([ii/period,ii/period],[1e-8,1e3],'Color','r')
% 
%     end
% end

%% Filtering

% 25 min -> 700 points = 1/(0.5182 days period) * 365
% 5 min -> 5000 points

case3.aFilt = MAF(case3.a,700);
case3.eFilt = MAF(case3.e,700);
case3.iFilt = MAF(case3.i,700);
case3.OmegaFilt = MAF(case3.Omega,700);
case3.omegaFilt = MAF(case3.omega,700);
case3.thetaFilt = MAF(case3.theta,700);

case4.aFilt = MAF(case4.a,700);
case4.eFilt = MAF(case4.e,700);
case4.iFilt = MAF(case4.i,700);
case4.OmegaFilt = MAF(case4.Omega,700);
case4.omegaFilt = MAF(case4.omega,700);
case4.thetaFilt = MAF(case4.theta,700);

%% Comparing filtered and unfiltered in frequency domain

figure
sgtitle('Filtered and unfiltered data in frequency domain')
subplot(3,2,1)
[case3.frequencyRange,case3.aAmplitude] = spectralAnalysis(case3.a,case3.time); hold on
[case4.frequencyRange,case4.aAmplitude] = spectralAnalysis(case4.a,case4.time);
[case3.frequencyRangeFilt,case3.aAmplitudeFilt] = spectralAnalysis(case3.aFilt,case3.time,2.5);
[case4.frequencyRangeFilt,case4.aAmplitudeFilt] = spectralAnalysis(case4.aFilt,case4.time,1.5);
title('Semi-major axis')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
ylabel('Amplitude [km]')
%xlim([-6e-3,6e-3])

subplot(3,2,3)
[case3.frequencyRange,case3.eAmplitude] = spectralAnalysis(case3.e,case3.time); hold on
[case4.frequencyRange,case4.eAmplitude] = spectralAnalysis(case4.e,case4.time);
[case3.frequencyRangeFilt,case3.eAmplitudeFilt] = spectralAnalysis(case3.eFilt,case3.time,2.5);
[case4.frequencyRangeFilt,case4.eAmplitudeFilt] = spectralAnalysis(case4.eFilt,case4.time,1.5);
title('Eccentricity')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
ylabel('Amplitude [-]')
%xlim([-6e-3,6e-3])

subplot(3,2,5)
[case3.frequencyRange,case3.iAmplitude] = spectralAnalysis(rad2deg(case3.i),case3.time); hold on
[case4.frequencyRange,case4.iAmplitude] = spectralAnalysis(rad2deg(case4.i),case4.time);
[case3.frequencyRangeFilt,case3.iAmplitudeFilt] = spectralAnalysis(rad2deg(case3.iFilt),case3.time,2.5);
[case4.frequencyRangeFilt,case4.iAmplitudeFilt] = spectralAnalysis(rad2deg(case4.iFilt),case4.time,1.5);
title('Inclination')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

subplot(3,2,2)
[case3.frequencyRange,case3.OmegaAmplitude] = spectralAnalysis(rad2deg(case3.Omega),case3.time); hold on
[case4.frequencyRange,case4.OmegaAmplitude] = spectralAnalysis(rad2deg(case4.Omega),case4.time);
[case3.frequencyRangeFilt,case3.OmegaAmplitudeFilt] = spectralAnalysis(rad2deg(case3.OmegaFilt),case3.time,2.5);
[case4.frequencyRangeFilt,case4.OmegaAmplitudeFilt] = spectralAnalysis(rad2deg(case4.OmegaFilt),case4.time,1.5);
title('RAAN')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

subplot(3,2,4)
[case3.frequencyRange,case3.omegaAmplitude] = spectralAnalysis(rad2deg(case3.omega),case3.time); hold on
[case4.frequencyRange,case4.omegaAmplitude] = spectralAnalysis(rad2deg(case4.omega),case4.time);
[case3.frequencyRangeFilt,case3.omegaAmplitudeFilt] = spectralAnalysis(rad2deg(case3.omegaFilt),case3.time,2.5);
[case4.frequencyRangeFilt,case4.omegaAmplitudeFilt] = spectralAnalysis(rad2deg(case4.omegaFilt),case4.time,1.5);
title('Argument of perigee')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

subplot(3,2,6)
[case3.frequencyRange,case3.thetaAmplitude] = spectralAnalysis(rad2deg(case3.theta),case3.time,0.5); hold on
[case4.frequencyRange,case4.thetaAmplitude] = spectralAnalysis(rad2deg(case4.theta),case4.time,0.5);
[case3.frequencyRangeFilt,case3.thetaAmplitudeFilt] = spectralAnalysis(rad2deg(case3.thetaFilt),case3.time,2.5);
[case4.frequencyRangeFilt,case4.thetaAmplitudeFilt] = spectralAnalysis(rad2deg(case4.thetaFilt),case4.time,1.5);
title('True anomaly')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
ylabel('Amplitude [deg]')
%xlim([-6e-3,6e-3])

%% Comparing filtered and unfiltered in time domain

figure
sgtitle('Filtered and unfiltered data in time domain')
subplot(3,2,1),
hold on
plot(case3.time/86400,case3.a)
plot(case4.time/86400,case4.a)
plot(case3.time/86400,case3.aFilt,"LineWidth",2.5)
plot(case4.time/86400,case4.aFilt,"LineWidth",1.5)
xlabel('Time [days]')
ylabel('a [km]')
axis tight, grid on, title('Semi-major axis')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
xticks(0:30:360)

subplot(3,2,3), hold on
plot(case3.time/86400,case3.e)
plot(case4.time/86400,case4.e)
plot(case3.time/86400,case3.eFilt,"LineWidth",2.5)
plot(case4.time/86400,case4.eFilt,"LineWidth",1.5)
xlabel('Time [days]')
ylabel('e [-]')
axis tight, grid on, title('Eccentricity')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
xticks(0:30:360)

subplot(3,2,5), hold on
plot(case3.time/86400,rad2deg(case3.i))
plot(case4.time/86400,rad2deg(case4.i))
plot(case3.time/86400,rad2deg(case3.iFilt),"LineWidth",2.5)
plot(case4.time/86400,rad2deg(case4.iFilt),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('i [deg]')
axis tight, grid on, title('Inclination')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
xticks(0:30:360)

subplot(3,2,2), hold on
plot(case3.time/86400,rad2deg(case3.Omega))
plot(case4.time/86400,rad2deg(case4.Omega))
plot(case3.time/86400,rad2deg(case3.OmegaFilt),"LineWidth",2.5)
plot(case4.time/86400,rad2deg(case4.OmegaFilt),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('\Omega [deg]')
axis tight, grid on, title('RAAN')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
xticks(0:30:360)

subplot(3,2,4), hold on
plot(case3.time/86400,rad2deg(case3.omega))
plot(case4.time/86400,rad2deg(case4.omega))
plot(case3.time/86400,rad2deg(case3.omegaFilt),"LineWidth",2.5)
plot(case4.time/86400,rad2deg(case4.omegaFilt),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('\omega [deg]')
axis tight, grid on, title('Argument of perigee')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
xticks(0:30:360)

subplot(3,2,6), hold on
plot(case3.time/86400,rad2deg(case3.theta))
plot(case4.time/86400,rad2deg(case4.theta))
plot(case3.time/86400,rad2deg(case3.thetaFilt),"LineWidth",2.5)
plot(case4.time/86400,rad2deg(case4.thetaFilt),"LineWidth",1.5)
xlabel('Time [days]')
ylabel('\theta [deg]')
axis tight, grid on, title('True anomaly')
legend('Cartesian','Gaussian','Filtered Cartesian','Filtered Gaussian')
yticks([0 90 180 270 360])
xticks(0:30:360)

fprintf('\nComparing results for J2 + SRP Cartesian and Gauss VOP: ')
toc

%% Orbit representation

% figure, hold on
% Earth3d
% colorPlot3(case3.time(1:end),case3.cart(1:end,:),jet)
% c = colorbar;
% c.Label.String = 'Time [days]';
% 
% c.Ticks = ([3 258*1/10 258*2/10 258*3/10 258*4/10 258*5/10 258*6/10 258*7/10 258*8/10 258*9/10 255]);
% c.TickLabels = ([0 36 36*2 36*3 36*4 36*5 36*6 36*7 36*8 36*9 360]);
% 
% grid on
% pbaspect([1 1 1])
% title('Orbit evolution over time, with J2 and SRP perturbations')
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% 
% view(40,20)

%% Comparison with real data

tic

opts = odeset('Reltol',1e-12,'AbsTol',1e-15);
mu = 398600;
C_R = 1;
A2m = 5;
perturbations = {1 [C_R A2m]};

filename10 = 'horizons_results_10d_2min.txt';

headerlinesOut10 = 1285;
delimiterOut10 = ',';

X10 = importdata(filename10,delimiterOut10,headerlinesOut10);

x10 = X10.data;
a10 = x10(:,10);
e10 = x10(:,1);
i10 = x10(:,3);
Om10 = x10(:,4);
om10 = x10(:,5);
th10 = x10(:,9);

time10 = linspace(0,seconds(days(10)),length(a10));

kep10 = [a10(1),e10(1),deg2rad(i10(1)),deg2rad(Om10(1)),deg2rad(om10(1)),deg2rad(th10(1))];
[rr10,vv10] = par2car(kep10);

[case10.time,case10.cart] = ode113(@(t,y) ode_2bodyPerturb(t, y, mu, perturbations, 'cart'), time10, [rr10,vv10], opts); % SRP + J2, cartesian

for jj = 1:length(case10.time)
    
    [case10.a(jj),case10.e(jj),case10.i(jj),case10.Omega(jj),case10.omega(jj),case10.theta(jj)] = car2par(case10.cart(jj,1:3),case10.cart(jj,4:6));

end

figure
sgtitle('Comparison with real data (10 days)')
subplot(3,2,1),
hold on
plot(case10.time/86400,case10.a)
plot(case10.time/86400,a10)
xlabel('Time [days]')
ylabel('a [km]')
axis tight, grid on, title('Semi-major axis')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')

subplot(3,2,3), hold on
plot(case10.time/86400,case10.e)
plot(case10.time/86400,e10)
xlabel('Time [days]')
ylabel('e [-]')
axis tight, grid on, title('Eccentricity')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')

subplot(3,2,5), hold on
plot(case10.time/86400,rad2deg(case10.i))
plot(case10.time/86400,i10)
xlabel('Time [days]')
ylabel('i [deg]')
axis tight, grid on, title('Inclination')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')

subplot(3,2,2), hold on
plot(case10.time/86400,rad2deg(case10.Omega))
plot(case10.time/86400,Om10)
xlabel('Time [days]')
ylabel('\Omega [deg]')
axis tight, grid on, title('RAAN')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')

subplot(3,2,4), hold on
plot(case10.time/86400,rad2deg(case10.omega))
plot(case10.time/86400,om10)
xlabel('Time [days]')
ylabel('\omega [deg]')
axis tight, grid on, title('Argument of perigee')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')

subplot(3,2,6), hold on
plot(case10.time/86400,rad2deg(case10.theta))
plot(case10.time/86400,th10)
xlabel('Time [days]')
ylabel('\theta [deg]')
yticks([0 90 180 270 360])
axis tight, grid on, title('True anomaly'), ylim([0,360]) 
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')





filename365 = 'horizons_results_1y_25min.txt';

headerlinesOut365 = 1285;
delimiterOut365 = ',';

X365 = importdata(filename365,delimiterOut365,headerlinesOut365);

x365 = X365.data;
a365 = x365(:,10);
e365 = x365(:,1);
i365 = x365(:,3);
Om365 = x365(:,4);
om365 = x365(:,5);
th365 = x365(:,9);

time365 = linspace(0,seconds(days(365)),length(a365));

kep365 = [a365(1),e365(1),deg2rad(i365(1)),deg2rad(Om365(1)),deg2rad(om365(1)),deg2rad(th365(1))];
[rr365,vv365] = par2car(kep365);

[case365.time,case365.cart] = ode113(@(t,y) ode_2bodyPerturb(t, y, mu, perturbations, 'cart'), time365, [rr365,vv365], opts); % SRP + J2, cartesian

for jj = 1:length(case365.time)
    
    [case365.a(jj),case365.e(jj),case365.i(jj),case365.Omega(jj),case365.omega(jj),case365.theta(jj)] = car2par(case365.cart(jj,1:3),case365.cart(jj,4:6));

end

figure
sgtitle('Comparison with real data (365 days)')
subplot(3,2,1),
hold on
plot(case365.time/86400,case365.a)
plot(case365.time/86400,a365)
xlabel('Time [days]')
ylabel('a [km]')
axis tight, grid on, title('Semi-major axis')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')
xticks(0:30:360)

subplot(3,2,3), hold on
plot(case365.time/86400,case365.e)
plot(case365.time/86400,e365)
xlabel('Time [days]')
ylabel('e [-]')
axis tight, grid on, title('Eccentricity')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')
xticks(0:30:360)

subplot(3,2,5), hold on
plot(case365.time/86400,rad2deg(case365.i))
plot(case365.time/86400,i365)
xlabel('Time [days]')
ylabel('i [deg]')
axis tight, grid on, title('Inclination')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')
xticks(0:30:360)

subplot(3,2,2), hold on
Om = rad2deg(case365.Omega);
for j = 1:length(Om)
    if Om(j) > 180
        Om(j) = Om(j) - 360;
    end
end
for j = 1:length(Om365)
    if Om365(j) > 180
        Om365(j) = Om365(j) - 360;
    end
end

plot(case365.time/86400,Om)
plot(case365.time/86400,Om365)
xlabel('Time [days]')
ylabel('\Omega [deg]')
axis tight, grid on, title('RAAN')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')
xticks(0:30:360)

subplot(3,2,4), hold on
plot(case365.time/86400,rad2deg(case365.omega))
plot(case365.time/86400,om365)
xlabel('Time [days]')
ylabel('\omega [deg]')
axis tight, grid on, title('Argument of perigee')
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')
xticks(0:30:360)

subplot(3,2,6), hold on
plot(case365.time/86400,rad2deg(case365.theta))
plot(case365.time/86400,th365)
xlabel('Time [days]')
ylabel('\theta [deg]')
axis tight, grid on, title('True anomaly'), ylim([0,360]) 
yticks([0 90 180 270 360])
legend('Numerical propagation','NASA-HORIZONS ephemerides generator')
xticks(0:30:360)

fprintf('\nComparison with real data: ')
toc

%% Orbit evolution video

OrbitEvolutionVideo

fprintf('\nTotal ')
toc(totalT)