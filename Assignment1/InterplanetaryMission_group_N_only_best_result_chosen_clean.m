%% INTERPLANETARY MISSION
%--------------------------------------------------------------------------
% GROUP 2304
% commento da parte degli amici su cosa fa

%--------------------------------------------------------------------------

addpath(genpath('MATLAB functions\')), addpath('textures\')
clc, clear, close all, format long g
%% SETUP for plots and initial data

N = 30;             % N+1 units of Dvtot will be shown in Porkchop (Default:  20)
STEP = 2;           % Steps between Porkchop lines in Dvtot units  (Default:  4)

DD = 48;            % Sinodic period discretization to refine time grid 
D = 8;              % Chosen window to show only one result (from 1 to 48)

% True anomaly step for discretization of time windows [deg]
% dth1 --> first arc || dth2 --> second arc
dth1 = [2,2];          
dth2 = [dth1(2),2]; 

tmax = date2mjd2000([2058,1,1,0,0,0]);  
mu_S = astroConstants(4);   % Sun gravitational parameter   [km^3/s^2]
ID = [1,2,81];              % IDs of bodies to use in uplanet.m | ephNEO.m
G = astroConstants(1); % Universal gavity constant
mS = mu_S/G;
mu_V = astroConstants(12);
mV = mu_V/G;
R_V = astroConstants(22);

tic

% Retrieve keplerian parameters for all bodies with their orbital period
for s = 1:3
    if ID(s) < 12
        kep(s,:) = uplanet(0,ID(s));
    else
        kep(s,:) = ephNEO(0,ID(s));
    end
    a = kep(s,1);                           
    e = kep(s,2);
    T(s) = 2*pi*sqrt(a^3/mu_S)/(3600*24);   % 
end     

%% Lambert arcs for minimum dt (parabolic)
% from a - dt plot

% orbitType = 0; % (0: prograde, 1: retrograde)
% for ARC = 1:2
%     a1 = kep(ARC,1);a2 = kep(ARC+1,1);
%     e1 = kep(ARC,2);e2 = kep(ARC+1,2);
%     i1 = kep(ARC,3);i2 = kep(ARC+1,3);
%     Om1 = kep(ARC,4);Om2 = kep(ARC+1,4);
%     om1 = kep(ARC,5);om2 = kep(ARC+1,5);
% 
%     jj = 0:5:360; %jj = 0; % every n degs 
%     kk = 0:5:360;% kk = 50;
%     mm = 0:.1:25; %T(ARC+1); % every n days
% 
% figure
%     J = 0;
%     for j = jj
%         J = J+1;
%         K = 0;
%         for k = kk
%             K = K+1;
%             M = 0;
%             for m = mm 
%                 M = M+1;
%                 for s = 1:2
%                     switch s
%                         case 1
%                             a = a1; e = e1; i = i1; Om = Om1; om = om1;
%                             theta = j;
%                         case 2
%                             a = a2; e = e2; i = i2; Om = Om2; om = om2;
%                             theta = k;
%                     end
%                 
%                     rvect(:,s) = kep2car(a,e,i,Om,om,theta*pi/180,mu_S);
%                 end
%                     RI = rvect(:,1); % [km]     position 1
%                     RF = rvect(:,2); % [km]     position 2
%                     TOF = m*3600*24;
% 
%                     [A(J,K,M),~,~,~,~,~,TPAR(J,K,M),~] = lambertMR(RI,RF,TOF,mu_S,orbitType);
%             end
% %             TOF = mm;
% %             AA = squeeze(A(J,K,:));
% %             plot(AA,TOF), grid on, hold on
% %             xlim([0,6e8])
% %             ylim([0,25])
%         end
%     end
%     TPARmin(ARC) = min(TPAR,[],'all');
% end
% toc
%
%--------
% Result: dtparabolic from plot of horizontal asymptote
% absolute minimum

dtp1 = 8; % days
dtp2 = 5;

%% SYNODIC PERIOD 3 planets
% by Borislav Borisov - 22.09.2011
% "Three-planet resonances in the Solar system"

a1 = kep(1,1); a2 = kep(2,1); a3 = kep(3,1);
Tsyn12 = 1/abs(1/T(2)-1/T(1)); % first arc
Tsyn23 = 1/abs(1/T(3)-1/T(2)); % second arc
Tsyn13 = 1/abs(1/T(3)-1/T(1)); % second arc

%
Tsynmax = max([Tsyn12,Tsyn23,Tsyn13]); % [days]
% P1 < P2 < P3
P1 = T(1);
P2 = T(2);
P3 = T(3);
P13 = Tsyn13;

n1 = sqrt(mu_S/a1^3); n2 = sqrt(mu_S/a2^3); n3 = sqrt(mu_S/a3^3); 
ratio = (n1-n2)/(n2-n3);

i = 0; j = 0; 
accuracy = 2; accuracy_ = 1;
fprintf('\nSYNODIC PERIOD 3 planets\nby Borislav Borisov:\naccuracy =          ') 
err = 0.01;
tol = 0.005;

while accuracy - err > tol || accuracy > accuracy_
    if ratio > 1
        i = i+1;
        j = round(i*ratio);
    elseif ratio <= 1
        j = j+1;
        i = round(j/ratio);
    end

    Pij = (i+j)*P13/365.2417;
    P1_ = i/((i+j)/P2-j/P3);
    dP1 = P1_-P1;
    accuracy = abs(dP1)/P1;

    fprintf('\b\b\b\b\b\b\b\b\b')
    fprintf('%.6f\n',accuracy)
    
    if accuracy - err <= tol
        if ratio > 1
            i_ = i+1;
            j_ = round(i_*ratio);
        elseif ratio <= 1
            j_ = j+1;
            i_ = round(j_/ratio);
        end
        Pij_ = (i+j)*P13/365.2417;
        P1_ = i_/((i_+j_)/P2-j_/P3);
        dP1_ = P1_-P1;
        accuracy_ = abs(dP1_)/P1;
    end
end
Tsyn123 = Pij*365.2417;
fprintf('\b < %.2f\ni:j = %d:%d\nTsyn123 = %.4fy = %.2fd\n',err,i,j,Pij,Tsyn123)

%% GRID SEARCH
DVTOT = [];
    % close all
    % clearvars -except D DD DVTOT t2_ N STEP TH T tmax Tsyn123 dtp1 dtp2 mu_S tic0 dth1 dth2
    fprintf('\n# %d\n',D)
for arc = 1:2
    switch arc 
        case 1
% MERCURY - VENUS
            dth = dth1;                                              
            ID = [1,2];
            t1 = [2028,1,1,0,0,0]; 
            t1 = mjd20002date(date2mjd2000(t1)+(D-1)*Tsyn123/DD);
            t2 = mjd20002date(date2mjd2000(t1)+Tsyn123/DD);

            tw1 = [t1;t2];
            tw1_ = tw1;

            ToF = [1.2*dtp1/365.2417;T(2)/365.2417];
        case 2
% VENUS - NEO n.81 2005EL70
            dth = dth2;
            ID = [2,81];
            t1 = mjd20002date(date2mjd2000(t1)+ToFG(1)*365.2417);
            t2 = mjd20002date(date2mjd2000(t2)+ToFG(2)*365.2417);
            tw1 = [t1; t2];
            ToF = [1.2*dtp2/365.2417;T(3)/365.2417];
    end
% days: time step for tw1 & tw2 & ToF corresponding to "dth" average true anomaly
    dt = 1/360*[dth(1)*T(arc),dth(2)*T(arc+1),20*360];
    t3 = mjd20002date(date2mjd2000(t1)+ToF(1)*365.2417);
    t4 = mjd20002date(date2mjd2000(t2)+ToF(2)*365.2417);
    tw2 = [t3;t4];


% NEO n.81 2005EL70

% NUMBER OF ELEMENTS
    for s = 1:3
        switch s
            case 1
                ti = date2mjd2000(tw1(1,:));
                tf = date2mjd2000(tw1(2,:));
            case 2
                ti = date2mjd2000(tw2(1,:));
                tf = date2mjd2000(tw2(2,:));
            case 3
                ti = ToF(1)*365.2417; 
                tf = ToF(2)*365.2417;
        end
  
        n(s) = ceil((tf-ti)/dt(s));
    end

% PORKCHOP PLOT
TOFMAX = T(arc+1);
[Dv1,Dv2,Dv1_,Dv2_,DV1,DV2] = porkchop(tw1,ToF,n,ID,arc,TOFMAX);
Dvmin1 = min(Dv1(:));  % C3 departure / vinf_plus , tw2 mode
Dvmin1_ = min(Dv1_(:));  % ToF mode
Dvmin2 = min(Dv2(:));  % vinf_minus / C3 arrival, tw2 mode
Dvmin2_ = min(Dv2_(:));  % ToF mode

[row1,column1] = find(Dv1 == Dvmin1);
[row1_,column1_] = find(Dv1_ == Dvmin1_);
[row2,column2] = find(Dv2 == Dvmin2);
[row2_,column2_] = find(Dv2_ == Dvmin2_);

% porkchop domain setup

tw1G = tw1; % Gregorian
for j = 1:2
    t2max = tw2(j,1); % [day]
end
tw2G = tw2;
ToFG = ToF;

i = [5,3,2];
for K = 1:2
for j = 1:3
    if j < 3
        if K == 1
            Dvmin = Dvmin1;
        elseif K == 2
            Dvmin = Dvmin2;
        end
        figure(5*arc-4/K)
        nexttile(j)
        h = findobj(gcf);
        h1 = h(10-j);
        
    else
        if K == 1
            Dvmin = Dvmin1_;
        elseif K == 2
            Dvmin = Dvmin2_;
        end
        figure(5*arc-floor(3/K))
        h = findobj(gcf);
        h1 = h(4);
    end

    clim([floor(Dvmin)-1 floor(Dvmin)+N]), 
    h1.LevelList = floor(Dvmin)-1:STEP:floor(Dvmin)+N;
    c = h(i(j)); c.Label.Position = [-1.173322737216949,floor(Dvmin)+N+1,0];
    ax = h(i(j)+1);
    ax.XTick = ax.XTick(1:ceil(length(ax.XTick)/10):end);
    ax.XTickLabel = ax.XTickLabel(1:ceil(length(ax.XTickLabel)/10):end);
    ax.YTick = ax.YTick(1:ceil(length(ax.YTick)/10):end);
    ax.YTickLabel = ax.YTickLabel(1:ceil(length(ax.YTickLabel)/10):end);
    if t2max > tmax && j < 3 
        ax.YLim = ax.XLim;
        ax.YTick = ax.XTick;
        ax.YTickLabel = ax.XTickLabel;
    end
end

if K == 1
    column = column1;
    row = row1;
    column_ = column1_;
    row_ = row1_;
elseif K == 2
    column = column2;
    row = row2;
    column_ = column2_;
    row_ = row2_;
end

tw1 = linspace(datetime(tw1G(1,:)),datetime(tw1G(2,:)),n(1))'; 
tw2 = linspace(datetime(tw2G(1,:)),datetime(tw2G(2,:)),n(2))';
ToF = linspace(ToFG(1),ToFG(2),n(3))';
x01 = str2num(char(datetime([tw1(column);tw2(row)],'Format','yyyy,MM,dd,hh,mm,ss.SS')));
x02 = [str2num(char(datetime([tw1(column_)],'Format','yyyy,MM,dd,hh,mm,ss.SS')));ToF(row_),zeros(1,5)];

x1(1,K) = date2mjd2000(x01(1,:)); x1(2,K) = date2mjd2000(x01(2,:)); 
x2(1) = date2mjd2000(x02(1,:)); x2(2) = x02(2,1); %
if K == 1
    Dvmin = Dvmin1;
    Dvmin_ = Dvmin1_;
elseif K == 2
    Dvmin = Dvmin2;
    Dvmin_ = Dvmin2_;
end
for j = 1:2
    figure(5*arc-4/K)
    nexttile(j)
    p1 = plot3(x1(1,K),x1(2,K),Dvmin,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k','DisplayName','\Deltav initial guess unconstrained');
    clim([floor(Dvmin)-1 floor(Dvmin)+N])
end
figure(5*arc-floor(3/K)),p2 = plot3(x2(1),x2(2),Dvmin_,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k','DisplayName','\Deltav initial guess unconstrained');
clim([floor(Dvmin_)-1 floor(Dvmin_)+N])

figure(5*arc-4/K)
nexttile(2), h = findobj(gcf); c = h(8); c.LevelListMode = 'auto';

end

figure(5*arc),
plot3(x1(1,1),x1(2,1),Dvmin1,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k','DisplayName','\Deltav initial guess unconstrained')
zlim([Dvmin1-3 Dvmin1+N])
clim([floor(Dvmin1)-1 floor(Dvmin1)+N])

% Extract velocities --- quadratic arrays
    switch arc
        case 1  % ( n x m )
            C3_dep = Dv1;
            Vinf_minus = DV2;  % [3] vector per each cell
            nn = n(1);
            mm = n(2);
            DT1 = dt;
            DTH1 = dth;
        case 2  % ( m x l )
            Vinf_plus = DV1;
            C3_arr = Dv2;
            ll = n(2);
            DT2 = dt;
            DTH2 = dth;
    end
end

%-------------------------------------------------------------------------------
% Dates
tw_dep = tw1_; % [date]
tw_fb = tw1G;
tw_arr = tw2G;

for j = 1:2
    TW_dep(j) = date2mjd2000(tw_dep(j,:));
    TW_fb(j) = date2mjd2000(tw_fb(j,:));
    TW_arr(j) = date2mjd2000(tw_arr(j,:));
end

TW_dep = linspace(TW_dep(1),TW_dep(2),nn)'; % [days]
TW_fb = linspace(TW_fb(1),TW_fb(2),mm)';
TW_arr = linspace(TW_arr(1),TW_arr(2),ll)';

%------------------------------------------------------------------------------
% Recombinations of arrays
% m x ( n + l ) --- quadratic + linear / nonlinear operation

Dv_tot = zeros(ll,nn,mm); Dvflyby_tot = Dv_tot; 

for m = 1:mm
    C3_dep_fbdate = C3_dep(m,:);
    Vinf_minus_fbdate = Vinf_minus(m,:);
    Vinf_plus_fbdate = Vinf_plus(:,m);
    C3_arr_fbdate = C3_arr(:,m);

%----------
% linear operations

    Dv_tot(:,:,m) = C3_dep_fbdate + C3_arr_fbdate;
        for n = 1:nn
            for l = 1:ll
                Vinf_minus_FB = Vinf_minus_fbdate{n};
                Vinf_plus_FB = Vinf_plus_fbdate{l};

                if anynan(Vinf_minus_FB) || anynan(Vinf_plus_FB)
                    Dvflyby_tot(l,n,m) = NaN;
                    continue
                end
            Dvflyby_tot(l,n,m) = norm(Vinf_plus_FB - Vinf_minus_FB);
            end
        end
end

Dv_tot = Dv_tot + Dvflyby_tot;

%-----------------------------------------------------------------------------
% choose date M with minimum Dv_tot

h_atm_venus = 500; % [km]
rp = R_V+h_atm_venus; % inizial guess for fzero
KK = 0;
% REMOVE collision prevention COMPUTING SO TO SPEED UP THE PROCESS
while rp <= R_V + h_atm_venus   %--> a method to exclude the colliding solutions
KK = KK+1;
[Dv_totmin,index] = min(Dv_tot,[],'all');
[row,column_] = find(Dv_tot(index) == Dv_tot);

M = ceil(column_/nn);
column = mod(column_,nn);
if column == 0
    column = nn;
end

Dv_inf = Dvflyby_tot(row,column,M);


%----------------------------------------------------------------------------
% Flyby results after minimization
l = row;
m = M;
n = column;

% nonlinear operation: 
% Powered gravity assist fly-by algorithm --> rp --> dvp

% get rV for rSoi (sphere of influence)
time = TW_fb(m);   % mjd2000 
kep = uplanet(time,2);
a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
rV = kep2car(a,e,i,Omega,omega,theta,mu_S);
rV = norm(rV);

rSoi = rV*(mV/mS)^(2/5);

Vinf_minus_fbdate = Vinf_minus(M,:);
Vinf_plus_fbdate = Vinf_plus(:,M);
Vinf_minus_FB = Vinf_minus_fbdate{n};
Vinf_plus_FB = Vinf_plus_fbdate{l};

if anynan(Vinf_minus_FB) || anynan(Vinf_plus_FB)
    dvp = NaN;
    dtSoi = NaN;
    Dv_tot(index) = NaN; % crucial to restart minimization without this solution
    continue
end

a_minus = -mu_V/norm(Vinf_minus_FB)^2;
a_plus = -mu_V/norm(Vinf_plus_FB)^2;
rp = R_V+h_atm_venus; % initial guess close enough to solution (Venus radius)

options = optimset('TolX',1e-6,'Display','off');
rp = fzero(@fd,rp,options,Vinf_minus_FB,Vinf_plus_FB,mu_V);

if rp <= R_V + h_atm_venus
    dvp = NaN;
    rp = R_V;
    dtSoi = NaN;
    Dv_tot(index) = NaN;  % crucial to restart minimization without this solution
    continue
end
end
[~,~,e_minus,e_plus,~,~] = fd(rp,Vinf_minus_FB,Vinf_plus_FB,mu_V);

vp_minus = sqrt(norm(Vinf_minus_FB)^2+2*mu_V/rp);
vp_plus = sqrt(norm(Vinf_plus_FB)^2+2*mu_V/rp);
dvp = abs(vp_minus-vp_plus);

E_minus = acosh((-rSoi/a_minus+1)/e_minus);
E_plus = acosh((-rSoi/a_plus+1)/e_plus);

dt_minus = (e_minus*sinh(E_minus)-E_minus)/sqrt(-mu_V/a_minus^3);
dt_plus = (e_plus*sinh(E_plus)-E_plus)/sqrt(-mu_V/a_plus^3);

dtSoi = dt_minus + dt_plus;
dtSOI = seconds(dtSoi);
dtSOI.Format = 'hh:mm:ss.SS';

% Final sum
Dv_totmin = Dv_totmin - Dv_inf + dvp;


TW_DEP = datetime(mjd20002date(TW_dep(n)),'Format','yyyy MMM dd hh:mm:ss'); 
TW_FB = datetime(mjd20002date(TW_fb(m)),'Format','yyyy MMM dd hh:mm:ss'); 
TW_ARR = datetime(mjd20002date(TW_arr(l)),'Format','yyyy MMM dd hh:mm:ss'); 

if rp <= R_V + h_atm_venus
    str = '<';
    STR = 'Constraint NOT satisfied';
elseif rp > R_V+h_atm_venus
    str = '>';
    STR = 'Constraint satisfied';
else
    str = '?';
    STR = 'rp = undefined NaN';
end

%---------------------------------------------------------------------------

fprintf('\nDates from time grid: Tsyn123/%d\n%c%c: %d°, %d°, %d° - %ct: %d d, %d d, %d d\nDeparture: %s - %.4f days (MJD2000) \nFly-by: %s - %.4f days\nArrival: %s - %.4f days\n', ...
    DD,char(916),char(952),DTH1(1),DTH1(2),DTH2(2),char(916),DT1(1),DT1(2),DT2(2),TW_DEP,TW_dep(n),TW_FB,TW_fb(m),TW_ARR,TW_arr(l))
fprintf('\n%cv_tot: %.3f km/s\n%cv%c: %.3f km/s\n%cv_powered_GA: %.3f km/s\n%cv_GA/%cv%c: %.3f\n%ct_soi: %.3f hours (%s)\nrp: %.2f km %s %.2f km (RVenus+hatm)\n', ...
    char(916),Dv_totmin,char(916),char(8734),Dv_inf,char(916),dvp,char(916),char(916),char(8734),dvp/Dv_inf,char(916),hours(seconds(dtSoi)),dtSOI,rp,str,R_V+h_atm_venus)
disp(STR)
if anynan(rp)
    DVTOT(D).num = D;
    DVTOT(D).Dvtot = NaN;
    DVTOT(D).dvfb_ratio = NaN;
    toc(tic1)
end
    
%% FMINCON
X0 = [TW_dep(n);TW_fb(m);TW_arr(l)]; % initial guess from time grid

% Dvtotmin = cost_tot(X0,h_atm_venus,T)
options = optimoptions("fmincon",'Display','off','Algorithm','active-set', ...                                                              
'ConstraintTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-6, ...
'MaxFunctionEvaluations',1e4,'MaxIterations',9e3);

A = []; b = []; Aeq = []; beq = []; % Linear equality & inequality constraints
% lower & upper bounds 
lb = [X0(1)-20 X0(2)-10 X0(3)-1300]; 
ub = [X0(1)+30 X0(2)+10 X0(3)+1500];
try
    try   % active-set is more accurate, but could be unstable (if error occurs, use default option for reliability)
        [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
    catch
        options.Algorithm = 'interior-point';
        [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
    end
catch
    DVTOT(D).num = D;
    DVTOT(D).Dvtot = NaN;
    DVTOT(D).dvfb_ratio = NaN;
    Xmin = NaN;
    fprintf('\nFmincon not executed\n')
    toc(tic1)
end

if anynan(Xmin)
    options.Algorithm = 'interior-point';
    [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
    if anynan(Xmin)
        DVTOT(D).num = D;
        DVTOT(D).Dvtot = NaN;
        DVTOT(D).dvfb_ratio = NaN;
        toc(tic1)
    end
end

TW_DEP = datetime(mjd20002date(Xmin(1)),'Format','yyyy MMM dd hh:mm:ss'); 
TW_FB = datetime(mjd20002date(Xmin(2)),'Format','yyyy MMM dd hh:mm:ss'); 
TW_ARR = datetime(mjd20002date(Xmin(3)),'Format','yyyy MMM dd hh:mm:ss'); 

%--------------------------------------------------------------------------
% Radius SOI

kep = uplanet(Xmin(2),2);

a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
rV = kep2car(a,e,i,Omega,omega,theta,mu_S);
rV = norm(rV);

rSoi = rV*(mV/mS)^(2/5);

%--------------------------------------------------------------------------
% D time SOI

[~,rp,a,Vinf,dvp_GA,C3] = cost_tot(Xmin,h_atm_venus,T);
a_minus = a(1);
a_plus = a(2);
Vinf_minus = Vinf(1,:);
Vinf_plus = Vinf(2,:);

[~,d,e_minus,e_plus,d_minus,d_plus] = fd(rp,Vinf_minus,Vinf_plus,mu_V);

E_minus = acosh((-rSoi/a_minus+1)/e_minus);
E_plus = acosh((-rSoi/a_plus+1)/e_plus);

dt_minus = (e_minus*sinh(E_minus)-E_minus)/sqrt(-mu_V/a_minus^3);
dt_plus = (e_plus*sinh(E_plus)-E_plus)/sqrt(-mu_V/a_plus^3);

dtSoi = dt_minus + dt_plus;
dtSOI = seconds(dtSoi);
dtSOI.Format = 'hh:mm:ss.SS';

Dvinf = norm(Vinf_plus-Vinf_minus);

%---------------------------------------------------------------------------
fprintf('\nDates from fmincon refinement:\nDeparture: %s - %.4f days (MJD2000) \nFly-by: %s - %.4f days\nArrival: %s - %.4f days\n', ...
    TW_DEP,Xmin(1),TW_FB,Xmin(2),TW_ARR,Xmin(3))
fprintf('\n%cv_tot: %.3f km/s\n%cv%c: %.3f km/s\n%cv_powered_GA: %.3e km/s\n%cv_GA/%cv%c: %.3e\n%ct_soi: %.3f hours (%s)\nrp: %.16f km > %.2f km (RVenus+hatm)\n\n', ...
    char(916),Dvtotmin,char(916),char(8734),Dvinf,char(916),dvp_GA,char(916),char(916),char(8734),dvp_GA/Dvinf,char(916),hours(seconds(dtSoi)),dtSOI,rp,R_V+h_atm_venus)

DVTOT(D).num = D;
DVTOT(D).Dvtot = Dvtotmin;
DVTOT(D).dvfb_ratio = dvp_GA/Dvinf;
%---------------------------------------------------------------------------

lb = [X0(1)-20 X0(2)-10 X0(3)-1300]; 
ub = [X0(1)+30 X0(2)+10 X0(3)+1500];

tw_dep = linspace(lb(1),ub(1),nn);
tw_arr = linspace(lb(3),ub(3),ll);

for n = 1:nn
    for l = 1:ll

        X = [tw_dep(n);Xmin(2);tw_arr(l)];
        [Dvtot(l,n),rp,~,Vinf] = cost_tot(X,h_atm_venus,T);
        if anynan(Vinf) %|| rp <= R_V + h_atm_venus
            Dvtot(l,n) = NaN;
            Dvflybytot(l,n) = NaN;
            continue
        end
        Dvflybytot(l,n) = norm(Vinf(1,:) - Vinf(2,:));
    end
end
%--------------------------------------------------------------------------
% plots

Dv_tot_fbdate = Dv_tot(:,:,M);
Dvflyby_tot_fbdate = Dvflyby_tot(:,:,M);

costplots(TW_dep,TW_arr,Dv_tot_fbdate,Dvflyby_tot_fbdate,X0,Dv_totmin,Dv_inf, ...
          Dvtot,Dvflybytot,Xmin,Dvtotmin,Dvinf,lb,ub,STEP,N)

% update C3 and Vinf plots
C3_dep = C3(1);
C3_arr = C3(2);

for s = 1:2
    figure(1), nexttile(s)
    plot3(Xmin(1),Xmin(2),C3_dep,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 2, legend('','ToF[days]','$\mathbf{\Delta v_{departure}}$ min','$\mathbf{\Delta v_{departure}}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(3), nexttile(s)
    plot3(Xmin(1),Xmin(2),Vinf_minus,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 2, legend('','ToF[days]','$\mathbf{v_{\infty}^-}$ min','$\mathbf{v_{\infty}^-}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(6), nexttile(s)
    plot3(Xmin(2),Xmin(3),Vinf_plus,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 2, legend('','ToF[days]','$\mathbf{v_{\infty}^+}$ min','$\mathbf{v_{\infty}^+}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(8), nexttile(s)
    plot3(Xmin(2),Xmin(3),C3_arr,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 2, legend('','ToF[days]','$\mathbf{\Delta v_{arrival}}$ min','$\mathbf{\Delta v_{arrival}}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(12), nexttile(s)
    if s == 1, legend('','ToF[days]','$\mathbf{\Delta v_{tot}}$ min from time grid','$\mathbf{\Delta v_{tot}}$ min from refinement','FontSize',15,'interpreter','latex','Location','southeast'), end
    if s == 2, legend('','ToF[days]','$\mathbf{\Delta v_{\infty}}$ mission from time grid','$\mathbf{\Delta v_{\infty}}$ mission from refinement','FontSize',15,'interpreter','latex','Location','southeast'), end
end

% debug plot 3d mesh
figure
mesh(tw_dep,tw_arr,Dvtot), hold on, grid on
plot3(Xmin(1),Xmin(3),Dvtotmin,'o','MarkerFaceColor','g','MarkerSize',10)
view(90,0)
xlim([Xmin(1)-10 Xmin(1)+10])
ylim([Xmin(3)-30 Xmin(3)+30])
zlim([Dvtotmin-1 Dvtotmin+4])
clim([Dvtotmin Dvtotmin+4])
xlabel('x')
ylabel('y')
zlabel('z')

close figure 5 figure 10 figure 13

dd = [DVTOT.Dvtot]';
ddd = [DVTOT.dvfb_ratio]';
[DVTOTMIN,INDEX] = min(dd);
fprintf('\n# %d\n%cv_tot: %.3f km/s\n%cv_GA/%cv%c: %.3e\n\n', ...
    D,char(916),DVTOTMIN,char(916),char(916),char(8734),ddd(INDEX))

toc

%% FUNCTIONS

function [t,x] = R2BP_propagation(x0,tspan,mu,J2)
% Restricted 2 body problem propagation
% INPUT:
%   x0[6x1] Initial State of the body ( rx, ry, rz, vx, vy, vz) [ L, L/T ]
%   tspan[N] Propagation time, boundaries[2] or vector[N] can be borth provided [s]
%   mu[1] Gravitational constant for celestial bodies [km^3/s^2]
%   J2[1]  boolean value, adds J2 perturbation
%
% OUTPUT
%   t[N] Propagation time vector [ s ]
%   x[N,6] Propagated state

if nargin < 4
    J2 = 0;
end
if nargin < 3
    mu = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
end

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[t,x] = ode113(@ode_2bp,tspan,x0,options,J2,mu);
end
%--------------------------------------------------------------------------
function y = ode_2bp(~,x,J2,mu)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y )
%
% INPUT:
%   t[1]   Time (can be omitted, as the system is autonomous) [T]
%   y[6x1] State of the body ( rx, ry, rz, vx, vy, vz) [ L, L/T ]
%   J2[1]  boolean value, adds J2 perturbation
%
% OUTPUT:
%   dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% 
% VERSIONS:
%
% mu_E[1] Gravitational parameter of the primary [L^3/T^2]

r = x(1:3);
v = x(4:6);
y = zeros(6,1);
y(1:3) = v;

rn = norm(r);
switch J2
    case 0
        aj2 = 0;
    case 1
        Re = astroConstants(23); % mean Earth radius
%       j2 = 0.00108263; % from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
        j2 = astroConstants(9);
        aj2 = 3/2*j2*mu*Re^2/(rn^4)*r/rn.*(5*(r(3)/rn)^2-[1;1;3]); % J2 perturbation
end
y(4:6) = -mu*r/rn^3+aj2; % dynamics
end

%--------------------------------------------------------------------------

function globe(name,R,r_input)

% plots globe (sphere) with specified texture
% INPUTS:
%       name[char/string]: 'name.extension' of texture file
%       R[1]: radius of globe [km]
%       r_input[3]: position of globe [km], optional

if nargin < 3
    r_input = zeros(3,1);
end
% R = astroConstants(23);
[X,Y,Z] = sphere(100);
I = imread(name); % I = flip(I,1);
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.CData = I;
surf(-R*X+r_input(1),-R*Y+r_input(2),-R*Z+r_input(3),props), hold on, grid on
axis equal padded, pbaspect([1,1,1])

end

%--------------------------------------------------------------------------

function [rvect,vvect,p,r,vpf] = kep2car(a,e,i,Omega,omega,theta,mu)
% Converts Keplerian parameters into Cartesian coordinates and velocities.
% 
% PROTOTYPE:
%   [rvect,vvect,p,r,vpf] = kep2car(a,e,i,Omega,omega,theta,mu)
% 
% INPUT:
%   a[1] = semi-major axis [km]
%   e[1] = eccentricity [-]
%   i[1] = inclination [rad]
%   Omega[1] = RAAN [rad]
%   omega[1] = argument of periapsis [rad]
%   theta[1] = true anomaly [rad]
%   mu[1] = gravitational parameter [km^3/s^2]
%   
% (alternatively, if nargin <3)
%   a[6] = Keplerian parameters set
%   e[1] = gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%   rvect[3] = position vector, can either be row or column [km]
%   vvect[3] = velocity vector, can either be row or column [km/s]
%   p[1] = semi-latus rectum [km]
%   r[1] = scalar distance from attractor [km]
%   vpf[3] = velocity vector (column) in perifocal reference frame [km/s]
% 
% FUNCTIONS CALLED:
%   (none)
%
% -------------------------------------------------------------------------

% If the user doesn't define mu in the inputs, default to Earth's value:

if nargin < 7 
    mu = 398600;
end

% If the user inputs Keplerian parameters as a vector, as opposed to single
% scalar values:

if nargin < 3
    if nargin ==1
        mu = 398600;
    else
        mu = e;
    end
    % a = a(1);
    e = a(2);
    i = a(3);
    Omega = a(4);
    omega = a(5);
    theta = a(6);
    a = a(1);
end

p = a*(1-e^2);
r = p/(1+e*cos(theta));
rpf = r*[cos(theta);sin(theta);0];
vpf = sqrt(mu/(p))*[-sin(theta);e+cos(theta);0];

ROmega = [cos(Omega) sin(Omega) 0
         -sin(Omega) cos(Omega) 0
         0 0 1];

Romega = [cos(omega) sin(omega) 0
         -sin(omega) cos(omega) 0
         0 0 1];

Ri = [1 0 0
      0 cos(i) sin(i)
      0 -sin(i) cos(i)];

T = Romega*Ri*ROmega;

rvect = T'*rpf;
vvect = T'*vpf;
end

%--------------------------------------------------------------------------

function [Dv1,Dv2,Dv1_,Dv2_,DV1,DV2] = porkchop(tw1,ToF,n,ID,arc,TOFMAX)

% INPUTs:
%       tw1[2 char]: departure time window, (Gregorian) [YYYY,MM,DD,hh,mm,ss]  
%       tw2[2 char]: arrival time window // ToF: Time of Flight[2] [years]
%       n[3]: number of discretization points for tw1 & tw2 (& ToF)
%       ID[2]: bodies identifier for uplanet / ephNEO
%       arc[1]: 1 or 2, first or second arc

% OUTPUT:
%       Dvtot[nxn]: grid of total cost of maneuver

Tw1 = datetime(tw1(1,:)):caldays:datetime(tw1(2,:));
tw2 = zeros(2,6);
for j = 1:2
    tw1(j,1) = date2mjd2000(tw1(j,:)); % [d]
    tw2(j,1) = tw1(j,1) + ToF(j) * 365.2417;
    tw2(j,:) = mjd20002date(tw2(j,1));
end

Tw2 = datetime(tw2(1,:)):calmonths:datetime(tw2(2,:));

for j = 1:2
    tw2(j,1) = date2mjd2000(tw2(j,:));
end
tw1 = tw1(:,1);
tw2 = tw2(:,1);
tw1 = linspace(tw1(1),tw1(2),n(1))';
tw2 = linspace(tw2(1),tw2(2),n(2))';
ToF = linspace(ToF(1),ToF(2),n(3))';

MU = astroConstants(4); % Sun's gravitational parameter

TOF = zeros(n(2),n(1)); Dvtot1 = TOF; Dvtot2 = TOF;
for j = 1:n(1)
    for k = 1:n(2)
        mod = 1;
            [Dv1(k,j),Dv2(k,j),TOF(k,j),~,DV1{k,j},DV2{k,j}] = cost([tw1(j),tw2(k)],ID,mod,TOFMAX);
        
    end
    for k = 1:n(3)
        mod = 2;
        [Dv1_(k,j),Dv2_(k,j)] = cost([tw1(j),ToF(k)],ID,mod,TOFMAX);
    end
end

ttw1 = linspace(tw1(1),tw1(end),length(Tw1)); 
ttw2 = linspace(tw2(1),tw2(end),length(Tw2));

Tw1 = datetime(Tw1,'Format','yyyy MMM dd');
Tw2 = datetime(Tw2,'Format','yyyy MMM dd');

for K = 1:2
    if K == 1
        Dv = Dv1;
        Dv_ = Dv1_;
    else
        Dv = Dv2;
        Dv_ = Dv2_;
    end
figure('Name',sprintf('PORKCHOP PLOT: TIME GRID'),'NumberTitle','off')
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
str = ["on","off"];
for j = 1:3
    if j < 3
        nexttile(j)
        [~,h] = contour(tw1,tw2,Dv,'ShowText',str(j),'LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
        ylabel('Arrival date','FontSize',15)
        title('Cost plot','FontSize',18)
        yticks(ttw2)
        yticklabels(char(Tw2)), ytickangle(45)
        if j == 2
            nexttile(2)
            I = find(TOF<0);
            TOF(I) = 0;
            contour(tw1,tw2,TOF/(3600*24),60:60:300,'ShowText','on','LineWidth',1,'EdgeColor','k')
            title('Cost plot','FontSize',18)
        end
    else
        figure('Name',sprintf('PORKCHOP PLOT: DDate & ToF'),'NumberTitle','off')
        [~,h] = contour(tw1,ToF,Dv_,'ShowText',str(1),'LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
        ylabel('ToF [years]','FontSize',15)
        title('Cost plot','FontSize',18)
    end

c = colorbar(gca,"eastoutside");
warning off MATLAB:handle_graphics:exceptions:SceneNode

if K == 1 
    if arc == 1
        c.Label.String = '$\mathbf{\Delta v_{departure}\ [km/s]}$';
    elseif arc == 2
        c.Label.String = '$\mathbf{v_{\infty}^{+}\ [km/s]}$';
    end
elseif K == 2
    if arc == 1
        c.Label.String = '$\mathbf{v_{\infty}^{-}\ [km/s]}$';
    elseif arc == 2
        c.Label.String = '$\mathbf{\Delta v_{arrival}\ [km/s]}$';
    end
end
c.Label.Interpreter = 'latex';
c.Label.FontSize = 15;
c.Label.Rotation = 0;
xlabel('Departure date','FontSize',15)
xticks(ttw1)
xticklabels(char(Tw1)), xtickangle(45)
end
end

figure('Name',sprintf('Debug mesh plot of \x394V')), 
mesh(tw1,tw2,Dv1,'DisplayName','\Deltav'), hold on, grid on, axis padded
xlabel('x'), ylabel('y'), zlabel('z')
view(0,0)

end

%--------------------------------------------------------------------------

function [Dv1,Dv2,TOF,cartesian,DV1,DV2] = cost(x,ID,mod,TOFMAX)
% Cost (objective) function for each arc, part of "time grid minimization", 
% and part of total cost function fro fmincon (cost_tot)
% Can be upgraded so to receive gregorian date as input, which can be 
% converted here to mjd2000, rather than outside
%   INPUTs: 
%       x[2]: mjd2000 "initial tr. time" & "final transfer time or ToF"
%       ID[2]
%       mod[1]: 1 for tw2 input, 2 for ToF input
%       TOFMAX[1]: days
%   OUTPUTs:
%       Dvtot[1]
%       TOF[1]
%       cartesian[3x6]: [RI,RF,vi,vf,VI',VF']

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

%--------------------------------------------------------------------------
function [Dvtot,rp,a,Vinf,dvp_GA,C3] = cost_tot(X,h_atm_venus,T)
% INPUT:
%       X[3] departure, fly-by, arrival mjd2000
%       h_atm_venus[1]: km
%       T[3]: days  period of revolution of Mercury, Venus, NEO81
% OUTPUT:
%       Dvtot[1] [km/s] total mission cost   
%       Vinf [2x3] [km/s]

ID = [1,2,81];
mod = 1;
[C3_dep,v_inf_minus,TOF,cartesian,DV1,V_inf_minus] = cost(X(1:2),ID(1:2),mod,T(2));
[v_inf_plus,C3_arr,TOF,cartesian,V_inf_plus,DV4] = cost(X(2:3),ID(2:3),mod,T(3));

if anynan(v_inf_minus) || anynan(v_inf_plus)
    Dvtot = NaN;
    Vinf = NaN;
    rp = NaN;
    a = NaN;
    dvp_GA = NaN;
    C3 = NaN;
    return
end

mu_V = astroConstants(12);
R_V = astroConstants(22);

a_minus = -mu_V/v_inf_minus^2;
a_plus = -mu_V/v_inf_plus^2;
rp0 = R_V+h_atm_venus; % initial guess close enough to solution (Venus radius)

options = optimset('TolX',1e-14,'Display','off');
rp = fzero(@fd,rp0,options,V_inf_minus,V_inf_plus,mu_V);

vp_minus = sqrt(v_inf_minus^2+2*mu_V/rp);
vp_plus = sqrt(v_inf_plus^2+2*mu_V/rp);
dvp_GA = abs(vp_minus - vp_plus);

Dvtot = C3_dep + C3_arr + dvp_GA;

a = [a_minus;a_plus];
Vinf = [V_inf_minus';V_inf_plus'];
C3 = [C3_dep;C3_arr];

end
%--------------------------------------------------------------------------

function [c,ceq] = nonlcon(X,h_atm_venus,T)
% non linear constraint function for fmincon
%   INPUTs: 
%       x[3]: mjd2000 departure, fly-by & and arrival transfer time mjd2000
%       h_atm_venus[1] km
%       T[3]: days  period of revolution of Mercury, Venus, NEO81
%   OUTPUTs:
%       c <= 0   constraint function
%       ceq == 0

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

%--------------------------------------------------------------------------

function [Fd,d,e_minus,e_plus,d_minus,d_plus] = fd(rp,vinf_minus,vinf_plus,mu)
% rp solver function for fzero or fsolve, for hyperbolic gravity assist
% INPUTs:
%         rp[1]: unknown pericenter radius [km]
%         vinf_minus,vinf_plus[3x1]: incoming & outcoming velocities 
%         mu[1]: planetary constant [km^3/s^2]
% OUTPUTs:
%         Fd: equation to be solved
%         d: total turning angle [deg] 
%         e_minus, e_plus: eccentricity of each arc of hyperbola
%         d_minus, d_plus: turning angle provided by each arc

e_minus = 1+rp*norm(vinf_minus)^2/mu;
d_minus = 2*asind(1/e_minus);

e_plus = 1+rp*norm(vinf_plus)^2/mu;
d_plus = 2*asind(1/e_plus);

d = acosd(dot(vinf_minus,vinf_plus)/(norm(vinf_minus)*norm(vinf_plus)));

Fd = d - (d_minus+d_plus)/2;

end

%--------------------------------------------------------------------------
function costplots(TW_dep,TW_arr,Dv_tot,Dvflyby_tot,X0,Dv_totmin,Dv_infmin, ...
    Dvtot,Dvflybytot,Xmin,Dvtotmin,Dvinfmin,lb,ub,STEP,N)



tw_dep(1,:) = mjd20002date(lb(1));
tw_dep(2,:) = mjd20002date(ub(1));
tw_arr(1,:) = mjd20002date(lb(3));
tw_arr(2,:) = mjd20002date(ub(3));

Tw1 = datetime(tw_dep(1,:)):10*caldays:datetime(tw_dep(2,:));
Tw2 = datetime(tw_arr(1,:)):calmonths:datetime(tw_arr(2,:));
Tw1 = datetime(Tw1,'Format','yyyy MMM dd');
Tw2 = datetime(Tw2,'Format','yyyy MMM dd');

ttw1 = linspace(lb(1),ub(1),length(Tw1)); 
ttw2 = linspace(lb(3),ub(3),length(Tw2));

for s = 1:2
switch s
    case 1
        X = X0;
        DV = Dv_tot;
        DVfb = Dvflyby_tot;
        DVtotmin = Dv_totmin;
        DVinfmin = Dv_infmin;
        TOF = zeros(length(TW_arr),length(TW_dep));
        for n = 1:length(TW_dep)
            for l = 1:length(TW_arr)
                TOF(l,n) = TW_arr(l) - TW_dep(n);
            end
        end 
    case 2
        X = Xmin;
        DV = Dvtot;
        DVfb = Dvflybytot;
        DVtotmin = Dvtotmin;
        DVinfmin = Dvinfmin;
        TW_dep = linspace(lb(1),ub(1),length(TW_dep));
        TW_arr = linspace(lb(3),ub(3),length(TW_arr));
        TOF = zeros(length(TW_arr),length(TW_dep));
        for n = 1:length(TW_dep)
            for l = 1:length(TW_arr)
                TOF(l,n) = TW_arr(l) - TW_dep(n);
            end
        end  
        I = find(TOF<0); J = find(TOF>1400);
        TOF(I) = 0; TOF(J) = NaN;
end

    figure
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
    nexttile
    [~,h] = contour(TW_dep,TW_arr,DV,'ShowText','off','LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
    contour(TW_dep,TW_arr,TOF,'LevelListMode','auto','ShowText','on','LineWidth',1,'EdgeColor','k')

    plot3(X0(1),X0(3),Dv_totmin,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k')
    if s == 2
        plot3(X(1),X(3),Dvtotmin,'o','MarkerFaceColor','g','MarkerSize',10,'MarkerEdgeColor','k');
    end
    xlabel('Departure date','FontSize',15)
    ylabel('Arrival date','FontSize',15)
    title('Cost plot','FontSize',18)
    h.LevelList = floor(DVtotmin)-1:STEP/2:floor(DVtotmin)+N;
    c = colorbar(gca,"eastoutside");
    c.Label.String = '$\mathbf{\Delta v_{tot}\ [km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 15;
    c.Label.Rotation = 0;
    c.Label.Position = [.564177274703979,floor(DVtotmin)+N+1,0];
    clim([floor(DVtotmin)-1 floor(DVtotmin)+N])
    xticks(ttw1)
    yticks(ttw2)
    xticklabels(char(Tw1)), xtickangle(45)
    yticklabels(char(Tw2)), ytickangle(45)
    ax = gca;
    ax.XTick = ax.XTick(1:ceil(length(ax.XTick)/10):end);
    ax.XTickLabel = ax.XTickLabel(1:ceil(length(ax.XTickLabel)/10):end);
    ax.YTick = ax.YTick(1:ceil(length(ax.YTick)/10):end);
    ax.YTickLabel = ax.YTickLabel(1:ceil(length(ax.YTickLabel)/10):end);
    xlim([lb(1) ub(1)])
    ylim([lb(3) ub(3)])

    nexttile
    [~,h] = contour(TW_dep,TW_arr,DVfb,'ShowText','off','LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
    contour(TW_dep,TW_arr,TOF,'LevelListMode','auto','ShowText','on','LineWidth',1,'EdgeColor','k')
    plot3(X0(1),X0(3),Dv_infmin,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k')
    if s == 2
        plot3(X(1),X(3),DVinfmin,'o','MarkerFaceColor','g','MarkerSize',10,'MarkerEdgeColor','k');
    end
    xlabel('Departure date','FontSize',15)
    ylabel('Arrival date','FontSize',15)
    title('Cost plot','FontSize',18)
    h.LevelList = floor(DVinfmin)-1:STEP/2:floor(DVinfmin)+N;
    c = colorbar(gca,"eastoutside");
    c.Label.String = '$\mathbf{\Delta v_{\infty}\ [km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 15;
    c.Label.Rotation = 0;
    c.Label.Position = [-.123322725296021,floor(DVinfmin)+N+1,0];
    clim([floor(DVinfmin)-1 floor(DVinfmin)+N])
    xticks(ttw1)
    yticks(ttw2)
    xticklabels(char(Tw1)), xtickangle(45)
    yticklabels(char(Tw2)), ytickangle(45)
    ax = gca;
    ax.XTick = ax.XTick(1:ceil(length(ax.XTick)/10):end);
    ax.XTickLabel = ax.XTickLabel(1:ceil(length(ax.XTickLabel)/10):end);
    ax.YTick = ax.YTick(1:ceil(length(ax.YTick)/10):end);
    ax.YTickLabel = ax.YTickLabel(1:ceil(length(ax.YTickLabel)/10):end);
    xlim([lb(1) ub(1)])
    ylim([lb(3) ub(3)])

end
end


