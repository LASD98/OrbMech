%% INTERPLANETARY MISSION
%--------------------------------------------------------------------------
% GROUP 2304
%
% Contributors:
% Aufiero Luca
% Caushi Andrea
% Luciardello Lecardi Matteo
%--------------------------------------------------------------------------

% FUNCTIONS CALLED:

% car2kep.m
% cost.m
% cost_tot.m
% costplots.m
% fd.m
% kep2car.m
% nonlcon.m
% ode_2bp.m
% planet3D.m
% porkchop.m
% R2BP_propagation.m
% rodrigues.m
% sec2time.m

addpath(genpath('MATLAB functions\')), addpath(genpath('functions\'))
clc, clear, close all, format long g

%% SETUP for plots and initial data

N = 30;             % N+1 units of Dvtot will be shown in contour plot (Default:  20)
STEP = 2;           % Steps between contour plot lines in Dvtot units  (Default:  4)

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
x01 = str2num(char(datetime([tw1(column);tw2(row)],'Format','yyyy,MM,dd,HH,mm,ss.SS')));
x02 = [str2num(char(datetime([tw1(column_)],'Format','yyyy,MM,dd,HH,mm,ss.SS')));ToF(row_),zeros(1,5)];

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
%             V_inf_dep = DV1;
            nn = n(1);
            mm = n(2);
            DT1 = dt;
            DTH1 = dth;
        case 2  % ( m x l )
            Vinf_plus = DV1;
%             V_inf_arr = DV2;
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


TW_DEP = datetime(mjd20002date(TW_dep(n)),'Format','yyyy MMM dd HH:mm:ss'); 
TW_FB = datetime(mjd20002date(TW_fb(m)),'Format','yyyy MMM dd HH:mm:ss'); 
TW_ARR = datetime(mjd20002date(TW_arr(l)),'Format','yyyy MMM dd HH:mm:ss'); 

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
fprintf('\n%cv_dep: %.3f km/s\n%cv_arr: %.3f km/s\n\n%cv_tot: %.3f km/s\n%cv%c: %.3f km/s\n%cv_powered_GA: %.3f km/s\n%cv_GA/%cv%c: %.3f\n%ct_soi: %.3f hours (%s)\nrp: %.2f km %s %.2f km (RVenus+hatm)\n', ...
    char(916),C3_dep(m,n),char(916),C3_arr(l,m),char(916),Dv_totmin,char(916),char(8734),Dv_inf,char(916),dvp,char(916),char(916),char(8734),dvp/Dv_inf,char(916),hours(seconds(dtSoi)),dtSOI,rp,str,R_V+h_atm_venus)
disp(STR)
if anynan(rp)
    DVTOT(D).num = D;
    DVTOT(D).Dvtot = NaN;
    DVTOT(D).dvfb_ratio = NaN;
%     toc(tic1)
end
    
%% FMINCON
X0 = [TW_dep(n);TW_fb(m);TW_arr(l)]; % initial guess from time grid

% Dvtotmin = cost_tot(X0,h_atm_venus,T)
options = optimoptions("fmincon",'Display','off','Algorithm','active-set', ...                                                              
'ConstraintTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-6, ...
'MaxFunctionEvaluations',1e4,'MaxIterations',9e3);

A = []; b = []; Aeq = []; beq = []; % Linear equality & inequality constraints
% lower & upper bounds 
lb = [X0(1)-30 X0(2)-10 X0(3)-30];  
ub = [X0(1)+30 X0(2)+10 X0(3)+30];

% lb = [];
% ub = [];

try
    try
        try   % active-set is more accurate, but could be unstable (if error occurs, use default option for reliability)
            [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
        catch
            disp('active-set failure')
            options.Algorithm = 'interior-point';
            [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
        end
    catch
        disp('interior-point failure')
        options.Algorithm = 'sqp';
        [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
    end
catch
    disp('sqp (sequential quadratic programming) failure')
    DVTOT(D).num = D;
    DVTOT(D).Dvtot = NaN;
    DVTOT(D).dvfb_ratio = NaN;
    Xmin = NaN;
    error('Fmincon not executed')
%     disp('Fmincon not executed')
%     toc(tic1)
end

if anynan(Xmin)
        options.Algorithm = 'interior-point';
        [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
        if anynan(Xmin)
            DVTOT(D).num = D;
            DVTOT(D).Dvtot = NaN;
            DVTOT(D).dvfb_ratio = NaN;
            %         toc(tic1)
        end
end

TW_DEP = datetime(mjd20002date(Xmin(1)),'Format','yyyy MMM dd HH:mm:ss'); 
TW_FB = datetime(mjd20002date(Xmin(2)),'Format','yyyy MMM dd HH:mm:ss'); 
TW_ARR = datetime(mjd20002date(Xmin(3)),'Format','yyyy MMM dd HH:mm:ss'); 

%--------------------------------------------------------------------------
% Radius SOI

kep = uplanet(Xmin(2),2);

a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
rV = kep2car(a,e,i,Omega,omega,theta,mu_S);
rV = norm(rV);

rSoi = rV*(mV/mS)^(2/5);

%--------------------------------------------------------------------------
% D time SOI

[~,rp,a,Vinf,dvp_GA,C3,DV,V,R] = cost_tot(Xmin,h_atm_venus,T);
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
% update C3 and Vinf plots
C3_dep = C3(1);
C3_arr = C3(2);

fprintf('\nDates from fmincon refinement:\nDeparture: %s - %.4f days (MJD2000) \nFly-by: %s - %.4f days\nArrival: %s - %.4f days\n', ...
    TW_DEP,Xmin(1),TW_FB,Xmin(2),TW_ARR,Xmin(3))
fprintf('\n%cv_dep: %.3f km/s\n%cv_arr: %.3f km/s\n\n%cv_tot: %.3f km/s\n%cv%c: %.3f km/s\n%cv_powered_GA: %.3e km/s\n%cv_GA/%cv%c: %.3e\n%ct_soi: %.3f hours (%s)\nrp: %.16f km > %.2f km (RVenus+hatm)\n\n', ...
    char(916),C3_dep,char(916),C3_arr,char(916),Dvtotmin,char(916),char(8734),Dvinf,char(916),dvp_GA,char(916),char(916),char(8734),dvp_GA/Dvinf,char(916),hours(seconds(dtSoi)),dtSOI,rp,R_V+h_atm_venus)
Vv = V(2,:) - DV(2,:);
dv = DV(3,:)- DV(2,:);
fprintf(['\n%cV = [%.4f,%.4f,%.4f] [km/s]\n[%.4f,%.4f,%.4f] [km/s]\n[%.4f,%.4f,%.4f] [km/s]\n[%.4f,%.4f,%.4f] [km/s]\n\n' ...
    'V = [%.4f,%.4f,%.4f] [km/s]\n[%.4f,%.4f,%.4f] [km/s]\n[%.4f,%.4f,%.4f] [km/s]\n[%.4f,%.4f,%.4f] [km/s]\n\n' ...
    'V_Venus = [%.4f,%.4f,%.4f] [km/s]\n%cv_flyby =  [%.4f,%.4f,%.4f] [km/s]\n'],char(916),DV',V',Vv,char(916),dv)
kep1 = car2kep(R(1,:),V(1,:),mu_S); kep2 = car2kep(R(3,:),V(3,:),mu_S);
kep1(3:6) = kep1(3:6)*180/pi; kep2(3:6) = kep2(3:6)*180/pi;
fprintf(['\na1 = %.4e km, e1 = %.4f, i1 = %.4f deg, Om1 = %.4f deg, om1 = %.4f deg, f1 = %.4f deg\n' ...
    'a2 = %.4e km, e2 = %.4f, i2 = %.4f deg, Om2 = %.4f deg, om2 = %.4f deg, f2 = %.4f deg\n' ...
    'r_soi = %.4f km\na_hyp_minus = %.4f km, e_hyp_minus = %.4f\na_hyp_plus = %.4f km, e_hyp_plus = %.4f\n'], ...
    kep1,kep2,rSoi,a_minus,e_minus,a_plus,e_plus)

DVTOT(D).num = D;
DVTOT(D).Dvtot = Dvtotmin;
DVTOT(D).dvfb_ratio = dvp_GA/Dvinf;
%---------------------------------------------------------------------------

lb = [X0(1)-50 X0(2)-30 X0(3)-1000]; 
ub = [X0(1)+30 X0(2)+30 X0(3)+200];

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

% % update C3 and Vinf plots
% C3_dep = C3(1);
% C3_arr = C3(2);

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
xlabel('$x [km]$','FontSize',15,'Interpreter','latex')
ylabel('$y [km]$','FontSize',15,'Interpreter','latex')
zlabel('$z [km]$','FontSize',15,'Interpreter','latex')

close figure 5 figure 10 figure 13

dd = [DVTOT.Dvtot]';
ddd = [DVTOT.dvfb_ratio]';
[DVTOTMIN,INDEX] = min(dd);
fprintf('\n# %d\n%cv_tot: %.3f km/s\n%cv_GA/%cv%c: %.3e\n\n', ...
    D,char(916),DVTOTMIN,char(916),char(916),char(8734),ddd(INDEX))

toc
%% Interplanetary trajectory characterization with plots

mu_Sun    = astroConstants(4);         % [km^3/s^2] Sun planetary constant
R_V = astroConstants(22);
mu_V = astroConstants(12);
MercuryID = 1;
VenusID   = 2;
NeoID     = 81;
km2AU     = 6.684587122e-9;
days2sec  = 24*60*60;
Np        = 50000;
h_atm     = 300;

% Mercury position at departure
[kep0_M] = uplanet(Xmin(1), MercuryID);
[r0_M, ~] = kep2car(kep0_M(1),kep0_M(2),kep0_M(3),kep0_M(4),kep0_M(5),kep0_M(6),mu_Sun);

% Mercury orbital period
T_M = 2*pi*sqrt(kep0_M(1)^3/mu_Sun);

% Set tspan of points
tspanM_orbit = linspace(0,T_M,Np);
[h, min, sec] = sec2time (tspanM_orbit(2)-tspanM_orbit(1));
fracDay = hms2fracday(h, min, sec);

% Mercury orbit propagation
rM_orbit = zeros(length(tspanM_orbit),3);
rM_orbit(1,1:3) = r0_M;
t = 0;
for i = 2 : length(tspanM_orbit)
    t = t + fracDay;
    kepM_orbit = uplanet(Xmin(1)+t, MercuryID);
    [rM,~] = kep2car (kepM_orbit(1),kepM_orbit(2),kepM_orbit(3),kepM_orbit(4),kepM_orbit(5),kepM_orbit(6), mu_Sun);
    rM_orbit(i,1:3) = rM';
end
rM_AU = rM_orbit*km2AU;

% Venus position at flyby
[kep0_V] = uplanet(Xmin(2), VenusID);
[r0_V, vV] = kep2car(kep0_V(1),kep0_V(2),kep0_V(3),kep0_V(4),kep0_V(5),kep0_V(6),mu_Sun);

% Venus orbital period
T_V = 2*pi*sqrt(kep0_V(1)^3/mu_Sun);

% Venus orbit propagation
tspanV_orbit = linspace(0,T_V,Np);
[h, min, sec] = sec2time (tspanV_orbit(2)-tspanV_orbit(1));
fracDay = hms2fracday(h, min, sec);
rV_orbit = zeros(length(tspanV_orbit),3);
rV_orbit(1,1:3) = r0_V;
t = 0; % initialize t
for i = 1 : length(tspanV_orbit)
    t = t + fracDay;
    kepV_orbit = uplanet(Xmin(2)+t, VenusID);
    [rJ,~] = kep2car (kepV_orbit(1),kepV_orbit(2),kepV_orbit(3),kepV_orbit(4),kepV_orbit(5),kepV_orbit(6), mu_Sun);
    rV_orbit(i,1:3) = rJ';
end
rV_AU = rV_orbit*km2AU;

% NEO position at arrival
[kep0_NEO, ~] = ephNEO(Xmin(3), NeoID);
[r0_NEO, ~] = kep2car(kep0_NEO(1),kep0_NEO(2),kep0_NEO(3),kep0_NEO(4),kep0_NEO(5),kep0_NEO(6),mu_Sun);
T_NEO = 2*pi*sqrt(kep0_NEO(1)^3/mu_Sun);

% NEO orbit propagation
tspanNEO_orbit = linspace(0,(21/20)*T_NEO,Np);
[h, min, sec] = sec2time(tspanNEO_orbit(2)-tspanNEO_orbit(1));
fracDay = hms2fracday(h, min, sec);

rNEO_orbit = zeros(length(tspanNEO_orbit),3);
rNEO_orbit(1,1:3) = r0_NEO;
t = 0; % initialize t
for i = 1 : length(tspanNEO_orbit)
    t = t + fracDay;
    kepNEO_orbit = ephNEO(Xmin(3)+t, NeoID);
    [rNEO,~] = kep2car (kepNEO_orbit(1),kepNEO_orbit(2),kepNEO_orbit(3),kepNEO_orbit(4),kepNEO_orbit(5),kepNEO_orbit(6), mu_Sun);
    rNEO_orbit(i,1:3) = rNEO';
end
rNEO_AU = rNEO_orbit*km2AU;

% FIRST ARC
% Time of flight of the first transfer arc and tspan
TOF_arc1 = (Xmin(2) - Xmin(1))*days2sec;
tspan_arc1 = linspace(Xmin(1)*days2sec,Xmin(2)*days2sec,Np);

% Obtain initial and final velocity of the arc from Lambert
[~,~,~,~,vi_arc1,vf_arc1,~,~] = lambertMR(r0_M,r0_V,TOF_arc1,mu_Sun,0,0,0,0);

% Propagate state in the given tspan
[~,y_arc1]= R2BP_propagation([r0_M,vi_arc1'],tspan_arc1,mu_Sun);
y_arc1 = y_arc1*km2AU;
[a1,e1,i1,OM1,om1,th1] = car2kep (r0_M,vi_arc1',mu_Sun);

% SECOND ARC
% Time of flight of the second transfer arc and tspan
TOF_arc2 = (Xmin(3) - Xmin(2))*days2sec;
% Obtain initial and final velocity of the arc from Lambert
[~,~,~,~,vi_arc2,vf_arc2,~,~] = lambertMR(r0_V,r0_NEO,TOF_arc2,mu_Sun,0,0,0,0);
tspan_arc2 = linspace(Xmin(2)*days2sec,Xmin(3)*days2sec,Np);

% Propagate state in the given tspan
[~,y_arc2]= R2BP_propagation([r0_V,vi_arc2'],tspan_arc2,mu_Sun);
y_arc2 = y_arc2*km2AU;
[a2,e2,i2,OM2,om2,th2] = car2kep(r0_V,vi_arc2',mu_Sun);


%-------SUN OPTIONS------------
optSun.Position = [0 0 0];
optSun.Units = '20AU';
%-------MERCURY OPTIONS----------
optMercury.Position = [rM_AU(1,1) rM_AU(1,2) rM_AU(1,3)];
optMercury.Units = '1kAU';
%-------VENUS OPTIONS----------
optVenus.Position = [rV_AU(1,1) rV_AU(1,2) rV_AU(1,3)];
optVenus.Units = '1kAU';

figure(5)
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
az = [0, 0, -90];
el = [90, 0, 0];
for j = 1:3
if j == 1
    J = [2,1];
elseif j == 2 
    J = 2;
else
    J = 4;
end
nexttile(J)
title("Mercury-Venus-NEO81 Transfer",'FontSize',18)
xlabel("$x [AU]$",'FontSize',18,'interpreter','latex'); ylabel("$y [AU]$",'FontSize',18,'interpreter','latex'); zlabel("$z [AU]$",'FontSize',18,'interpreter','latex')
hold on; grid on; box on;
axis equal;
planet3D('Sun',optSun);
plot3( rM_AU(:,1), rM_AU(:,2), rM_AU(:,3), '--','Color', '#4DBEEE',LineWidth=1);
planet3D('Mercury',optMercury);
plot3( rV_AU(:,1), rV_AU(:,2), rV_AU(:,3), '--','Color', '#c16d43',LineWidth=1);
planet3D('Venus',optVenus);
plot3( rNEO_AU(:,1), rNEO_AU(:,2), rNEO_AU(:,3), '--','Color', 'k',LineWidth=1);
plot3( rNEO_AU(1,1), rNEO_AU(1,2), rNEO_AU(1,3),'*','Color', 'k', LineWidth=5);

plot3( y_arc1(:,1), y_arc1(:,2), y_arc1(:,3), 'Color', '#00891d',LineWidth=2);
plot3( y_arc2(:,1), y_arc2(:,2), y_arc2(:,3), 'Color', '#ff0084',LineWidth=2);

legend('','Mercury''s orbit','', 'Venus''s orbit','','NEO''s orbit','',...
    'First Lambert Arc','Second Lambert Arc','FontSize',18,'Location','southoutside');

axis padded
if j == 1
%     pbaspect([1,1,1])
    axis tight
end
view(az(j),el(j))
end

%--------------------------Flyby plots-------------------------------------
% Relative velocities wrt Venus before and after flyby/gravity assist
v_inf_minus = vf_arc1' - vV;                            % [km/s]
v_inf_plus = vi_arc2' - vV;                             % [km/s]  

% Turning angle
delta = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));
delta_deg = rad2deg(delta);

% Radius of pericentre
fun = @(x) (asin(1./(1+((x.*dot(v_inf_minus,v_inf_minus))./mu_V)))+...
            asin(1./(1+((x.*dot(v_inf_plus,v_inf_plus))./mu_V)))  - ...
            delta);
opt = optimset('Display','off', 'TolFun',1e-12);
rp_norm = fsolve(fun, R_V+h_atm, opt);

if rp_norm < (R_V+h_atm) || rp_norm > R_V*101.7
    disp('Pericentre radius is NOT physically feasible');
end

% Altitude at which the impulsive manoeuvre is performed
h_ga = rp_norm-R_V;                                   % [km]

vp_minus_norm = sqrt( (norm(v_inf_minus))^2 + (2*mu_V)/rp_norm);
vp_plus_norm = sqrt( (norm(v_inf_plus))^2 + (2*mu_V)/rp_norm);

% Magnitude of the manoeuvre at pericentre
delta_vp = norm(vp_minus_norm-vp_plus_norm);          % [km/s]

% Magnitude of the flyby manoeuvre
delta_v_flyby_norm = norm(v_inf_plus-v_inf_minus);    % [km/s]
delta_v_flyby = v_inf_plus - v_inf_minus;             % [km/s]

% Hyperbola before gravity assist
e_hyp_minus = 1+((rp_norm*dot(v_inf_minus,v_inf_minus))/mu_V);
a_hyp_minus = rp_norm/(1-e_hyp_minus);

% Hyperbola after gravity assist
e_hyp_plus = 1+((rp_norm*dot(v_inf_plus,v_inf_plus))/mu_V);
a_hyp_plus = rp_norm/(1-e_hyp_plus);

figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

for j = 1:2
% reference_frame = 'Inertial Venus-centered frame';
% reference_frame = 'Perifocal frame';
    if j == 1
        reference_frame = 'Perifocal frame';
    else
        reference_frame = 'Inertial Venus-centered frame';
    end
    switch reference_frame
        case 'Perifocal frame'
            apse_line = [1; 0; 0];
            h_dir = [0; 0; 1];
        case 'Inertial Venus-centered frame'
            h_dir = cross(v_inf_minus,v_inf_plus)/norm(cross(v_inf_minus,v_inf_plus));
            half_delta_plus = asin(1/e_hyp_plus);
            x_vec = rodrigues(v_inf_plus,h_dir,pi/2-half_delta_plus);
            apse_line = - x_vec/norm(x_vec);
    end
    
    rp = rp_norm * apse_line;                             % [km]
    
    vp_minus = vp_minus_norm * (cross(h_dir,apse_line));  % [km/s]
    vp_plus = vp_plus_norm * (cross(h_dir,apse_line));    % [km/s]
    
    
    % Define radius of Venus's SOI
    G = astroConstants(1);
    mass_S = mu_Sun/G;
    mass_V = mu_V/G;
    rSOI_V = (1/km2AU)*5.2*((mass_V/mass_S)^(2/5));                %[km]
    
    % Flyby time in Venus's SOI - first arc
    h_hyp_minus = cross(rp,vp_minus);
    Th1 = acos(norm(h_hyp_minus)^2/(rSOI_V*mu_V*e_hyp_minus)-(1/e_hyp_minus));
    F_minus = 2*(atanh(tan(Th1/2)/sqrt((e_hyp_minus+1)/(e_hyp_minus-1))));
    deltaT_flyby_1 = (e_hyp_minus*sinh(F_minus)-F_minus)*sqrt((-a_hyp_minus)^3/mu_V);
    
    % Flyby time in Venus's SOI - second arc
    h_hyp_plus = cross(rp,vp_plus);
    Th2 = acos(norm(h_hyp_plus)^2/(rSOI_V*mu_V*e_hyp_plus)-1/e_hyp_plus);
    F_plus = 2*(atanh(tan(Th2/2)/sqrt((e_hyp_plus+1)/(e_hyp_plus-1))));
    deltaT_flyby_2 = (e_hyp_plus*sinh(F_plus)-F_plus)*sqrt((-a_hyp_plus)^3/mu_V);
    
    % Time that takes the spacecraft to go through the SOI of Venus
    flyby_time = deltaT_flyby_1+deltaT_flyby_2;
    
    
    %--------------------------Hyperbola propagation--------------------------
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    t_hyp_minus = linspace(0, -deltaT_flyby_1, Np);
    t_hyp_plus = linspace(0, deltaT_flyby_2, Np);
    
    [~,y_minus_SOI] = ode113(@(t,s) ode_2bp(t,s,0,mu_V),t_hyp_minus,[rp;vp_minus],options);
    [~,y_plus_SOI] = ode113(@(t,s) ode_2bp(t,s,0,mu_V),t_hyp_plus,[rp;vp_plus],options);
    y_minus_SOI_RV = y_minus_SOI/R_V;
    y_plus_SOI_RV = y_plus_SOI/R_V;
    
    %-------VENUS OPTIONS----------
    optVenus.Position = [0 0 0];
    optVenus.Units = 'RV';
    
    %-------SUN OPTIONS------------
    optSun.Position = [-r0_V(1)/R_V -r0_V(2)/R_V -r0_V(3)/R_V];
    optSun.Units = '250AU';
    
    % optVenus.Position = [0 0 0];
    % optVenus.Units = 'RV';
    
    % SOI components for plot
    [X_u,Y_u,Z_u] = sphere(50);
    X_SOI = X_u.*(rSOI_V/R_V);
    Y_SOI = Y_u.*(rSOI_V/R_V);
    Z_SOI = Z_u.*(rSOI_V/R_V);
    
    
    % PLOT WITH SPHERE OF INFLUENCE
    % figure()
    % hold on;
    % planet3D('Venus',optVenus);
    % plot3(y_minus_SOI_RV(:,1), y_minus_SOI_RV(:,2), y_minus_SOI_RV(:,3),'Color','#00891d',LineWidth=2);
    % plot3(y_plus_SOI_RV(:,1), y_plus_SOI_RV(:,2), y_plus_SOI_RV(:,3),'Color','#ff0084',LineWidth=2);
    % lightGrey = 0.9*[1 1 1];
    % surface(X_SOI,Y_SOI,Z_SOI,'FaceColor','none','EdgeColor',lightGrey)
    % % planet3D('Sun',optSun);
    % % quiver3(0,0,0,vJ(1)*10,vJ(2)*10,vJ(3)*10,'-k','LineWidth',1)
    %
    % title("Venus's Sphere of Influence",'FontSize',18)
    % legend('','Hyperbolic orbit before flyby','Hyperbolic orbit after flyby','Venus velocity',...
    %        'Location','northeast','FontSize',18)
    % xlabel("x $[R_V]$",'interpreter','Latex','FontSize',18); ylabel("y $[R_V]$",'interpreter','Latex','FontSize',18); zlabel("z $[R_V]$",'interpreter','Latex','FontSize',18)
    % axis equal
    % grid on
    % view(0, 90)
    
    %--------------------------Plots with asymptotes--------------------------
    % Asymptotes direction
    t_hyp_minus = linspace(0, -deltaT_flyby_1/150, Np);
    t_hyp_plus = linspace(0, deltaT_flyby_2/150, Np);
    
    [~,y_minus_SOI] = ode113(@(t,s) ode_2bp(t,s,0,mu_V),t_hyp_minus,[rp;vp_minus],options);
    [~,y_plus_SOI] = ode113(@(t,s) ode_2bp(t,s,0,mu_V),t_hyp_plus,[rp;vp_plus],options);
    y_minus_SOI_RV = y_minus_SOI/R_V;
    y_plus_SOI_RV = y_plus_SOI/R_V;
    
    c_minus_vec = (rp_norm+norm(a_hyp_minus)).*apse_line;
    c_minus_vec_RV = c_minus_vec/R_V;
    c_plus_vec = (rp_norm+norm(a_hyp_plus))*apse_line;
    c_plus_vec_RV = c_plus_vec/R_V;
    rp_RV = rp/R_V;
    
    switch reference_frame
        case 'Perifocal frame'
            x_apse = [-2 c_minus_vec_RV(1)];
            y_apse = [0 0];
            z_apse = [0 0];
        case 'Inertial Venus-centered frame'
            y_apse = [c_minus_vec_RV(2) 1.92];
    
            apse_line_xy = apse_line(1:2)/norm(apse_line(1:2));
            apse_line_yz = apse_line(2:3)/norm(apse_line(2:3));
            angleA = acosd(dot(-apse_line_xy,[0;1]));
            angleB = acosd(dot(-apse_line_yz,[1;0]));
            x_apse = c_minus_vec_RV(1) - tand(angleA)*(y_apse-c_minus_vec_RV(2));
            z_apse = c_minus_vec_RV(3) - tand(angleB)*(y_apse-c_minus_vec_RV(2)) ;
    end
    
    optVenus.Position = [0 0 0];
    optVenus.Units = 'RV';
    
    % figure()
    nexttile
    hold on;
    planet3D('Venus',optVenus);
    plot3(y_minus_SOI_RV(:,1), y_minus_SOI_RV(:,2), y_minus_SOI_RV(:,3), ...
        'Color','#0072bd','LineWidth',3);
    plot3(y_plus_SOI_RV(:,1), y_plus_SOI_RV(:,2), y_plus_SOI_RV(:,3), ...
        'Color','#ff6929','LineWidth',3);
    
    % vJ_norm = norm(vJ) * (cross(h_dir,apse_line));
    % quiver3(0,0,0,vJ_norm(1)*10,vJ_norm(2)*10,vJ_norm(3)*10,'-k','LineWidth',1)
    
    line([c_minus_vec_RV(1) y_minus_SOI_RV(end,1)],...
        [c_minus_vec_RV(2) y_minus_SOI_RV(end,2)],...
        [c_minus_vec_RV(3) y_minus_SOI_RV(end,3)], ...
        'color','#0072bd','linestyle','-.','LineWidth',3);
    hold on;
    line([c_plus_vec_RV(1) y_plus_SOI_RV(end,1)],...
        [c_plus_vec_RV(2) y_plus_SOI_RV(end,2)],...
        [c_plus_vec_RV(3) y_plus_SOI_RV(end,3)], ...
        'color','#ff6929','linestyle','-.','LineWidth',3);
    
    % line([-rp_RV(1)*5e+1 rp_RV(1)*0.5e+1],[-rp_RV(2)*5e+1 rp_RV(2)*0.5e+1],...
    %      [-rp_RV(3)*5e+1 rp_RV(3)*0.5e+1],'color','k','linestyle','--',LineWidth=1);
    
    % draw apse line
    plot3(x_apse,y_apse,z_apse,'color','k','linestyle','-.','LineWidth',3)
    
    axis equal
    switch reference_frame
        case 'Perifocal frame'
            view(0,90)
            xlim([x_apse(1) 2])
            title("Flyby in Perifocal frame",'FontSize',18)
        case 'Inertial Venus-centered frame'
            view(-69.7882418771531,17.6496756347232)
            title("Flyby in Inertial Venus-centered frame",'FontSize',18)
    end
    
    legend('','Incoming hyperbola','Outcoming hyperbola','Asymptote for incoming hyperbola','Asymptote for outcoming hyperbola','Apse line','Location','northeast','FontSize',12)
    xlabel("$x [R_V]$",'interpreter','Latex','FontSize',18)
    ylabel("$y [R_V]$",'interpreter','Latex','FontSize',18)
    zlabel("$z [R_V]$",'interpreter','Latex','FontSize',18)
    grid on

end

%-------------------------- Velocity Triangles-----------------------------

delta_v_tot_flyby = vi_arc2-vf_arc1;

figure()
hold on
grid on
axis equal
quiver3(0,0,0,vV(1),vV(2),vV(3),'r','LineWidth',2,...
       'MaxHeadSize',0.1,'AutoScaleFactor',1)
quiver3(vV(1),vV(2),vV(3),v_inf_minus(1),v_inf_minus(2),v_inf_minus(3),'g','LineWidth',2,...
       'MaxHeadSize',0.1,'AutoScaleFactor',1)
quiver3(0,0,0,vf_arc1(1),vf_arc1(2),vf_arc1(3),'b','LineWidth',2,...
       'MaxHeadSize',0.1,'AutoScaleFactor',1)
quiver3(vV(1),vV(2),vV(3),v_inf_plus(1),v_inf_plus(2),v_inf_plus(3),'m','LineWidth',2,...
       'MaxHeadSize',0.1,'AutoScaleFactor',1)
quiver3(0,0,0,vi_arc2(1),vi_arc2(2),vi_arc2(3),'c','LineWidth',2,...
       'MaxHeadSize',0.1,'AutoScaleFactor',1)
quiver3(vf_arc1(1),vf_arc1(2),vf_arc1(3),delta_v_tot_flyby(1),delta_v_tot_flyby(2),delta_v_tot_flyby(3),...
       'k','LineWidth',2,'MaxHeadSize',0.1,'AutoScaleFactor',1)
title('Velocity Triangles','FontSize',18)
xlabel('$V_x [km/s]$','interpreter','Latex','FontSize',18)
ylabel('$V_y [km/s]$','interpreter','Latex','FontSize',18)
zlabel('$V_z [km/s]$','interpreter','Latex','FontSize',18)
legend('$V_{\mathrm{Venus}}$','$v_\infty^-$','$V^-$','$v_\infty^+$','$V^+$','${\Delta}v$','Location','northeast','interpreter','Latex','FontSize',18)
