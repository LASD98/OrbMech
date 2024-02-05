%% INTERPLANETARY MISSION
% GROUP 2304
cd 'D:\Courses\4\1.Orbital Mechanics\2023_24\codes'
addpath(genpath('MATLAB functions\')), addpath('textures\')

%--------------------------------------------------------------------------
clc, clear, close all
tic0 = tic;

% SETUP

% n = 3e2; % discretization of time grid
N = 30;  % N+1 units of Dvtot will be shown in Porkchop (Default:  20)
STEP = 2; % Steps between Porkchop lines in Dvtot units  (Default:   4)
DD = 48; % Tsyn discretization for time grid

% tw1 = [2028,1,1,0,0,0;2057,1,1,0,0,0]; tw2 = [2029,4,1,0,0,0;2058,1,1,0,0,0]; % [date]

% tw = [2028,1,1,0,0,0;2058,1,1,0,0,0];
tmax = date2mjd2000([2058,1,1,0,0,0]);
mu_S = astroConstants(4); % [km^3/s^2]
ID = [1,2,81];
% d_trueanomaly = 2; % [deg]   % use for eccentric anomaly or true anomaly for circular orbit

for s = 1:3
%     if s == 1
%         tw = tw1;
%     elseif s == 2
%         tw = tw2;
%     elseif s == 3
%         tw = tw3;
%     end
    if ID(s) < 12
        kep(s,:) = uplanet(0,ID(s));
    else
        kep(s,:) = ephNEO(0,ID(s));
    end

    a = kep(s,1);
    e = kep(s,2);
    T(s) = 2*pi*sqrt(a^3/mu_S)/(3600*24); % [days]

%     E = 2*atand(sqrt((1-e)/(1+e))*tand(d_trueanomaly/2)); % [deg]
%     dt_step(s) = (E * pi/180 - e*sind(E))/sqrt(mu_S/a^3)/(3600*24); % [days]

%     dt_step(s) = d_trueanomaly * T(s)/360; % [days]  circular orbit 
%     ti = date2mjd2000(tw(1,:)); 
%     tf = date2mjd2000(tw(2,:));
%     n(s) = (tf(1)-ti(1))/dt_step(s);
end     

%--------------------
% Labert arcs for minimum dt (parabolic)
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

% FLY-BY preliminary --> real rSOI of Venus is function of rVenus at flyby
% Rvenus = astroConstants(22); % [km]
% ratio_Rsoi_Rvenus = 101.7; 
% Rsoi = ratio_Rsoi_Rvenus*Rvenus; % [km]

% rp = 0.1635758079559320;
% ra = 4.725549239901490;
% a = (ra+rp)/2
% e = (ra-rp)/(ra+rp)

%----------------------------------------------------------------------------
% SYNODIC PERIOD 3 planets
% by Borislav Borisov - 22.09.2011
% "Three-planet resonances in the Solar system"

a1 = kep(1,1); a2 = kep(2,1); a3 = kep(3,1);
Tsyn12 = 1/abs(1/T(2)-1/T(1)); % first arc
Tsyn23 = 1/abs(1/T(3)-1/T(2)); % second arc
Tsyn13 = 1/abs(1/T(3)-1/T(1)); % second arc

%
Tsynmax = max([Tsyn12,Tsyn23,Tsyn13]); % [days]
% P1 < P2 < P3
P1 = T(1); kepjup = uplanet(0,5);% a3 = kepjup(1);
P2 = T(2); P3 = 2*pi*sqrt(a3^3/mu_S)/(3600*24); P13 = 89.7924087;
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
%         j = i*round(ratio);
    elseif ratio <= 1
        j = j+1;
        i = round(j/ratio);

    end
    Pij = (i+j)*P13/365.2417;
%     Pij = 360*(i+j)/(n1-n3);
    P1_ = i/((i+j)/P2-j/P3);
%     syms P1_ real
%     P1_ = solve(i/P1_-(i+j)/P2+j/P3,P1_);
    dP1 = P1_-P1;
    accuracy = abs(dP1)/P1;

    fprintf('\b\b\b\b\b\b\b\b\b')
    fprintf('%.6f\n',accuracy)
%     pause(.5)
    if accuracy - err <= tol
        if ratio > 1
            i_ = i+1;
            j_ = round(i_*ratio);
        elseif ratio <= 1
            j_ = j+1;
            i_ = round(j_/ratio);
        end
        Pij_ = (i+j)*P13/365.2417;
%         Pij = 360*(i+j)/(n1-n3);
        P1_ = i_/((i_+j_)/P2-j_/P3);
%         syms P1_ real
%         P1_ = solve(i/P1_-(i+j)/P2+j/P3,P1_);
        dP1_ = P1_-P1;
        accuracy_ = abs(dP1_)/P1;
    end
end
Tsyn123 = Pij*365.2417;
fprintf('\b < %.2f\ni:j = %d:%d\nTsyn123 = %.4fy = %.2fd\n',err,i,j,Pij,Tsyn123)

%-------------------------------------------------------------------------------

% Hohmann transfer orbit (hto) Period & dv
for hto = 1:2
    a1 = kep(hto,1);
    e1 = kep(hto,2);
    i1 = kep(hto,3);
    a2 = kep(hto+1,1);
    e2 = kep(hto+1,2);
    rp = a1*(1-e1);
    ra = a2*(1+e2);
    aH = (rp+ra)/2;
    
    TH(hto) = 2*pi*sqrt(aH^3/mu_S) /(3600*24*365.2417); % [years]
    
    ra1 = a1*(1+e1);
    rp2 = a2*(1-e2);
    h1 = sqrt(2*mu_S*rp*ra1/(rp+ra1));
    h2 = sqrt(2*mu_S*rp*ra/(rp+ra));
    h3 = sqrt(2*mu_S*rp2*ra/(rp2+ra));
    v1 = h1/rp; v2 = h2/rp; v3 = h2/ra; v4 = h3/ra;
    dv1(hto) = abs(v2-v1);
    dv2(hto) = abs(v4-v3);
    
end
DvtotH = dv1(1)+dv2(2);


dth = 2; % deg  True Anomaly average step, doesn't consider eccentricity 
%               of orbital path, ideal for nearly circular orbits


DVTOT = [];
% t2_ = [2028,1,1,0,0,0];

 for D = 1:DD
% D = 8;
    close all
    clearvars -except D DD DVTOT t2_ N STEP TH T tmax Tsyn123 dtp1 dtp2 mu_S tic0
    tic1 = tic;
    fprintf('\n# %d\n',D)
for arc = 1:2
    switch arc 
        case 1
% MERCURY - VENUS
            
            dth = [10,10]; % deg true anomaly step
%             dt = [10,10,20]; % days   time step for tw1 & tw2 & ToF
%             dt = dth/360*[T(arc),T(arc+1),20*360/dth]; % days: time step for tw1 & tw2 & ToF
%                                                      corresponding to "dth" average true anomaly
                                                    
            ID = [1,2];
            t1 = [2028,1,1,0,0,0]; 
%             t1 = t2_;

%             t1 = mjd20002date(date2mjd2000(t1)+Tsyn123/48);
%             t1 = mjd20002date(date2mjd2000(t1)+Tsyn123/3);
%             t1 = mjd20002date(date2mjd2000(t1)+Tsyn123);
%             t1 = mjd20002date(date2mjd2000(t1)+Tsyn123/DD);
            t1 = mjd20002date(date2mjd2000(t1)+(D-1)*Tsyn123/DD);

%             t2 = [2058,1,1,0,0,0];
%             t2 = mjd20002date(date2mjd2000(t1)+Tsyn123);
%             t2 = mjd20002date(date2mjd2000(t1)+Tsyn123/48);
            t2 = mjd20002date(date2mjd2000(t1)+Tsyn123/DD);

%             t2 = [2038,1,1,0,0,0];
%             t1 = [2028,1,1,0,0,0];
%             t2 = [2029,4,1,0,0,0];

            t2_ = t2;

            tw1 = [t1;t2];
            tw1_ = tw1;
%             tw2 = [2028,1,1,0,0,0;2058,1,1,0,0,0];
%             ToF = [0;30];
%             ToF = [.5*TH(1)/2;1.5*TH(1)/2]; % [years]
%             ToF = [TH(1)/2-Tsynmax/2;TH(1)/2+Tsynmax/2];
            ToF = [1.2*dtp1/365.2417;T(2)/365.2417];
%             ToF = [1.1*dtp1/365.2417;TH(1)]; % TH is 2*(half Hohmann transfer orbit period) -> maximum for "zero-revolution" transfer
        case 2
% VENUS - NEO n.81 2005EL70
            dth = [dth(2),10]; % deg true anomaly step
%             dt = [dt(2),10,20]; % days
%             dt = dth/360*[T(arc),T(arc+1),20*360/dth]; % time step for tw1 & tw2 & ToF

            ID = [2,81];
%             tw1 = [2028,1,1,0,0,0;2058,1,1,0,0,0];
            t1 = mjd20002date(date2mjd2000(t1)+ToFG(1)*365.2417);
            t2 = mjd20002date(date2mjd2000(t2)+ToFG(2)*365.2417);
            tw1 = [t1; t2];
%             tw2 = [2028,1,1,0,0,0;2058,1,1,0,0,0];
%             ToF = [0;30]; 
%             ToF = [.7*TH(2);1.3*TH(2)]; % [years]
            ToF = [1.2*dtp2/365.2417;T(3)/365.2417];
%             ToF = [1.1*dtp2/365.2417;TH(2)];
    end

 dt = 1/360*[dth(1)*T(arc),dth(2)*T(arc+1),20*360]; % days: time step for tw1 & tw2 & ToF
%                                             corresponding to "dth" average true anomaly


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
%     dt_step(s) = d_trueanomaly * T(arc)/360; % [days]  circular orbit 
  
    n(s) = ceil((tf-ti)/dt(s));
end

% PORKCHOP PLOT

% figure('Name',sprintf('PORKCHOP PLOT: TIME GRID'),'NumberTitle','off')
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
% tw2 = zeros(2,6);
for j = 1:2
%     tw1(j,1) = date2mjd2000(tw1(j,:)); % [d]
%     tw2(j,1) = tw1(j,1) + ToF(j) * 365.2417;
    t2max = tw2(j,1); % [day]
%     tw2(j,:) = mjd20002date(tw2(j,1));
end
tw2G = tw2;
ToFG = ToF;

% figure(3*arc-2), 
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

    clim([floor(Dvmin)-1 floor(Dvmin)+N]), %h = findobj(gcf);
    h1.LevelList = floor(Dvmin)-1:STEP:floor(Dvmin)+N;
    c = h(i(j)); c.Label.Position = [.564177274703979,floor(Dvmin)+N+.9,0];
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


% tw1G = tw1; % Gregorian
% tw2 = zeros(2,6);
% for j = 1:2
%     tw1(j,1) = date2mjd2000(tw1(j,:)); % [d]
%     tw2(j,1) = tw1(j,1) + ToF(j) * 365.2417;
%     t2max = tw2(j,1); % [day]
%     tw2(j,:) = mjd20002date(tw2(j,1));
% end
% tw2G = tw2;
% ToFG = ToF;

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
% legend(p2,'Location','southeast','FontSize',15)

% x0(1) = date2mjd2000(x02(1,:)); x0(2) = date2mjd2000(x02(2,:));

% FMINCON setup

% options = optimoptions("fmincon",'Algorithm','active-set','Display','off', ...
% 'ConstraintTolerance',1e-10,'StepTolerance',1e-16,'FunctionTolerance',1e-16,'OptimalityTolerance',1e-6, ...
% 'MaxFunctionEvaluations',10e3,'MaxIterations',9e3);
% 
% A = []; b = []; Aeq = []; beq = []; % Linear equality & inequality constraints
% lower & upper bounds 
% lb = [date2mjd2000(tw1G(1,:)) date2mjd2000(tw2G(1,:))]; 
% ub = [date2mjd2000(tw1G(2,:)) date2mjd2000(tw2G(2,:))];
% [xmin,fmin] = fmincon(@cost,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options,ID,vinf);
% 
figure(5*arc-4/K)
% for j = 1:3
%     if j < 3
% nexttile(j)
% p1 = plot3(x0(1),x0(2),Dvmin,'o','MarkerFaceColor',"#D95319",'MarkerSize',8,'MarkerEdgeColor','k','DisplayName','Unconstrained \DeltaV minimum');
% p2 = plot3(xmin(1),xmin(2),fmin,'o','MarkerFaceColor',[0 0.4470 0.741],'MarkerSize',8,'MarkerEdgeColor','k','DisplayName','Constrained     \DeltaV minimum');
% end
% nexttile(1), legend([p1],'Location','southeast','FontSize',15)
% nexttile(1), legend([p1,p2],'Location','southeast','FontSize',15)
nexttile(2), h = findobj(gcf); c = h(8); c.LevelListMode = 'auto';
% legend(p1,'Location','southeast','FontSize',15)

end

figure(5*arc),
% plot3(xmin(1),xmin(2),fmin,'o','MarkerFaceColor','g','MarkerSize',10,'MarkerEdgeColor','k','DisplayName','\Deltav constrained refined')
plot3(x1(1,1),x1(2,1),Dvmin1,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k','DisplayName','\Deltav initial guess unconstrained')
zlim([Dvmin1-3 Dvmin1+N]), %legend('Location','northeast','FontSize',15)
clim([floor(Dvmin1)-1 floor(Dvmin1)+N])

% 
% xmind = [mjd20002date(xmin(1));mjd20002date(xmin(2))];
% % xmind = datetime(xmin,'Format','yyyy/MM/dd hh:mm:ss.SS a') % 12h format
% xmind = datetime(xmind,'Format','yyyy/MM/dd HH:mm:ss.SS');
% [~,TOF,cartesian,dv1] = cost(xmin,ID); 
% 
% % TRANSFER TRAJECTORY PROPAGATION
% Transfer(xmin,[tw1G;tw2G],ID)

%----------------------------------------------------------------------------
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


%-------
% constants - Masses of Sun and Venus and Radius of Venus
G = astroConstants(1); % Universal gavity constant
mS = mu_S/G;
mu_V = astroConstants(12);
mV = mu_V/G;
R_V = astroConstants(22);

%------------------------------------------------------------------------------
% Recombinations of arrays
% m x ( n + l ) --- quadratic + linear / nonlinear operation

Dv_tot = zeros(ll,nn,mm); Dvflyby_tot = Dv_tot; 
% dvp = Dv_tot; rp_ = Dv_tot; dtSoi = Dv_tot;

for m = 1:mm
    C3_dep_fbdate = C3_dep(m,:);
    Vinf_minus_fbdate = Vinf_minus(m,:);
    Vinf_plus_fbdate = Vinf_plus(:,m);
    C3_arr_fbdate = C3_arr(:,m);

%----------
% linear operations

    Dv_tot(:,:,m) = C3_dep_fbdate + C3_arr_fbdate;% + dvp;
%
%     for n = 1:nn
%         Vinf_minus_fbdate(n,:) = Vinf_minus_fbdate_{n}';
%     end
%     for n = 1:ll
%         Vinf_plus_fbdate(l,:) = Vinf_plus_fbdate_{l}';
%     end
%
        for n = 1:nn
            for l = 1:ll
                Vinf_minus_FB = Vinf_minus_fbdate{n};
                Vinf_plus_FB = Vinf_plus_fbdate{l};

                if anynan(Vinf_minus_FB) || anynan(Vinf_plus_FB)
                    Dvflyby_tot(l,n,m) = NaN;
                    continue
                end

            Dvflyby_tot(l,n,m) = norm(Vinf_plus_FB - Vinf_minus_FB);

%             Dvflyby_tot(:,:,m) = norm(Vinf_plus_fbdate - Vinf_minus_fbdate');

            end
        end
    

% nonlinear operation: COMPUTATIONALLY TOO EXPENSIVE HERE IN TIME GRID
% Powered gravity assist fly-by algorithm --> rp --> dvp

% get rV for rSoi (sphere of influence)
% time = TW_fb(m);   % mjd2000 
% kep = uplanet(time,2);
% a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
% rV = kep2car(a,e,i,Omega,omega,theta,mu_S);
% rV = norm(rV);
% 
% rSoi = rV*(mV/mS)^(2/5);
% 
% for n = 1 : nn
%     for l = 1 : ll
% 
%         Vinf_minus_FB = Vinf_minus_fbdate{n};
%         Vinf_plus_FB = Vinf_plus_fbdate{l};
% 
%         if anynan(Vinf_minus_FB) || anynan(Vinf_plus_FB)
%             dvp(l,n,m) = NaN;
%             rp_(l,n,m) = NaN;
%             dtSoi(l,n,m) = NaN;
%             continue
%         end
% 
%         a_minus = -mu_V/norm(Vinf_minus_FB)^2;
%         a_plus = -mu_V/norm(Vinf_plus_FB)^2;
%         rp = R_V; % initial guess close enough to solution (Venus radius)
%         
%         options = optimset('TolX',1e-14,'Display','off');
%         rp = fzero(@fd,rp,options,Vinf_minus_FB,Vinf_plus_FB,mu_V);
% 
%         if rp <= R_V
%             dvp(l,n,m) = NaN;
%             rp_(l,n,m) = NaN;
%             dtSoi(l,n,m) = NaN;
%             continue
%         end
% 
%         [~,d,e_minus,e_plus,d_minus,d_plus] = fd(rp,Vinf_minus_FB,Vinf_plus_FB,mu_V);
% 
%         vp_minus = sqrt(norm(Vinf_minus_FB)^2+2*mu_V/rp);
%         vp_plus = sqrt(norm(Vinf_plus_FB)^2+2*mu_V/rp);
%         dvp(l,n,m) = abs(vp_minus-vp_plus);
%         rp_(l,n,m) = rp;
% 
%         E_minus = acosh((-rSoi/a_minus+1)/e_minus);
%         E_plus = acosh((-rSoi/a_plus+1)/e_plus);
%         
%         dt_minus = (e_minus*sinh(E_minus)-E_minus)/sqrt(-mu_V/a_minus^3);
%         dt_plus = (e_plus*sinh(E_plus)-E_plus)/sqrt(-mu_V/a_plus^3);
% 
%         dtSoi(l,n,m) = dt_minus + dt_plus;
%     end
% end
% 
% % Final sum
% 
% Dv_tot(:,:,m) = Dv_tot(:,:,m) + dvp(:,:,m);

end

% experimental: try adding the full Dv, even though it includes natural fly
% by, non powered

Dv_tot = Dv_tot + Dvflyby_tot;

%-----------------------------------------------------------------------------
% choose date M with minimum Dv_tot

h_atm_venus = 500; % [km]
rp = R_V+h_atm_venus; % inizial guess for fzero
KK = 0;
% REMOVE collision prevention COMPUTING SO TO SPEED UP THE PROCESS
while rp <= R_V + h_atm_venus   %--> a method to exclude the colliding solutions
% while 1   % minimal cycle to avoid NaN solutions
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
%     rp = NaN;
    dtSoi = NaN;
    Dv_tot(index) = NaN; % crucial to restart minimization without this solution
    continue
end

a_minus = -mu_V/norm(Vinf_minus_FB)^2;
a_plus = -mu_V/norm(Vinf_plus_FB)^2;
rp = R_V+h_atm_venus; % initial guess close enough to solution (Venus radius)

options = optimset('TolX',1e-6,'Display','off');
rp = fzero(@fd,[rp],options,Vinf_minus_FB,Vinf_plus_FB,mu_V);

if rp <= R_V + h_atm_venus
    dvp = NaN;
    rp = R_V;
    dtSoi = NaN;
    Dv_tot(index) = NaN;  % crucial to restart minimization without this solution
    continue
end
break
end
[~,d,e_minus,e_plus,d_minus,d_plus] = fd(rp,Vinf_minus_FB,Vinf_plus_FB,mu_V);

vp_minus = sqrt(norm(Vinf_minus_FB)^2+2*mu_V/rp);
vp_plus = sqrt(norm(Vinf_plus_FB)^2+2*mu_V/rp);
dvp = abs(vp_minus-vp_plus);
%         rp_(l,n,m) = rp;

E_minus = acosh((-rSoi/a_minus+1)/e_minus);
E_plus = acosh((-rSoi/a_plus+1)/e_plus);

dt_minus = (e_minus*sinh(E_minus)-E_minus)/sqrt(-mu_V/a_minus^3);
dt_plus = (e_plus*sinh(E_plus)-E_plus)/sqrt(-mu_V/a_plus^3);

dtSoi = dt_minus + dt_plus;
dtSOI = seconds(dtSoi);
dtSOI.Format = 'hh:mm:ss.SS';


% Final sum

Dv_totmin = Dv_totmin - Dv_inf + dvp;
%

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
    continue
end
    
%---------------------------------------------------------------------------
% FMINCON
X0 = [TW_dep(n);TW_fb(m);TW_arr(l)]; % initial guess from time grid
% X0 = [10837.4268;11716.0648;12981.0216];
% X0 = [10226.5000;10253.5853;11717.0570];
% X0 = [11283.0761;11494.5950;11800.7779];
% X0 = [11543.4532;11706.9236;12948.2478];
% X0 = [10409.5423;10608.4813;11442.8291];
% X0 = [10237.2961;10386.4152;11419.0214];
% X0 = [12610.2430;12833.3477;13044.6101];

% Dvtotmin = cost_tot(X0,h_atm_venus,T)
options = optimoptions("fmincon",'Display','off','Algorithm','active-set', ...
                                                               ...'sqp', ...
                                                               ...'interior-point',...
'ConstraintTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-6, ...
'MaxFunctionEvaluations',1e4,'MaxIterations',9e3);

A = []; b = []; Aeq = []; beq = []; % Linear equality & inequality constraints
% lower & upper bounds 
% lb = [date2mjd2000(tw1G(1,:)) date2mjd2000(tw2G(1,:))]; 
% ub = [date2mjd2000(tw1G(2,:)) date2mjd2000(tw2G(2,:))];
% lb = [TW_dep(1) TW_fb(1) TW_arr(1)];
% ub = [TW_dep(end) TW_fb(end) TW_arr(end)];
lb = [X0(1)-20 X0(2)-10 X0(3)-1300]; %lb = [];   % experimental w/o boundaries
ub = [X0(1)+30 X0(2)+10 X0(3)+1500]; %ub = [];
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
    fprinft('Fmincon not executed')
    toc(tic1)
    continue
end

if anynan(Xmin)
    options.Algorithm = 'interior-point';
    [Xmin,Dvtotmin] = fmincon(@cost_tot,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options,h_atm_venus,T);
    if anynan(Xmin)
        DVTOT(D).num = D;
        DVTOT(D).Dvtot = NaN;
        DVTOT(D).dvfb_ratio = NaN;
        toc(tic1)
        continue
    end
end

% Xmin = X0;

TW_DEP = datetime(mjd20002date(Xmin(1)),'Format','yyyy MMM dd hh:mm:ss'); 
TW_FB = datetime(mjd20002date(Xmin(2)),'Format','yyyy MMM dd hh:mm:ss'); 
TW_ARR = datetime(mjd20002date(Xmin(3)),'Format','yyyy MMM dd hh:mm:ss'); 

%--------------------------------------------------------------------------
% Radius SOI

% mu_S = astroConstants(4);
G = astroConstants(1); % Universal gavity constant
mS = mu_S/G;
mu_V = astroConstants(12);
mV = mu_V/G;
R_V = astroConstants(22);
kep = uplanet(Xmin(2),2);

a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
rV = kep2car(a,e,i,Omega,omega,theta,mu_S);
rV = norm(rV);

rSoi = rV*(mV/mS)^(2/5);

%--------------------------------------------------------------------------
% D time SOI

% a_minus = -mu_V/norm(Vinf_minus_FB)^2;
% a_plus = -mu_V/norm(Vinf_plus_FB)^2;

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
%  fb date refined cost 
% tic

lb = [Xmin(1)-20 Xmin(2)-10 Xmin(3)-1300]; 
ub = [Xmin(1)+30 Xmin(2)+10 Xmin(3)+1500];

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
% vinf_minus = norm(Vinf(1,:));
% vinf_plus = norm(Vinf(2,:));
%
for s = 1:2
    figure(1), nexttile(s)
    plot3(Xmin(1),Xmin(2),C3_dep,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 1, legend('','$\mathbf{C_{3 departure}}$ min','$\mathbf{C_{3 departure}}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(3), nexttile(s)
    plot3(Xmin(1),Xmin(2),Vinf_minus,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 1, legend('','$\mathbf{v_{\infty}^-}$ min','$\mathbf{v_{\infty}^-}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(6), nexttile(s)
    plot3(Xmin(2),Xmin(3),Vinf_plus,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 1, legend('','$\mathbf{v_{\infty}^+}$ min','$\mathbf{v_{\infty}^+}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(8), nexttile(s)
    plot3(Xmin(2),Xmin(3),C3_arr,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10)
    if s == 1, legend('','$\mathbf{C_{3 arrival}}$ min','$\mathbf{C_{3 arrival}}$ mission (min $\mathbf{\Delta v_{tot}}$)','FontSize',15,'interpreter','latex','Location','southeast'), end
    figure(12), nexttile(s)
    if s == 1, legend('','$\mathbf{\Delta v_{tot}}$ min from time grid','$\mathbf{\Delta v_{tot}}$ min from refinement','FontSize',15,'interpreter','latex','Location','southeast'), end
    if s == 2, legend('','$\mathbf{\Delta v_{\infty}}$ mission from time grid','$\mathbf{\Delta v_{\infty}}$ mission from refinement','FontSize',15,'interpreter','latex','Location','southeast'), end
end

%
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

toc(tic1)
end

dd = [DVTOT.Dvtot]';
ddd = [DVTOT.dvfb_ratio]';
[DVTOTMIN,INDEX] = min(dd);
fprintf('\n# %d\n%cv_tot: %.3f km/s\n%cv_GA/%cv%c: %.3e\n\n', ...
    INDEX,char(916),DVTOTMIN,char(916),char(916),char(8734),ddd(INDEX))



toc(tic0)

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

Tw1 = datetime(tw1(1,:)):calmonths:datetime(tw1(2,:));
tw2 = zeros(2,6);
for j = 1:2
    tw1(j,1) = date2mjd2000(tw1(j,:)); % [d]
    tw2(j,1) = tw1(j,1) + ToF(j) * 365.2417;
    tw2(j,:) = mjd20002date(tw2(j,1));
end

Tw2 = datetime(tw2(1,:)):calmonths:datetime(tw2(2,:));

for j = 1:2
%     tw1(j,1) = date2mjd2000(tw1(j,:)); % [d]
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
        title('Porkchop plot: D & A Dates','FontSize',18)
        yticks(ttw2)
        yticklabels(char(Tw2)), ytickangle(45)
        if j == 2
            nexttile(2)
            I = find(TOF<0);
            TOF(I) = 0;
            contour(tw1,tw2,TOF/(3600*24),60:60:300,'ShowText','on','LineWidth',1,'EdgeColor','k')
            title(sprintf('Porkchop plot (%st lines [days])',char(916)),'FontSize',18,'Interpreter','latex')
        end
    else
        figure('Name',sprintf('PORKCHOP PLOT: DDate & ToF'),'NumberTitle','off')
        [~,h] = contour(tw1,ToF,Dv_,'ShowText',str(1),'LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
        ylabel('ToF [years]','FontSize',15)
        title('Porkchop plot: Departure Date & ToF','FontSize',18)
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
% c.Label.Position = [0.564177274703979,10.239066930161309,0];
c.Label.Rotation = 0;
% h.LevelStep = 1;
% h.LevelList = 5:10;
% clim([5 10])  % CRUCIAL!!! for colorbar domain of Dv
xlabel('Departure date','FontSize',15)
% ylabel('Arrival date','FontSize',15)
% title('Porkchop plot: D & A Dates','FontSize',18)
xticks(ttw1)
% yticks(ttw2)
xticklabels(char(Tw1)), xtickangle(45)
% yticklabels(char(Tw2)), ytickangle(45)
end
end

% nexttile(2)
% I = find(TOF<0);
% TOF(I) = 0;
% contour(tw1,tw2,TOF/(3600*24),60:60:300,'ShowText','on','LineWidth',1,'EdgeColor','k')
% title(sprintf('Porkchop plot (%st lines [days])',char(916)),'FontSize',18,'Interpreter','latex')

figure('Name',sprintf('Debug mesh plot of \x394V')), 
% tiledlayout(1,1,"TileSpacing","tight","Padding","tight")
mesh(tw1,tw2,Dv1,'DisplayName','\Deltav'), hold on, grid on, axis padded
xlabel('x'), ylabel('y'), zlabel('z')
view(0,0)

% figure('Name',sprintf('PORKCHOP PLOT: DDate & ToF'),'NumberTitle','off')

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

% mu_S = astroConstants(4);
% G = astroConstants(1); % Universal gavity constant
% mS = mu_S/G;
mu_V = astroConstants(12);
% mV = mu_V/G;
R_V = astroConstants(22);
% kep = uplanet(X(2),2);
% 
% a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); theta = kep(6);
% rV = kep2car(a,e,i,Omega,omega,theta,mu_S);
% rV = norm(rV);
% 
% rSoi = rV*(mV/mS)^(2/5);

a_minus = -mu_V/v_inf_minus^2;
a_plus = -mu_V/v_inf_plus^2;
rp0 = R_V+h_atm_venus; % initial guess close enough to solution (Venus radius)

options = optimset('TolX',1e-14,'Display','off');
rp = fzero(@fd,rp0,options,V_inf_minus,V_inf_plus,mu_V);

% if rp <= rp0  % experimental
%     Dvtot = NaN;
%     return
% end
% [~,d,e_minus,e_plus,d_minus,d_plus] = fd(rp,Vinf_minus,Vinf_plus,mu_V);

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

function Transfer(xmin,tw,ID)

% NOT ACTIVE!!!!

% TRANSFER TRAJECTORY PROPAGATION + PLOT



switch ID(2)
    case 1
        name = 'Mercury'; TEX = '2k_mercury.jpg';
    case 2
        name = 'Venus';   TEX = '2k_venus_surface.jpg';
    case 81
        name = '2005EL70'; TEX = '2k_moon.jpg';
end

tw1G = tw(1:2,:); tw2G = tw(3:4,:);
[~,TOF,cartesian] = cost(xmin,ID);

km2AU = 1/astroConstants(2); % [1AU/km]
MU = astroConstants(4); % Sun gravitational parameter
RI = cartesian(:,1); RF = cartesian(:,2); 
vi = cartesian(:,3); vf = cartesian(:,4);
VI = cartesian(:,5); VF = cartesian(:,6);

% Earth
[t,xi] = R2BP_propagation([RI;vi],[0 TOF],MU); % Earth orbit during Transfer
ai = 1/(2/norm(RI) - dot(vi,vi)/MU);
Ti = 2*pi*sqrt(ai^3/MU);
[t,xi_] = R2BP_propagation(xi(end,:),[0 Ti-TOF],MU); % Earth orbit
% Departure window
tw1s(1) = date2mjd2000(tw1G(1,:))*3600*24; % [s] beginning of initial time window
tw1s(2) = date2mjd2000(tw1G(2,:))*3600*24; % [s] ending of initial time window
[t,dwin] = R2BP_propagation([RI;vi],[0 -(xmin(1)*3600*24-tw1s(1))],MU); % backward propagation to get beginning state
[t,dwin] = R2BP_propagation(dwin(end,:),[0 tw1s(2)-tw1s(1)],MU); 
% Transfer
[t,xt] = R2BP_propagation([RI;VI],[0 TOF],MU); % Transfer arc
at = 1/(2/norm(RI) - dot(VI,VI)/MU);
Tt = 2*pi*sqrt(at^3/MU);
[t,xt_] = R2BP_propagation(xt(end,:),[0 Tt-TOF],MU); % unused transfer arc
% SECOND PLANET
[t,xf] = R2BP_propagation([RF;vf],[0 -TOF],MU); % SECOND PLANET orbit during Transfer
af = 1/(2/norm(RF) - dot(vf,vf)/MU);
Tf = 2*pi*sqrt(af^3/MU);
[t,xf_] = R2BP_propagation(xf(1,:),[0 Tf-TOF],MU); % Mars orbit
% Arrival window
tw2s(1) = date2mjd2000(tw2G(1,:))*3600*24; % [s] beginning of final time window
tw2s(2) = date2mjd2000(tw2G(2,:))*3600*24; % [s] ending of final time window
[t,awin] = R2BP_propagation([RF;vf],[0 -(xmin(2)*3600*24-tw2s(1))],MU); % backward propagation to get beginning state
[t,awin] = R2BP_propagation(awin(end,:),[0 tw2s(2)-tw2s(1)],MU); 

xi = km2AU*xi; xi_ = km2AU*xi_; dwin = km2AU*dwin;
xt = km2AU*xt; xt_ = km2AU*xt_; awin = km2AU*awin;
xf = km2AU*xf; xf_ = km2AU*xf_;

az = [-37.5,0];
el = [30,90];

% Rescaling for visibility, possibly maintain relative dimension between planets
Re = astroConstants(20+ID(1));    Re = 3e3 * Re; % 3e3 default
Rp = astroConstants(20+ID(2));    Rp = 3e3 * Rp; % 3e3 default
Rs = astroConstants(3);           Rs = 1e2 * Rs; % 4e1 default

% Another customized rescale specific per some missions
switch ID(2)
    case 1
        Re = 2e3 * astroConstants(23);
        Rp = 2e3 * astroConstants(21);
        Rs = 5e1 * astroConstants(3);
    case 2
        Re = 2.5e3 * astroConstants(23);
        Rp = 2.5e3 * astroConstants(22);
        Rs = 8e1 * astroConstants(3);
    
end


% PLOT
f = figure('Name',sprintf('EX4 2: Transfer trajectory - %s Express',name),'NumberTitle','off');
T1 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for j = 1:2
    nexttile(T1,j)
    % EARTH orbit
    globe("2k_earth_daymap.jpg",Re*km2AU,xi(1,1:3)); hold on, grid on, axis equal padded, pbaspect([1,1,1])
    globe("2k_earth_daymap.jpg",Re*km2AU,xi(end,1:3));
    p1 = plot3(xi_(:,1),xi_(:,2),xi_(:,3),'color',[0 0.4470 0.7410],'LineWidth',3,'LineStyle','--');
    p2 = plot3(xi(:,1),xi(:,2),xi(:,3),'color',[0 0.4470 0.7410],'LineWidth',3);
    p3 = plot3(dwin(:,1),dwin(:,2),dwin(:,3),'color',[0 .4470 .7410 .3],'LineWidth',15); % use {'Color',[R G B alpha]} for transparency | [.7882 .9176 1]
    % SECOND PLANET orbit
    globe(TEX,Rp*km2AU,xf(1,1:3));
    globe(TEX,Rp*km2AU,xf(end,1:3));
    p4 = plot3(xf_(:,1),xf_(:,2),xf_(:,3),'color',[0.8500 0.3250 0.0980],'LineWidth',3,'LineStyle','--');
    p5 = plot3(xf(:,1),xf(:,2),xf(:,3),'color',[0.8500 0.3250 0.0980],'LineWidth',3);
    p6 = plot3(awin(:,1),awin(:,2),awin(:,3),'color',[.8500 .3250 .098 .3],'LineWidth',15); % [.9804 .8863 .7451]
    % SUN
    globe("2k_sun.jpg",Rs*km2AU);
    % TRANSFER
    p7 = plot3(xt_(:,1),xt_(:,2),xt_(:,3),'color',[0.4660 0.6740 0.1880],'LineWidth',3,'LineStyle','--');
    p8 = plot3(xt(:,1),xt(:,2),xt(:,3),'color',[0.4660 0.6740 0.1880],'LineWidth',3);
    xlabel('$x$ [AU]','FontSize',15,'interpreter','latex')
    ylabel('$y$ [AU]','FontSize',15,'interpreter','latex')
    zlabel('$z$ [AU]','FontSize',15,'interpreter','latex')
    title(sprintf('%s Express',name),'FontSize',12)
%     lgn = legend([p1 p2 p3 p4 p5 p6 p7 p8],'Earth orbit','Earth motion during transfer','Departure window', ...
%           sprintf('%s orbit',name),sprintf('%s motion during transfer',name),'Arrival window','Unused arc of transfer orbit','Transfer arc', ...
%           'FontSize',15,'interpreter','latex','Location','northeast');
    view(az(j),el(j))
end
delete(lgn)

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

% Tw1 = datetime(tw_dep(1,:)):calmonths:datetime(tw_dep(end,:));
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
    case 2
        X = Xmin;
        DV = Dvtot;
        DVfb = Dvflybytot;
        DVtotmin = Dvtotmin;
        DVinfmin = Dvinfmin;
        TW_dep = linspace(lb(1),ub(1),length(TW_dep));
        TW_arr = linspace(lb(3),ub(3),length(TW_arr));
end

    figure
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
    nexttile
    [~,h] = contour(TW_dep,TW_arr,DV,'ShowText','on','LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
    plot3(X0(1),X0(3),Dv_totmin,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k')
    if s == 2
        plot3(X(1),X(3),Dvtotmin,'o','MarkerFaceColor','g','MarkerSize',10,'MarkerEdgeColor','k');
    end
    xlabel('Departure date','FontSize',15)
    ylabel('Arrival date','FontSize',15)
    title('Cost plot: D & A Dates','FontSize',18)
    % h.LevelStep = 1;
    h.LevelList = floor(DVtotmin)-1:STEP/2:floor(DVtotmin)+N;
    c = colorbar(gca,"eastoutside");
    c.Label.String = '$\mathbf{\Delta v_{tot}\ [km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 15;
    c.Label.Rotation = 0;
    c.Label.Position = [.564177274703979,floor(DVtotmin)+N+.9,0];
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
    [~,h] = contour(TW_dep,TW_arr,DVfb,'ShowText','on','LineWidth',1); grid on, hold on, axis tight, pbaspect([1,1,1])
    plot3(X0(1),X0(3),Dv_infmin,'o','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k')
    if s == 2
        plot3(X(1),X(3),DVinfmin,'o','MarkerFaceColor','g','MarkerSize',10,'MarkerEdgeColor','k');
    end
    xlabel('Departure date','FontSize',15)
    ylabel('Arrival date','FontSize',15)
    title('Cost plot: D & A Dates','FontSize',18)
    % h.LevelStep = 1;
    h.LevelList = floor(DVinfmin)-1:STEP/2:floor(DVinfmin)+N;
    c = colorbar(gca,"eastoutside");
    c.Label.String = '$\mathbf{\Delta v_{\infty}\ [km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 15;
    c.Label.Rotation = 0;
    c.Label.Position = [.564177274703979,floor(DVinfmin)+N+.9,0];
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


