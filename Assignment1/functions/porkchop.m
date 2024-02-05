function [Dv1,Dv2,Dv1_,Dv2_,DV1,DV2] = porkchop(tw1,ToF,n,ID,arc,TOFMAX)

% INPUTs:
%       tw1[2 char]: departure time window, (Gregorian) [YYYY,MM,DD,hh,mm,ss]  
%       tw2[2 char]: arrival time window // ToF: Time of Flight[2] [years]
%       n[3]: number of discretization points for tw1 & tw2 (& ToF)
%       ID[2]: bodies identifier for uplanet / ephNEO
%       arc[1]: 1 or 2, first or second arc
%	    TOFMAX[1]: maximum ToF constraint [days]

% OUTPUTs:
%       Dv1[nxm]: grid of cost of first maneuver [km/s]
%	    Dv2[nxm]: grid of cost of second maneuver [km/s]
%       DV1[3xnxm]: cell grid of vectorial cost of first maneuver [km/s]
%	    DV2[3xnxm]: cell grid of vectorial cost of second maneuver [km/s]
%       Dv1_[nxm]: [unused] cost of first maneuver [km/s]
%	    Dv2_[nxm]: [unused] cost of second maneuver [km/s]

% FUNCTIONS CALLED:
% cost.m

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
