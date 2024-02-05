function costplots(TW_dep,TW_arr,Dv_tot,Dvflyby_tot,X0,Dv_totmin,Dv_infmin, ...
    Dvtot,Dvflybytot,Xmin,Dvtotmin,Dvinfmin,lb,ub,STEP,N)

% this function plots the deltav cost for recombined doublets 
% of departure and arrival dates at same flyby date

% INPUTs: 
%	TW_dep[n],TW_arr[l]: grid-search used dates, 
%			     for departure and arrival windows [MJD2000, days]
%	Dv_tot,Dvflyby_tot,X0,Dv_totmin,Dv_infmin: [unused] as below, 
%			       but from grid search optimal flyby date 
%			       (included in X0, fmincon first guess) [km/s]
%	Dvtot,Dvflybytot[nxl]: grid of costs at solution flyby date from fmincon, 
%       		       total dv and flyby dv [km/s]
%	Xmin[3]: solutions of optimization problem from fmincon, dates [MJD2000, days]
%	Dvtotmin, Dvinfmin[1]: solutions of optimization problem from fmincon,
%			       total cost & flyby cost [km/s]
%	lb,ub [3]: lower & upper boundaries, limits for plotting [days]
%	STEP[1]: step between contour plot lines in Dvtot units
%	N[1]: N+1 units of Dvtot will be shown in contour plot

% OUTPUTs:

% FUNCTIONS CALLED:


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
