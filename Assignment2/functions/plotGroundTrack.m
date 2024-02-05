function p = plotGroundTrack(lon,lat,color,i)

% Plots the ground track corresponding to the given latitudes and
% longitudes. It can optionally check whether the orbit is retrograde or
% prograde to display the ground track correctly.
% 
% PROTOTYPE:
%   p = plotGroundTrack(lon,lat,color,i)
% 
% INPUT:
%   lon[n] = longitude vector [deg]
%   lat[n] = latitude vector [deg]
%   color = ground track color, can either be a character array or a RGB
%       decimal triplet.
%   i[1] = inclination of the orbit [rad]
% 
% OUTPUT:
%   p[3] = vector of initial point, final point and ground track proper
%       (line element)
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

hold on

if nargin < 4
    i = 1; % if not specified, assume prograde orbit
    if nargin < 3
        color = 'g'; % if not specified, color is green
    end
end

ps = plot(lon(1),lat(1),'Marker','o','Color',color,'linestyle','none','LineWidth',4);
pe = plot(lon(end),lat(end),'Marker','d','Color',color,'linestyle','none','LineWidth',4);

for ii = 1:length(lon)-1
    if i < pi/2
        if lon(ii+1)*lon(ii) < 0 && lon(ii+1) < 0
%             plot(rad2deg([lon(ii-1),pi]),rad2deg(lat([ii-1,ii])),'Color','g')
%             plot(rad2deg([-pi,lon(ii+1)]),rad2deg(lat([ii,ii+1])),'Color','g')
            lon(ii) = NaN;
        end
    else
        if lon(ii+1)*lon(ii) < 0 && lon(ii+1) > 0
            lon(ii) = NaN;
        end
    end
end

pp = plot(lon,lat,'.','Color',color,'MarkerSize',2);
% plot(lon(a:end),lat(a:end),'LineWidth',1,'Color','g');

hold off

p = [ps pp pe];