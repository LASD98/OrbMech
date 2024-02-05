function colorPlot3(T,Y,chosenColormap)

% Plots the color-coded time evolution of a given state vector.
% 
% PROTOTYPE:
%   colorPlot3(T,Y,chosenColormap)
% 
% INPUT:
%   T[n] = time vector, can either be row or column [s]
%   Y[n,6] = Cartesian state vector, positions first
%   chosenColormap = chosen colormap (i.e. parula, hsv, hot, cold just to
%   name a few)
% 
% OUTPUT:
%   (none)
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%    Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

if nargin < 3
    chosenColormap = hsv;
end

c = T;

cmap = colormap(chosenColormap);

c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));

x = Y(:,1);
y = Y(:,2);
z = Y(:,3);

%%
hold on

for ii = 1:length(T)-1
    plot3([x(ii),x(ii+1)],[y(ii),y(ii+1)],[z(ii),z(ii+1)],'color',cmap(c(ii),:))
end