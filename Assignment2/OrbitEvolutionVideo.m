% Generates the orbit evolution video
% Called at the end of the main script
%
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

warning('off')

tic
fprintf('\nOrbit evolution video\n')

figure
Earth3d
hold on

c = case3.time;

chosenColormap = jet;
cmap = colormap(chosenColormap);

c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));

cbar = colorbar;
cbar.Label.String = 'Time [days]';

cbar.Ticks = ([3 258*1/10 258*2/10 258*3/10 258*4/10 258*5/10 258*6/10 258*7/10 258*8/10 258*9/10 255]);
cbar.TickLabels = ([0 36 36*2 36*3 36*4 36*5 36*6 36*7 36*8 36*9 360]);

grid on
pbaspect([1 1 1])

xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

view(40,20)

nn = length(case3.time);

%%
F = [];

for jj = 1:100:nn
    
    o1 = plotOrbit(case3.a(jj),case3.e(jj),case3.i(jj),case3.Omega(jj),case3.omega(jj),0,2*pi,2*pi/100);
    set(o1,'Color',cmap(c(jj),:))
    F = [F getframe(gcf)];
    pause(0.000001)
end

for kk = 1:40
F = [F getframe(gcf)];
end

%% Video writing
fprintf('Encoding initiated\n')

video = VideoWriter('OrbitEvolution','MPEG-4');
set(video,'FrameRate',30);
open(video)
writeVideo(video,F);
close(video)
toc
fprintf('Encoding complete\n')