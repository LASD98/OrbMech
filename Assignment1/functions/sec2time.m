function [h, min, sec] = sec2time (s)

%--------------------------------------------------------------------------
% 
% This function converts time expressed in seconds into hours, minutes, 
% seconds 
%
% PROTOTYPE:
% [h, min, sec] = sec2time (s)
%
% INPUT:
% Time expressed in seconds  [s]
%
% OUTPUT:
% Time expressed in hours, minutes, seconds  [h, min, sec]
%--------------------------------------------------------------------------

s = floor(s);

minutes = floor(s/60);
sec = floor( (s/60 - floor(minutes))*60 );

hours = floor(minutes/60);
min = floor( (minutes/60-floor(hours))*60 );

h = floor( hours );

end