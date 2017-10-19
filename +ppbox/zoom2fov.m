function [ fovx, fovy ] = zoom2fov(zoom, micID, objID)

%   INPUTS: 

%   zoom:  the query zoom 
%   objID:  the objective used, either 'x16' (default) or 'x20'.
%   micID:  the microscope used, either 'bscope' (default) or 'mom' or ' b2'. 

%   zoom2fov(zoom) returns the width in micrometers of the FOV of the 
%   microscope at a query zoom.
%   zoom2foc(zoom, micID) allows you to specify which microscope you used
%   zoom2foc(zoom, micID, objID) allows you to specify which scope and
%   which objective you used for imaging.

%   OUTPUTS: 
%   fovx:  the size of the FOV along scanning lines
%   fovy:  the size of the FOV across scanning lines

%   2016-01-24 function and 2p measurements by L. Federico Rossi
%   2016-02-05 adapted from infoPixUm in Suite2P folder (SS)

if nargin < 2
    micID = 'bscope';
elseif ~any(strcmp(micID, {'b', 'bscope', 'b2', 'mom'}))
    display('WARNING: unknown type of microscope.')
end

% we assume a 16x objective in all cases so far
if nargin < 3
    objID = 'x16';
end

switch micID
    case {'bscope', 'b'}
        zooms = [1.3, 1.4, 1.5, 1.6, 1.7, 2, 2.2, 2.3, 2.5, 3, ...
            3.5, 4, 5, 6.1, 7.1, 9.1];
        measuredHoriz = [855, 740, 680, 635, 615.5, 524, 491, 460, 431.5, 371, ...
            318, 287, 237.5, 199, 172.5, 140];
        measuredVert = [790, 715, 664.5, 615, 593, 503.5, 452.5, 430, 401, 337.5, ...
            284.5, 249, 179.5, 161, 141, 111.5];
    case 'b2'
        zooms = [1.6 1.9 2 2.2];
        measuredHoriz = [772 653 619.5 565.5];
        measuredVert = [755 640 605.5 547];
    case 'mom'
        zooms = [3 4];
        measuredHoriz = [155 117];
        measuredVert = [155 117];
end

curveX = fit(zooms', measuredHoriz', 'exp2');
fovx = curveX(zoom);
curveY = fit(zooms', measuredVert', 'exp2');
fovy = curveY(zoom);

%% Visualize fit, measured data, and queried zoom
% figure
% hold on
% x = min(1,zoom):0.05:max(10,zoom);
% plot(x, curveX(x))
% plot(zooms, measuredHoriz, 'ko')
% plot(zoom, fovx, 'ko', 'MarkerFaceColor', 'k')
% title('Field of view in X')
% figure
% hold on
% x = min(1,zoom):0.05:max(10,zoom);
% plot(x, curveY(x))
% plot(zooms, measuredVert, 'ko')
% plot(zoom, fovy, 'ko', 'MarkerFaceColor', 'k')
% title('Field of view in Y')

end

