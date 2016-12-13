function [aligned, offset, scenes] = vaStimuli(varargin)
% Create the pair of vernier stimuli
%
% The scene is made at 2x the FOV of the planned cone mosaic size.
%
% The offset of 1 pixel for a test with 210,210 pixels and a field of view
% of 0.35 deg is 6 sec.  That is a useful value. If you scale the number of
% pixels and the field of view together, then the offset remains the same.
%
% TODO:
%   * Control line length 
%   * Control the gap between the lines
%
% Notes
%   Westheimer and McKee, 1977
%   Luminance 2 log millilamberts
%   Line width 1 minute of arc.
%

%%
p = inputParser;
p.addParameter('barWidth',1,@isscalar);
p.addParameter('barOffset',1,@isscalar);

p.parse(varargin{:});

barWidth  = p.Results.barWidth;
barOffset = p.Results.barOffset;

%% Display model

display = displayCreate('LCD-Apple');

%% Gaussian time series.

% We should make some functions 
std = 0.2;
tseries = ieScale(gausswin(150,1/std),0,1);
% vcNewGraphWin; plot(tseries);
% tseries = ieScale(fspecial('gaussian',[1,150],40),0,1);

clear sparams;  % Scene parameters in general

% The size of the integration for the cone mosaic is 5 minutes for a line.
% There are two lines, so the stimulus should be more than 10 minutes.  We
% set up the cone mosaic to be about 15 minutes, and we will try shrinking
% it later.  And enlarging it.
%
% We create the scene and thus the oi to be a bit bigger to allow for eye
% movements.  We set the scene to be .35*60 ~ 21 minutes.  The oi is even a
% little bigger to allow for blur.
sparams.fov      = (0.35); 
sparams.distance = (0.3);    % Meters

% Basic vernier parameters for the oiSequence
clear vparams;
for ii = 3:-1:1
    vparams(ii) = vernierP;
    vparams(ii).display = display;
    % For a fov of 0.35 and this size, 1 pixel offset is 6 sec of arc
    vparams(ii).sceneSz = ([210 210]);
end

% Uniform field, no line
vparams(1).name     = 'uniform'; 
vparams(1).bgColor  = 0.5; 
vparams(1).barWidth = 0;

% Offset Line on a zero background
vparams(2).name     = 'offset';  
vparams(2).bgColor  = 0; 
vparams(2).offset   = barOffset;
vparams(2).barWidth = barWidth;   % The bar width is 1 minute, 10 * 6sec

% Aligned lines on a zero background
vparams(3).name     = 'aligned'; 
vparams(3).bgColor  = 0; 
vparams(3).offset   = 0;
vparams(3).barWidth = barWidth;

[offset, scenes] = oisCreate('vernier','add', tseries, ...
    'testParameters',vparams([1 2]),...
    'sceneParameters',sparams);
% offset.visualize;
% ieAddObject(offset.oiFixed); ieAddObject(offset.oiModulated); oiWindow;
% ieAddObject(scenesO{2}); sceneWindow;

% Offset lines
offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

aligned = oisCreate('vernier','add', tseries,...
    'testParameters',vparams([1 3]),...
    'sceneParameters',sparams);
% aligned.visualize;

end
