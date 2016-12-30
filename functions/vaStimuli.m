function [aligned, offset, scenes, tseries] = vaStimuli(varargin)
% Create the pair of vernier stimuli (aligned and offset)
%
% There is usually one input argument that is a struct with these
% parameters
%
%  vernier   - Parameters for the vernier stimuli; Default is vernierP
%  tsamples  - Time samples (sec)
%  timesd    - Time standard deviation
%  display   - Display struct, default displayCreate('LCD-Apple')
%  sceneFOV  - Degrees, default 0.35
%  distance  - Meters, viewing distance to display, default is 0.3 m
%
% When the LCD-Apple display is 210, 210 pixels and this is 0.35, then 1
% pixel is 6 arc sec

%%
% The default display is an Apple LCD monitor.
%
% The offset of 1 pixel for a test with 210,210 pixels and a field of view
% of 0.35 deg is 6 sec.  That is a useful value. If you scale the number of
% pixels and the field of view together, then the offset remains the same.
%
% Examples:
%    clear p; p.barOffset = 3; p.barWidth  = 3;
%    p.tsamples = [-60:70]; p.timesd = 20;  
%    [~,offset,~,tseries] = vaStimuli(p);
%    offset.visualize;
%    vcNewGraphWin; plot(tseries)
%
% TODO:
%   * Control line length 
%   * Control the gap between the lines
%   * Control the background level
%   * Control the color and contrast of the lines
%   
%
% See also:
%   imageVernier, sceneVernier, sceneCreate, s_vaAbsorptions.m
%
% Notes
%   Westheimer and McKee, 1977
%   Luminance 2 log millilamberts
%   Line width 1 minute of arc.
%
% BW, ISETBIO Team, 2016

%%
p = inputParser;

p.KeepUnmatched = true;   % Sometimes we overload p for SVM and cMosaic

p.addParameter('vernier',vernierP,@isstruct);

p.addParameter('tsamples',(-50:100),@isvector);  % Time samples (ms)
p.addParameter('timesd',20,@isscalar);           % Time standard deviation
p.addParameter('display',displayCreate('LCD-Apple'),@isstruct);

% When the LCD-Apple display is 210, 210 pixels and this is 0.35, then 1
% pixel is 6 arc sec
p.addParameter('sceneFOV',0.35,@isscalar);  % Degrees. 
p.addParameter('distance',0.3,@isscalar);   % Meters
p.parse(varargin{:});

vernier   = p.Results.vernier;
barWidth  = vernier.barWidth;
barOffset = vernier.offset;
bgColor   = vernier.bgColor;

tsamples  = p.Results.tsamples;
timesd    = p.Results.timesd;
display   = p.Results.display;

sceneFOV = p.Results.sceneFOV;
distance = p.Results.distance;

%% Build Gaussian time series.

tseries = exp(-(tsamples/timesd).^2);
tseries = ieScale(tseries,0,1);

%%  Scene parameters in general
clear sparams; 

% The size of the integration for the cone mosaic is 5 minutes for a line.
% There are two lines, so the stimulus should be more than 10 minutes.  We
% set up the cone mosaic to be about 15 minutes, and we will try shrinking
% it later.  And enlarging it.
%
% We create the scene and thus the oi to be a bit bigger to allow for eye
% movements.  We set the scene to be .35*60 ~ 21 minutes.  The oi is even a
% little bigger to allow for blur.
sparams.fov      = sceneFOV;
sparams.distance = distance;    % Meters

% Basic vernier parameters for the oiSequence.  Reverse order forces the
% allocation first so the array does not grow over the loop.
clear vparams;
for ii = 3:-1:1
    vparams(ii) = vernier;
    vparams(ii).display = display;
end

% Uniform field, no line, just the background color
vparams(1).name     = 'uniform'; 
vparams(1).bgColor  = bgColor; 
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

% Offset lines
[offset, scenes] = oisCreate('vernier','add', tseries, ...
    'sampleTimes',tsamples, ...
    'testParameters',vparams([1 2]),...
    'sceneParameters',sparams);
% offset.visualize;
% ieAddObject(offset.oiFixed); ieAddObject(offset.oiModulated); oiWindow;
% ieAddObject(scenesO{2}); sceneWindow;

% Aligned lines
aligned = oisCreate('vernier','add', tseries,...
    'sampleTimes',tsamples, ...
    'testParameters',vparams([1 3]),...
    'sceneParameters',sparams);
% aligned.visualize;

% Print out the offset in degrees of arc sec 
offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

end
