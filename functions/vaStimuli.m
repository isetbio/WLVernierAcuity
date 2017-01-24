function [aligned, offset, scenes, tseries, fname] = vaStimuli(varargin)
% VASTIMULI - Create the pair of vernier stimuli (aligned and offset)
%
%  [aligned, offset, scenes, tseries, fname] = vaStimuli(varargin)
%
% There is one input argument that is a struct with these parameters
%
%  vernier   - Parameters for the vernier stimuli; Default is vernierP
%  tsamples  - Time samples (sec)
%  timesd    - Time standard deviation
%  display   - Display struct, default displayCreate('LCD-Apple')
%  sceneFOV  - Degrees, default 0.35
%  distance  - Meters, viewing distance to display, default is 0.3 m
%
%  When the LCD-Apple display is [210, 210] pixels and the sceneFOV is 0.35,
%  then 1 pixel is 6 arc sec.
%
%  See also s_EIParameters, s_EIMosaicSize
%
% The default display is an Apple LCD monitor.
%
% The offset of 1 pixel for a test with 210,210 pixels and a field of view
% of 0.35 deg is 6 sec.  That is a useful value. If you scale the number of
% pixels and the field of view together, then the offset remains the same.
%
% Examples:
%    clear p; p.barOffset = 3; p.barWidth  = 3;
%    p.tsamples = [-60:70]; p.timesd = 20;  
%    [~,offset,~,tseries, fname] = vaStimuli(p);
%    offset.visualize;
%    vcNewGraphWin; plot(tseries)
%
% Notes
%   Westheimer and McKee, 1977
%   Luminance 2 log millilamberts
%   Line width 1 minute of arc.
%
% BW, ISETBIO Team, 2016
%
% See also imageVernier, sceneVernier, sceneCreate, s_vaAbsorptions.m

%%
params = varargin{1};

% Stored imageBasis filename with the same parameters as this
[~,fname] = vaFname(params);
if exist(fname,'file'), delete(fname); end  % Forcing create
if exist(fname,'file')
    disp('Loading stimulus from file - parameters match')
    try
        load(fname,'aligned','offset','scenes','tseries');
        return;
    catch
        disp('File found, but not the variables.  Creating.')
    end
else 
    disp('Creating and saving stimulus file - no match found')
end


%%
p = inputParser;

p.KeepUnmatched = true;   % Sometimes we overload p for SVM and cMosaic

p.addParameter('vernier',vernierP,@isstruct);

p.addParameter('tsamples',(-50:100),@isvector);  % Time samples (ms)
p.addParameter('timesd',20,@isscalar);           % Time standard deviation

% When the LCD-Apple display is 210, 210 pixels and sceneFOV is 0.35, then 1
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
P.sampleTimes = tsamples;
P.testParameters = vparams([1 2]);
P.sceneParameters = sparams;
if isfield(params,'oi'), P.oi = params.oi; end
[offset, scenes] = oisCreate('vernier','add', tseries, P);
% offset.visualize;
% ieAddObject(offset.oiFixed); ieAddObject(offset.oiModulated); oiWindow;
% ieAddObject(scenesO{2}); sceneWindow;

% Aligned lines
P.testParameters = vparams([1 3]);
aligned = oisCreate('vernier','add', tseries, P);
% aligned.visualize;

%%
save(fname,'aligned','offset','scenes','tseries','P');

% Print out the offset in degrees of arc sec 
% offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
% fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

end
