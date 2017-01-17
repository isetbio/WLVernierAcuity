function [stimUniform, stimHarmonic, scenes, tseries] = csfStimuli(varargin)
% Create the pair of harmonic stimuli (uniform and harmonic)
%
% There is usually one input argument that is a struct with these
% parameters
%
%  harmonic  - Parameters for the harmonic stimuli
%  tsamples  - Time samples (sec)
%  timesd    - Time standard deviation
%  display   - Display struct, default displayCreate('LCD-Apple')
%  sceneFOV  - Degrees, default 0.35
%  distance  - Meters, viewing distance to display, default is 0.3 m
%
% When the LCD-Apple display is 210, 210 pixels and this is 0.35, then 1
% pixel is 6 arc sec
%
% BW, ISETBIO Team, 2016

%%
params = varargin{1};

% Stored imageBasis filename with the same parameters as this
[~,fname] = csfFname(params);

if exist(fname,'file')
    disp('Loading stimulus from file - parameters match')
    try
        load(fname,'stimUniform','stimHarmonic','scenes','tseries');
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

p.addParameter('harmonic',harmonicP,@isstruct);

p.addParameter('tsamples',(-50:100),@isvector);  % Time samples (ms)
p.addParameter('timesd',20,@isscalar);           % Time standard deviation

% When the LCD-Apple display is 210, 210 pixels and sceneFOV is 0.35, then 1
% pixel is 6 arc sec
p.addParameter('sceneFOV',0.35,@isscalar);  % Degrees. 
p.addParameter('distance',0.3,@isscalar);   % Meters
p.parse(varargin{:});

harmonic   = p.Results.harmonic;

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
for ii=2:-1:1
    hparams(ii) = harmonic;
end

% The first parameter has zero contrast, the second parameter has some contrast
hparams(1).contrast = 0;

% Offset lines
P.sampleTimes     = tsamples;
P.testParameters  = hparams;
P.sceneParameters = sparams;
P.oi = params.oi;
[stimHarmonic, scenes] = oisCreate('harmonic','blend', tseries, P);

% vcNewGraphWin; stimHarmonic.visualize;
% ieAddObject(offset.oiFixed); ieAddObject(offset.oiModulated); oiWindow;
% ieAddObject(scenesO{2}); sceneWindow;

% Zero contrast for both the modulated and fixed.
hparams(2).contrast = 0;
P.testParameters  = hparams;
stimUniform = oisCreate('harmonic','blend', tseries, P);
% stimUniform.visualize;

%%
save(fname,'stimUniform','stimHarmonic','scenes','tseries','P');

% Print out the offset in degrees of arc sec 
% offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
% fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

end
