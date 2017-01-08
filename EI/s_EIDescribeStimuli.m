%% Description of the stimuli for plotting
%
% Show the dependence on bar length for the computational observer. Use the
% match on bar length with behavior as an indicator of the spatial summation
% region of the human eye

nTrials = 600;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 0.35;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.2;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = (sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%% Summarize spatial parameters
% Each pixel size is 6 arc sec per pixel when sc = 1.  Finer resolution when sc
% is higher.
secPerPixel = (6 / sc);
minPerPixel = (6 / sc) / 60;
degPerPixel = minPerPixel/60;
barLength = params.vernier.barLength*minPerPixel;
barWidth   = params.vernier.barWidth*minPerPixel;
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%%  Build the stimuli

[~, offset,scenes,tseries] = vaStimuli(params);

ieAddObject(scenes{2}); sceneWindow;
ieAddObject(offset.oiModulated); oiWindow;

oiGet(offset.oiModulated,'angular resolution')*3600
