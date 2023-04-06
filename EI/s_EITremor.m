%% Impact of tremor
%
% Show the dependence on the tremor parameter for the computational observer.
%
% We also do this matching the Westheimer and McKee parameters a bit more
% closely.

%%
disp('**** EI Tremor')

nTrials = 1000;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 0.4;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = 4*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%% More like the McKee-Westheimer condition

% Helps find the ones that are like this 
params.matchHuman = true;     % Black background, matching human

params.vernier.bgColor = 0;   % Bright bar on a zero background

% Spatial summation is about 6 arcmin.  
% The field of view is 0.25*60 min = 15 min.  
% We want the bar to be about 6 arcmin.
bLengthScale = 0.5;         
params.vernier.barLength = round(params.vernier.sceneSz(1)*bLengthScale);  % 
params.timesd    = 200e-3;    %  Reduce effect of eye movements
params.em.emFlag = [1 1 0];   % Microsaccades are suppressed for hyperacuity

% Reduce the tremor amplitude
% params.em = emSet(params.em,'tremor interval', 0.0120/2);

%%  Build the stimuli if you want to check stuff
%
[~, offset,scenes,tseries] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Summarize spatial parameters
% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
barLength  = params.vernier.barLength*minPerPixel;
barWidth   = params.vernier.barWidth*minPerPixel;
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%% Initialize offsets and lengths

barOffset  = [0 1 2 3 4];           % Pixels on the display
tremorAmplitude = [0, 0.0024, 0.0037, 0.0049, 0.0073]; 

PC = zeros(length(barOffset),length(barLengths));
fprintf('Max bar length %.2f\n',degPerPixel*max(barLengths));
fprintf('Mosaic size %.2f\n',coneMosaicFOV);

%% Run for all the bar lengths
if isempty(gcp), parpool('local'); end

tic;
svmMdl = cell(1, length(tremorAmplitude));
parfor pp=1:length(tremorAmplitude)
    fprintf('Starting %d of %d ...\n',pp, length(tremorAmplitude));
    thisParam = params;
    thisParam.em = emSet(thisParam.em,'tremor amplitude',tremorAmplitude(pp)); 
    [P,thisMdl] = vaAbsorptions(barOffset,thisParam);
    svmMdl{pp} = thisMdl;
    PC(:,pp) = P(:);
    fprintf('Finished %d\n',pp);
end
toc

%% Plot

% Legend
lStrings = cell(1,length(tremorAmplitude));
for pp=1:length(tremorAmplitude)
    lStrings{pp} = sprintf('%g amp',tremorAmplitude(pp));
end

title(sprintf('Scene FOV %.1f',sceneGet(scenes{1},'fov')),'FontSize',14)

barOffsetSec = secPerPixel*barOffset;

h = vcNewGraphWin;
plot(barOffsetSec,PC,'-o');
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%% Save
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures','tremor',['tremor-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset', 'tremorAmplitude','scenes');

%%

