%% Match the Westheimer bar length curve 
%
% It is likely that 

%%
disp('**** EI Westheimer')

nTrials = 1000;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  
sceneFOV = 0.5;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.12;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = 2*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%% More like the McKee-Westheimer condition

% Helps find the ones that are like this 
params.matchHuman = true;     % Black background, matching human

params.vernier.bgColor = 0;   % Bright bar on a zero background

params.timesd    = 200e-3;    %  Reduce effect of eye movements
params.em.emFlag = [1 1 0];   % Microsaccades are suppressed for hyperacuity
% params.em = emSet(params.em,'tremor amplitude',0.004);  % Reduce the amplitude

%%  Build the stimuli if you want to check stuff
%
[~, offset,scenes,tseries] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Initialize offsets and lengths

% barLength Min is 60 * barLength/2
barLengths = params.vernier.sceneSz(1)*2*[0.2, 0.3, 0.5, 1 1.5]/3; 
barOffset  = [0 1 2 3 4 5];           % Pixels on the display
PC = zeros(length(barOffset),length(barLengths));
fprintf('Max bar length %.2f min\n',minPerPixel*max(barLengths)/2);
fprintf('Mosaic size %.2f min\n',coneMosaicFOV*60);

% Confirm barlengths
barLengthsMin = barLengths*minPerPixel/2
barOffset*secPerPixel
coneMosaicFOV*60    % Half of this will cap bar length performance

%% Run for all the bar lengths
if isempty(gcp), parpool('local'); end

tic;
svmMdl = cell(1, length(barLengths));
parfor pp=1:length(barLengths)
    fprintf('Starting %d of %d ...\n',pp, length(barLengths));
    thisParam = params;
    thisParam.vernier.barLength = barLengths(pp); 
    [P,thisMdl] = vaAbsorptions(barOffset,thisParam);
    svmMdl{pp} = thisMdl;
    PC(:,pp) = P(:);
    fprintf('Finished %d\n',pp);
end
toc

%% Plot

% Legend
lStrings = cell(1,length(barLengths));
for pp=1:length(barLengths)
    lStrings{pp} = sprintf('%g min',barLengthsMin(pp));
end

barOffsetSec = secPerPixel*barOffset;

h = vcNewGraphWin;
plot(barOffsetSec,PC,'-o');
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%% Save
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures','westheimer',['tremor-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset', 'barLengths','scenes');

%%

