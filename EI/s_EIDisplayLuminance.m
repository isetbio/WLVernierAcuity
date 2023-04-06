%% Display luminance variation
%
% Parameters are chosen to match those in the Westheimer case.
% 
%%
disp('**** EI Display Luminance')

nTrials = 1000;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 0.3;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.15;   % First bunch run with 0.12 like W/M experiment

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = 2*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

% Helps find the ones that are like this 
params.matchHuman = true;     % Black background, matching human

params.vernier.bgColor = 0;     % Bright bar on a zero background
params.vernier.barLength = 200; % Matches Westheimer (6 min)
params.timesd    = 200e-3;      %  Reduce effect of eye movements
params.em.emFlag = [1 1 0];     % Microsaccades are suppressed for hyperacuity

%%  Build the stimuli if you want to check stuff
%
[~, offset,scenes,tseries] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ill = oiCalculateIlluminance(oi); oi = oiSet(oi,'illuminance',ill); ieAddObject(oi); oiWindow;

degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

fprintf('Bar length %.2f min\n',minPerPixel*params.vernier.barLength/2);
fprintf('Mosaic size %.2f min\n',coneMosaicFOV*60);

%% Darken the scene or not
sceneLumFactor = logspace(-2,1,7);
spd  = displayGet(display,'spd');

barOffset  = [0:2:8];                  % Pixels on the display
PC = zeros(length(barOffset),length(sceneLumFactor));
parfor pp = 1:length(sceneLumFactor)
    thisParams = params;
    display = thisParams.vernier.display;
    thisParams.vernier.display = displaySet(display,'spd',spd*sceneLumFactor(pp));
    pLum = displayGet(thisParams.vernier.display,'peak luminance');
    fprintf('%d:  Peak luminance %.2f cd/m2\n',pp,pLum);
    [PC(:,pp),svmMdl] = vaAbsorptions(barOffset,thisParams);
    fprintf('Finished %d\n',pp);
end

%% Plot

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC,'-o','LineWidth',2);
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on;
set(gca,'FontSIze',16)

%%
thresh = zeros(1,length(sceneLumFactor));
barOffsetInterp = min(barOffset(:)):0.2:max(barOffset(:));
for ii=1:length(sceneLumFactor)
    p = interp1(barOffset,PC(:,ii),barOffsetInterp);
    [v,idx] = min(abs(p - 75));
    thresh(ii) = barOffsetInterp(idx);
end

%%
vcNewGraphWin;
tmp = displayGet(display,'white point');
Y = tmp(2);
semilogx(sceneLumFactor*Y, thresh,'-o','LineWidth',2);
grid on
xlabel('Log white point luminance (cd/m2)');
ylabel('Offset threshold (arcsec)');
set(gca,'FontSize',16);

xlabel(''); ylabel('');


%% Save
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures','displayPeakLum',['luminance-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset','scenes','sceneLumFactor');

%%

