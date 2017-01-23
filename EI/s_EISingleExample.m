%% Run a single example for figures and little tests
% 
%%
disp('**** EI Single Example')

nTrials = 1000;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 0.4;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.35;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = 2*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);


%%  Build the stimuli if you want to check stuff
%
[~, offset,scenes,tseries] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ill = oiCalculateIlluminance(oi); oi = oiSet(oi,'illuminance',ill); ieAddObject(oi); oiWindow;

degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%%
%   barColor: [1 1 1]
%     barLength: 229
%      barWidth: 11
%       bgColor: 0.5000
%       display: [1x1 struct]
%           gap: 0
%          name: 'unknown'
%        offset: 1
%       pattern: []
%       sceneSz: [240 240]

%% Darken the scene or not
sceneLumFactor = 0.7;
display = params.vernier.display;
spd  = displayGet(display,'spd');
params.vernier.display = displaySet(display,'spd',spd*sceneLumFactor);

%%
pLum = displayGet(params.vernier.display,'peak luminance');
fprintf('Peak luminance %.2f\n',pLum);

%% Run it

barOffset  = [0 1 2 3 4];                  % Pixels on the display
tic;
[PC,svmMdl] = vaAbsorptions(barOffset,params);
fprintf('Finished %d\n',pp);
toc

%% Plot

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC,'-o');
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on;

%% Save
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['singleExample-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset','scenes');

%%

