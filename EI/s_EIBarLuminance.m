%% Impact of background luminance


disp('**** EI Bar Luminance')

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
sc = (sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);


%%  Build the stimuli if you want to check stuff
%
[~, offset,scenes,tseries] = vaStimuli(params);
scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;


%% Initialize offsets and lengths

barOffset  = [0 1 2 3 4];                  % Pixels on the display
bgColors   = [0 .25 .5 .8 1];
params.vernier.barColor = ones(1,3)*1;   % Low intensity bar

PC = zeros(length(barOffset),length(bgColors));
fprintf('Background level %.2f\n',bgColors);
fprintf('Mosaic size %.2f deg\n',coneMosaicFOV);

%% Run for all the background colors
if isempty(gcp), parpool('local'); end

% This reduces the contrast.
tic;
svmMdl = cell(1, length(bgColors));
for pp=1:length(bgColors)
    fprintf('Starting %d of %d ...\n',pp, length(bgColors));
    thisParam = params;
    thisParam.vernier.bgColor = ones(1,3)*bgColors(pp);
    [P,thisMdl] = vaAbsorptions(barOffset,thisParam);
    svmMdl{pp} = thisMdl;
    PC(:,pp) = P(:);
    fprintf('Finished %d\n',pp);
end
toc

%% Plot

% Legend
lStrings = cell(1,length(bgColors));
for pp=1:length(bgColors)
    lStrings{pp} = sprintf('%.2f level',bgColors(pp));
end

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC,'-o');
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%% Save
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['luminance-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset', 'bgColors','scenes');

%%

