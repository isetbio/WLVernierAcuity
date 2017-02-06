%% Figure 4  Luminance panel
%
% Where the parameters came from
%  load('example-20170121T140456.mat','params');
%  fname = fullfile(wlvRootPath,'EI','figures','displayPeakLum','displayPeakLum-Params');
%  save(fname,'params');

%% Load up the parameters

fname = fullfile(wlvRootPath,'EI','figures','displayPeakLum','displayPeakLum-Params');
load(fname,'params');

[~, ~,scenes] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

% Get these now just to check/confirm that we are OK with parameter selections
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Set up luminance
sceneLumFactor = logspace(-2,0.3,7);
display = params.vernier.display;
spd  = displayGet(display,'spd');
barOffset  = 0:2:8;                  % Pixels on the display
PC = zeros(length(barOffset),length(sceneLumFactor));

parfor pp = 1:length(sceneLumFactor)
    thisParams = params;
    thisParams.vernier.display = ...
        displaySet(thisParams.vernier.display,'spd',spd*sceneLumFactor(pp));
    pLum = displayGet(thisParams.vernier.display,'peak luminance');
    fprintf('%d:  Peak luminance %.2f cd/m2\n',pp,pLum);
    [PC(:,pp),svmMdl] = vaAbsorptions(barOffset,thisParams);
    fprintf('Finished %d\n',pp);
end

%%
fname = fullfile(wlvRootPath,'EI','figures','displayPeakLum','Figure4-Luminance.mat');
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset','scenes','sceneLumFactor');

%%
fname = fullfile(wlvRootPath,'EI','figures','displayPeakLum','Figure4-Luminance.mat');
load(fname);

barOffsetSec = barOffset*secPerPixel;

vcNewGraphWin;
plot(barOffsetSec,PC,'LineWidth',2)
set(gca,'FontName','Georgia','FontSize',14);

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure4-LuminanceCurves.png');
saveas(gcf,fname,'png')

%%
thresh = zeros(1,size(PC,2));
barSteps = barOffsetSec(1):barOffsetSec(end);
for ii=1:(size(PC,2))
    p = interp1(barOffsetSec,PC(:,ii),barSteps);
    [v,idx] = min(abs(p - 80));
    thresh(ii) = barSteps(idx);
end

%
vcNewGraphWin;
set(gca,'FontName','Georgia','FontSize',14)
whiteP = displayGet(display,'white point');
semilogx(whiteP(2)*sceneLumFactor,thresh,'-o','LineWidth',2);
xlabel('Log_{10} luminance (cd/m^2)')
ylabel('Offset threshold (arc sec)')
grid on

%%
title(''); xlabel(''); ylabel('')
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure4-LuminanceThresholds.png');
saveas(gcf,fname,'png')

%%