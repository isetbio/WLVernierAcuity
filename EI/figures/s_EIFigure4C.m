%% Figure 4  Bar length
%
% Where the parameters came from
%  load('spatialBarLength-20170109T003935.mat','params');
%  fname = fullfile(wlvRootPath,'EI','figures','barLength','barLength-Params');
%  save(fname,'params');

%% Load up the parameters
fname = fullfile(wlvRootPath,'EI','figures','barLength','barLength-Params');
load(fname,'params');

[~, ~,scenes] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

% Get these now just to check/confirm that we are OK with parameter selections
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Run for all the bar lengths

barOffset  = [0 1 2 3 4];           % Pixels on the display
barLengths = [70 90 120 240 350];   % Bar length is the top and bottom

if isempty(gcp), parpool('local'); end

PC = zeros(length(barOffset),length(barLengths));
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

%%
fname = fullfile(wlvRootPath,'EI','figures','barLength','Figure4-barLength.mat');
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset', 'barLengths');

%%
fname = fullfile(wlvRootPath,'EI','figures','barLength','Figure4-barLength.mat');
barOffsetSec = barOffset*secPerPixel;
load(fname);
vcNewGraphWin;
plot(barOffsetSec,PC,'-o','LineWidth',2);
grid on
title('Bar length'); xlabel('Bar offset (arcsec)'); ylabel('Percent correct classification');

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure4-BarLengthCurves.png');
title(''); xlabel(''); ylabel('');
saveas(gcf,fname,'png')

%%
thresh = zeros(1,size(PC,2));
barSteps = barOffsetSec(1):barOffsetSec(end);
for ii=1:(size(PC,2))
    p = interp1(barOffsetSec,PC(:,ii),barSteps);
    [v,idx] = min(abs(p - 80));
    thresh(ii) = barSteps(idx);
end

barLengthMin = barLengths*minPerPixel/2;
vcNewGraphWin;
plot(barLengthMin,thresh,'-o','LineWidth',2);
xlabel('Bar length (arcmin)')
ylabel('Offset threshold (arcsec)');
grid on

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure4-BarLengthThresholds.png');
title(''); xlabel(''); ylabel('');
saveas(gcf,fname,'png')

