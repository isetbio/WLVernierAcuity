%% Figure 4 - Defocus panel
%
% Where the parameters came from
%  load('vaDefocus-20170123T223201.mat','params');
%  fname = fullfile(wlvRootPath,'EI','figures','vaDefocus','vaDefocus-Params');
%  save(fname,'params');

%%
fname = fullfile(wlvRootPath,'EI','figures','vaDefocus','vaDefocus-Params');
load(fname,'params');

[~, ~,scenes] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

% Get these now just to check/confirm that we are OK with parameter selections
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%%
barOffset = [0 1 2 3 4];
defocus = [0 0.5 1 1.5 2];   % Microns
PC = zeros(length(barOffset),length(defocus));
svmMdl = cell(1, length(defocus));
parfor pp = 1:length(defocus)
    fprintf('Starting %d ...\n',pp);

    thisParams = params;
    thisParams.defocus = defocus(pp);
    
    % Make the oi
    thisParams.oi = oiDefocus(defocus(pp));
    
    % Compute classification accuracy
    [P, thisMdl] = vaAbsorptions(barOffset, thisParams);
    PC(:,pp) = P;
    svmMdl{pp} = thisMdl;
    fprintf('Finished defocus level %d\n',defocus(pp));
end

%%
fname = fullfile(wlvRootPath,'EI','figures','vaDefocus','Figure4-Defocus.mat');
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl','barOffset','secPerPixel','defocus');

%%
fname = fullfile(wlvRootPath,'EI','figures','vaDefocus','Figure4-Defocus.mat');
load(fname)
barOffsetSec = barOffset*secPerPixel;

%%
vcNewGraphWin;
set(gca,'FontName','Georgia','FontSize',14)
plot(barOffsetSec,PC,'LineWidth',2);
xlabel('Vernier offset (arcsec)');
ylabel('Percent correct classification');
grid on

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure4-DefocusCurves.png');
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

defocusDiopters = wvfDefocusMicronsToDiopters(defocus,3);
vcNewGraphWin;
set(gca,'FontName','Georgia','FontSize',14)
plot(defocusDiopters,thresh,'-o','LineWidth',2);
xlabel('Defocus (diopters)');
ylabel('Vernier offset (arcsec)');
grid on

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure4-DefocusThreshold.png');
title(''); xlabel(''); ylabel('');
saveas(gcf,fname,'png')

%%