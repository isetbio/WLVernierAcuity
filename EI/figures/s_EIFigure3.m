%%  Figure 3 - Eye Movements
%
% This is how I selected the initial parameters.
%
%  load('FineEyeMovements-20170118T143607.mat','params');
%  params.cmFOV = 0.35;
%  save(fname,'params');

%% Load default parameters
fname = fullfile(wlvRootPath,'EI','figures','EyeMovements','Figure3-Parameters.mat');
load(fname,'params');

%%  This is the long bar calculation

[~,~,scenes] = vaStimuli(params);
secPerPixel = sceneGet(scenes{1},'degrees per sample')*3600;
minPerPixel = secPerPixel/60;
fprintf('Bar length %.2f\n',params.vernier.barLength*minPerPixel/2);

%%
barOffset = [0 1 2 4 6];
PC = zeros(numel(barOffset),5);

c  = gcp; if isempty(c), parpool('local'); end
emTypes = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 1 1 1]';

svmMdl = cell(1,size(emTypes,2));
parfor pp=1:size(emTypes,2)
    fprintf('Starting %d of %d ...\n',pp,size(emTypes,2));
    thisParam = params;
    thisParam.em.emFlag = emTypes(:,pp);
    [P, thisMdl] = vaAbsorptions(barOffset,thisParam);
    PC(:,pp) = P(:);
    svmMdl{pp} = thisMdl;
    fprintf('Finished %d\n',pp);
end

%% Save the file for plotting
longFile = fullfile(wlvRootPath,'EI','figures','EyeMovements','Figure3-Long.mat');
fprintf('Saving %s\n',longFile);
save(longFile, 'PC', 'params', 'svmMdl', 'emTypes','barOffset');

%%  Show a graph to see that things are generally OK
load(longFile)
barOffsetSec = barOffset*secPerPixel;
vcNewGraphWin;
plot(barOffsetSec,PC,'LineWidth',2);
grid on
set(gca,'FontName','Georgia','FontSize',14)
set(gca,'ylim',[40 110]);
xlabel('Vernier offset (arcsec)');
ylabel('Percent correct classification');

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure3-Long.png');
title(''); xlabel(''); ylabel('');
saveas(gcf,fname,'png')

%%  This is the short bar calculation - Matches except the bar length.

% Change the bar length
params.vernier.barLength = params.vernier.barLength*(3/10);

[~,~,scenes] = vaStimuli(params);
secPerPixel = sceneGet(scenes{1},'degrees per sample')*3600;
minPerPixel = secPerPixel/60;
fprintf('Bar length %.2f\n',params.vernier.barLength*minPerPixel/2);

svmMdl = cell(1,size(emTypes,2));
parfor pp=1:size(emTypes,2)
    fprintf('Starting %d of %d ...\n',pp,size(emTypes,2));
    thisParam = params;
    thisParam.em.emFlag = emTypes(:,pp);
    [P, thisMdl] = vaAbsorptions(barOffset,thisParam);
    PC(:,pp) = P(:);
    svmMdl{pp} = thisMdl;
    fprintf('Finished %d\n',pp);
end

%% Save the file for the short barlength panel

shortFile = fullfile(wlvRootPath,'EI','figures','EyeMovements','Figure3-Short.mat');
fprintf('Saving %s\n',shortFile);
save(shortFile, 'PC', 'params', 'svmMdl', 'emTypes','barOffset');

%%  Show a graph to see that things are generally OK
load(shortFile);
vcNewGraphWin;
plot(barOffsetSec,PC,'LineWidth',2);
grid on
set(gca,'FontName','Georgia','FontSize',14)
set(gca,'ylim',[40 110]);
xlabel('Vernier offset (arcsec)');
ylabel('Percent correct classification');

%%
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure3-Short.png');
title(''); xlabel(''); ylabel('');
saveas(gcf,fname,'png')

