%% Figure 5 Westheimer calculation
%
%% Find the files and store the parameters
%
%   ddir = fullfile(wlvRootPath,'EI','figures','westheimer');
%   chdir(ddir); dfiles = dir('tremor*'); nFiles = length(dfiles);
%   ii = 6; load(dfiles(ii).name,'params');
%   params.vernier
%
% Where the parameters came from
%  load('tremor-20170129T123448.mat','params');
%  fname = fullfile(wlvRootPath,'EI','figures','westheimer','westheimer-Params');
%  save(fname,'params');

%% Load up the parameters

fname = fullfile(wlvRootPath,'EI','figures','westheimer','westheimer-Params');
load(fname,'params');
params.nTrials = 2000;

[~, offset,scenes,tseries] = vaStimuli(params);
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
% oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

% Get these now just to check/confirm that we are OK with parameter selections
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Set up model parameters

% Choose bar lengths to run out to about 15 min of arc
barLengths = params.vernier.sceneSz(1)*2*[0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5]/3; 
barOffset  = [0 1 2 3 4 5];           % Pixels on the display
fprintf('Max bar length %.2f min\n',minPerPixel*max(barLengths)/2);
fprintf('Mosaic size %.2f min\n',coneMosaicFOV*60);

% Confirm barlengths
barLengthsMin = barLengths*minPerPixel/2
barOffSetSec = barOffset*secPerPixel

%% Run for all the bar lengths
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
fname = fullfile(wlvRootPath,'EI','figures','westheimer','Figure5-Westheimer.mat');
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl', 'barOffset','scenes','barLengths');

%%
fname = fullfile(wlvRootPath,'EI','figures','westheimer','Figure5-Westheimer.mat');
load(fname);
vcNewGraphWin;
plot(barOffSetSec(:),PC,'LineWidth',2)
grid on
set(gca,'FontName','Georgia');
set(gca,'FontSize',14);
xlabel('Offset (arc sec)')
ylabel('Percent correct');
title('Westheimer approximation')

%%
xlabel(''); ylabel(''); title('');
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure5-WestheimerCurves.png');
saveas(gcf,fname,'png')

%%
thresh = zeros(1,size(PC,2));
barSteps = barOffsetSec(1):barOffsetSec(end);
for ii=1:size(PC,2)
    p = interp1(barOffSetSec,PC(:,ii),barSteps);
    [v,idx] = min(abs(p - 80));
    thresh(ii) = barSteps(idx);
end

fname = fullfile(wlvRootPath,'EI','figures','westheimer','WestheimerCR');
load(fname);
fname = fullfile(wlvRootPath,'EI','figures','westheimer','WestheimerSM');
load(fname);

vcNewGraphWin;
plot(barLengthsMin(:),thresh,'-o','LineWidth',2)
set(gca,'FontName','Georgia','FontSize',14);
hold on; plot(WestheimerCR(:,1),WestheimerCR(:,2),'ro','markersize',10,'markerface','r');
hold on; plot(WestheimerSM(:,1),WestheimerSM(:,2),'rs','markersize',10,'markerface','r');
set(gca,'xlim',[0 16]);set(gca,'ylim',[0 16])
grid on
xlabel('Bar length (arcmin)');
ylabel('Vernier offset (arcsec)');
%%
xlabel(''); ylabel(''); title('');
fname = fullfile(wlvRootPath,'EI','figures','PNG','Figure5-WestheimerThreshold.png');
saveas(gcf,fname,'png')

%%