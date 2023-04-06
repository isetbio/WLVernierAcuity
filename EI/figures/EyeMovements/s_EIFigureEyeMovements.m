%% Figure 3 - Eye movement impact
%
%  To create the data file we only need to run
% 
% shortBarFile = 'FineEyeMovements-20170126T092551.mat';
% load(shortBarFile,'params');
%
% You can see all of the parameters in the params file, and particularly
% params.vernier
%
% This code  executes to create the classification accuracy
%
% barOffset = [0 1 2 4 6];
% PC = zeros(numel(barOffset),5);
% c  = gcp; if isempty(c), parpool('local'); end
% emTypes = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ; 1 1 1]';
% 
% svmMdl = cell(1,size(emTypes,2));
% parfor pp=1:size(emTypes,2)
%     fprintf('Starting %d of %d ...\n',pp,size(emTypes,2));
%     thisParam = params;
%     thisParam.em.emFlag = emTypes(:,pp);
%     [P, thisMdl] = vaAbsorptions(barOffset,thisParam);
%     PC(:,pp) = P(:);
%     svmMdl = thisMdl;
%     fprintf('Finished %d\n',pp);
% end

%%  Eye movement results are stored here

shortBarFile = fullfile(wlvRootPath,'EI','figures','EyeMovements','Figure3-Short.mat');
%% Panel A (short bar)

load(shortBarFile,'params');  % PC, barOffset, params, scenes
[~,~,scenes] = vaStimuli(params);

secPerPixel = sceneGet(scenes{1},'degrees per sample')*3600;
minPerPixel = secPerPixel/60;
barOffsetSec = barOffset*secPerPixel;

vcNewGraphWin;
% We skip the microsaccade data
plot(barOffsetSec,PC,'-o','LineWidth',2);
lStrings = cell({'None','Tremor only','Drift only','Tremor and Drift'});
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on;
l = legend(lStrings,'Location','SouthEast');
set(l,'FontSize',14)
set(gca,'FontSize',14);
set(gca,'ylim',[40 110]);
set(gca,'xlim',[0 12]);
title(sprintf('Bar length %.1f arcmin',params.vernier.barLength*minPerPixel/2));

%% Panel B (long bar)
longBarFile = fullfile(wlvRootPath,'EI','figures','EyeMovements','Figure3-Long.mat');
[aligned,offset,scenes] = vaStimuli(params);

secPerPixel = sceneGet(scenes{1},'degrees per sample')*3600;
minPerPixel = secPerPixel/60;
barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;

vcNewGraphWin;
plot(barOffsetSec,PC,'-o','LineWidth',2);
lStrings = cell({'None','Tremor only','Drift only','Tremor and Drift'});
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on;
l = legend(lStrings,'Location','SouthEast');
set(l,'FontSize',14)
set(gca,'FontSize',14);
set(gca,'ylim',[40 110]);
title(sprintf('Bar length %.1f arcmin',params.vernier.barLength*minPerPixel/2));


%% Extra

% Run through some of files for group averages.  These are similar enough to the
% single example that I don't bother averaging.
PCall = zeros(5,5);
for ii = 1:3   % Good one is 8,9 Bad 3
    load(dfiles(ii+7).name);  % PC, barOffset, params, scenes
    barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;
    PCall = PCall + PC;
end
PCall = PCall/3;

%
vcNewGraphWin;
% plot(barOffsetSec,PCall,'-o','LineWidth',2);
% lStrings = cell({'None','Tremor only','Drift only','Microsaccade only','All'});
plot(barOffsetSec,PCall(:,[1:3,5]),'-o','LineWidth',2);
lStrings = cell({'None','Tremor only','Drift only','Tremor and Drift'});

xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; 
l = legend(lStrings,'Location','SouthEast');
set(l,'FontSize',14)
set(gca,'FontSize',14);
set(gca,'ylim',[40 110]);
title(sprintf('Bar length is %.2f arcmin\n',barLengthMin));

%% Scroll throgh each of the individual files

ddir = fullfile(wlvRootPath,'EI','figures','EyeMovements');
chdir(ddir);
dfiles = dir('FineEye*');
nFiles = length(dfiles);
for ii=1:nFiles
    
    load(dfiles(ii).name);  % PC, barOffset, params, scenes
    barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;
    PCall = PC;
    
    vcNewGraphWin;
    plot(barOffsetSec,PCall,'-o','LineWidth',2);
    lStrings = cell({'None','Tremor only','Drift only','Microsaccade only','All'});
    xlabel('Offset arc sec'); ylabel('Percent correct')
    grid on;
    l = legend(lStrings,'Location','SouthEast');
    set(l,'FontSize',14)
    set(gca,'FontSize',14);
    set(gca,'ylim',[40 110]);
    title(sprintf('%s',dfiles(ii).name)); drawnow;
    pause(1);
end

%%

