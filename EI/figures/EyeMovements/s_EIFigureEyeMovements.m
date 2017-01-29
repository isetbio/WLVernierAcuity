%% Eye movement impact
%
% TODO:
%   Make video of eye movements and absorptions?
%   
% 
%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/EyeMovements';
ddir = fullfile(wlvRootPath,'EI','figures','EyeMovements');
chdir(ddir);
dfiles = dir('FineEye*');
% nFiles = length(dfiles);

%%  The bar length is about half the scene, which is 0.35 deg.
% scene = sceneAdd(scenes{1},scenes{2}); ieAddObject(scene); sceneWindow;
barLengthMin = (params.vernier.barLength*sceneGet(scenes{1},'degrees per sample')*60)/2;
fprintf('Bar length is %.2f arcmin\n',barLengthMin);

%%
PCall = zeros(5,5);
for ii = 1:3   % Good one is 8,9 Bad 3
    load(dfiles(ii+7).name);  % PC, barOffset, params, scenes
    barOffsetSec = barOffset*sceneGet(foo.scenes{2},'degrees per sample')*3600;
    PCall = PCall + PC;
end
PCall = PCall/3;

%%
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

%% Maybe for the talk?
%
% vaImageBasis(params);

%% Show each of the individual files

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

%% Show the short bar example

shortBarFile = 'FineEyeMovements-20170126T092551.mat';
load(shortBarFile);  % PC, barOffset, params, scenes

secPerPixel = sceneGet(scenes{2},'degrees per sample')*3600;
minPerPixel = secPerPixel/60;
barOffsetSec = barOffset*secPerPixel;
PCall = PC;

vcNewGraphWin;
plot(barOffsetSec,PCall(:,[1:3,5]),'-o','LineWidth',2);
lStrings = cell({'None','Tremor only','Drift only','Tremor and Drift'});
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on;
l = legend(lStrings,'Location','SouthEast');
set(l,'FontSize',14)
set(gca,'FontSize',14);
set(gca,'ylim',[40 110]);
set(gca,'xlim',[0 12]);
title(sprintf('Bar length %.1f arcmin',params.vernier.barLength*minPerPixel/2));

%%
longBarFile  = 'FineEyeMovements-20170118T143607.mat';
load(longBarFile);  % PC, barOffset, params, scenes
barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;
PCall = PC;

vcNewGraphWin;
plot(barOffsetSec,PCall(:,[1:3,5]),'-o','LineWidth',2);
lStrings = cell({'None','Tremor only','Drift only','Tremor and Drift'});
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on;
l = legend(lStrings,'Location','SouthEast');
set(l,'FontSize',14)
set(gca,'FontSize',14);
set(gca,'ylim',[40 110]);
title(sprintf('Bar length %.1f arcmin',params.vernier.barLength*minPerPixel/2));

%%

