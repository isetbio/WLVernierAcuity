%% Impact of cone mosaic size
%
% This complements the bar length analysis.  It should tell a similar story.
%

% Show the dependence on spatial size of the cone mosaic for the computational
% observer.
nTrials = 600;
nBasis  = 40;

% Integration time 
tStep   = 10;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Original scene
sceneFOV = 0.8;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = (sceneFOV/0.35); 

s_EIParameters;

% Make the bar length a little less than the scene size
params.vernier.barLength = round(params.vernier.sceneSz(1))-1;

%% Summarize

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
secPerPixel = (6 / sc);
minPerPixel = (6 / sc) / 60;
degPerPixel = minPerPixel/60;
barLength = params.vernier.barLength*minPerPixel;
barWidth   = params.vernier.barWidth*minPerPixel;
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
%%
fprintf('Bar offset per pixel is %.1f sec\n',secPerPixel);
barOffset = [0 1 2 3 4];

cmFOV = [0.15 0.3]; % [0.2 0.4 0.80];    % Degress of visual angle
PC = zeros(length(barOffset),length(cmFOV));

%% Compute classification accuracy
tic;
c = gcp; if isempty(c), parpool('local'); end
parfor pp = 1:length(cmFOV)
    thisParam = params;
    thisParam.cmFOV = cmFOV(pp);
    P = vaAbsorptions(barOffset, thisParam);
    PC(:, pp) = P(:);
end
toc
% mesh(PC)

%% Make summary graph
% Legend
lStrings = cell(1,length(cmFOV));
for pp=1:length(cmFOV)
    lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
end

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC,'o-');
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['mosaicSize-',str,'.mat']);
save(fname, 'PC','params', 'barOffset','secPerPixel','cmFOV');

%%
% ddir = fullfile(wlvRootPath,'EI','figures');
% dfiles = dir(fullfile(ddir,'mosaicSize*'));

% To return the scene and such you could read params and run
% d = load(dfiles(1).name,'params');
% [aligned, offset, scenes, tseries] = vaStimuli(d.params);

%%
% vcNewGraphWin;
% for pp=1:length(dfiles)
%     d = load(dfiles(pp).name);
%     secPerPixel = sceneGet(scenes{1},'degrees per sample')*3600;
%     offsets = d.barOffset*secPerPixel;
%     plot(offsets,PC,'o-'); hold on
% end
% 
% % Legend
% lStrings = cell(1,length(cmFOV));
% for pp=1:length(cmFOV)
%     lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
% end
% xlabel('Offset arc sec'); ylabel('Percent correct')
% grid on; l = legend(lStrings);
% set(l,'FontSize',12)

%%
% h = vcNewGraphWin;
% cnt = 0;
% for ii=1:length(dfiles)
%     load(dfiles(ii).name);
%     ii, size(PC)
%     PC
%     if ii == 1
%         PCall = PC; cnt = 1;
%     else
%         if size(PC) == size(PCall)
%             PCall = PCall + PC;
%             cnt = cnt + 1;
%         end
%     end
% end
% PCall = PCall/cnt;
% plot(secPerPixel*barOffset,PCall,'-o');
% xlabel('Offset arc sec'); ylabel('Percent correct')
% grid on; l = legend(lStrings);
% set(l,'FontSize',12)

%%
