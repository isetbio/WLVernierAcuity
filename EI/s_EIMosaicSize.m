%% Impact of cone mosaic size
%
% This complements the bar length analysis.  It should tell a similar story.
%

% Show the dependence on spatial size of the cone mosaic for the computational
% observer.
nTrials = 1000;
nBasis  = 40;

% Integration time 
tStep   = 10;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Original scene
sceneFOV = 0.6;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 2*(sceneFOV/0.35);   % If you do not multiply by a scalar, offset is 6 arc sec

s_EIParameters;

% Make the bar length a little less than the scene size
params.vernier.barLength = params.vernier.sceneSz(1)-1;

%%  Build the stimuli if you want to check stuff
% 
[~, offset,scenes,tseries] = vaStimuli(params);
% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;

degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;


%% Summarize

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (params.vernier.barLength*minPerPixel),...
    (params.vernier.barLength*degPerPixel),...
    (params.vernier.barWidth*minPerPixel));
fprintf('Bar offset per pixel is %.1f sec\n',secPerPixel);
barOffset = [0 1 2 3 4];
fprintf('Offsets in seconds %.1f\n',barOffset*secPerPixel);

cmFOV = [0.15 0.25 0.5]; % [0.2 0.4 0.80];    % Degress of visual angle
PC = zeros(length(barOffset),length(cmFOV));

%% Compute classification accuracy
tic;
c = gcp; if isempty(c), parpool('local'); end
parfor pp = 1:length(cmFOV)
    fprintf('Starting %d of %d ...\n',pp,length(cmFOV));
    thisParam = params;
    thisParam.cmFOV = cmFOV(pp);
    P = vaAbsorptions(barOffset, thisParam);
    PC(:, pp) = P(:);
    fprintf('Finished %d\n',pp);
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
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'barOffset','secPerPixel','cmFOV');
%%
