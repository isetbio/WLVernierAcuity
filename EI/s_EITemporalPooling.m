%% Impact of temporal pooling
%
% BW, ISETBIO Team, 2017

%%
disp('**** EI Temporal Pooling')

% Initialize parameters for s_EIParameters call
nTrials = 1000;
nBasis  = 40;

% Integration time 
tStep   = 10;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Original scene
sceneFOV = 0.35;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = (sceneFOV/0.35);  % 6 arc sec steps.

maxOffset = 5;         % Out to 30 arc sec for vaPCA

s_EIParameters;

params.vernier.barWidth = 10;
params.tsamples  = (-500:tStep:500)*1e-3;   % In second M/W was 200 ms

%% Summnarize

[~, offset,scenes,tseries] = vaStimuli(params);

% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Summarize spatial parameters
% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
barLength = params.vernier.barLength*minPerPixel;
barWidth   = params.vernier.barWidth*minPerPixel;
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%%  Sweep out different durations

% Total of 1 sec duration
barOffset = 2;
sd = [50 100 200 400 600]*1e-3;
PC = zeros(1,length(sd));
svmMdl = cell(1, length(sd));
tic
parfor pp=1:length(sd)
    fprintf('Starting %d of %d ...\n',pp,length(sd));
    thisParams = params;
    thisParams.timesd  = sd(pp);                  % In seconds                 
    [P, thisMdl] = vaAbsorptions(barOffset,thisParams);
    PC(:,pp) = P(:);
    svmMdl{pp} = thisMdl;
    fprintf('Finished %d\n',pp);
end
toc

%%
vcNewGraphWin;
plot(sd,squeeze(PC),'-o');
xlabel('Duration (ms)'); ylabel('Percent correct')
grid on; set(l,'FontSize',12);

% disp(params)
% disp(params.vernier)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['temporalPooling-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC', 'params','svmMdl', 'barOffset', 'sd');

%%
