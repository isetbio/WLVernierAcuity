%% Impact of temporal pooling
%
% BW, ISETBIO Team, 2017

%% Initialize parameters for s_EIParameters call
nTrials = 300;
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

params.vernier.barWidth = 20;
params.tsamples  = (-500:tStep:500)*1e-3;   % In second M/W was 200 ms

%% Summarize spatial parameters
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
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%%  Sweep out different durations

% Total of 1 sec duration
barOffset = 2;
sd = [50 100 200 400 600]*1e-3;
PC = zeros(1,length(sd));

for pp=1:length(sd)
    fprintf('Loop %d of %d\n',pp,length(sd));
    params.timesd  = sd(pp);                  % In seconds                 
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end

%%
vcNewGraphWin;
plot(sd,squeeze(PC),'-o');
xlabel('Duratison (ms)'); ylabel('Percent correct')
grid on; set(l,'FontSize',12);

% disp(params)
% disp(params.vernier)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['temporalPooling-',str,'.mat']);
save(fname, 'PC', 'params', 'barOffset', 'sd');

%%
fname = fullfile(wlvRootPath,'EI','figures','temporalPooling.mat');
load(fname)
