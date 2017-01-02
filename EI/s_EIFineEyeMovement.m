%% Fine eye movements and different types
%  Note that this is done at high spatial (2 sec) and temoral (10 ms)
%  resolution. So it takes a while to run.
% 
% It would be very nice to figure out a way to run this using parfor
%

nTrials = 500;

% Integration time 
tStep   = 10;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Scene field of view
sceneFOV = 0.35;

% Spatial scale that controls visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 3*(sceneFOV/0.35);  

s_EIParameters;

%%

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
minPerPixel = (6 / sc) / 60;

% So the bar length in arc minutes should be this number times minPerPixel
params.vernier.barLength = 360;   % This is 0.2 deg

% [aligned, offset, scenes] = vaStimuli(params);
% ieAddObject(scenes{2}); sceneWindow;

% Maybe compare the prob. correct at 6 sec when there are no eye movements to
% the standard eye movement parameter
barOffset = [0 1 2 4 6 7];
PC = zeros(numel(barOffset),5);

params.em.emFlag = [0 0 0]';
s_vaAbsorptions;
PC(:,1) = P(:);

params.em.emFlag = [1 0 0]';
s_vaAbsorptions;
PC(:,2) = P(:);

params.em.emFlag = [0 1 0]';
s_vaAbsorptions;
PC(:,3) = P(:);

params.em.emFlag = [0 0 1]';
s_vaAbsorptions;
PC(:,4) = P(:);

params.em.emFlag = [1 1 1]';
s_vaAbsorptions;
PC(:,5) = P(:);


%% Plot

lStrings = cell({'No em','tremor only','drift only','msaccade only','All'});

vcNewGraphWin;
plot(minPerPixel*barOffset*60,PC,'o-');
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)
set(gca,'ylim',[40 110]);

%%S
fname = fullfile(wlvRootPath,'EI','figures','FineEyeMovement.mat');
save(fname, 'PC', 'params', 'barOffset', 'scenes');

%%

load(fname)
