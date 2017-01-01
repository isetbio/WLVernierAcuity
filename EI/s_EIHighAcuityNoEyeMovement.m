%% FIne eye movements and different types

nTrials = 500;

% Integration time 
tStep   = 30;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Spatial scale to control visual angle of each display pixel
sc = 3;

s_EIParameters;

params.vernier.barLength = 120;

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

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetMin = offsetDeg*60;

lStrings = cell({'No em','tremor only','drift only','msaccade only','default'});

vcNewGraphWin;
plot(offsetMin*barOffset*60,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

fname = fullfile(wlvRootPath,'EI','figures','fineEyeMovements.mat');
save(fname, 'PC', 'params', 'barOffset');

