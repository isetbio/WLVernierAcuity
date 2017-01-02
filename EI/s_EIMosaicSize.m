%% Impact of cone mosaic size
%
% This complements the bar length analysis.  It should tell a similar story.
%

% Show the dependence on spatial size of the cone mosaic for the computational
% observer.
nTrials = 300;

% Integration time 
tStep   = 30;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Original scene
sceneFOV = 0.35;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = (sceneFOV/0.35);  % 4 arc sec

s_EIParameters;

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
secPerPixel = (6 / sc);
minPerPixel = (6 / sc) / 60;
degPerPixel = minPerPixel/60;
fprintf('Bar length is %.1f deg, %.1f min\n',...
    (params.vernier.barLength*degPerPixel),...
    (params.vernier.barLength*minPerPixel));

%%
fprintf('Bar offset per pixel is %.1f sec\n',secPerPixel);
barOffset = [0  1 2 3  4 ];

%%
cmFOV = [.1 .15 .25 .35];    % Degress of visual angle
PC = zeros(length(barOffset),length(cmFOV));

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
minPerPixel = (6 / sc) / 60;

localParams = params;

for pp=1:length(vals)
    params.cmFOV = cmFOV(pp);
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end
% mesh(PC)

%% 

% Offset per sample on the display
offsetDeg = sceneGet(scenes{1},'degrees per sample');
offsetSec = offsetDeg*3600;

% Legend
lStrings = cell(1,length(vals));
for pp=1:length(vals)
    lStrings{pp} = sprintf('%.1f arc min',offsetSec*vals(pp)/60);
end

h = vcNewGraphWin;
plot(offsetSec*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
fname = fullfile(wlvRootPath,'EI','figures','spatialBarLength.mat','scenes');
save(fname, 'PC','params', 'barOffset', 'vals');

%%
load(fname)
