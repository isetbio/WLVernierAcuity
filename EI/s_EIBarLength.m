%% Impact of bar length

% Show the dependence on bar length for the computational observer. Use the
% match on bar length as an indicator of the spatial summation region of the
% human eye
nTrials = 300;

% Integration time 
% Captures eye movements up to 50HZ
% Adequate for absorptions (ms)
tStep   = 20;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 1;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.5;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 1.5*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%%
barOffset = [0 1 2 3 ];       % Pixels on the display
vals = [30 60 120 240 360];   % Bar length is half the FOV (6/1.5)*max(vals)/3600
PC = zeros(length(barOffset),length(vals));

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
minPerPixel = (6 / sc) / 60;

%%
for pp=1:length(vals)
    params.vernier.barLength = vals(pp);
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
    lStrings{pp} = sprintf('%.2f deg',offsetSec*vals(pp)/3600);
end

title(sprintf('Scene FOV %.1f',sceneGet(scenes{1},'fov')),'FontSize',14)
h = vcNewGraphWin;
plot(offsetSec*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%

fname = fullfile(wlvRootPath,'EI','figures','spatialBarLength.mat');
save(fname, 'PC','params', 'barOffset', 'vals','scenes');

%%

load(fname)
