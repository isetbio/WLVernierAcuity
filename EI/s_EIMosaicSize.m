%% Impact of bar length

% Show the dependence on bar length for the computational observer.
% Use the match on bar length as an indicator of the spatial summation region of
% the human eye
nTrials = 300;

% Integration time 
tStep   = 30;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Spatial scale to control visual angle of each display pixel
sc = 1;

s_EIParameters;
barOffset = [0  1  3  5 ];
cmSize = [.1 .15 .20 .25 .30];    % Degress of visual angle
PC = zeros(length(barOffset),length(vals));

for pp=1:length(vals)
    params. = vals(pp);
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end
% mesh(PC)

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

fname = fullfile(wlvRootPath,'EI','figures','spatialBarLength.mat','scenes');
save(fname, 'PC','params', 'barOffset', 'vals');