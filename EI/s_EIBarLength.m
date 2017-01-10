%% Impact of bar length
% Show the dependence on bar length for the computational observer. Use the
% match on bar length with behavior as an indicator of the spatial summation
% region of the human eye

nTrials = 1000;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 0.6;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.5;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = (sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);


%%  Build the stimuli if you want to check stuff

[~, offset,scenes,tseries] = vaStimuli(params);

ieAddObject(scenes{2}); sceneWindow;
ieAddObject(offset.oiModulated); oiWindow;
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Summarize spatial parameters
% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
barLength  = params.vernier.barLength*minPerPixel;
barWidth   = params.vernier.barWidth*minPerPixel;
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%% Initialize offsets and lengths

barOffset  = [0 1 2 3 4];           % Pixels on the display
% Make this less than
% params.vernier.sceneSz(1)
barLengths = [30 60 120 240 350];   % Bar length is the top and bottom
if max(barLengths) > params.vernier.sceneSz(1)
    error('Bar length is too long');
end
PC = zeros(length(barOffset),length(barLengths));
fprintf('Max bar length %.2f\n',degPerPixel*max(barLengths));
fprintf('Mosaic size %.2f\n',coneMosaicFOV);

%% Run for all the bar lengths

tic;
c = gcp; if isempty(c), parpool('local'); end
for pp=1:length(barLengths)
    fprintf('Starting %d of %d ...\n',pp, length(barLengths));
    thisParam = params;
    thisParam.vernier.barLength = barLengths(pp);
    P = vaAbsorptions(barOffset,thisParam);
    PC(:,pp) = P(:);
    fprintf('Finished %d\n',pp);
end
toc

%% Plot

% Legend
lStrings = cell(1,length(barLengths));
for pp=1:length(barLengths)
    lStrings{pp} = sprintf('%.2f deg',degPerPixel*barLengths(pp));
end

title(sprintf('Scene FOV %.1f',sceneGet(scenes{1},'fov')),'FontSize',14)

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC,'-o');
xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%% Save
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['spatialBarLength-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'barOffset', 'barLengths','scenes');

%%

