%% Impact of bar length
% Show the dependence on bar length for the computational observer. Use the
% match on bar length with behavior as an indicator of the spatial summation
% region of the human eye

nTrials = 600;
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
sc = 1.5*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

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

%%  Build the stimuli if you want to check stuff
% 
% [~, offset,scenes,tseries] = vaStimuli(params);
% 
% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;
% 
% oiGet(offset.oiModulated,'angular resolution')*3600

%%
barOffset = [0 1 2 3 ];           % Pixels on the display
vals = [30 60 120 240 300 360];   % Bar length is half the cmFOV degPerPixel*max(vals)
PC = zeros(length(barOffset),length(vals));
fprintf('Max bar length %.2f\n',degPerPixel*max(vals));
fprintf('Half mosaic size %.2f\n',coneMosaicFOV/2);

%%
for pp=1:length(vals)
    params.vernier.barLength = vals(pp);
    s_vaAbsorptions;
    PC(:,pp) = P(:);
end
% mesh(PC)

%% Plot

% Legend
lStrings = cell(1,length(vals));
for pp=1:length(vals)
    lStrings{pp} = sprintf('%.2f deg',degPerPixel*vals(pp));
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
save(fname, 'PC','params', 'barOffset', 'vals','scenes');

%% Something like this.
ddir = fullfile(wlvRootPath,'EI','figures');
dfiles = dir(fullfile(ddir,'spatialBarLength*'));

h = vcNewGraphWin;
cnt = 0;
for ii=1:length(dfiles)
    load(dfiles(ii).name);
    ii, size(PC)
    PC
    if ii == 1
        PCall = PC; cnt = 1;
    else
        if size(PC) == size(PCall)
            PCall = PCall + PC;
            cnt = cnt + 1;
        end
    end
end
PCall = PCall/cnt;
plot(secPerPixel*barOffset,PCall,'-o');

% Legend
lStrings = cell(1,length(vals));
for pp=1:length(vals)
    lStrings{pp} = sprintf('%.2f deg',degPerPixel*vals(pp));
end

xlabel('Offset arc sec'); ylabel('Percent correct')
set(gca,'ylim',[45 100]);
grid on; l = legend(lStrings);
set(l,'FontSize',12)

