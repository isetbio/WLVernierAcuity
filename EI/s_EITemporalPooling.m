%% Impact of temporal pooling
%
% BW, ISETBIO Team, 2017

%%
disp('**** EI Temporal Pooling')

nTrials = 1000;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  
sceneFOV = 0.5;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.12;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = 2*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%% More like the McKee-Westheimer condition

% Helps find the ones that are like this 
params.matchHuman = true;     % Black background, matching human

params.vernier.bgColor = 0;   % Bright bar on a zero background

params.timesd    = 200e-3;    %  Reduce effect of eye movements
params.em.emFlag = [1 1 0];   % Microsaccades are suppressed for hyperacuity

%% Summnarize

[~, offset,scenes,tseries] = vaStimuli(params);

% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Summarize spatial parameters
% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher
barOffset = [0 1 2 4];
barLength = params.vernier.barLength*minPerPixel;
barWidth   = params.vernier.barWidth*minPerPixel;
fprintf('\nOffset %.2f\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barOffset*secPerPixel),...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%%  Sweep out different durations

% Total of 1 sec duration
sd = [50 100 200 400 600]*1e-3;
PC = zeros(length(barOffset),length(sd));
svmMdl = cell(1, length(sd));
tic
if isempty(gcp), parpool('local'); end
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
plot(barOffset*secPerPixel,PC,'-o','LineWidth',2);
xlabel('Offset (arcsec)'); ylabel('Percent correct classification')
grid on; 
set(gca,'FontSize',16,'FontName','Georgia');
barOffset*secPerPixel
xlabel(''); ylabel(''); title('')

% disp(params)
% disp(params.vernier)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures','temporal',['temporalPooling-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC', 'params','svmMdl', 'barOffset', 'sd');

%%
