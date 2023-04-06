%% Impact of defocus on the VA threshold using photocurrent
%
%  Just using photocurrent itself is useful.

%%
disp('**** EI Defocus VA Current')

%%
% Show the dependence on the cone mosaic size for the computational
% observer.
nTrials = 1000;
nBasis  = 40;

% Integration time 
tStep   = 10;         % Adequate for photocurrent (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.35;

% Original scene
sceneFOV = 0.4;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 3*(sceneFOV/0.35);   % If you do not multiply by a scalar, offset is 6 arc sec

s_EIParameters;

% Make the bar length a little less than the scene size
params.vernier.barLength = params.vernier.sceneSz(1)-1;

%%  Build the stimuli if you want to check stuff

[~, offset,scenes,tseries] = vaStimuli(params);

% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%% Summarize

% Each pixel size is 6 arc sec per pixel when sc =  1.  Finer resolution when sc
% is higher.
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (params.vernier.barLength*minPerPixel),...
    (params.vernier.barLength*degPerPixel),...
    (params.vernier.barWidth*minPerPixel));
fprintf('Bar offset per pixel is %.1f sec\n',secPerPixel);
barOffset = [0 2 4 6 8];
fprintf('Offsets in seconds %.1f\n',barOffset*secPerPixel);


%% Compute for all defocus and store
% if isempty(gcp), parpool('local'); end

tic
defocus = [0 0.5 1 1.5]; % [0 0.5 1 1.5];
PC = zeros(length(barOffset),length(defocus));
svmMdl = cell(1,length(defocus));
parfor pp = 1:length(defocus)
    fprintf('Starting %d ...\n',pp);
    thisParams = params;
    thisParams.defocus = defocus(pp);
    
    % Some problem with the files, so deleting.
    %     [~,fname] = vaFname(params);
    %     delete(fname);
    
    % Make the oi
    thisParams.oi = oiDefocus(defocus(pp));
    
    % Compute classification accuracy
    [P, thisMdl] = vaCurrent(barOffset, thisParams);
    PC(:,pp) = P;
    svmMdl{pp} = thisMdl;
    fprintf('Finished defocus level %d\n',defocus(pp));
end
toc

%%
vcNewGraphWin;
plot(barOffset*secPerPixel,PC,'-o');
xlabel('Offset (arc sec)');
ylabel('Percent correct');
grid on
title('Current')

lStrings = cell(1,length(length(defocus)));
for pp=1:length(defocus)
    lStrings{pp} = sprintf('%.2f D',defocus(pp));
end
legend(lStrings);

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['vaDefocus-Current-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC','params', 'svmMdl','barOffset','secPerPixel','defocus');

%%
