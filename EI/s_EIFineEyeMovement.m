%% Fine eye movements and different types
%  Note that this is done at high spatial (2 sec) and temoral (10 ms)
%  resolution. So it takes a while to run.
% 
% It would be very nice to figure out a way to run this using parfor
%
disp('**** EI Fine Eye Movement')

% Can be changed without recomputing everything
nTrials = 1000;
nBasis  = 40;

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

% params.vernier.barLength = 360;   % This is 0.2 deg
%% Summarize

[~, offset,scenes,tseries] = vaStimuli(params);

% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;

%%
barLength = params.vernier.barLength*minPerPixel;
barWidth  = params.vernier.barWidth*minPerPixel;
fprintf('\nBar length %.1f min (%3.1f deg)\nBar width  %3.1f min\n',...
    (barLength),...
    (params.vernier.barLength*degPerPixel),...
    (barWidth));
fprintf('Bar offset %3.1f sec/pixel\n',secPerPixel);

%%

% Maybe compare the prob. correct at 6 sec when there are no eye movements to
% the standard eye movement parameter
barOffset = [0 1 2 4 6];
PC = zeros(numel(barOffset),5);
c  = gcp; if isempty(c), parpool('local'); end
emTypes = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ; 1 1 1]';
tic;
svmMdl = cell(1,size(emTypes,2));
parfor pp=1:size(emTypes,2)
    fprintf('Starting %d of %d ...\n',pp,size(emTypes,2));
    thisParam = params;
    thisParam.em.emFlag = emTypes(:,pp);
    [P, thisMdl] = vaAbsorptions(barOffset,thisParam);
    PC(:,pp) = P(:);
    svmMdl = thisMdl;
    fprintf('Finished %d\n',pp);
end
toc

%% Plot

lStrings = cell({'No em','tremor only','drift only','msaccade only','All'});

vcNewGraphWin;
plot(minPerPixel*barOffset*60,PC,'o-');
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)
set(gca,'ylim',[40 110]);

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['FineEyeMovements-',str,'.mat']);
fprintf('Saving %s\n',fname);
save(fname, 'PC', 'params', 'svmMdl', 'emTypes','barOffset', 'scenes');

%%

