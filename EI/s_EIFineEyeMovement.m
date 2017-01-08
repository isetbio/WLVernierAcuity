%% Fine eye movements and different types
%  Note that this is done at high spatial (2 sec) and temoral (10 ms)
%  resolution. So it takes a while to run.
% 
% It would be very nice to figure out a way to run this using parfor
%

% Can be changed without recomputing everything
nTrials = 500;
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

params.vernier.barLength = 360;   % This is 0.2 deg

%%
secPerPixel = (6 / sc);
minPerPixel = (6 / sc) / 60;
degPerPixel = minPerPixel / 60;
barLength = params.vernier.barLength*minPerPixel;
barWidth  = params.vernier.barWidth*minPerPixel;
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

%%
% str = datestr(now,30);
% fname = fullfile(wlvRootPath,'EI','figures',['FineEyeMovements-',str,'.mat']);
% save(fname, 'PC', 'params', 'barOffset', 'scenes');

%%
% ddir = fullfile(wlvRootPath,'EI','figures');
% dfiles = dir(fullfile(ddir,'FineEyeMovements*'));
% 
% h = vcNewGraphWin;
% cnt = 0;
% for ii=1:length(dfiles)
%     load(dfiles(ii).name);
%     ii, size(PC)
%     PC
%     if ii == 1
%         PCall = PC; cnt = 1;
%     else
%         if size(PC) == size(PCall)
%             PCall = PCall + PC;
%             cnt = cnt + 1;
%         end
%     end
% end
% PCall = PCall/cnt;
% plot(secPerPixel*barOffset,PCall,'-o');
% xlabel('Offset arc sec'); ylabel('Percent correct')
% lStrings = cell({'No em','tremor only','drift only','micro saccade only','All'});
% grid on; l = legend(lStrings);
% set(l,'FontSize',12)
% set(gca,'ylim',[40 110]);

