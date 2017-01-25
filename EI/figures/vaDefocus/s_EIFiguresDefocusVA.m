%% s_EIFiguresDefocusVA
% 

%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/IntroFigure';
ddir = fullfile(wlvRootPath,'EI','figures','vaDefocus');
chdir(ddir);
dfiles = dir('vaDefocus*');
nFiles = length(dfiles);

% [basisFile, stimFile] = vaFname(params);
% delete(stimFile); delte(basisFile);

% [~, offset,scenes,tseries] = vaStimuli(params);
% 
% ieAddObject(scenes{2}); sceneWindow;
% ieAddObject(offset.oiModulated); oiWindow;
% degPerPixel = sceneGet(scenes{2},'degrees per sample');
% minPerPixel = degPerPixel * 60;
% secPerPixel = minPerPixel * 60;

%%

nD = length(defocus);
PCall = zeros(5,4);
% They all have the same defocus parameters
for ii = 1:nFiles  
   load(dfiles(ii).name,'params','defocus','secPerPixel','barOffset','PC');
   % params.vernier
   % defocus
   % barOffset
   PCall = PCall + PC;
end
PCall = PCall/nFiles;

%% Now load up the special case with 2 diopters

ddir = fullfile(wlvRootPath,'EI','figures','vaDefocus');
chdir(ddir);
dfiles = dir('2vaDefocus*');
nFiles = length(dfiles);
load(dfiles(1).name,'params','defocus','secPerPixel','barOffset','PC');
PCall(:,end+1) = PC*ones(3,1)*(1/3);


%%
vcNewGraphWin;
plot(barOffset*secPerPixel,PCall,'-o');
xlabel('Offset (arc sec)');
ylabel('Percent correct');
grid on

defocus = [0 0.5 1.0 1.5 2];
lStrings = cell(1,length(length(defocus)));
for pp=1:length(defocus)
    lStrings{pp} = sprintf('%.2f D',defocus(pp));
end
legend(lStrings);

%%