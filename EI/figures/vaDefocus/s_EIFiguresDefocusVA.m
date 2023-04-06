%% s_EIFiguresDefocusVA
% 

%%  Defocus for analyzing vernier acuity

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
PCall = zeros(5,4);
tmp = zeros(5,1); kk = 0;
% They all have the same defocus parameters
for ii = 1:nFiles % nFiles  
   load(dfiles(ii).name,'params','defocus','secPerPixel','barOffset','PC');
   % params.vernier
   % sii, defocus
   % barOffsets
   PCall = PCall + PC(:,1:4);   
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
plot(barOffset*secPerPixel,PCall,'-o','LineWidth',2);
xlabel('Offset (arc sec)');
ylabel('Percent correct');
grid on

defocus = [0 0.5 1.0 1.5 2];
lStrings = cell(1,length(length(defocus)));
for pp=1:length(defocus)
    lStrings{pp} = sprintf('%.2f D',defocus(pp));
end
legend(lStrings,'Location','Southeast');
set(gca,'fontsize',14)

%%
saveas(gcf,'defocusVA.png','png')

%%
barOffsetSec = barOffset*secPerPixel;
barOffsetInterp = barOffsetSec(1):0.5:barOffsetSec(end);
thresh = zeros(1,length(defocus));
for ii=1:length(defocus)
    p = interp1(barOffsetSec,PCall(:,ii), barOffsetInterp);
    [v,idx] = min(abs(p - 75));
    thresh(ii) = barOffsetInterp(idx);
end
vcNewGraphWin;
plot(defocus,thresh,'-o','LineWidth',2);
grid on
xlabel('Zernike defocus parameter');
ylabel('Threshold offset (arc sec)');
grid on
set(gca,'FontSize',16);
set(gca,'FontName','Georgia');
set(gca,'ylim',[0 10])
set(gca,'xtick',[0 1 2],'ytick',[0:2:10]);

xlabel(''), ylabel('')

%%
saveas(gcf,'defocusThreshVA.png','png');

%%

%%
