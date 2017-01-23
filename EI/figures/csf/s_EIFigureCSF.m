%% Eye movement impact
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/IntroFigure';
ddir = fullfile(wlvRootPath,'EI','figures','csf');
chdir(ddir);
dfiles = dir('csf*');
nFiles = length(dfiles);

%%
defocusList = [0 .5 1 1.5 2];
nD = length(defocusList);
f = cell(nD,nFiles);
fileDefocus = zeros(1,nFiles);
hContrast = cell(nD,nFiles);
PCall = cell(nD,nFiles);

for ii = 1:nFiles  
    load(dfiles(ii).name,'params','PC','contrasts');  % PC, barOffset, params, scenes
    d = find(params.defocus == defocusList);
    fileDefocus(ii) = params.defocus;
    f{d,ii} = params.freqSamples;
    PCall{d,ii} = PC;
    hContrast{d,ii} = contrasts;
end

%%
vcNewGraphWin;
for ii=1:nFiles
    fprintf('%d - %f \n',ii, params.defocus);
    semilogx(hContrast{ii},PCall{ii}');
    pause;
end

defocus = cell2mat(defocus);
for ii=1:length(defocusList)
    for jj=1:nFiles
    lst = ~isempty(f{ii,jj})
    

%% Image parameters

foo = load(dfiles(20).name);  % PC, barOffset, params, scenes
barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;

params.harmonic.freq = 10;
% The scene is not saved correctly
% ieAddObject(scenes{2}); sceneWindow;

% This is the retinal irradiance
[~,harmonic] = csfStimuli(params);
oi = harmonic.frameAtIndex(20);
ieAddObject(oi); oiWindow; 

% oi2 = oiCrop(oi); ieAddObject(oi2); oiWindow;

%%
cm = coneMosaic;
cm.setSizeToFOV(params.cmFOV);
cm.integrationTime = params.tStep*1e-3;
cm.emGenSequence(100);
cm.compute(oi);
cm.window;

%%
vcNewGraphWin;
plot(barOffsetSec,PC,'-o','LineWidth',2);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; 
set(l,'FontSize',18)
set(gca,'ylim',[40 110]);

%% Maybe for the talk?
%
% vaImageBasis(params);

imageBasis = vaPCA(params);   % There is probably a way to load this from VCS
img = vaImageSVM(svmMdl,imageBasis,params);
vcNewGraphWin;
imagesc(flipud(sum(img,3)));   % Not sure why the flipud
colormap('default'); axis image
colorbar


%%
