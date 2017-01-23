%% Eye movement impact
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/IntroFigure';
chdir(ddir);
dfiles = dir('single*');
nFiles = length(dfiles);

%% Image parameters

load(dfiles(3).name);  % PC, barOffset, params, scenes
barOffsetSec = barOffset*sceneGet(foo.scenes{2},'degrees per sample')*3600;

% Make the offset big so people can see it
params.vernier.offset = 4;

% This image is read into the displayWindow
I = imageVernier(params.vernier);
imwrite(I,'vernierStimulus.png','png');

% This is the retinal irradiance
[~,offset] = vaStimuli(params);
oi = offset.frameAtIndex(20); % save('offsetOI','oi');
ieAddObject(oi); oiWindow; 

oi2 = oiCrop(oi); ieAddObject(oi2); oiWindow;

% I loaded up the oi 
% load('offsetOI'); ieAddObject(oi); oi2 = oiCrop; oiAddObject(oi2); oiWindow;

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
