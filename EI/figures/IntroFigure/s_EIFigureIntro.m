%% Eye movement impact
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%%  Peak Luminance
ddir = fullfile(wlvRootPath,'EI','figures','IntroFigure');
chdir(ddir);
dfiles = dir('single*');
nFiles = length(dfiles);

%% Image parameters

load(dfiles(3).name);  % PC, barOffset, params, scenes
barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;

% Make the offset big so people can see it
params.vernier.offset = 4;

% This image is read into the displayWindow
I = imageVernier(params.vernier);
imwrite(I,'vernierStimulus.png','png');

d = displayCreate('LCD-Apple');
d = displaySet(d,'main image',I);
ieAddObject(d); displayWindow;

%% This is the retinal irradiance
[aligned,offset,scenes] = vaStimuli(params);
oi = offset.frameAtIndex(20); % save('offsetOI','oi');
ieAddObject(oi); oiWindow; 
oi2 = oiCrop(oi); ieAddObject(oi2); oiWindow;

oi = aligned.frameAtIndex(24); % save('offsetOI','oi');
ieAddObject(oi); oiWindow; 
oi2 = oiCrop(oi); ieAddObject(oi2); oiWindow;

scene = sceneAdd(scenes{1},scenes{2});
ieAddObject(scene); sceneWindow;

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
