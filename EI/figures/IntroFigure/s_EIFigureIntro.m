%% Eye movement impact
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?
%
disp('**** EI Intro')

nTrials = 1000;
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
sc = 3*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%% Make the vernier stimuli
params.vernier.offset = 8;
[~, offset,scenes,tseries] = vaStimuli(params);

% Show it in a window.
scene = sceneAdd(scenes{1}, scenes{2});
scene = sceneSet(scene,'name','offset');
ieAddObject(scene); sceneWindow;

%% 
oi = offset.frameAtIndex(20);
ieAddObject(oi); oiWindow;
oi = oiCrop(oi);
oi = oiSet(oi,'name','offset');
ieAddObject(oi); oiWindow;

%%
degPerPixel = sceneGet(scenes{2},'degrees per sample');
minPerPixel = degPerPixel * 60;
secPerPixel = minPerPixel * 60;


%%
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
