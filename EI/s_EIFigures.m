%% Initial figure to illustrate the basic measurement
%
% Components would be the stimulus on the display, the OI, and the cone
% mosaic, and the data figure.
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%
ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures';
% ddir = fullfile(wlvRootPath,'EI','figures');
% dfiles = dir(fullfile(ddir,'spatialBarLength*'));
dfiles = dir(fullfile(ddir,'csf*'));

% To return the scene and such you could read params and run
chdir(ddir);
d = load('csf-20170118T155052.mat');
d = load(dfiles(1).name);

%%

% Either this
imageBasis = csfPCA(d.params);
% Or this
imageBasis = vaPCA(d.params);

% Rename this to svmImage()??
img = vaImageSVM(d.svmMdl{3},imageBasis,params);
vcNewGraphWin;
colormap('default');
imagesc(sum(img,3));

%% To verify some of the parameters, you could do this
[aligned, offset, scenes, tseries] = vaStimuli(d.params);
degPerSample = sceneGet(scenes{1},'degrees per sample');
minPerSample = degPerSample*60;
secPerSample = minPerSample*60;
sceneGet(scenes{1},'fov')

scene = sceneAdd(scenes{1},scenes{2});
scene = sceneSet(scene,'name','Vernier stimulus');
ieAddObject(scene); sceneWindow;

% Visualize the stimulus
% offset.visualize;

% Show the OI
% oi = oiAdd(offset.oiFixed,offset.oiModulated,[1 1]);
% ieAddObject(oi); oiWindow;
%
% [oi,rect] = oiCrop(vcGetObject('oi'),[]);
% ieAddObject(oi); oiWindow;
 
%%
vcNewGraphWin;
plot(d.barOffset*secPerSample,d.PC(:,3),'-o')
xlabel('Offset (arc sec)'); ylabel('Percent correct');
set(gca,'ylim',[40 110])
grid on
d.barLengths*degPerSample

% Visualize the cone mosaic
% P.resampling = 8; P.vDensity = false; P.customLambda = [];
% cm = coneMosaicHex(8,false,[]);
cm = coneMosaic;
cm.setSizeToFOV(0.3);
cm.integrationTime = d.params.tStep*1e-3;
cm.emGenSequence(numel(d.params.tsamples));
cm.compute(oi);
cm.window;


%%
vcNewGraphWin;
plot(d.barOffset*secPerSample,d.PC,'-o')

grid on; xlabel('Offset (arc sec)');
ylabel('Percent correct')

% Legend
lStrings = cell(1,length(cmFOV));
for pp=1:length(cmFOV)
    lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
end
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)