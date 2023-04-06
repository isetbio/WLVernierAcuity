%% Figure for JEF Gaussian line spreads
%

%% Create the scene and initiate the optics

vParams = vernierP;
scene = sceneCreate('vernier','display',vParams);
ieAddObject(scene); sceneWindow;
scene = sceneSet(scene,'fov',1.2);

%%
pupilMM = 3;
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
zCoefs(5) = 0;
oi = oiCreate('wvf human',pupilMM,zCoefs);

oi = oiCompute(oi,scene);
ieAddObject(oi); oiWindow;
% axis image; set(gca,'xlim',[-25 25],'ylim',[-25 25]);
% title(''); xlabel(''); ylabel(''); colormap(gray(256));
% set(gca,'xtick',[],'ytick',[]);
%% For different purposes we change the integration time and eye movement steps.
cm = coneMosaic;
cm.integrationTime = 0.005;
cm.emGenSequence(100);
cm.setSizeToFOV(1);
cm.compute(oi);
cm.computeCurrent;
cm.window;

%%
sharpConesTop = get(gca,'userdata');
sharpConesBot = get(gca,'userdata');
%%
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
zCoefs(5) = 2;
oi = oiCreate('wvf human',pupilMM,zCoefs);

oi = oiCompute(oi,scene);
ieAddObject(oi); oiWindow;
%%
cmBlur = coneMosaic;
cmBlur.integrationTime = 0.05;
cmBlur.setSizeToFOV(1);
% cmBlur.emGenSequence(100);
cmBlur.compute(oi);
cmBlur.window;

%%
defocusConesTop = get(gca,'userdata');
defocusConesBot = get(gca,'userdata');

save('defocusIllustration','defocusConesTop','defocusConesBot','sharpConesTop','sharpConesBot')

%%
pos = 50:100;
defocusT = interp1(defocusConesTop.pos{1},defocusConesTop.data{1},pos);
defocusB = interp1(defocusConesBot.pos{1},defocusConesBot.data{1},pos);
vcNewGraphWin;
plot(pos, defocusT, 'r-',pos,defocusB,'r--','LineWidth',2);
grid on
xlabel('Cone position')
ylabel('Photon absorptions (50 ms)');
set(gca,'ylim',[0 350]);


sharpT = interp1(sharpConesTop.pos{1},sharpConesTop.data{1},pos);
sharpB = interp1(sharpConesBot.pos{1},sharpConesBot.data{1},pos);
vcNewGraphWin;
plot(pos, sharpT, 'r-',pos,sharpB,'r--','LineWidth',2);
grid on
xlabel('Cone position')
ylabel('Photon absorptions (50 ms)');
set(gca,'ylim',[0 350]);

%%
scene = sceneCreate('point array');
scene = sceneSet(scene,'fov',0.3);
ieAddObject(scene); sceneWindow;


%% This is the way to do it by hand, and we need to do better with controlling the coefficients 

%% How to make the oi with different zernike coefficients
defocus = 0;
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wave = 400:10:700; wave = wave(:);

% Create wavefront parameters
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
wvfP = wvfSet(wvfP,'zcoeffs',defocus,{'defocus'});
wvfP = wvfComputePSF(wvfP);
% [u,p,f] = wvfPlot(wvfP,'2d psf space','um',550);
% set(gca,'xlim',[-20 20],'ylim',[-20 20]);

oi = wvf2oi(wvfP);
oi = oiSet(oi,'optics lens',Lens);
oi = oiSet(oi,'optics model', 'shiftInvariant');
oi = oiSet(oi, 'oi name', sprintf('human-wvf-%d',round(10*defocus)));

%%
