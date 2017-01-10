%% Initialize the base parameters for the EI paper
%
% Typically, we run this before the segments.  Then we over-ride these values to
% produce different curves.
%
% We also set different parameters, such as the barOffset in pixels.
%
% See ...

%% These are oisequence and other parameters
clear params

params.tsamples  = (-200:tStep:200)*1e-3;   % (sec) McKee/Westheimer was 200 ms
params.timesd    = 100e-3;                  % (sec) +/- 1 std is 200 ms               
params.nTrials   = nTrials;
params.tStep     = tStep;  % In milliseconds?  That's odd.  Let's fix everywhere
params.sc        = sc;
params.nBasis    = nBasis;
params.cmFOV     = coneMosaicFOV;           % Cone mosaic field of view (deg)
params.sceneFOV    = sceneFOV;              % Scene field of view (deg)
params.distance    = 0.3;
params.em          = emCreate;
params.em.emFlag   = [1 1 1]';
params.defocus     = defocus;               % Defocus in diopters

%% Set basic parameters for the harmonic stimulus

h = harmonicP;
h.row = 210*params.sc;
h.col = 210*params.sc;
h.GaborFlag = 0.2;

% Attach the vernier parameters
params.harmonic = h;

%% Set basic parameters for the oi
pupilMM = 3;

zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wave = 400:10:700; wave = wave(:);

% Create wavefront parameters
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
wvfP = wvfSet(wvfP,'zcoeffs',0,{'defocus'});
wvfP = wvfComputePSF(wvfP);
% [u,p,f] = wvfPlot(wvfP,'2d psf space','um',550);
% set(gca,'xlim',[-20 20],'ylim',[-20 20]);
oi = wvf2oi(wvfP);
oi = oiSet(oi,'optics lens',Lens);
oi = oiSet(oi,'optics model', 'shiftInvariant');
oi = oiSet(oi, 'oi name', 'human-wvf');

params.oi = oi;

%%

