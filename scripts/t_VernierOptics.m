%% t_VernierOptics
%
% Illustrate the effects of human optics parameters with vernier scene
%
%  HJ/BW, ISETBIO TEAM, 2015

%% Initialize a new ISETBIO session
ieInit;

%% Create a verneir scene
%  See t_VernierScene for options in vernier scene create. Here, we use a
%  basic version, composing of two misaligned lines
p.bgColor = 0.5;
p.offset  = 1;
scene = sceneCreate('vernier', 'display', p);
vcAddObject(scene); sceneWindow;

%% Create standard human optics
%  create with Marimont and Wandell algorithm
oi = oiCreate('human');

% compute irradiance image
oi = oiCompute(oi, scene);
vcAddObject(oi); oiWindow;

% create with wavefront measurements
oi = oiCreate('wvf human');

% compute again and see
oi = oiCompute(oi, scene);
vcAddObject(oi); oiWindow;

%% Plot properties of optics
%  plot optical transfer function
oiPlot(oi, 'otf', [], 450); % otf at 450 nm
oiPlot(oi, 'otf', [], 550); % otf at 550 nm

%  plot point spread function
oiPlot(oi, 'psf', [], 450); % psf at 450 nm
oiPlot(oi, 'psf', [], 550); % psf at 550 nm

% line spread
nPoints = 50;
oiPlot(oi, 'ls wavelength', [], nPoints);

%  plot irradiance slice through one horizontal line
oiPlot(oi, 'irradiance hline', [1 20]);

%% Create customized human optics - pupil diameter, defocus, etc.
%  Load Zernike coefficient
pupilSize = 3; % pupil size in mm
zCoefs = wvfLoadThibosVirtualEyes(pupilSize);

%  Create wavefront structure
wvf = wvfCreate('wave', wave, 'zcoeffs', zCoefs, 'name', 'human optics');
wvf = wvfSet(wvf, 'calc pupil size', pupilSize); 

% Adjust for individuals
% Here, we use defocus as an example. For more adjustable entries, see
% wvfOSAIndexToVectorIndex

% ajust zernike coefficient for defocus 
% if we need to use defocus in diopters, use wvfDefocusDioptersToMicrons to
% do the conversion
defocus = -2.0; % diopters
zDefocus = wvfDefocusDioptersToMicrons(defocus); 
wvf = wvfSet(wvf, 'zcoeffs', zDefocus, {'defocus'});

% compute psf and convert to optical image structure
wvf = wvfComputePSF(wvf);
oi = wvf2oi(wvf, 'human');
oi = oiSet(oi, 'name', sprintf('Human WVF %.1f mm', pupilSize));

% let's have a look
vcAddObject(oi); oiWindow;
oiPlot(oi, 'otf', [], 450); % otf at 450 nm
oiPlot(oi, 'otf', [], 550); % otf at 550 nm
oiPlot(oi, 'psf', [], 450); % psf at 450 nm
oiPlot(oi, 'psf', [], 550); % psf at 550 nm

%% Other notes
%    1. Shift-invariant is assumed for the optics. 
%    2. lens transmittance is not handled in optics. Instead, lens
%       transmittance is handled together with macular transmittance in
%       sensor structure

%% END