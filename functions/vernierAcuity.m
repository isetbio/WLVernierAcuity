function [s, params] = vernierAcuity(params)
%% vernierAcuity
%    Compute vernier acuity for test case defined by params
%
%  Inputs:
%    params.scene.d            - display structure
%    params.scene.fov          - scene field of view
%    params.scene.vDist        - viewing distance of the scene
%    params.scene.barWidth     - bar width in number of pixels
%    params.scene.offset       - offset of vernier bar in pixels
%    params.scene.bgColor      - background color for the scene
%    params.scene.barColor     - vernier bar color
%    params.scene.meanLum      - mean luminance of the scene
%
%    params.oi.defocus         - defocus of huamn optics in diopters
%    params.oi.pupil_d         - pupil diameter in mm
%    
%    params.sensor.density     - spatial density of K,L,M,S cones
%    params.sensor.fov         - sensor field of view
%    params.sensor.adapt_noise - bool, indicate to add cone noise or not
%    params.sensor.nFrames     - number of frames to compute for each scene
%
%    params.svm.opts           - svm options
%    params.svm.nFolds         - number of folds for cross validation
%    params.svm.method         - svm method, usually use 'linear'
%
%  Outputs:
%    s.sensor                  - sensor used in simulation, with no data
%
%    s.absorption.acc          - svm accuracy based on photon absorptions
%    s.absorption.err          - svm cross validation standard variation
%    s.absorption.weights      - svm weights
%    s.absorption.data         - cell array of cone photon absorption data
%
%    s.adaptation.acc          - svm accuracy based on cone adaptation
%    s.adaptation.err          - svm cross validation standard variation
%    s.adaptation.weights      - svm weights
%    s.adaptation.data         - cell array of cone current
%
% HJ/BW, ISETBIO TEAM, 2015

%% Create scene
%  load parameters
try d = params.scene.d; catch, d = displayCreate('LCD-Apple'); end
try vDist = params.scene.vDist; catch, vDist = 1.0; end

%  compute image size
ppi = displayGet(d, 'dpi');
try imgFov = params.scene.fov; catch, imgFov = [0.5 0.5]; end
if isscalar(imgFov), imgFov = [imgFov, imgFov]; end
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/39.37/vDist);     % Actual fov

%  set up params for sceneCreate('vernier')
p.display = d;
p.sceneSz = imgSz;
p.offset  = 0;  % for reference scene, no offset

try p.barWidth = params.scene.barWidth; catch, p.barWidth = 1; end
try p.bgColor = params.scene.bgColor; catch, p.bgColor = 0.5; end
try p.barColor = params.scene.barColor; catch, p.barColor = 0.9; end
try p.meanLum = params.scene.meanLum; catch, p.meanLum = []; end

%  create vernier scene
scene = cell(2, 1);
scene{1} = sceneCreate('vernier', 'display', p);

try p.offset = params.scene.offset; catch, p.offset = 1; end
scene{2} = sceneCreate('vernier', 'display', p);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

%% Create Human Lens
%  Load Zernike Coefficient
try pupilSize = params.oi.pupil_d; catch, pupilSize = 3; end
zCoefs = wvfLoadThibosVirtualEyes(pupilSize);

%  Create wavefront structure
wave = sceneGet(scene{1}, 'wave');
wvf = wvfCreate('wave', wave, 'zcoeffs', zCoefs, 'name', 'human optics');
wvf = wvfSet(wvf, 'calc pupil size', pupilSize);

%  Adjust for defocus
try defocus = params.oi.defocus; catch, defocus = 0; end
zDefocus = wvfDefocusDioptersToMicrons(defocus, pupilSize);
wvf = wvfSet(wvf, 'zcoeffs', zDefocus, {'defocus'});

%  Convert wvf to isetbio oi structure
wvf = wvfComputePSF(wvf);
oi = wvf2oi(wvf, 'human');

% Compute optical image
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);

%% Create Sensor
%  Load sensor parameters
try sensorFov = params.sensor.fov; catch, sensorFov = [.2 .2]; end
if isscalar(sensorFov), sensorFov = [sensorFov, sensorFov]; end
try nFrames = params.sensor.nFrames; catch, nFrames = 3000; end
try density = params.sensor.density; catch, density = [0 .6 .3 .1]; end
cone = coneCreate('human', 'spatial density', density);
sensor = sensorCreate('human', [], cone);

%  Set exposure time to 1 ms
expTime = sensorGet(sensor, 'exp time');
emDuration = 0.001;
emPerExposure = expTime / emDuration;
sensor = sensorSet(sensor, 'exp time', emDuration);
sensor = sensorSetSizeToFOV(sensor, sensorFov(1), scene{1}, OIs{1});
sz = sensorGet(sensor, 'size');

% Generate eyemovement
p.nSamples = nFrames * emPerExposure;
sensor = eyemoveInit(sensor, p);

% store sensor to output
s.sensor = sensor;

% Compute the cone absopritons
sensor1 = coneAbsorptions(sensor, OIs{1});

% Store the photon samples and add photons in one exposure time
pSamples1 = sensorGet(sensor1, 'photons');
pSamples1 = sum(reshape(pSamples1, [sz nFrames emPerExposure]), 4);
s.absorption.data{1} = pSamples1;

% Compute cone absorptions for the second stimulus and store photon
% absorptions
sensor2 = coneAbsorptions(sensor, OIs{2});
pSamples2 = sensorGet(sensor2, 'photons');
pSamples2 = sum(reshape(pSamples2, [sz nFrames emPerExposure]), 4);
s.absorption.data{2} = pSamples2;

%% SVM using cone absorptions
% Classification
try nFolds = params.svm.nFolds; catch, nFolds = 5; end
try opts = params.svm.opts; catch, opts = '-s 2 -q'; end
try method = params.svm.method; catch, method = 'linear'; end

labels = [ones(nFrames, 1); -1*ones(nFrames, 1)];
data = cat(1, RGB2XWFormat(pSamples1)', RGB2XWFormat(pSamples2)');

[acc, w] = svmClassifyAcc(data, labels, nFolds, method, opts);

% store classification results
s.absorption.acc = acc(1);
s.absorption.err = acc(2);
s.absorption.weights = reshape(mean(w, 2), sz);

%% SVM using the adapted cone current
% With cone noise added
try p.addNoise = params.sensor.adapt_noise; catch, p.addNoise = true; end

[~, adapted1] = coneAdapt(sensor1, 'rieke', p);
adapted1 = adapted1(:, :, 1:emPerExposure:end);
s.adaptation.data{1} = adapted1;

[~,adapted2] = coneAdapt(sensor2, 'rieke', params);
adapted2 = adapted2(:, :, 1:emPerExposure:end);
s.adaptation.data{2} = adapted2;

% Spatial time by position
data = cat(1, RGB2XWFormat(adapted1)', RGB2XWFormat(adapted2)');

% Perform the classification
[acc, w] = svmClassifyAcc(data, labels, nFolds, method, opts);

% store classification results
s.adaptation.acc = acc(1);
s.adaptation.err = acc(2);
s.adaptation.weights = reshape(mean(w, 2), sz);

%% Record experiment params
if nargout > 1
    params.scene.d = d;
    params.scene.fov = imgFov;
    params.scene.vDist = vDist;
    params.scene.barWidth = p.barWidth;
    params.scene.offset = p.offset;
    params.scene.bgColor = p.bgColor;
    params.scene.barColor = p.barColor;
    params.scene.meanLum = p.meanLum;
    
    params.oi.defocus = defocus;
    params.oi.pupil_d = pupilSize;
    
    params.sensor.density = density;
    params.sensor.fov = sensorFov;
    params.sensor.adapt_noise = p.addNoise;
    params.sensor.nFrames = nFrames;
    
    params.svm.opts = opts;
    params.svm.nFolds = nFolds;
    params.svm.method = method;
end

end