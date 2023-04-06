function [acc, w] = coVernier(varargin)
% compute discrimination probability for vernier acuity
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) ISETBIO TEAM, 2016

% Parse inputs
p = parseInputs(varargin{:});

d = p.Results.display;            % display structure
imgFov = p.Results.imgFov;        % image field of view
vDist = p.Results.vDist;          % viewing distance
ppi = displayGet(d, 'ppi');       % display resolution

% Convert to pixels and then round based on actual number of integer pixels
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);  % number of pixels in image
imgFov = atand(imgSz/ppi/39.37/vDist);         % Actual fov

% Eye movement structure
em      = p.Results.em;         % eye movement structure
nFrames = p.Results.nFrames;    % number of frames in training

% Cone mosaic exposure time.
expTime = p.Results.expTime;      % exposure time

%% Create Scene
params.display = d;
params.sceneSz = [imgSz imgSz];
params.barWidth = 1;
params.barColor = p.Results.barColor;
params.bgColor = 0.5;

% Number of pixels for the offset
scene = cell(2, 1);
params.offset = 0; scene{1} = sceneCreate('vernier', 'display', params);
params.offset = 1; scene{2} = sceneCreate('vernier', 'display', params);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end
% fprintf('Scene size %d\n',sceneGet(scene{1},'size'));

%% Create Human optics
oi = oiCreate('wvf human');

% Compute optical image
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);

%% Create cone mosaic and eye movement sequence
cm = p.Results.cm;
cm.integrationTime = expTime;
cm.setSizeToFOV(p.Results.spatialInt,'sceneDist', ...
    vDist,'focalLength',oiGet(oi,'optics focal length'));
sz = cm.mosaicSize;

% Generate eye movements
emPerExposure = expTime/cm.sampleTime;
cm.emGenSequence(nFrames*emPerExposure,'em',em);
cm.integrationTime = cm.sampleTime;

%% Compute response
resp = cell(2, 1);
for ii = 1 : 2
    cm.compute(OIs{ii});
    switch ieParamFormat(p.Results.stage)
        case 'absorptions'
            resp{ii} = sum(reshape(cm.absorptions, ...
                [sz nFrames emPerExposure]), 4);
        case 'photocurrent'
            resp{ii} = reshape(cm.current, [sz nFrames emPerExposure]);
            resp{ii} = resp{ii}(:,:,:,10);
        case 'bipolar'
            bp = bipolar(cm.os);
            bp.set('sRFcenter', 1);
            bp.set(sRFsurround', 1);
            bp.compute(cm.os);
            resp{ii} = bp.responseCenter - bp.responseSurround;
        case 'rgc'
            error('NYI');
    end
end

%% prepare data for SVM linear classification
nFolds = 5;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
data = cat(1, RGB2XWFormat(resp{1})', RGB2XWFormat(resp{2})');

% Normalize data
data = bsxfun(@rdivide, bsxfun(@minus, data, mean(data)), std(data));

% Choose parameter C in svm and compute cross-validation error
opts = sprintf('-s 2 -q -C -v %d', nFolds);
res = train(labels, sparse(data), opts);
acc = res(2);

% Get weights of svm linear classifier
if nargout > 1
    opts = sprintf('-s 2 -q -c %e', res(1));
    svmStruct = train(labels, sparse(data), opts);
    w = reshape(svmStruct.w, sz);
end

end

function p = parseInputs(varargin)
% Parse input parameters in name-value pairs
%
% Support parameters:
%   'display'    - display structure
%   'imgFov'     - image field of view (degrees)
%   'spatialInt' - spatial integration size (degrees)
%   'vDist'      - viewing distance (meters)
%   'expTime'    - exposure time of human cone
%   'barColor'   - color of the bar, can be scalar or 3 element vector
%
% HJ, ISETBIO TEAM, 2016

p = inputParser;

p.addParameter('display', displayCreate('LCD-Apple', 'dpi', 400));
p.addParameter('imgFov', 0.5, @isnumeric);      % Degrees
p.addParameter('spatialInt', 0.2, @isnumeric);  % Degrees
p.addParameter('expTime', 0.05, @isnumeric);
p.addParameter('nFrames', 1000, @isnumeric);
p.addParameter('vDist', 1.0, @isnumeric);
p.addParameter('barColor', 0.99, @isnumeric);
p.addParameter('stage', 'absorptions', @ischar);
p.addParameter('cm', coneMosaic, @(x) isa(x, 'coneMosaic'));
p.addParameter('em', emCreate, @isstruct);

p.parse(varargin{:});

end