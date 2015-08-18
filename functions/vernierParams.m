function params = vernierParams(varargin)
% Initialize parameters for vernier acuity function
%
%    params = vernierParams;
%
% Parameter list
%
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
% HJ/BW ISETBIO Team, Copyright 2015

% Maybe we should set up a set/get thing

params.scene.d           = displayCreate('LCD-Apple');
params.scene.fov         = [0.5 0.5];
params.scene.vDist       = 1.0;
params.scene.barWidth    = 1;
params.scene.bgColor     = 0.5;
params.scene.barColor    = 0.9;
params.scene.meanLum     = [];

params.oi.defocus         = 0;
params.oi.pupil_d         = 3;  % mm

params.sensor.density     = [0 .6 .3 .1];
params.sensor.fov         = [.2 .2];
params.sensor.adapt_noise = true;
params.sensor.nFrames     = 3000;
params.svm.opts           = '-s 2 -q';
params.svm.nFolds         = 5;
params.svm.method         = 'linear';

% Some day, if we make set/get
% for ii=1:2:length(varargin)
% end


end
