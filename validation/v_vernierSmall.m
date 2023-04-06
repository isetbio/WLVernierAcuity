%% v_vernierSmall
%
% Validation for the vernier function on a small field of view.
% We want something that is relatively quick.
%
% This is a coarse display and a small field of view.
% It produces a high accuracy and the mean weight is close to zero.
%
% HJ/BW ISETBIO Team Copyright 2015

%% Init Parameters
ieInit;

% Set up all the default vernier params
params = vernierParams;

d = displayCreate('LCD-Apple');
d = displaySet(d, 'dpi', 200);
params.scene.d = d;
params.scene.vDist = 0.5;

params.sensor.density = [0 1.0 0 0]; % monochrome with only L cones
params.sensor.fov = [0.1 0.1];
params.sensor.nFrames = 500;

%% Execute the function
[s, params] = vernierAcuity(params);

%% Check

% Accuracy should be very high
assert(abs(s.absorption.acc - 1) <0.01)

% Check that the weights are stable?
w = mean(s.absorption.weights(:));
assert(abs(w) < 0.05 , sprintf('Mean of weights is %f.  Hoping for 0\n',w));

% imagesc(s.absorption.weights)

%% END
