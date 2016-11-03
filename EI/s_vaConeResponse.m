%% Starting with two oiSequences, create cone mosaic response
%

%% This script creates the two sequences, oiSeqAligned and oiSeqOffset
s_vaStimulus;

%%
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.6 * imgFov);
cMosaic.integrationTime = 0.001;
cMosaic.emGenSequence(tSamples);

%%

% I put a check in compute() to call computeForOISequence()
% The direct call to computeForOISequence() should also work.
cMosaic.compute(oiSeqOffset);

cMosaic.compute(oiSeqOffset,'currentFlag',true);

% Could do these calls instead
%   cMosaic.computeForOISequence(oiSeqOffset);
%   cMosaic.computeCurrent;

cMosaic.window;

%%