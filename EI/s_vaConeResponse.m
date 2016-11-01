%% Starting with two oiSequences, create cone mosaic response
%

%% This script creates the two sequences, oiSeqAligned and oiSeqOffset
s_vaStimulus;
who

%%
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.4 * imgFov);
cMosaic.integrationTime = 0.001;
cMosaic.emGenSequence(tSamples);

%%
cMosaic.computeForOISequence(oiSeqOffset);
cMosaic.computeCurrent;
cMosaic.window;

%%