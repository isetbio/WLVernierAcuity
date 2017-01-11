function [fname,fnameStimulus] = csfFname(params)
% Return the filename for storing basis and stimulus in the va analysis.
%
% BW, ISETBIO Team, 2017

% Stored imageBasis filename with the same parameters as this.  
% We do not care about certain parameters, so we null them out here before
% calculating the file name.
params.nBasis  = [];
params.nTrials = [];
params.oi = [];

fname = fullfile(wlvRootPath,'local',['csf-',md5(savejson([],params)),'.mat']);

fnameStimulus = fullfile(wlvRootPath,'local',['csf-',md5(savejson([],params)),'-Stimulus.mat']);

end
