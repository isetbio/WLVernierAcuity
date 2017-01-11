function [fname,fnameStimulus] = vaFname(params)
% Return the filename for storing basis and stimulus in the va analysis.
%
% BW, ISETBIO Team, 2017

% Stored imageBasis filename with the same parameters as this.  Originally, I
% used these parameters and counted them.  But in fact, they should not count.
% We can use the same params to create Stimulus and imageBasis with different
% numbers of used basis or trials.  So, in this function we dummy up the number
% to what they were originally, 
params.nBasis  = 40;
params.nTrials = 300;

fname = fullfile(wlvRootPath,'local',[md5(savejson([],params)),'.mat']);

fnameStimulus = fullfile(wlvRootPath,'local',[md5(savejson([],params)),'-Stimulus.mat']);

end
