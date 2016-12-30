function [imageBasisA,imageBasisC] = vaPCA(varargin)
% Make the principal component images for the absorptions and photocurrent
%
% The input argument is a struct that must include all the parameters
% needed to produce a vernier stimulus.
%
% We used to do this separately for absorptions and current, but no
% longer. Now, it is done for both, here.
%
% Steps:
%   Create a set of stimuli  with different offsets.
%   Calculate the responses over time
%   Combine the absorptions and currents into their large matrices
%   Compute the PCA across all the stimuli
%   Save two files with the parameters and the image bases.
%
% BW, ISETBIO Team, Copyright 2016

%% Basic parameters

% Figure out a sensible way to set this.  We average some number of trials
% to reduce noise, I think.
nTrials = 20;


%% Parse the inputs
p = inputParser;

p.KeepUnmatched = true;

% These parameters influence the image.  Though I wonder whether we might
% have the same basis images for different time steps?
p.addParameter('vernier',vernierP,@isstruct);

p.addParameter('tStep',5,@isscalar);
p.addParameter('tsamples',[],@isvector);
p.addParameter('timesd',20*1e-3,@isscalar);
p.addParameter('fov',0.25,@isscalar);

p.parse(varargin{:});

% Stimulus parameters for the vernier target
vernier = p.Results.vernier;

% We want plenty of basis functions stored.  We may not use all of them.
nBasis = 100;

% Time samples for the oiSequence
if isempty(p.Results.tsamples), tsamples = (-60:tStep:70)*1e-3;
else                            tsamples = p.Results.tsamples;
end

% Saved in the output files
basisParameters = varargin{1}; %#ok<NASGU>

%% Create the matched vernier stimuli

cMosaic = coneMosaic;

% Set the mosaic size to 15 minutes (.25 deg) because that is the spatial
% pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(p.Results.fov);

% Not sure why these have to match, but there is a bug and they do.
cMosaic.integrationTime = tsamples(2) - tsamples(1);

% Should we have noise?  Or not?
cMosaic.os.noiseFlag = 'random';
cMosaic.noiseFlag    = 'random';

%%  Store up the absorptions and current across multiple trials

absorptions = [];
current = [];
for offset = 0:1:7      % A large range of pixel offsets
    
    basisParameters.vernier.offset = offset;     % Pixels on the display
    [~,thisStim] = vaStimuli(basisParameters);
    tSamples = thisStim.length;

    emPaths = cMosaic.emGenSequence(tSamples,'nTrials',nTrials);
    
    [thisAbsorptions,thisCurrent] = cMosaic.compute(thisStim, ...
        'emPaths',emPaths, ...
        'currentFlag',true);
    % cMosaic.current = squeeze(mean(alignedC,1));
    % cMosaic.window;
    
    absorptions = cat(1,absorptions,thisAbsorptions);
    current = cat(1,current,thisCurrent);
    % Plot for one of the trials
    %  tmp = squeeze(thisAbsorptions(10,:,:,:));
    %  ieMovie(tmp);
    %
    % Show we did the cat the right way
    %  tmp = squeeze(absorptions(50,:,:,:));
    %  ieMovie(tmp);

end

%% Check that the movies make sense

% vcNewGraphWin;
% for ii=1:10:size(absorptions,1)
%     tmp = squeeze(absorptions(ii,:,:,:));
%     ieMovie(tmp);
%     pause(0.5);
% end

%% Convert the shape of the absorptions so we can perform the svd

tAbsorptions = trial2Matrix(absorptions,cMosaic);

% We make bases for each type of stimulus and then concatenate them.
[~,~,V] = svd(tAbsorptions,'econ');
imageBasis = V(:,1:nBasis); %#ok<*NASGU>

% Save the result and the parameters used to create the result
% When we load in s_vaAbsorptions we decide whether or not to recompute
% based on the match of basisParameters
save('imageBasisAbsorptions','imageBasis','basisParameters')
imageBasisA = imageBasis;

%% Convert the shape of the absorptions so we can perform the svd

tCurrent = trial2Matrix(current,cMosaic);

% We make bases for each type of stimulus and then concatenate them.
[~,~,V] = svd(tCurrent,'econ');
imageBasis = V(:,1:nBasis);

% Save the result and the parameters used to create the result When we load
% in s_vaCurrent we decide whether or not to recompute based on the match
% of basisParameters
save('imageBasisCurrent','imageBasis','basisParameters')
imageBasisC = imageBasis;

%% Visualize
% mx = max(imageBasis(:));
% mn = min(imageBasis(:));
% 
% % Have a look at the resulting bases
% vcNewGraphWin; 
% colormap(gray(256))
% for ii=1:nBasis
%     imagesc(reshape(imageBasis(:,ii),cMosaic.rows,cMosaic.cols),[1.2*mn 1.2*mx]);
%     title(sprintf('Basis %d',ii));
%     pause(1);
% end

end
