function [imageBasisAbsorptions,imageBasisCurrent] = csfPCA(varargin)
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

%% Check if the PCA has already been computed

params = varargin{1};
fname = vaFname(params);

if exist(fname,'file')
    disp('Loading image basis from file - parameters match')
    load(fname,'imageBasisAbsorptions','imageBasisCurrent');
    return;
else 
    disp('Creating and saving image basis file - parameters do not match')
end

%% Basic parameters

% Figure out a sensible way to set this.  We average some trials to reduce
% noise, I think.
nTrials = 20;

% We want plenty of basis functions stored.  We may not use all of them.
nBasis = 100;

%% Parse the inputs - probably not needed any more.
% Just accept that params will be passed in

p = inputParser;

p.KeepUnmatched = true;

% These parameters influence the image.  Though I wonder whether we might
% have the same basis images for different time steps?
p.addParameter('harmonic',[]);

p.addParameter('tStep',5,@isscalar);
p.addParameter('tsamples',[],@isvector);
p.addParameter('timesd',20*1e-3,@isscalar);
p.addParameter('cmFOV',0.25,@isscalar);
p.addParameter('freqSamples',1,@isvector);

p.parse(varargin{:});

% Stimulus parameters for the vernier target
harmonic   = p.Results.harmonic;

% These are the frequency we train for in the PCA
freqSamples = p.Results.freqSamples;

% Time samples for the oiSequence
if isempty(p.Results.tsamples), tsamples = (-60:tStep:70)*1e-3;
else                            tsamples = p.Results.tsamples;
end

%% Create the matched vernier stimuli

cMosaic = coneMosaic;

% Set the mosaic size to 15 minutes (.25 deg) because that is the spatial
% pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(p.Results.cmFOV);

% Not sure why these have to match, but there is a bug and they do.
cMosaic.integrationTime = tsamples(2) - tsamples(1);

% Should we have noise?  Or not?
cMosaic.os.noiseFlag = 'random';
cMosaic.noiseFlag    = 'random';

%%  Store up the absorptions and current across multiple trials

absorptions = [];
current = [];
for ff = 1:length(freqSamples)      % A large range of pixel offsets
    
    params.harmonic.freq = freqSamples(ff);
    [~,thisStim] = csfStimuli(params);
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
imageBasisAbsorptions = V(:,1:nBasis); %#ok<*NASGU>

%% Convert the shape of the absorptions so we can perform the svd

tCurrent = trial2Matrix(current,cMosaic);

% We make bases for each type of stimulus and then concatenate them.
[~,~,V] = svd(tCurrent,'econ');
imageBasisCurrent = V(:,1:nBasis);

%% Save

save(fname,'imageBasisAbsorptions','imageBasisCurrent','params');

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
