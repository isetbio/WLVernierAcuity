%% s_vaAbsorptionsPCA
%
% Make the principal components for the absorption images.  
%
%   Create a set of stimuli  with different offsets.  
%   Calculate the absorption patterns over time
%   Combine the absorptions into one large matrix of absorptions
%   Compute the PCA across all the stimuli
%
% BW

%%
ieInit

nTrials = 20;
tStep   = 5;   % ms

% Set basic parameters for the vernier stimulus
clear p; 
p.barOffset = 0;     % Pixels on the display
p.barWidth  = 3;     % Pixels on the display
p.tsamples = (-60:tStep:70)*1e-3;   % In second
p.timesd = 20*1e-3;                 % In seconds

%% Create the matched vernier stimuli

cMosaic = coneMosaic;

% Set the mosaic size to 15 minutes (.25 deg) because that is the spatial
% pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(0.25);

% Not sure why these have to match, but there is a bug and they do.
cMosaic.integrationTime = p.tsamples(2) - p.tsamples(1);

nFrames = length(p.tsamples);

cMosaic.os.noiseFlag = 'random';
cMosaic.noiseFlag    = 'random';

%%
absorptions = [];
for offset = 0:1:7      % A large range of pixel offsets
    
    p.barOffset = offset;     % Pixels on the display
    [~,thisStim] = vaStimuli(p);
    tSamples = thisStim.length;

    emPaths = cMosaic.emGenSequence(tSamples,'nTrials',nTrials);
    
    thisAbsorptions = cMosaic.compute(thisStim, ...
        'emPaths',emPaths, ...
        'currentFlag',false);
    % cMosaic.current = squeeze(mean(alignedC,1));
    % cMosaic.window;
    
    absorptions = cat(1,absorptions,thisAbsorptions);
    % Plot for one of the trials
    %  tmp = squeeze(thisAbsorptions(10,:,:,:));
    %  ieMovie(tmp);
    %
    % Show we did the cat the right way
    %  tmp = squeeze(absorptions(50,:,:,:));
    %  ieMovie(tmp);

end

%% Check that the movies make sense

vcNewGraphWin;
for ii=1:10:size(absorptions,1)
    tmp = squeeze(absorptions(ii,:,:,:));
    ieMovie(tmp);
    pause(0.5);
end

%% Convert the shape of the absorptions so we can perform the svd

tAbsorptions = trial2Matrix(absorptions,cMosaic);

% We make bases for each type of stimulus and then concatenate them.
nComponents = 20;
[~,~,V] = svd(tAbsorptions,'econ');
imageBasis = V(:,1:nComponents);
mx = max(imageBasis(:));
mn = min(imageBasis(:));

% Have a look at the resulting bases
vcNewGraphWin; 
colormap(gray(256))
for ii=1:size(imageBasis,2)
    imagesc(reshape(imageBasis(:,ii),cMosaic.rows,cMosaic.cols),[1.2*mn 1.2*mx]);
    title(sprintf('Basis %d',ii));
    pause(1);
end

%%  We will do better, but for now

save('imageBasisAbsorptions','imageBasis')
