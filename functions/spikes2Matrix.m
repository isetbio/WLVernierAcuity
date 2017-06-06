function trialSortedSpikes = spikes2Matrix(nTrialsSpikes)
% Convert a cell array of nTrialsSpikes from N trials of RGC spike
% responses to a set of N trials x (K Cells * T time bins)
%
%
%%

nTrials = size(nTrialsSpikes{1},1);

for iTrial = 1:nTrials
    
    % spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
    % spikesout2 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{2},'spikes'));
    % spikesout3 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{3},'spikes'));
    % spikesout4 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{4},'spikes'));
    
    spikesout  = RGB2XWFormat(squeeze(nTrialsSpikes{1}(iTrial,:,:,:)));
    spikesout2 = RGB2XWFormat(squeeze(nTrialsSpikes{2}(iTrial,:,:,:)));
    spikesout3 = RGB2XWFormat(squeeze(nTrialsSpikes{3}(iTrial,:,:,:)));
    spikesout4 = RGB2XWFormat(squeeze(nTrialsSpikes{4}(iTrial,:,:,:)));
    
    timeBins = max([size(spikesout,2) size(spikesout2,2) size(spikesout3,2) size(spikesout4,2)]);
    
    spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins,'uint8');
    spikesoutsm(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
    spikesoutsm(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;
    
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;
    
    clear  spikesout1 spikesout2 spikesout3 spikesout4
    
    %%
    
    spikesout = double(spikesoutsm);
    pointer = 0;%(blockNum-1)*blocklength;
    spikeResp = zeros(size(spikesoutsm,1),size(spikesoutsm,2)/10);
    for i = 1:size(spikesoutsm,2)/10
        blocksize = 10;
        endval = i*blocksize;
        if endval > size(spikesout,2)
            endval = size(spikesout,2);
        end
        startval = (i-1)*blocksize + 1;
        spikeResp(:,pointer+i) = sum(spikesout(:,startval:endval),2);
    end
    
    % trialSortedSpikes(iTrial,:) = spikeResp(:);
    szCol = size(spikeResp,2);
    trialSortedSpikes((iTrial-1)*szCol+1:iTrial*szCol,:) = spikeResp';
end