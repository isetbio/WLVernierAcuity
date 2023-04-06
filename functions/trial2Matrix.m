function imgList = trial2Matrix(coneData,cMosaic)
% Convert the 4D multitrial coneMosaic data to an image list
%
%
% BW, ISETBIO Team, 2016

%% Set up input parser

% 4D array input
tSamples = size(coneData,4);
nTrials  = size(coneData,1);

rows = cMosaic.rows;
cols = cMosaic.cols;

%%  Build the matrix

imgList = zeros(nTrials*tSamples,rows*cols);
for tt = 1:nTrials
    lst = (1:tSamples) + (tSamples)*(tt-1);
    thisTrial = squeeze(coneData(tt,:,:,:));
    thisTrial = permute(thisTrial,[3 1 2]);  
    thisTrial = reshape(thisTrial,tSamples,[]);
    imgList(lst,:) = thisTrial;
end


end
