function vaImageSVM(mdl,imageBasis,params)
% vaImageSVM Visualize the beta weights of the SVM
%
% See 'Support Vector Machines for Binary Classification' in fitcsvm
%
%  mdl.Beta appears to be the term that we multiply as an inner product
%  with the data to determine whether the stimulus is type (A) or (not A).
%  The sign of this value determines
%
%     x'*mdl.Beta + Bias
%
%  The data in mdl.Beta have size of tSamples*nBasis. For example, if there
%  are 150 tSamples in a trial, and 10 spatial image components, then Beta
%  is 1500.
%
%  If we want to express the meaning of these values we need to create a
%  movie in image space by multiplying the imageBasis by the weights at
%  each of the time points. So the classifier is a movie.  Find the vector
%  of nBasis weights, multiply that vector times the image basis, and that
%  provides the classifier weights at that moment in time.  Then move to
%  the next moment in time and multiply again.
%
% BW ISETBIO Team, 2017

nBasis = params.nBasis; 
nSamples = numel(params.tsamples);

% These are the beta weights, nBasis of them for each time point.
beta = mdl.Beta;
rows = sqrt(size(imageBasis,1));
cols = rows;

% This is the classifier image over time
img = zeros(rows,cols,nSamples);


for ii=1:nSamples
    % Get the nBasis weights for each time sample
    lst = (1:nBasis) + (ii-1)*nBasis;
    
    % Multiply the beta coefficients times the imageBasis to form the iamge
    % at this moment in time
    tmp = imageBasis*beta(lst);
    
    % Put the image into a movie and show it
    img(:,:,ii) = reshape(tmp,rows,cols);
    
end

vcNewGraphWin;
ieMovie(img);

end

