function data = vaImageBasis(params)
% VAIMAGEBASIS - Compute image basis for vernier acuity task
%
%
% Probably should just return the data and plot with ieMovie on the return

%%
data = vaPCA(params);
rows = sqrt(size(data,1));
cols = rows;
data = reshape(data,rows,cols,[]);

%%  Maybe this should be an ieMovie call?
ieMovie(data,'FrameRate',4);

% vcNewGraphWin; 
% colormap(gray(256))
% mx = max(data(:)); mn = min(data(:));
% for ii=1:params.nBasis
%     imagesc(reshape(data(:,ii),rows,cols),[mn mx]);
%     title(sprintf('Basis %d',ii));
%     pause(0.5);
% end

end
