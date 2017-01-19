function vaImageBasis(params)
%
%
%
%

data = vaPCA(params);
rows = sqrt(size(data,1));
cols = rows;
vcNewGraphWin; 
colormap(gray(256))
mx = max(data(:)); mn = min(data(:));
for ii=1:params.nBasis
    imagesc(reshape(data(:,ii),rows,cols),[mn mx]);
    title(sprintf('Basis %d',ii));
    pause(0.5);
end

end
