%% Initial figure to illustrate the basic measurement
%
% Components would be the stimulus on the display, the OI, and the cone
% mosaic, and the data figure.
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/displayPeakLum';
ddir = fullfile(wlvRootPath,'EI','figures','displayPeakLum');
chdir(ddir);
dfiles = dir('example*');
nFiles = length(dfiles);

%%
peakLum = zeros(1,nFiles);
PCall = zeros(numel(barOffset),nFiles);
for ii=1:nFiles
    load(dfiles(ii).name);
    display = params.vernier.display;
    peakLum(ii) = displayGet(display,'peak luminance');
    PCall(:,ii) = PC(:);
end

% mesh(PCall)

%%
[p,idx] = sort(log10(peakLum));
PCsorted = PCall(:,idx);

p(3) = p(3) + 0.001;  % Include the duplicate.  Could average, but ...

barOffsetSec = barOffset*sceneGet(scenes{2},'degrees per sample')*3600;
vcNewGraphWin; mesh(p,barOffsetSec,PCsorted)

% Interpolate the data
[X,Y] = meshgrid(p(1):0.2:p(end),barOffsetSec(1):2:barOffsetSec(end));
PCInterp = interp2(p,barOffsetSec,PCsorted,X,Y); 

%%
vcNewGraphWin;
s = surf(X,Y,PCInterp);
set(gca,'fontsize',14);
xlabel('Max luminance (log 10 cd/m^2)');
ylabel('Offset (arc sec)')
zlabel('Percent correct')

% view(-23.5, 22);   % Reasonable view
% view(90,-90); colorbar;  % Good view for curve

%%
