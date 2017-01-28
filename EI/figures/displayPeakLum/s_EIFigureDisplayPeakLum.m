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

%% Five bar offsets

PCall = zeros(5,nFiles);
peakLum = zeros(1,nFiles);
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

%% Interpolate the data
peakInterp = p(1):0.1:p(end);
[X,Y] = meshgrid(peakInterp,barOffsetSec(1):2:barOffsetSec(end));
PCInterp = interp2(p,barOffsetSec,PCsorted,X,Y); 

%%
vcNewGraphWin;
s = surf(X,Y,PCInterp);
set(gca,'fontsize',14);
xlabel('Max luminance (log_{10} cd/m^2)');
ylabel('Offset (arc sec)')
zlabel('Percent correct')

% view(-23.5, 22);   % Reasonable view
% view(90,-90); colorbar;  % Good view for curve

%% 
thresh = zeros(1,size(PCInterp,2));
barSteps = barOffsetSec(1):barOffsetSec(end);
for ii=1:size(PCInterp,2)
    p = interp1(Y(:,1),PCInterp(:,ii),barSteps);
    [v,idx] = min(abs(p - 75));
    thresh(ii) = barSteps(idx);
end

vcNewGraphWin;
set(gca,'FontName','Georgia','FontSize',14)

plot(peakInterp,thresh,'-o','LineWidth',2)
xlabel('Log_{10} luminance (cd/m^2)')
ylabel('Offset threshold (arc sec)')
grid on
saveas(gcf,'peakLuminance.png','png')

% Why doesn't this work?  It does work in the CSF script, I think.
% [fitProbC,thresh(ii),params] = PALweibullFit(barOffsetSec(:),PCInterp(:,3)/100,0.81,1000,barSteps);


%%
vcNewGraphWin;
for ii=1:size(Y,2)
    plot(Y(:,ii),PCInterp(:,ii),'-o')
    hold on
end
xlabel('Offset (arc sec)')
ylabel('Percent correct')
set(gca,'FontSize',14)
grid on

%%

