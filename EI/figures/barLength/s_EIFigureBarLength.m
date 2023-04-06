%% Show the dependency on bar length
%
%  This is from the first Results figure.  It includes an image of the display
%  as well as the curves.
%

%%  Find the files  
ddir = fullfile(wlvRootPath,'EI','figures','barLength');
chdir(ddir);
dfiles = dir('spatial*');
nFiles = length(dfiles);

%% Load the probability correct matrices 
% Each is PC vs. bar offset

% Checking that the parameters match OK.  Chose the last bunch.
% for ii=1:nFiles
%     load(dfiles(ii).name);
%     barLengths
%     barOffset
% end

PCall = zeros(5,5);
for ii=5:nFiles
    load(dfiles(ii).name);
    PCall = PCall + PC;
end
PCall = PCall/4;

%%   Make images illustrating the scene

% This is just show time
params.vernier.offset = 8;
[aligned, offset, scenes] = vaStimuli(params);
scene = sceneAdd(scenes{1},scenes{2});
ieAddObject(scene); sceneWindow;
oi = offset.frameAtIndex(20); ieAddObject(oi); oiWindow;

% This is the offset used for the demonstration
secPerSample*params.vernier.offset

%% Make the image showing the cures

secPerSample = sceneGet(scenes{2},'degrees per sample')*3600;
barOffsetSec = barOffset*secPerSample;
thresh = zeros(1,size(PC,2));
barSteps = barOffsetSec(1):barOffsetSec(end);

vcNewGraphWin;
for ii=1:size(PC,2)
    plot(barOffsetSec,PC(:,ii),'-o','LineWidth',2);
    % [fitProbC,thresh(ii)] = PALweibullFit(barOffsetSec,PC(:,ii)/100,0.81,1000,barSteps);
    p = interp1(barOffsetSec,PC(:,ii),barSteps);
    [v,idx] = min(abs(p - 75));
    thresh(ii) = barSteps(idx);
    hold on
end

grid on
set(gca,'FontSize',14);
xlabel('Offset (arc sec)')
ylabel('Percent correct')
lstr = cell(1,numel(barLengths));
for ii=1:numel(barLengths)
    lstr{ii} = sprintf('%2.f arc min',secPerSample*barLengths(ii)/60);
end
legend(lstr);

%% Now show the threshold curve.  Not used.

vcNewGraphWin;
plot(barLengths*secPerSample/60,thresh,'-o','LineWidth',2);
grid on
xlabel('Bar length (min)');
ylabel('Offset threshold (arc sec)')

%mesh(bgColors,barOffsetSec,PC)

%% Screwing around with other representations

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
