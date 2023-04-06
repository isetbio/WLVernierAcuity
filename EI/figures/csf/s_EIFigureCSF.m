%% Eye movement impact
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/IntroFigure';
ddir = fullfile(wlvRootPath,'EI','figures','csf');
chdir(ddir);
dfiles = dir('csf*');
nFiles = length(dfiles);

%%  Cleaned up the files
%
% defocusList = [0 .5 1 1.5 2];
% nD = length(defocusList);
% f = cell(nD,nFiles);
% fileDefocus = zeros(1,nFiles);
% hContrast = cell(nD,nFiles);
% PCall = cell(nD,nFiles);
% 
% for ii = 1:nFiles  
%     load(dfiles(ii).name,'params','PC','contrasts');  % PC, barOffset, params, scenes
%     d = find(params.defocus == defocusList);
%     fileDefocus(ii) = params.defocus;
%     f{d,ii} = params.freqSamples;
%     PCall{d,ii} = PC;
%     hContrast{d,ii} = contrasts;
% end

%%

% Set up the cell arrays
defocusList = [0 .5 1 1.5 2];
nD = length(defocusList);
maxFiles = 7;  % Number of 0 defocus is longest
% f = cell(nD,maxFiles);
% hContrast = cell(nD,maxFiles);
% PCall = cell(nD,maxFiles);

fileDefocus = zeros(1,nFiles);
for ii = 1:nFiles  
    load(dfiles(ii).name,'params');  % PC, barOffset, params, scenes
    d = find(params.defocus == defocusList);
    fileDefocus(ii) = params.defocus;
end

% for ii=1:length(lst)            % For each of the defocus groups
%     thisLst = lst{ii};
%     for jj=1:length(thisLst)    % This the data from a particular defocus
%         load(dfiles(thisLst(jj)).name,'params','PC','contrasts','freqSamples');  % PC, barOffset, params, scenes
%         fileDefocus = find(params.defocus == defocusList);
%         % f{ii,jj} = freqSamples;
%         % PCall{ii,jj} = PC;
%         % hContrast{ii,jj} = contrasts;
%     end
% end

% These are the files for each defocus
lst{5} = find(fileDefocus == 2);
lst{4} = find(fileDefocus == 1.5);
lst{3} = find(fileDefocus == 1);
lst{2} = find(fileDefocus == 0.5);   
lst{1} = find(fileDefocus == 0); 
    
%%

% A few minor deletions to make my life simpler
lst{3} = [8 10 17 18];   % 
lst{1} = [7  9 15 16];   %[2     3     4     7     9    15    16]
lst{2} = 5;

%%
PCall = cell(1,length(lst));
dFreqSamples = cell(1,length(lst));
dContrasts = cell(1,length(lst));
for dd = 1:5;   % Which defocus
    for ii=1:length(lst{dd})
        thisFile = lst{dd}(ii);
        load(dfiles(thisFile).name,'params','PC','contrasts','freqSamples');  % PC, barOffset, params, scenes
        dContrasts{dd} = contrasts;
        dFreqSamples{dd} = freqSamples;
        if ii==1, PCall{dd} = PC;
        else      PCall{dd} = PCall{dd} + PC;
        end
        % Should check that all the contrasts and freqSamples are the same
    end
end
for dd=1:5, PCall{dd} = PCall{dd}/length(lst{dd}); end

%%
for dd = 1:nD;
    vcNewGraphWin;
    p = semilogx(dContrasts{dd},PCall{dd},'-o','LineWidth',2);
    set(gca,'FontSize',14);
    grid on;
    xlabel('Log contrast'); ylabel('Percent correct');
end

%%
dd = 1;
fitLevels = logspace(log10(dContrasts{dd}(1)),log10(dContrasts{dd}(end)),100);
[fitProbC,fitThresh,fitParams] = PALweibullFit(dContrasts{dd},PCall{dd}(:,3)/100,0.81,1000,fitLevels);
vcNewGraphWin;
semilogx(fitLevels,fitProbC,'-',dContrasts{dd},PCall{dd}(:,3)/100,'o')

%%
nD = length(lst);
fThresh = zeros(nD,length(freqSamples));
tContrasts = logspace(-2.5,0,50);
for dd=1:nD
    probC = PCall{dd}/100;
    vcNewGraphWin;
    for ii=1:size(PCall{dd},2)
        % defocusList(dd), dFreqSamples{dd}(ii)
        [fitProbC,fThresh(dd,ii)] = PALweibullFit(dContrasts{dd},probC(:,ii),0.81,1000,tContrasts);
        p = interp1(dContrasts{dd},probC(:,ii),tContrasts,'pchip');
        semilogx(dContrasts{dd},probC(:,ii),'o',tContrasts,p,'LineWidth',2); hold on;
    end
    grid on; xlabel('Frequency (cpd)'); ylabel('Probability correct');
    set(gca,'FontSize',14); title(sprintf('Defocus %1.1f',defocusList(dd)));
end

% Patch up the NaNs as if they are infinite threshold
fThresh(isnan(fThresh)) = 1;

%%
% nD = length(lst);
% fThresh = zeros(nD,length(freqSamples));
% tContrasts = logspace(-2.5,0,50);
% for dd=1:5
%     for ii=1:size(PCall{dd},2)
%         % defocusList(dd), dFreqSamples{dd}(ii)
%         p = interp1(dContrasts{dd},PCall{dd}(:,ii),tContrasts,'pchip');
%         [v,idx] = min(abs(p - 75));
%         if v > 15, fThresh(dd,ii) = 1;
%         else fThresh(dd,ii) = tContrasts(idx);
%         end
%     end
% end
%%
vcNewGraphWin;
for dd=1:nD
    semilogy(freqSamples,1./fThresh(dd,:),'-o','LineWidth',2);
    hold on;
end
grid on;
xlabel('Frequency (cpd)'); ylabel('Contrast sensitivity')
set(gca,'FontSize',14);
legend({'0.0','0.5','1.0','1.5','2.0'});
%%
