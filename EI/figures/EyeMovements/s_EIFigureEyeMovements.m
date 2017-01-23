%% Eye movement impact
%
% Should we also try to show the SVM classifier?
%
% Movie of the eye movements?

%%  Peak Luminance
% ddir  = '/Volumes/users/wandell/github/WL/WLVernierAcuity/EI/figures/EyeMovements';
ddir = fullfile(wlvRootPath,'EI','figures','EyeMovements');
chdir(ddir);
dfiles = dir('FineEye*');
nFiles = length(dfiles);

%%
PCall = zeros(5,5);
for ii = 1:3   % Good one is 8,9 Bad 3
    load(dfiles(ii+7).name);  % PC, barOffset, params, scenes
    barOffsetSec = barOffset*sceneGet(foo.scenes{2},'degrees per sample')*3600;
    PCall = PCall + PC;
end
PCall = PCall/3;

%%
vcNewGraphWin;
plot(barOffsetSec,PCall,'-o','LineWidth',2);
lStrings = cell({'None','Tremor only','Drift only','Microsaccade only','All'});
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; 
l = legend(lStrings,'Location','SouthEast');
set(l,'FontSize',14)
set(gca,'FontSize',14);
set(gca,'ylim',[40 110]);

%% Maybe for the talk?
%
% vaImageBasis(params);
