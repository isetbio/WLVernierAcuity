%% Combine data
%
ddir = fullfile(wlvRootPath,'EI','figures');
dfiles = dir(fullfile(ddir,'mosaicSize*'));

% To return the scene and such you could read params and run
d = load(dfiles(1).name);

% To verify some of the parameters, you could do this
[aligned, offset, scenes, tseries] = vaStimuli(d.params);
degPerSample = sceneGet(scenes{1},'degrees per sample');
minPerSample = degPerSample*60;
secPerSample = minPerSample*60;
sceneGet(scenes{1},'fov')

%%
vcNewGraphWin;
plot(d.barOffset*secPerSample,d.PC,'-o')

grid on; xlabel('Offset (arc sec)');
ylabel('Percent correct')

% Legend
lStrings = cell(1,length(cmFOV));
for pp=1:length(cmFOV)
    lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
end
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
% h = vcNewGraphWin;
% cnt = 0;
% for ii=1:length(dfiles)
%     load(dfiles(ii).name);
%     ii, size(PC)
%     PC
%     if ii == 1
%         PCall = PC; cnt = 1;
%     else
%         if size(PC) == size(PCall)
%             PCall = PCall + PC;
%             cnt = cnt + 1;
%         end
%     end
% end
% PCall = PCall/cnt;
% plot(secPerPixel*barOffset,PCall,'-o');
% xlabel('Offset arc sec'); ylabel('Percent correct')
% grid on; l = legend(lStrings);
% set(l,'FontSize',12)

%%