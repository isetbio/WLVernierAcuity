%% Impact of blurring the optics on CSF and Vernier
%

% Show the dependence on spatial size of the cone mosaic for the computational
% observer.
nTrials = 100;
nBasis  = 40;

% Integration time 
tStep   = 50;  % Adequate for absorptions (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.25;

% Original scene
sceneFOV = 0.35;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg then 0.5/0.35
sc = 1*(sceneFOV/0.35);  

s_EIParametersCSF;

%% Summarize

%% Set up for the CSF calculation

freqSamples = [10 30 50];  % CPD
PC = zeros(length(freqSamples),length(2));

%% 
for pp=1:length(freqSamples)
    s_csfAbsorptions;
    PC(:,pp) = P(:);
end
% mesh(PC)

%% Make summary graphs

% Legend
lStrings = cell(1,length(cmFOV));
for pp=1:length(cmFOV)
    lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
end

h = vcNewGraphWin;
plot(secPerPixel*barOffset,PC);
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
str = datestr(now,30);
fname = fullfile(wlvRootPath,'EI','figures',['mosaicSize-',str,'.mat']);
save(fname, 'PC','params', 'barOffset', 'cmFOV','scenes');

%%
ddir = fullfile(wlvRootPath,'EI','figures');
dfiles = dir(fullfile(ddir,'mosaicSize*'));

% Legend
lStrings = cell(1,length(cmFOV));
for pp=1:length(cmFOV)
    lStrings{pp} = sprintf('%.2f deg',cmFOV(pp));
end

h = vcNewGraphWin;
cnt = 0;
for ii=1:length(dfiles)
    load(dfiles(ii).name);
    ii, size(PC)
    PC
    if ii == 1
        PCall = PC; cnt = 1;
    else
        if size(PC) == size(PCall)
            PCall = PCall + PC;
            cnt = cnt + 1;
        end
    end
end
PCall = PCall/cnt;
plot(secPerPixel*barOffset,PCall,'-o');
xlabel('Offset arc sec'); ylabel('Percent correct')
grid on; l = legend(lStrings);
set(l,'FontSize',12)

%%
