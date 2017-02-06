%% Figure 6 - L-cone curves for default human and defocus
%
% Illustrate the 

disp('**** EI Figure 6')

nTrials = 100;
nBasis  = 40;

% Integration time 
% Captures eye movements up to 100HZ
% Adequate for absorptions (ms)
tStep   = 10;       

% Scene FOV.  Larger than usual to show the continuing pooling
sceneFOV = 0.6;

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.5;

% Spatial scale to control visual angle of each display pixel The rule is 6/sc
% arc sec for a 0.35 deg scene. The scale factor in the front divides the 6, so
% 6/1.5 is the secPerPixel here.
sc = 3*(sceneFOV/0.35);  
                           
s_EIParameters;
% If you want to initiate imageBasis by hand, do this
% vaPCA(params);

%%
integrationTime = 0.020;
params.vernier.offset = 16;
whichCone = 1;
[~,fname] = vaFname(params);
delete(fname);
[aligned, offset, scenes] = vaStimuli(params);

%%
cm = coneMosaic;
cm.integrationTime = integrationTime;
cm.compute(offset.frameAtIndex(20));
cm.window;

%% Plot the LMS values at the top, click ...

topSharp    = get(gca,'userdata');

bottomSharp = get(gca,'userdata');

vcNewGraphWin;
pos = 2*topSharp.pos{whichCone};
pos = pos - mean(pos);
plot(pos,topSharp.data{whichCone},'r-o'); hold on

pos = 2*bottomSharp.pos{whichCone};
pos = pos - mean(pos);
plot(pos,bottomSharp.data{whichCone},'r--o')
set(gca,'ylim',[0 200])
xlabel('Position (um)');
ylabel('Absorptions (20 ms)');
 off

%% Blur the oi and do it again
[~,fname] = vaFname(params);
delete(fname);
params.oi = oiDefocus(2);
[aligned, offset, scenes] = vaStimuli(params);

%%
cm = coneMosaic;
cm.integrationTime = integrationTime;
cm.compute(offset.frameAtIndex(20));
cm.window;

%% Plot the LMS values at the top, click ...
topBlur    = get(gca,'userdata');
bottomBlur = get(gca,'userdata');

vcNewGraphWin;
pos = 2*topBlur.pos{whichCone};
pos = pos - mean(pos);
plot(pos,topBlur.data{whichCone},'r-o'); hold on

pos = 2*bottomBlur.pos{whichCone};
pos = pos - mean(pos);
plot(pos,bottomBlur.data{whichCone},'r--o')
set(gca,'ylim',[0 200])

%%
