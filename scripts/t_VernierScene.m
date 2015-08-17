%% t_VernierScene
%
% Illustrate creating different vernier display scenes using an image and a
% display.
%
%  HJ/BW, ISETBIO TEAM, 2015

% Initialize a new ISETBIO session
ieInit;

%% General observation about scenes
rd = ieRdata('create');
val = ieRdata('load data',rd,'benchHDR.mat','scene');
ieAddObject(val.scene); sceneWindow;

% Show the depth map

%% Create the RGB images for testing

[~, p] = imageVernier();   % Mainly to get the parameters

p.pattern = 0.5*ones(1,33); p.pattern(17) = 1;
% p.bgColor = 0.5;

% Aligned
p.offset = 0;
imgA = imageVernier(p);

% Misaligned
p.offset = 2;
imgM = imageVernier(p);
        
vcNewGraphWin([], 'tall'); 
subplot(2,1,1); imshow(imgA); axis image; title('Aligned Image');
subplot(2,1,2); imshow(imgM); axis image; title('Misaligned Image');

%% Choose a display for simulation
% We start with an Sony OLED display that we once calibrated.
dpi = 110;
d = displayCreate('OLED-Sony','dpi',dpi);

viewDist = 0.3;  % Set the subject's viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);

% In general, RGB images are treated as DAC values, not linear RGB. They
% will be converted into linear RGB by display gamma in vcReadImage line
% 180. But we want to calculate here with linear RGB.  So we set the table
% to linear and thus there is no difference between DAC and linear RGB.
d = displaySet(d, 'gtable','linear');

% BW2HJ: Let's fix this by adding the field to the display
vcSESSION.imgData = imgM;
ieAddObject(d);
displayWindow;

%% Create Vernier Scene (full display radiance representation)
%
%  We create a vernier scene radiance image by specifying a image on a
%  calibrated display. This method makes each of the parameters explicit.
%  This provides flexibility for scripting.
%

% Create a scene with the image using the display parameters
% The scene spectral radiance is created using the RGB image and the
% properties of the display.
sceneA = sceneFromFile(imgA, 'rgb', [], d); % aligned
sceneM = sceneFromFile(imgM, 'rgb', [], d); % mis-aligned

fov = size(imgA,2)/displayGet(d,'dots per deg');
sceneA = sceneSet(sceneA,'fov',fov);
sceneM = sceneSet(sceneM,'fov',fov);

% Add the scenes to the database
vcAddObject(sceneA); 
vcAddObject(sceneM); 

% Bring up the window and have a look
sceneWindow;

%% Examine some of the scene properties
% This is the scene luminance at the different sample points on the display
sz = sceneGet(sceneM,'size');
scenePlot(sceneM,'luminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance of the same line shown as a
% spectrogram
scenePlot(sceneM,'radiance hline image',[1,round(sz(1)/2)]);

% Here the same data are shown as a surface plot
scenePlot(sceneM,'radiance hline',[1,round(sz(1)/2)]);

% vcNewGraphWin; 
% ph = sceneGet(sceneM, 'photons'); ph = squeeze(ph(:, 1, :));
% colormap('hot'); imagesc(ph);

%% Use a different display? Let's suppose a CRT
dpi = 85;
d = displayCreate('crt','dpi',dpi);
viewDist = 0.3;  % Set the subject's viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);

sceneM = sceneFromFile(imgM, 'rgb', [], d); % mis-aligned

% Bring up the window and have a look
% This is the scene luminance at the different sample points on the display
sz = sceneGet(sceneM,'size');
scenePlot(sceneM,'luminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance on the same line
scenePlot(sceneM,'radiance hline',[1,round(sz(1)/2)]);

ieAddObject(sceneM); sceneWindow;

%% A pair of L-cone lines
%
% We find the S-cone isolating direction and use the RGB in that direction
%

p.pattern = 0.5*ones(1,33); p.pattern(15:17) = 1;

% This converts the linear RGB values into LMS.
% This maps [r,g,b]* rgb2lms = [L,M,S]
% rgb2lms = displayGet(d,'rgb2lms');
lConeIsolating = unitLength([1 0 0]*displayGet(d,'lms2rgb'));
img = image1d(p.pattern,'rgb',lConeIsolating,'mean',0.5);

% This converts the linear RGB values into LMS.
% This maps [r,g,b]* rgb2lms = [L,M,S]
% rgb2lms = displayGet(d,'rgb2lms');
lConeIsolating = unitLength([1 0 0]*displayGet(d,'lms2rgb'));

% We are calculating as if we have a linear gamma table, for simplicity
dLinear = displaySet(d,'gtable','linear');
scene = sceneFromFile(img, 'rgb', [], dLinear); % mis-aligned
scenePlot(scene,'radiance hline',[1,round(sz(1)/2)]);
ieAddObject(scene); sceneWindow;


%% Render the scene through the human optics 

oi = oiCreate('wvf human');
oi = oiCompute(oi,scene);
ieAddObject(oi); oiWindow;

cones = sensorCreate('human');
cones = sensorSetSizeToFOV(cones,sceneGet(scene,'fov'),scene,oi);
cones = sensorSet(cones,'noiseflag',0);
cones = sensorCompute(cones,oi);
ieAddObject(cones); sensorWindow('scale',1);

row = round(sensorGet(cones,'rows')/2);

% Looks good.
sensorPlot(cones,'photons hline',[row,1]);

%% A pair of mainly S-cone lines
%
% We find the S-cone isolating direction and use the RGB in that direction
%

p.pattern = 0.5*ones(1,33); p.pattern(15:17) = 1;

% This converts the linear RGB values into LMS.
% This maps [r,g,b]* rgb2lms = [L,M,S]
% rgb2lms = displayGet(d,'rgb2lms');
sConeIsolating = unitLength([0 0 1]*displayGet(d,'lms2rgb'));
img = image1d(p.pattern,'rgb',sConeIsolating,'mean',0.5);

% This converts the linear RGB values into LMS.
% This maps [r,g,b]* rgb2lms = [L,M,S]
% rgb2lms = displayGet(d,'rgb2lms');
sConeIsolating = unitLength([0 0 1]*displayGet(d,'lms2rgb'));

% We are calculating as if we have a linear gamma table, for simplicity
dLinear = displaySet(d,'gtable','linear');
scene = sceneFromFile(img, 'rgb', [], dLinear); % mis-aligned
scenePlot(scene,'radiance hline',[1,round(sz(1)/2)]);
ieAddObject(scene); sceneWindow;


%% Notice that when we render the scene through the human optics 
%  onto a cone mosaic, a thin line is not that great at cone-isolation
%  This is because of chromatic aberration

%
oi = oiCreate('wvf human');
oi = oiCompute(oi,scene);
ieAddObject(oi); oiWindow;

%
cones = sensorCreate('human');
cones = sensorSetSizeToFOV(cones,sceneGet(scene,'fov'),scene,oi);
cones = sensorSet(cones,'noiseflag',0);
cones = sensorCompute(cones,oi);
ieAddObject(cones); sensorWindow('scale',1);

row = round(sensorGet(cones,'rows')/2);
sensorPlot(cones,'photons hline',[row,1]);

%% Does not work well when the mosaic is big.
% 
% Looks good for up to, say, 100 x 100
% We should write sensorCrop()
tmp = coneImageActivity(cones,[]);


%%