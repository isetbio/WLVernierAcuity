%% t_VernierScene
%
% Illustrate creating different vernier display scenes using an image and a
% display.
%
%  HJ/BW, ISETBIO TEAM, 2015


% Deal with the cone isolation stuff.
% Make a scenePlot case where it looks like a spectrogram with position on
% the x-axis and wavelength on the y-axis.

%  Initialize a new ISETBIO session
ieInit;

%% Create the RGB images for testing

% The img values are in DAC or linear RGB ????


[~,p] = imageVernier(p);   % Mainly to get the parameters

p.bgColor = 0.5;
p.offset = 0;
imgA = imageVernier(p);

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

% Add the scenes to the database
vcAddObject(sceneA); 
vcAddObject(sceneM); 

% Bring up the window and have a look
sceneWindow;

%% Examine some of the scene properties

% This is the scene luminance at the different sample points on the display
sz = sceneGet(sceneM,'size');
scenePlot(sceneM,'luminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance on the same line
scenePlot(sceneM,'radiance hline',[1,round(sz(1)/2)]);

%% Suppose use a different display?  Let's suppose a CRT
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

vcAddObject(sceneM); sceneWindow;

%% A pair of blue lines

% We could use this one to understand what happens with S-cones
% I would like to have the function that provides an S-cone direction for a
% display

% This converts the linear RGB values into LMS.
% This maps [r,g,b]* rgb2lms = [L,M,S]
rgb2lms = displayGet(d,'rgb2lms');

% To find the S cone isolating direction we calculate
% rgb*rgb2lms = [0 0 1] 
sconeRGB = [0,0,1]*inv(rgb2lms);
sconeRGB = sconeRGB/(2*max(sconeRGB));
p.barColor = [0.5 0.5 .5] + sconeRGB ;
p.barColor = ieLUTLinear(p.barColor,displayGet(d,'inverse gamma'));

imgM = imageVernier(p);
sceneM = sceneFromFile(imgM, 'rgb', [], d); % mis-aligned
scenePlot(sceneM,'radiance hline',[1,round(sz(1)/2)]);
vcAddObject(sceneM); sceneWindow;


%%