%% Calibrate offset
%
% These notes summarize how we control offset with the Vernier parameters
%
% The main idea is:  offset ~ FOV/size
%
%  If size = 100 and FOV = 1.0, offset is 36 sec
%  If size = 100 and FOV = 0.5, offset is 18 sec
%  If size = 200 and FOV = 1.0, offset is 18 sec
%
% So,  offset = k * FOV / size, where k = (3600)
%
% To set the size to achieve an offset, use
%
%    size = k * FOV / offset 
%
% For example, if you want to fix the FOV and adjust different offsets by
% controlling size, then
%
%    size =  (3600 * FOV) / offset
%
% To get a 10 arc sec offset, for 0.5 deg FOV, size = (3600 * 0.5)/10 = 180
%
% The scene.distance does not matter because the controlling variable is
% fov and number of samples per fov
%
% BW, Vistasoft Team, 2016

%%
ieInit

%% Init Parameters

% Gaussian time series
tseries = ieScale(fspecial('gaussian',[1,150],30),0,1);

display = displayCreate('LCD-Apple');

clear fov
clear offsetSec
clear sparams;  % Scene parameters in general

sparams.fov      = 0.5;  % Deg  - Changes offset

% Basic vernier parameters for the oiSequence
clear vparams;
for ii = 3:-1:1
    vparams(ii) = vernierP;
    vparams(ii).display = display;
    vparams(ii).sceneSz =[180 180];  % Changes offset
end

% Uniform field
vparams(1).name = 'uniform'; vparams(1).bgColor = 0.5; vparams(1).barWidth = 0;

% Offset Line
vparams(2).name = 'offset';  vparams(2).bgColor = 0; vparams(2).offset = 1;

% Aligned lines
vparams(3).name = 'aligned'; vparams(3).bgColor = 0; vparams(3).offset = 0;

[offset, scenes] = oisCreate('vernier','add', tseries,'tparams',vparams([1 2]),'sparams',sparams);
offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
offsetSec = offsetDeg*3600;
fprintf('Offset in arc secs %.2f\n',offsetSec);

%% This plots:  Offset ~ FOV

fov = (.2:.2:1);
offsetSec = zeros(size(fov));

for ff = 1:length(fov)
    clear sparams;  % Scene parameters in general
    sparams.fov      = fov(ff);  % Deg  - Changes offset
    
    % Basic vernier parameters for the oiSequence
    clear vparams;
    for ii = 3:-1:1
        vparams(ii) = vernierP;
        vparams(ii).display = display;
        vparams(ii).sceneSz =[56 56];  % Changes offset
    end
    
    % Uniform field
    vparams(1).name = 'uniform'; vparams(1).bgColor = 0.5; vparams(1).barWidth = 0;
    
    % Offset Line
    vparams(2).name = 'offset';  vparams(2).bgColor = 0; vparams(2).offset = 1;
    
    % Aligned lines
    vparams(3).name = 'aligned'; vparams(3).bgColor = 0; vparams(3).offset = 0;
    
    [offset, scenes] = oisCreate('vernier','add', tseries,'tparams',vparams([1 2]),'sparams',sparams);
    offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
    offsetSec(ff) = offsetDeg*3600;
    fprintf('FOV: %.1f\n,Offset in arc secs %.2f\n',fov(ff), offsetSec(ff));
end

%% 

vcNewGraphWin;
plot(fov,offsetSec,'-o');
xlabel('FOV (deg)'); ylabel('Offset (arc sec)');

%% This plots offset ~ 1 / SIZE

sz = 50:50:250;
offsetSec = zeros(size(sz));

for ss = 1:length(sz)
    clear sparams;  % Scene parameters in general
    sparams.fov      = 1;  % Deg  - Changes offset
    
    % Basic vernier parameters for the oiSequence
    clear vparams;
    for ii = 3:-1:1
        vparams(ii) = vernierP;
        vparams(ii).display = display;
        vparams(ii).sceneSz = [sz(ss), sz(ss)];  % Changes offset
    end
    
    % Uniform field
    vparams(1).name = 'uniform'; vparams(1).bgColor = 0.5; vparams(1).barWidth = 0;
    
    % Offset Line
    vparams(2).name = 'offset';  vparams(2).bgColor = 0; vparams(2).offset = 1;
    
    % Aligned lines
    vparams(3).name = 'aligned'; vparams(3).bgColor = 0; vparams(3).offset = 0;
    
    [offset, scenes] = oisCreate('vernier','add', tseries,'tparams',vparams([1 2]),'sparams',sparams);
    offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
    offsetSec(ss) = offsetDeg*3600;
    fprintf('Size: %.1f\nOffset in arc secs %.2f\n',sz(ss), offsetSec(ss));
end

%% Offset ~ 1/Size

vcNewGraphWin;
loglog(sz,offsetSec,'-o');
xlabel('Size (pixels)'); ylabel('Offset (arc sec)');

%%
