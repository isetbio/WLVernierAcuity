%% s_VernierAcuity_Params
%
%    Running vernier acuity under different params
%
% HJ/BW, ISETBIO TEAM, 2015


%% Init Parameters
ieInit;
d = displayCreate('LCD-Apple');
d = displaySet(d, 'dpi', 1000);
mpd = displayGet(d, 'meters per dot');

params.scene.d = d;
params.sensor.density = [0 1.0 0 0]; % monochrome with only L cones

command = '[s, params] = vernierAcuity(params);';

%% Test viewing distance
cprintf('*Keyword', 'Viewing Distance Effect\n');
params.scene.offset = 1;
for vDist = 0.2:0.1:0.6
    % print progress info
    offset_sec = atand(params.scene.offset * mpd / vDist) * 3600;
    fprintf('\t%.1f m (%.2f arcsec): ',vDist, offset_sec)
    
    % run simulation
    params.scene.vDist = vDist;
    eval(command)
    
    % print results
    fprintf('Absorption Acc: %.2f%%', s.absorption.acc*100)
    fprintf('\tAdaptation Acc: %.2f%%\n', s.adaptation.acc*100)
    
    % save data
    fname = sprintf('~/VernierResults/vDist/vd_%.1f.mat', vDist);
    save(fname, 's', 'params', 'command');
end

%% Test defocus
cprintf('*Keyword', 'Defocus Effect\n')
params.scene.vDist = 0.5;

for defocus = -3:0.5:1
    % run simulation
    params.oi.defocus = defocus;
    eval(command);
    
    % print progress info
    fprintf('\t(%.1f diopter) ', defocus)
    fprintf('Absorption Acc: %.2f%%', s.absorption.acc*100)
    fprintf('\tAdaptation Acc: %.2f%%\n', s.adaptation.acc*100)
    
    % save data
    fname = sprintf('~/VernierResults/defocus/defocus_%.1f.mat', defocus);
    save(fname, 's', 'params', 'command');
end

%% Test SNR
cprintf('*Keyword', 'Test SNR\n')
params.oi.defocus = 0;
mean_lum = [1, 5, 10, 20, 40, 100, 200, 1000];

for ii = 1 : length(mean_lum)
    % run simulation
    params.scene.meanLum = mean_lum(ii);
    eval(command)
    
    % print progress info
    fprintf('\t(%.1f cd/m2) ', mean_lum(ii))
    fprintf('Absorption Acc: %.2f%%', s.absorption.acc*100)
    fprintf('\tAdaptation Acc: %.2f%%\n', s.adaptation.acc*100)
    
    % save data
    fname = sprintf('~/VernierResults/snr/snr_%d.mat', mean_lum(ii));
    save(fname, 's', 'params', 'command');
end

%% Test line Contrast
cprintf('*Keyword', 'Line Contrast\n')
params.scene.meanLum = 100;
bar_color = 0.6:0.1:0.9;
for color = bar_color
    % run simulation
    params.scene.barColor = color;
    eval(command)
    
    % print progress info
    fprintf('\t(Bar color:%.1f) ', color)
    fprintf('Absorption Acc:%.2f%%', s.absorption.acc*100)
    fprintf('\tAdaptation Acc: %.2f%%\n', s.adaptation.acc*100)
    
    % save data
    fname = sprintf('~/VernierResults/contrast/contrast%.1f.mat', color);
    save(fname, 's', 'params', 'command');
end