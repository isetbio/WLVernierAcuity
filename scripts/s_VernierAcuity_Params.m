%% s_VernierAcuity_Params
%
%    Running vernier acuity under different params
%
% HJ/BW, ISETBIO TEAM, 2015


%% Init Parameters
ieInit;
d = displayCreate('LCD-Apple');
d = displaySet(d, 'ppi', 1000);

params.scene.d = d;
params.sensor.density = [0 1.0 0 0]; % monochrome with only L cones

command = '[s, params] = vernierAcuity(params);';

%% Test viewing distance
cprintf('*Keyword', 'Viewing Distance Effect');
for vDist = 0.2:0.1:0.6
    % run simulation
    params.scene.vDist = vDist;
    eval(command)
    
    % print progress info
    fprintf('\t(%.1f m) ',vDist)
    fprintf('Absorption Acc: %.2f%%', s.absorption.acc*100)
    fprintf('\tAdaptation Acc: %.2f%%\n', s.adaptation.acc*100)
    
    % save data
    fname = sprintf('~/VernierResults/vDist/vd_%.1f.mat', vDist);
    save(fname, 's', 'params', 'command');
end

%% Test defocus

%% Test SNR

%% Test line Contrast