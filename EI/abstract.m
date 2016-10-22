% Figure 3(a) shows the probability of an accurate discrimination, based on the
% cone absorptions, for a range of viewing distances and a particular choice of
% simulation parameters. Discrimination accuracy decreases as viewing distance
% increases. Figure 3(b) shows the learned weights for a linear SVM classifier
% at a viewing distance of one meter. The weights spread and the computational
% observer uses both positive and negative weights across the cone mosaic to
% maximize performance.

% In the final paper, we propose to report the position information in 
%
%   different stages of the retinal simulation 
%     absorptions, photocurrent, bipolar, and specific or all RGC
%   
%   retinal parameters
%      L vs. M spatial density, overall spatial density
%      Optics properties (defocus, astigmatism)
%      Effective eccentricity
%
% A variety of psychophysical experiments have been conducted to investigate the
% underlying mechanism and the marginal impact of parameters of stimuli on the
% positional acuity [21,22]. For example, Westheimer et al. measured Vernier
% acuity with different
%
%   line segment colors  - luminance, L cone, S cone, L-M isoluminant
%   lengths - 0.2 0.5 and 1 deg (should we use a small Gabor and phase)?
%   Temporal variations - Brief duration.  Flicker.  Bright line/ dark lsine
%
% in their seminal work [23,24].  
%
% We will simulate the effects of these stimulus manipulations and compare the
% simulations to the psychophysical and physiological literature.