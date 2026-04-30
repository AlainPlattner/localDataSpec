function rs=findDepthMinDiffSig_NZspec(spec,lrng,robs,rstart,Ltap,Lmax,sig,mistype)
% rs=findDepthMinDiff_NZspec(spec,lrng,robs,rstart,Ltap,Lmax,sig,mistype)
%
% Calculate source radius for which an upward-continued and regionalized
% Nonzonal spectrum (Langlais et al. 2014) has the minimal difference with 
% the given Mauersberger-Lowes spectrum. 
% Here, difference means rms(log(spec1) - log(spec2));
%  
% INPUT:
%
% spec     power spectrum at robs, including degree zero
% lrng     degrees for which to fit the spectrum
% robs     observation radius [km]
% rstart   start value for source radius testing
% Ltap     tapering bandwidth
% Lmax     maximum spherical-harmonic degree for the McLeod spectrum
% sig      localization uncertainty
% mistype  how to calculate spectral misfit: diff of log, or just diff
%
% OUTPUT:
%
% rs       source radius [rm]
%
% Last modified by plattner-at-alumni.ethz.ch, 4/30/2026

defval('mistype','log')
  
opts = optimset('MaxFunEvals',10000);

rs = fminsearch(@(x) mindiffSig_NZspec(spec,x,lrng,Ltap,robs,Lmax,sig,mistype), rstart, opts);
