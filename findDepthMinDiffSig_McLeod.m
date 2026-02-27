function rs=findDepthMinDiffSig_McLeod(spec,lrng,robs,rplanet,rstart,Ltap,Lmax,sig)
% rs=findDepthMinDiff_McLeod(spec,lrng,robs,rstart,Ltap,Lmax,bias)
%
% Calculate source radius for which an upward-continued and regionalized
% McLeod spectrum has the minimal difference with the given 
% Mauersberger-Lowes spectrum. 
% Here, difference means rms(log(spec1) - log(spec2));
%  
% INPUT:
%
% spec     power spectrum at robs, including degree zero
% lrng     degrees for which to fit the spectrum
% robs     observation radius [km]
% rplanet  planet radius [km] (required for bias)
% rstart   start value for source radius testing
% Ltap     tapering bandwidth
% Lmax     maximum spherical-harmonic degree for the McLeod spectrum
% sig      localization uncertainty
%
% OUTPUT:
%
% rs       source radius [rm]
%
% Last modified by plattner-at-alumni.ethz.ch, 02/27/2026
  
opts = optimset('MaxFunEvals',10000);

rs = fminsearch(@(x) mindiffSig_McLeod(spec,x,lrng,Ltap,robs,rplanet,Lmax,sig), rstart, opts);
