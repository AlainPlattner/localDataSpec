function rs=findDepthMinDiffSig_QFspec(spec,lrng,robs,rstart,Ltap,Lmax)
% rs=findDepthMinDiffSig_QFspec(spec,lrng,robs,rstart,Ltap,Lmax)
%
% Calculate source radius for which an upward-continued and regionalized
% Nonzonal spectrum (Langlais et al. 2014) has the minimal difference with 
% the given Mauersberger-Lowes spectrum. 
% Here, difference means rms(log(spec1) - log(spec2));
%
% QF spec can be obtained from a regular spec by
%
% [~,~,~,~,~,~,bigm,bigl]=addmon(Lmax);
% qf = (bigm==0) | (mod(bigm+bigl,2)==0);
% coef(~qf)=0;
%
% See Langlais et al. (2014), http://dx.doi.org/10.1016/j.epsl.2014.05.013
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
%
% OUTPUT:
%
% rs       source radius [rm]
%
% Last modified by plattner-at-alumni.ethz.ch, 2/27/2026
  
opts = optimset('MaxFunEvals',10000);

rs = fminsearch(@(x) mindiff_QFspec(spec,x,lrng,Ltap,robs,Lmax,sig), rstart, opts);
