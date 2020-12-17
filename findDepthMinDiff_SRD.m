function rs=findDepthMinDiff_SRD(spec,lrng,robs,rref,rstart,Ltap,Lmax,bias)
% rs=findDepthMinDiff_SRD(spec,lrng,robs,rref,rstart,Ltap,bias)
%
% Calculate source radius for which an upward-continued and regionalized
% Shell Randomly Oriented Dipole spectrum has the minimal difference with
%  the given Mauersberger-Lowes spectrum. 
% Here, difference means rms(log(spec1) - log(spec2));
%  
% INPUT:
%
% spec     power spectrum at robs, including degree zero
% lrng     degrees for which to fit the spectrum
% robs     observation radius [km]
% rref     reference radius (e.g. planetary radius)
% rstart   start value for source radius testing
% Ltap     tapering bandwidth
% Lmax     maximum spherical-harmonic degree for the McLeod spectrum
% bias     If local spec was obtained from tapering radial derivative
%          and then backtransforming, set this to "true" to take the
%          resulting bias into account
%  
% OUTPUT:
%
% rs       source radius [rm]
%
% Last modified by plattner-at-alumni.ethz.ch, 07/10/2020

  
opts = optimset('MaxFunEvals',10000);

rs = fminsearch(@(x) mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,bias), rstart, opts);
