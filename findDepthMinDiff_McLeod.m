function rs=findDepthMinDiff_McLeod(spec,lrng,robs,rplanet,rstart,Ltap,Lmax,bias)
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
%
% DEPRECATED INPUT:  
% bias     If local spec was obtained from tapering radial derivative
%          and then backtransforming, set this to "true" to take the
%          resulting bias into account. "True" here is deprecated.
%          This option is only kept here for testing.
%
% OUTPUT:
%
% rs       source radius [rm]
%
% Last modified by plattner-at-alumni.ethz.ch, 07/09/2020
defval('bias','false')
  
opts = optimset('MaxFunEvals',10000);

rs = fminsearch(@(x) mindiff_McLeod(spec,x,lrng,Ltap,robs,rplanet,Lmax,bias), rstart, opts);
