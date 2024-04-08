function varargout = findDepthMinDiff_SRD(spec,lrng,robs,rref,rstart,Ltap,Lmax,sig)
% [rs,A] = findDepthMinDiff_SRD(spec,lrng,robs,rref,rstart,Ltap,Lmax,sig)
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
% sig     standard deviation if want to use Wieczorek's minimization
%  
% OUTPUT:
%
% rs       source radius [rm]
% A        amplitude (multiply with SRD spec for the given source radius)
%
% Last modified by plattner-at-alumni.ethz.ch, 07/10/2020


defval('sig',[])
  
%opts = optimset('MaxFunEvals',10000);
%rs = fminsearch(@(x) mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,sig), rstart, opts);

% Trying to solve for rs and amp
xstart = [rstart,rms(spec)];
xopt = fminsearch(@(x) mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,sig), xstart);%, opts);

if isempty(sig)
    xopt = xopt(1);
    % Because without the sig, we just use bestA instead of solved-for A.
end

if nargout < 2
    varargout = {xopt(1)};
else
    varargout = {xopt(1),xopt(2)};
end