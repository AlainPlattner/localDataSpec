function varargout = findDepthMinDiff_SRD(spec,lrng,robs,rref,rstart,Ltap,Lmax,sig)
% [rs,A,chisq] = findDepthMinDiff_SRD(spec,lrng,robs,rref,rstart,Ltap,Lmax,sig)
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
% sig      standard deviation if want to use Wieczorek's minimization
%  
% OUTPUT:
%
% rs       source radius [rm]
% A        amplitude (multiply with SRD spec for the given source radius)
% chisq    chi-squared value of the solution
%
% Last modified by plattner-at-alumni.ethz.ch, 05/01/2024


defval('sig',[])
  
%opts = optimset('MaxFunEvals',10000);
%rs = fminsearch(@(x) mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,sig), rstart, opts);

% To speed up localization:
try
  M = mcouplings(Ltap,Lmax,0);
catch
  % In case the folder with the precalculated Wigner symbols
  % is empty, create one.
  wignercycle(1,0,0);
  M = mcouplings(Ltap,Lmax,0);
end

% Trying to solve for rs and amp
% First need to set good starting value for A
SRD_loc = SRD(rstart,robs,rref,Lmax,Ltap,M);
lsA = (min(lrng)+1) : (max(lrng)+1);
Astart = bestA(SRD_loc(lsA),spec(lsA));
xstart = [rstart,Astart];
[xopt,chisq] = fminsearch(@(x) mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,M,sig), xstart);%, opts);

if isempty(sig)
    xopt = xopt(1);
    % Because without the sig, we just use bestA instead of solved-for A.
end

if nargout < 2
    varargout = {xopt(1)};
else
    varargout = {xopt(1),xopt(2),chisq};
end
