function Smt=localizeSpec(S,Ltap)
  % Smt=localizeSpec(S,Ltap)
  %
  % Given a global power spectrum, calculate a localized power spectrum
  % for a given region and taper bandwidth. This can be done by using
  % the coupling matrix in Dahlen & Simons 2008, eq 145.
  % This is the same as using pyshtools approach
  % pyshtools.SHWindow.biased_spectrum
  %
  % INPUT:
  %
  % S          global spectrum. Must be from l=0 to Lmax=length(S)-1
  % Ltap       taper bandwidth
  %
  % OUTPUT:
  %
  % Smt        regional power spectrum
  %
  % Last modified by aplattner-at-alumni.ethz.ch, 05/19/2019

  Lmax = length(S)-1;
  try
    M = mcouplings(Ltap,Lmax,0);
  catch
    % In case the folder with the precalculated Wigner symbols
    % is empty, create one.
    wignercycle(1,0,0);
    M = mcouplings(Ltap,Lmax,0);
  end
    
  S=S(:)';
  Smt = S*M;
  



