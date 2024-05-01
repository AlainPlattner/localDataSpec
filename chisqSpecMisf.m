function csq = chisqSpecMisf(spec1,spec2,sig,nparam,lrng)
  % csq = chisqSpecMisf(spec1,spec2,sig,nparam,lrng)
  %
  % Calculates the Chi-squared misfit of two spectra (hopefully) similar
  % to Wieczorek (2018).
  % NOTE THAT THE PROVIDED SPECTRA MUST INCLUDE ALL DEGREES FROM 0 TO LMAX.
  % Otherwise the lrng subselection doesn't work properly
  %
  % INPUT:
  %
  % spec1, spec2   The two spectra (e.g. the analytical one and the data one)     
  % sig            Standard deviation of the spectrum (e.g. from the single taper
  %                spectra of the data one)
  % nparam         On how many independent parameters does the analytical spectrum
  %                depend
  % lrng           For which spherical harmonic degrees do you want to calculate
  %                the chi-squared misfit?
  %
  % OUTPUT:
  %
  % csq            The chi-squared value
  %
  % last modified by plattner-at-alumni.ethz.ch, 5/1/2024
  
  
  lrng = min(lrng):max(lrng);
  nparam = 2;
  csq = 1/(length(lrng)-nparam)  *  sum(  ( (spec1(lrng+1) - spec2(lrng+1))./sig(lrng+1) ).^2  );
