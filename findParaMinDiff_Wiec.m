function para = findParaMinDiff_Wiec(spec,lrng,rplanet,startPara,Ltap,Lmax,sig,optA)
  % para = findParaMinDiff_Wiec(spec,lrng,rplanet,startPara,Ltap,Lmax,sig)
  %
  % Calculate the parameters for Wieczorek 2018 equation 27 
  % (same as Gong & Wieczorek 2021 eq 2).
  %
  % INPUT:
  %
  % spec         spectrum to be fitted
  % lrng           degrees for which to fit spectrum
  % rplanet      planet radius
  % startPara  starting values for the parameters [rs,cTH,d]
  %                 rs is center of sills, cTH is sill radis, d is sill
  %                 thickness
  % Ltap         tapering bandwidth
  % Lmax       maximum spherical harmonic degree
  % sig        standard deviation per degree of the multitaper spectrum
  % optA       want to also optimize magnitude? Default: false
  %
  % OUTPUT:
  %
  % para        optimal parameters [rs,cTH,d,magnitude]
  %             or [rs,cTH,d], if optA is false
  % 
  % Last modified by plattner-at-alumni.ethz.ch  4/8/2024

  defval('sig',[])
  defval('optA',false)
  
  % To speed up localization:
  try
    M = mcouplings(Ltap,Lmax,0);
  catch
    % In case the folder with the precalculated Wigner symbols
    % is empty, create one.
    wignercycle(1,0,0);
    M = mcouplings(Ltap,Lmax,0);
  end
  
  %opts = optimset('MaxFunEvals',10000);
  if optA
    Sw_loc = specWiec(startPara(1),startPara(2),startPara(3),1,rplanet,Lmax,Ltap,M);
    A = bestA(Sw_loc,spec);
    xstart = [startPara(:)',A];
  else
    xstart = startPara;
  end
  para = fminsearch(@(x) mindiff_Wiec(spec, x, lrng, Ltap, rplanet, Lmax, M, sig, optA) , xstart);%, opts);
