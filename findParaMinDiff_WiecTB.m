function para = findParaMinDiff_WiecTB(spec,lrng,rplanet,startPara,Ltap,Lmax,sig,optA)
  % para = findParaMinDiff_WiecTB(spec,lrng,rplanet,startPara,Ltap,Lmax,sig)
  %
  % Calculate the parameters for Wieczorek 2018 equation 32
  %
  % INPUT:
  %
  % spec         spectrum to be fitted
  % lrng           degrees for which to fit spectrum
  % rplanet      planet radius
  % startPara  starting values for the parameters [rtop,rbot,cTH]
  %                 rtop is top radial position of sill, rbot is bottom
  %                 radial position of sill
  %                 cTH is sill diameter
  % Ltap         tapering bandwidth
  % Lmax       maximum spherical harmonic degree
  % sig        standard deviation per degree of the multitaper spectrum
  %            needed for chi-square calculation
  % optA       want to also optimize magnitude? Default: false
  %
  % OUTPUT:
  %
  % para     optimal parameters [rtop,rbot,cTH,magnitude]
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
    xstart = [startPara(:)',rms(spec)];
  else
    xstart = startPara;
  end
  para = fminsearch(@(x) mindiff_WiecTB(spec, x, lrng, Ltap, rplanet, Lmax, M, sig) ,xstart);%, opts);
  
