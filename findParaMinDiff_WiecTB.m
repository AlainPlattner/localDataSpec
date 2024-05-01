function varargout = findParaMinDiff_WiecTB(spec,lrng,rplanet,startPara,Ltap,Lmax,sig,optA)
  % [para,chisq] = findParaMinDiff_WiecTB(spec,lrng,rplanet,startPara,Ltap,Lmax,sig,optA)
  %
  % Calculate the parameters for Wieczorek 2018 equation 32
  %
  % INPUT:
  %
  % spec        spectrum to be fitted
  % lrng        degrees for which to fit spectrum
  % rplanet     planet radius
  % startPara   starting values for the parameters [rtop,rbot,cTH]
  %                 rtop is top radial position of sill,
  %                 rbot is bottom radial position of sill
  %                 cTH is sill diameter
  % Ltap        tapering bandwidth
  % Lmax        maximum spherical harmonic degree
  % sig         standard deviation per degree of the multitaper spectrum
  %             needed for chi-square calculation
  % optA        want to also optimize magnitude? Default: false
  %
  % OUTPUT:
  %
  % para        optimal parameters [rtop,rbot,cTH,magnitude]
  %             or [rtop,rbot,cTH], if optA is false
  % chisq       chi-squared value of the solution
  % 
  % Last modified by plattner-at-alumni.ethz.ch  5/1/2024

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
  %opts = optimset('Algorith','sqp');
  
  if optA
    Sw_loc = specWiecTB(startPara(1),startPara(2),startPara(3),1,rplanet,Lmax,Ltap,M);
    lsA = (min(lrng)+1) : (max(lrng)+1);
    Astart = bestA(Sw_loc(lsA),spec(lsA));
    %xstart = [startPara(:)',rms(spec)];
    xstart = [startPara(:)',Astart];
    %%%% Tried fmincon but never got good results
    %conmat = [-1,1,0,0]; % Meaning bot-rtop<=0 and -cTH<=0 
    %convec = [0];
  else
    xstart = startPara;
    %conmat = [-1,1,0]; % Meaning bot-rtop<=0 and -cTH<=0 
    %convec = [0];
  end
  [para,chisq] = fminsearch(@(x) mindiff_WiecTB(spec, x, lrng, Ltap, rplanet, Lmax, M, sig, optA) ,xstart);%, opts);

  if nargout < 2
    varargout = {para};
  else
    varargout = {para,chisq};
  end

%%% Haven't managed to get reasonable solutions with fmincon. Unfortunately, fminsearch can return rtop < rbot.
%%% maybe it's better to use findParMinDiff_Wiec instead ( using the thickness instead of top and bottom)
  %para = fmincon(@(x) mindiff_WiecTB(spec, x, lrng, Ltap, rplanet, Lmax, M, sig, optA), xstart, conmat, convec, [],[],[],[],[], opts);
