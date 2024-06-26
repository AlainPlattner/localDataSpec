function err = mindiff_WiecTB(spec, x, lrng, Ltap, rplanet, Lmax, M, sig, optA)
  % x = [rtop,rbot,cTH,(magnitude)]
  % spec needs to include the zero degree value

  defval('sig',[])
  defval('optA',false)

  Sw_loc = specWiecTB(x(1),x(2),x(3),1,rplanet,Lmax,Ltap,M);

  % Use only degrees in given range
  ls=min(lrng):max(lrng);
  %Sw_loc = Sw_loc(ls+1);
  %spec = spec(ls+1);

  % Normalize
  if optA
    Sw_loc = x(4)*Sw_loc;
  else
    if sig
      A = bestAsig(Sw_loc(ls+1),spec(ls+1),sig(ls+1));
    else
      A=bestA(Sw_loc(ls+1),spec(ls+1));
    end
    Sw_loc=A*Sw_loc;
  end
  
  %% Calc error
  if isempty(sig)
    err = rms(log(Sw_loc(ls+1)) - log(spec(ls+1)));
  else
    nparam = 4;
    err = chisqSpecMisf(Sw_loc,spec,sig,nparam,ls);
    %sig = sig(ls+1);
    %err = 1/length(spec)  *  sum(  ( (Sw_loc - spec)./sig ).^2  );
  end
  
  
