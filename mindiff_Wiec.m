function err = mindiff_Wiec(spec, x, lrng, Ltap, rplanet, Lmax, M, sig, optA)
  % x = [rs,cTH,(magnitude)]
  % spec needs to include the zero degree value

  defval('sig',[]);
  defval('optA',false)
  
  Sw_loc = specWiec(x(1),x(2),1,rplanet,Lmax,Ltap,M);

  % Use only degrees in given range
  ls=min(lrng):max(lrng);
  

  % Normalize
  if optA
    Sw_loc = x(3)*Sw_loc;
  else
    A=bestA(Sw_loc(ls+1),spec(ls+1));
    Sw_loc=A*Sw_loc;
  end
  
  %% Calc eerror
  if isempty(sig)
    err = rms(log(Sw_loc(ls+1)) - log(spec(ls+2)));
  else
    nparam = 3;
    err = chisqSpecMisf(Sw_loc,spec,sig,nparam,ls);
    %sig = sig(ls+1);
    %err = 1/length(spec)  *  sum(  ( (Sw_loc - spec)./sig ).^2  );
  end
