function [coefsO,coefsS] = createRandCoefNZ(spc, Lfull, rplanet)

  [m,l,mzero]=addmon(Lfull);
  c=randn(addmup(Lfull),2);
  c(mzero,2)=0;
  % Here is setting the zonal part to zero:
  c(mzero,1)=0;
  lmcosiO=[l m c];
  
  [sdl,el,bto]=plm2spec(lmcosiO,1);
  srep=addmin(sdl);
  for i=1:size(lmcosiO,1)
    lmcosiO(i,3)=lmcosiO(i,3)./sqrt(srep(i)).*sqrt(spc(l(i)+1));
    lmcosiO(i,4)=lmcosiO(i,4)./sqrt(srep(i)).*sqrt(spc(l(i)+1));
  end
  % Make L=0 coefficient zero. They are currently NaN
  lmcosiO(1,[3:4]) = [0,0]; 
    
  coefsO = lmcosi2coef(lmcosiO,0);
  coefsS = Olsen2Simons(coefsO)*rplanet;
