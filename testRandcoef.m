function testRandcoef(nz)

defval('nz',false)
  
Lfull = 20;
rplanet = 2440;
Ltap = 0;
rwish = 2000;


if nz
  spc = NZspec(rwish,rplanet,Lfull,Ltap);
  [coefsO,coefsS] = createRandcoefNZ(spc, Lfull, rplanet);
else
  spc = McLeod(rwish,rplanet,rplanet,Lfull,Ltap);
  [coefsO,coefsS] = createRandcoef(spc, Lfull, rplanet);
end
  
% Now get spec and compare
lmcosiO = coef2lmcosi(coefsO,0);
specML = plm2spec(lmcosiO,1);

semilogy(0:Lfull,specML,'+')
hold on
semilogy(0:Lfull,spc,'--x')
hold off
