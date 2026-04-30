function err=mindiffSig_NZspec(spec,rcore,lrng,Ltap,robs,Lmax,sig,mistype)
% err=mindiff_NZspec(spec,rcore,lrng,Ltap,robs,Lmax,sig)
%
% Calculates the root mean square error of the misfit of the logariths of
% the provided regional spectrum and the regionalized Nonzonal spectrum for
% the degrees given in lrng.

  
defval('mistype','log')

Snz_reg=NZspec(rcore,robs,Lmax,Ltap);

%%% Take only the degrees within the given range
ls=min(lrng):max(lrng);
Snz_reg = Snz_reg(ls+1);
spec = spec(ls+1);
sig = sig(ls+1);


%%% Find best-fitting factor for regional McLeoud
%A=bestA(spec,Smc_reg);
A=bestAsig(Snz_reg,spec,sig);
Snz_reg=A*Snz_reg;

switch mistype
  case 'log'
%%% Error is the difference of the log
    err = rms(log(Snz_reg) - log(spec));
  otherwise
    err = rms(Smc_reg - spec);
end
