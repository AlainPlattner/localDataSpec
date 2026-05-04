function err=mindiffSig_QFspec(spec,rcore,lrng,Ltap,robs,Lmax,sig,mistype)
% err=mindiffSig_QFspec(spec,rcore,lrng,Ltap,robs,Lmax,sig,mistype)
%
% Calculates the root mean square error of the misfit of the logariths of
% the provided regional spectrum and the regionalized Nonzonal spectrum for
% the degrees given in lrng.

defval('mistype','log')

Sqf_reg=QFspec(rcore,robs,Lmax,Ltap);

%%% Take only the degrees within the given range
ls=min(lrng):max(lrng);
Sqf_reg = Sqf_reg(ls+1);
spec = spec(ls+1);
sig = sig(ls+1);

%%% Find best-fitting factor for regional McLeoud
%A=bestA(spec,Smc_reg);
A=bestAsig(Sqf_reg,spec,sig);
Sqf_reg=A*Sqf_reg;

switch mistype
  case 'log'
%%% Error is the difference of the log
err = rms(log(Sqf_reg) - log(spec));
    case 'csq'
        nparam=2;
        err=1/(length(ls)-nparam)  *  sum(  ( (Sqf_reg - spec )./sig ).^2  );
  otherwise
    err = rms(Smc_reg - spec);
end
