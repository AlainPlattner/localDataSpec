function err=mindiffSig_QFspec(spec,rcore,lrng,Ltap,robs,Lmax,sig)
% err=mindiffSig_QFspec(spec,rcore,lrng,Ltap,robs,Lmax,sig)
%
% Calculates the root mean square error of the misfit of the logariths of
% the provided regional spectrum and the regionalized Nonzonal spectrum for
% the degrees given in lrng.


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

%%% Error is the difference of the log
err = rms(log(Sqf_reg) - log(spec));
