function err=mindiffSig_McLeod(spec,rcore,lrng,Ltap,robs,rplanet,Lmax,sig)
% err=mindiffSig_McLeod(spec,rcore,lrng,Ltap,robs,Lmax,bias)
%
% Calculates the root mean square error of the misfit of the logariths of
% the provided regional spectrum and the regionalized McLeod spectrum for
% the degrees given in lrng.
%
% To incorporate the bias created by tapering the radial derivative,
% then get the potential-field coefficients, then power spec, then
% "undoing" the (l+1)^2 factor, set bias=true  

Smc_reg=McLeod(rcore,robs,rplanet,Lmax,Ltap);

%%% Take only the degrees within the given range
ls=min(lrng):max(lrng);
Smc_reg = Smc_reg(ls+1);
spec = spec(ls+1);
sig = sig(ls+1);


%%% Find best-fitting factor for regional McLeoud
%A=bestA(spec,Smc_reg);
A=bestAsig(Smc_reg,spec,sig);
Smc_reg=A*Smc_reg;

%%% Error is the difference of the log
err = rms(log(Smc_reg) - log(spec));
