function err=mindiff_SRD(spec,rcore,lrng,Ltap,robs,rref,Lmax,bias)
% err=mindiff_SRD(spec,rcore,lrng,Ltap,robs,rref,Lmax)
%
% Calculates the root mean square error of the misfit of the
% logariths of the provided regional spectrum and the regionalized
% Shell Randomly Oriented Dipole (SRD) spectrum 
% (Voorhies, Sabaka, Purucker, 2002) for the degrees given in lrng.


Ssrd_reg=SRD(rcore,robs,rref,Lmax,Ltap,bias);

%%% Take only the degrees within the given range
ls=min(lrng):max(lrng);
Ssrd_reg = Ssrd_reg(ls+1);
spec = spec(ls+1);


%%% Find best-fitting factor for regional Randomly Oriented Dipoles
%A=bestA(spec,Smc_reg);
A=bestA(Ssrd_reg,spec);
Ssrd_reg=A*Ssrd_reg;

%%% Error is the difference of the log
err = rms(log(Ssrd_reg) - log(spec));
