function err=mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,M,sig)
% err=mindiff_SRD(spec,x,lrng,Ltap,robs,rref,Lmax,M,sig)
%
% Calculates the root mean square error of the misfit of the
% logariths of the provided regional spectrum and the regionalized
% Shell Randomly Oriented Dipole (SRD) spectrum 
% (Voorhies, Sabaka, Purucker, 2002) for the degrees given in lrng.
  defval('sig',[])
  defval('M',[])

rcore = x(1);
A = x(2);

Ssrd_reg=SRD(rcore,robs,rref,Lmax,Ltap,M);

%%% Take only the degrees within the given range
ls=min(lrng):max(lrng);
Ssrd_reg = Ssrd_reg(ls+1);
spec = spec(ls+1);


if isempty(sig)
    %%% Error is the difference of the log
    %%% In this case, provide A
    A=bestA(Ssrd_reg,spec);
    Ssrd_reg=A*Ssrd_reg;
    err = rms(log(Ssrd_reg) - log(spec));
else
    % This is the Wieczorek style. This is better
    Ssrd_reg=A*Ssrd_reg;
    sig = sig(ls+1);
    err = 1/length(spec)  *  sum(  ( (Ssrd_reg - spec)./sig ).^2  );
end
