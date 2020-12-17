function err=specflat(spec,rcore,robs,lrng)
% x(1) is radius
% x(2) is amplitude
% spec must be given only for lrng

% Idea: Downward continue the spectrum to a depth and then see how flat it
% is. The optimal solution is the one with a minimal variance

ls=min(lrng):max(lrng);

% Calculate downward continuation factor
fac = (robs/rcore).^(2*ls(:)+4);

% Downward continue spectrum to assumed core radius
spec = spec(:).*fac;

% Calculate variance
err = var(log(spec));

