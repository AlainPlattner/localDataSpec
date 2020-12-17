function spec=Olsen2Simons_spec(spec)
% Needs to have L=0 in there

warning('should put 4*pi in here')
  
Lmax=length(spec)-1;
ls=0:Lmax;

% Here we omit the (l+1) factor, because that is part of the radial
% derivative. This is solely to renormalize the potential-field
% coefficients
fac = 1./(2*ls(:)+1);

spec = spec(:).*fac;
