function coef=Olsen2Simons(coef)
% Needs to have L=0 in there, and be in ADDMON format

Lmax=sqrt(length(coef))-1;
[~,~,~,~,~,~,bigm,bigl]=addmon(Lmax);

% Here we omit the (l+1) factor, because that is part of the radial
% derivative. This is solely to renormalize the potential-field
% coefficients
fac = (-1).^(bigm(:))./(sqrt((2*bigl(:)+1)/(4*pi)));

coef = coef(:).*fac;
