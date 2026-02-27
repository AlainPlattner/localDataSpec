function Sqf_reg = QFspec(rsource,robs,Lmax,Ltap)
% QF spec can be obtained from a regular spec by
%
% [~,~,~,~,~,~,bigm,bigl]=addmon(Lmax);
% qf = (bigm==0) | (mod(bigm+bigl,2)==0);
% coef(~qf)=0;
%
% See Langlais et al. (2014), http://dx.doi.org/10.1016/j.epsl.2014.05.013
  
% For a given rsource, calculate spec, normalized to start at 1
ltot = (0:(Lmax))';
Sqf = (rsource/robs).^(2*ltot+4);
% For mag fields, the L=0 entry should be 0
Sqf(1) = 0;
% To regionalize, transform the non-zonal spectrum into power density
Sqf_pd = Sqf./(ltot(:)+1)./(2*ltot(:)+1);
Sqf_pd_reg = localizeSpec(Sqf_pd,Ltap)';
% Transform regionalized power density back to Mauersberger-Lowes
Sqf_reg = Sqf_pd_reg.*(ltot+1).*(2*ltot+1);

