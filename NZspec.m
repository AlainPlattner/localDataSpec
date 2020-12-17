function Snz_reg = NZspec(rsource,robs,Lmax,Ltap)

% For a given rsource, calculate spec, normalized to start at 1
ltot = (0:(Lmax))';
Snz = (rsource/robs).^(2*ltot+4);
% For mag fields, the L=0 entry should be 0
Snz(1) = 0;
% To regionalize, transform the non-zonal spectrum into power density
Snz_pd = Snz./(ltot(:)+1)./(2*ltot(:)+1);
Snz_pd_reg = localizeSpec(Snz_pd,Ltap)';
% Transform regionalized power density back to Mauersberger-Lowes
Snz_reg = Snz_pd_reg.*(ltot+1).*(2*ltot+1);

