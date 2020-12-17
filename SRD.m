function Ssrd_reg=SRD(rsource,robs,rref,Lmax,Ltap,bias)
% SRD means Shell of Random Dipoles  
% This is equation (5b) from Voorhies, Sabaka, and Purucker (2002) 
% a=rref, the reference radius
% To incorporate the bias created by tapering the radial derivative,
% then get the potential-field coefficients, then power spec, then
% "undoing" the (l+1)^2 factor, set bias=true

  
%%% For a given rcore, calculate spec, normalized to start at 1
ltot = (0:(Lmax))';
Ssrd = [0;(ltot(2:end)).*(ltot(2:end)+0.5).*(ltot(2:end)+1).*(rsource/rref).^(2*ltot(2:end)-2).*(rref/robs).^(2*ltot(2:end)+4)];

% To regionalize, transform the SRD spectrum into power density
Ssrd_pd = Ssrd./(ltot(:)+1)./(2*ltot(:)+1);
if bias
  biasfac = (ltot+1).^2.*(robs/rref).^(-2*ltot-4);
  Ssrd_pd = Ssrd_pd.*biasfac(:);
else
  biasfac = 1;
end
% Regionalize power density spectrum
Ssrd_pd_reg = localizeSpec(Ssrd_pd,Ltap)';
% Undo the bias factor
Ssrd_pd_reg = Ssrd_pd_reg./biasfac;
% Transform regionalized power density back to Mauersberger-Lowes
Ssrd_reg = Ssrd_pd_reg.*(ltot+1).*(2*ltot+1);
