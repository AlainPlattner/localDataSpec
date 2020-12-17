function Smc_reg=McLeod(rsource,robs,rplanet,Lmax,Ltap,bias)
% To incorporate the bias created by tapering the radial derivative,
% then get the potential-field coefficients, then power spec, then
% "undoing" the (l+1)^2 factor, set bias=true
defval('bias',false)
  
%%% For a given rcore, calculate spec, normalized to start at 1
ltot = (0:(Lmax))';
Smc = [0;(ltot(2:end)+0.5)./(ltot(2:end).*(ltot(2:end)+1)).*((rsource/robs).^(2*ltot(2:end)+4))];
% To regionalize, transform the McLeod spectrum into power density
Smc_pd = Smc./(ltot(:)+1)./(2*ltot(:)+1);
if bias
  biasfac = (ltot+1).^2.*(robs/rplanet).^(-2*ltot-4);
  Smc_pd = Smc_pd.*biasfac(:);
else
  biasfac = 1;
end
% Regionalize power density spectrum
%Smc_pd_reg = localizeSpec([0;Smc_pd],Ltap)';
Smc_pd_reg = localizeSpec(Smc_pd,Ltap)';
% Undo the bias factor
Smc_pd_reg = Smc_pd_reg./biasfac;
% Transform regionalized power density back to Mauersberger-Lowes
Smc_reg = Smc_pd_reg.*(ltot+1).*(2*ltot+1); 





