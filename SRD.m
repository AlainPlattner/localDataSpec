function Ssrd_reg=SRD(rsource,robs,rref,Lmax,Ltap,M)
% SRD means Shell of Random Dipoles  
% This is equation (5b) from Voorhies, Sabaka, and Purucker (2002) 
% a=rref, the reference radius

defval('M', [])
  
%%% For a given rcore, calculate spec, normalized to start at 1
ltot = (0:(Lmax))';
Ssrd = [0;(ltot(2:end)).*(ltot(2:end)+0.5).*(ltot(2:end)+1).*(rsource/rref).^(2*ltot(2:end)-2).*(rref/robs).^(2*ltot(2:end)+4)];

% To regionalize, transform the SRD spectrum into power density
Ssrd_pd = Ssrd./(ltot(:)+1)./(2*ltot(:)+1);

% Regionalize power density spectrum
if isempty(M)
  Ssrd_pd_reg = localizeSpec(Ssrd_pd,Ltap)';
else
  Ssrd_pd_reg = (Ssrd_pd'*M)';
end

% Transform regionalized power density back to Mauersberger-Lowes
Ssrd_reg = Ssrd_pd_reg.*(ltot+1).*(2*ltot+1);
