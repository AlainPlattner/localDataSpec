function Sl=specWiec(rs,cTH,d,Mag,rplanet,Lmax,Ltap,M)
  %  Sl=specWiec(rs,cTH,d,Mag,rplanet,Lmax,Ltap)
  %
  % Based on Wieczorek (2018) eq (27)
  % as well as Gong & Wieczorek (2021) eq (2)

  defval('M',[])
  defval('Ltap',0)

  cTH = cTH*pi/180;
  ls=(1:Lmax)';
  mu0 = 4*pi*1e7;

  Sl = mu0^2*Mag*(ls+1).*(ls+1).^2/12 .* (d/rplanet)^2 .* (rs/rplanet).^(2*ls+2) .* ...
       (0.5*(  intfun1(ls,cTH) + pval1(ls,cTH)  ).^2 + ...
         (  intfun0(ls,cTH) - pval1(ls,cTH) ).^2 );
       
  % Now localize
  ltot = (0:Lmax)';
  Sl = [0;Sl];

  if Ltap>0
    % Transform to density spectrum
     Sl_pd = Sl./(ltot+1)./(2*ltot+1);
    % Localize
    if isempty(M)
        Sl_pd_loc = localizeSpec(Sl_pd,Ltap)';
    else
        Sl_pd_loc = (Sl_pd'*M)';
    end
    % Transform back to ML spec
    Sl = Sl_pd_loc.*(ltot+1).*(2*ltot+1);
  end
  
end
  
function intvals=intfun0(ls,cTH)  
  intvals=nan(size(ls));
  x0=cos(cTH);
  for i=1:length(ls)
    intvals(i) = 2/(2*ls(i)+1)*legendreprodint(ls(i),0,1,0,x0,'gl');
  end
end

function intvals=intfun1(ls,cTH)  
  intvals=nan(size(ls));
  x0=cos(cTH);
  for i=1:length(ls)
    if ls(i)==1
      intvals(i) = 2/3 - x0 + x0^3/3;
    else
      intvals(i) = 4/(2*ls(i)-1)*legendreprodint(ls(i),1,1,1,x0,'gl');
    end
  end
end

function val=pval0(ls,cTH)
  Pl0 = sqrt(2./(2*ls+1)).*xlm(ls,0,cTH);
  val = Pl0*sin(cTH)^2./(ls+2);
end

function val=pval1(ls,cTH)
  % We have the condon-shortley phase. Remove it
  Pl1 = - sqrt(4./(2*ls+1)).*xlm(ls,1,cTH);
  val = Pl1*sin(cTH)*cos(cTH)./(ls+2);
end