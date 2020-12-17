function varargout=localspectrumDataVarAlt2(data,lon,lat,rad,rplanet,rsat,Lmax,Ltap,dom,Jcof,Jtap,rotcoord,Nspec,optn)
% [spec,specvar,spectap]=localspectrumDataVarAlt2(data,lon,lat,rad,rplanet,rsat,Lmax,Ltap,dom,Jcof,Jtap,rotcoord,Nspec,optn)
%
% Calculates the local multitaper spectrum using the formula described by
% Dahlen & Simons (2008), eq. 130, FOR ARBITRARILY SPACED DATA
%
% This one is different from ...VarAlt(..) as it directly inverts for the
% potential power spec and not the derivative power spec 
%
% INPUT:
%
% data      scalar values at locations lon, lat
% lon       longitudinal values of the data locations [in degrees, 0 to 360]
% lat       latitudinal values of the data locations [in degrees, -90 to 90]
% rad       radial position of the data
% rplanet   radius of the planet
% rsat      average radius of satellite data
% Lmax      Maximum sherical-harmonic degree for which to calculate the spectrum
% Ltap      Maximum spherical-harmonic degree for the multitapers with
%           which you want to calculate the local spectrum.
% dom       Named area or spherical cap semiopening angle
% Jcof      How many Slepian functions to invert for coefficients of
%           taper-data-product.
% Jtap      How many multitapers do you want to use for the spectral
%           estimation? [] for all (default)
% rotcoord  Center of spherical cap region if not northpole in degrees
%           [longitude colatitude], 0<=longitude<360, 0<=colatitude<=180
% Nspec     Global spectrum of the noise
% optn      0 (default) use division by (2*l+1) spectral normalization and 
%             divide by area*4*pi
%           1 use multiplication by (l+1) spectral normalization, no
%             division
%           2 normalize to obtain Mauersberger-Lowes spectrum
%           3 no normalization, just squaring and summation
%  
% OUTPUT:
%
% spec      Local power spectrum for provided spherical-harmonic degrees L
% specvar 	Error bars for the local spectrum
%
% Last modified by plattner-at-alumni.ethz.ch, 06/24/2020

defval('Jtap',[])
defval('rotcoord',[])
defval('onorout',[])
defval('Nspec',[])
defval('optn',0)
defval('mindeg',0)
%Lmin=0;
    
Lwid=Ltap;%2*Ltap+1;    

if ischar(dom)
   domarea=spharea(dom,0);
elseif length(dom)==1
   domarea=spharea(dom,1);
elseif length(dom)==2
   domarea=spharea(max(dom),1)-spharea(min(dom),1);
else
   error('Something is wrong with the domain')
end


if optn==1
    specnorm=1;
else
    specnorm=2;
end


% Make sure data is a column vector
data=data(:)';

% Get the multitaper coefficients
if ischar(dom) | isempty(rotcoord)
%  if length(dom)==2
%    [G,V]=glmalpharing(dom,Ltap);
%    [G2,V2]=glmalpharing(dom,Lmax);
  %  else
  [G,V]=glmalpha(dom,Ltap);
  % Do we need the upward-continued (but not derivative) slepian functions?
  % Or perhaps the normal scalar Slepian functions, because we undo the altitude with the alt-cog
  %[G,V]=glmalphapotup(dom,Ltap,rsat,rplanet);
  
  % Here we need the altitude-cognizant Slepian functions!
  %[G2,V2]=glmalpha(dom,[Lmin,Lmax]);
  %[G2,V2]=glmalphapotup(dom,[Lmin,Lmax],rsat,rplanet);
  %[G2,V2]=glmalphaup(dom,[Lmin,Lmax],rsat,rplanet);
  
    %  end
% Do the same below!!!
    
else
  [G,V]=glmalphapto(dom,Ltap,rotcoord(1),rotcoord(2),[]);
  %[G2,V2]=glmalphapto(dom,Lmax,rotcoord(1),rotcoord(2),[]);
  %[G2,V2]=glmalphapotuptoJp(dom,[Lmin,Lmax],rsat,rplanet,rotcoord(1),rotcoord(2),[],Jcof);
end

% If you don't want to sum all Slepian functions,
% then you need to order them
if ~ischar(dom)
  [V,isrt]=sort(V,'descend');
  G=G(:,isrt);
  %[V2,isrt2]=sort(V2,'descend');
  %G2=G2(:,isrt2);
end


% Define the number of Slepian functions for the tapering
if isempty(Jtap)
  Jtap=(Ltap+1)^2;
else
  Jtap=min(Jtap,(Ltap+1)^2);            
end

% For the g*d coefficient inversion, use the Shannon number
%defval('Jcof',round(domarea*((Lmax+1)^2-Lmin^2)));
% Make sure not larger than max number
%Jcof = min(Jcof,(Lmax+1)^2-Lmin^2);
Jcof = min(Jcof,(Lmax+1)^2);

% Evaluate the multitaper coefficients at the the data:
% First evaluate the spherical harmonics
% Remember: ylm doesn't have the (-1)^m phase, so just shift lon by 180
% The tapers don't start at Lmin!

% if Ltap==0
%   Y=ylm(0,0,(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1);   
%   Y=Y(:)';
%   G=full(G);
% else
%   Y=ylm([0 Ltap],[],(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1);        
% end

GJ=G(:,1:Jtap);
clear G

%%% Taper evaluation
%tapers=GJ'*Y;
tapers = evalManyCoefPotUpInt(GJ,pi*(90-lat)/180,pi/180*lon,ones(size(lat)),1,1);
% Evaluate tapers at data location including radial component
%tapers = evalManyCoefPotUpInt(GJ,pi*(90-lat)/180,pi/180*lon,rad,rplanet,1);
% Experimental: Try rad = rp/r to counteract the altitud variation
%tapers = evalManyCoefPotUpInt(GJ,pi*(90-lat)/180,pi/180*lon,rsat./rad,1,1);


% This is the product of the tapers times the data
gtimesd=tapers.*repmat(data,Jtap,1);

% Get the spectrum for the sum of tapers*data
spec=zeros(Lmax+1,1);
sumV=sum(V(1:Jtap));                            

%if Jcof>0
  % Evaluate Slepian functions for inversion:
%  Geval = evalManyCoefPotUpInt(G2(:,1:Jcof),pi*(90-lat)/180,pi/180*lon,rad,rplanet,1);  
%else
%  Y2=ylm([0 Lmax],[],(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1); 
%end

niter = 10; %0
coefs = LocalIntField_many(gtimesd',rad,(90-lat)*pi/180,lon*pi/180,...
                           dom,Lmax,Jcof,rplanet,rsat,rotcoord,[],[],niter);

for alpha=1:Jtap
  % Need to get spherical-harmonic coefficients of g*d. Use
  % Slepian function inversion
  if Jcof > 0
    %coef = LocalIntField(gtimesd(alpha,:)',rad,(90-lat)*pi/180,lon*pi/180,...
    %                     dom,Lmax,Jcof,rplanet,rsat,rotcoord);
    %lmcsgd = coef2lmcosi(coef,0);
    lmcsgd = coef2lmcosi(coefs{alpha},0);
  else
    % Sum
    %lmcsgd= coef2lmcosi(gtimesd(alpha,:)*Y2'/length(data)*4*pi*domarea,1);
    error('not implemented in this function')
  end  
  
  % This is the spectrum of each taper-data-product
  specprod = plm2spec(lmcsgd,specnorm);
  
  % If you want to see the individual tapered spectra: 
  if nargout>2
    spectap{alpha}=specprod;
  end
  % And normalize and sum them up
  spec=spec+V(alpha)/sumV*specprod;            
end        
           
if nargout<=2
  spectap=[];
end        


% Now get the error bars if you ask for them
if nargout>1
  specvar=mtvar(spec,(0:Lmax)',Lwid,dom);
else
  specvar=0;
end

switch optn
  case 2
    %% Normalize for Mauersberger-Lowes spectrum
    ls=(0:Lmax)';
    spec=spec.*(ls+1).*(2*ls+1);
    specvar=specvar.*(ls+1).*(2*ls+1);
  case 0
    %spec=spec*(rsat/rplanet)^2;
    %specvar=specvar*(rsat/rplanet)^2;
    %spec=spec/domarea*(rsat/rplanet)^2;%*(rsat/rplanet)^4;
    %specvar=specvar/domarea*(rsat/rplanet)^2;%*(rsat/rplanet)^4;
    spec=spec/domarea;%*(rsat/rplanet)^4;
    specvar=specvar/domarea;%*(rsat/rplanet)^4;
  case 3
    ls=(0:Lmax)';
    spec = spec.*(2*ls+1);
    specvar = specvar.*(2*ls+1);
end

% Now subtract the noise spectrum if provided
if ~isempty(Nspec)
  noiseLmax=length(Nspec)-1;
  M=mcouplings(Ltap,noiseLmax);
  disp('For the noise spectrum we just do eigenvalue weighted sum like in D.S. 2008 eq. (145)')
  spec=spec-M*Nspec(:);
end



  

varns={spec,specvar,spectap};
varargout=varns(1:nargout);

