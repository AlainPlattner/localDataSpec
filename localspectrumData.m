function varargout=localspectrumData(data,lon,lat,Lmax,Ltap,dom,Jcof,Jtap,rotcoord,Nspec,optn,Lmin)
% [spec,specvar,spectap]=localspectrumData(data,lon,lat,Lmax,Ltap,dom,Jcof,Jtap,rotcoord,Nspec,optn,Lmin)
%
% Calculates the local multitaper spectrum using the formula described by
% Dahlen & Simons (2008), eq. 130, FOR ARBITRARILY SPACED DATA
%
% INPUT:
%
% data      scalar values at locations lon, lat
% lon       longitudinal values of the data locations [in degrees, 0 to 360]
% lat       latitudinal values of the data locations [in degrees, -90 to 90]
% Lmax      Maximum sherical-harmonic degree for which to calculate the spectrum
% Ltap      Maximum spherical-harmonic degree for the multitapers with
%           which you want to calculate the local spectrum.
% dom       Named area or spherical cap semiopening angle
% Jcof      How many Slepian functions to invert for coefficients of
%           taper-data-product.
%           SET TO 0 FOR SUM INSTEAD OF INVERSION 
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
% Lmin      minimum degree (e.g. 1 for mag, 0 for grav)
%  
% OUTPUT:
%
% spec      Local power spectrum for provided spherical-harmonic degrees L
% specvar 	Error bars for the local spectrum
%
% Last modified by plattner-at-alumni.ethz.ch, 08/26/2019

defval('Jtap',[])
defval('rotcoord',[])
defval('onorout',[])
defval('Nspec',[])
defval('optn',0)
defval('mindeg',0)

    
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


% Evaluate the input coefficients on a grid
%[data,lon,lat]=plm2xyz(lmcosi);
% Already evaluated
data=data(:)';

% Get the multitaper coefficients
if ischar(dom) | isempty(rotcoord)
  if length(dom)==2
    [G,V]=glmalpharing(dom,Ltap);
    [G2,V2]=glmalpharing(dom,Lmax);
  else
    [G,V]=glmalpha(dom,Ltap);
    %[G2,V2]=glmalpha(dom,Lmax);
    %[G2,V2]=glmalpha(180,Lmax);
    [G2,V2]=glmalpha(dom,[Lmin,Lmax]);
  end
else
  [G,V]=glmalphapto(dom,Ltap,rotcoord(1),rotcoord(2),[]);
  [G2,V2]=glmalphapto(dom,Lmax,rotcoord(1),rotcoord(2),[]);
end

% If you don't want to sum all Slepian functions,
% then you need to order them
if ~ischar(dom)
  [V,isrt]=sort(V,'descend');
  G=G(:,isrt);
  [V2,isrt2]=sort(V2,'descend');
  G2=G2(:,isrt2);
end


% Define the number of Slepian functions for the tapering
if isempty(Jtap)
  Jtap=(Ltap+1)^2;
else
  Jtap=min(Jtap,(Ltap+1)^2);            
end

% For the g*d coefficient inversion, use the Shannon number
defval('Jcof',round(domarea*((Lmax+1)^2-Lmin^2)));
% Make sure not larger than max number
Jcof = min(Jcof,(Lmax+1)^2-Lmin^2);

% Evaluate the multitaper coefficients at the the data:
% First evaluate the spherical harmonics
% Remember: ylm doesn't have the (-1)^m phase, so just shift lon by
% 180
% The tapers don't start at Lmin!
if Ltap==0
  Y=ylm(0,0,(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1);   
  Y=Y(:)';
  G=full(G);
else
  Y=ylm([0 Ltap],[],(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1);
end
% evaluate the tapers
GJ=G(:,1:Jtap);
clear G
tapers=GJ'*Y;
% This is the product of the tapers times the data
gtimesd=tapers.*repmat(data,Jtap,1);

% Get the spectrum for the sum of tapers*data
spec=zeros(Lmax+1,1);
sumV=sum(V(1:Jtap));                            

if Jcof>0
  % Evaluate Slepian functions for inversion:
  Geval = evalManyCoef(G2(:,1:Jcof),pi*(90-lat)/180,pi/180*lon,1);
else
  Y2=ylm([0 Lmax],[],(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1); 
end


for alpha=1:Jtap
  % Need to get spherical-harmonic coefficients of g*d. Use
  % classical Slepian function inversion
  if Jcof > 0
    slepc=(Geval*Geval')\(Geval*gtimesd(alpha,:)');
    % Perhaps iteratively reweighted residuals to approximate
    % least absolute differences instead of least squares...
    slepc = itweighres(Geval,gtimesd(alpha,:)',slepc,10);
    % Project back to spherical harmonics
    lmcsgd = coef2lmcosi(G2(:,1:Jcof)*slepc,1);
  else
    % Sum
    lmcsgd= coef2lmcosi(gtimesd(alpha,:)*Y2'/length(data)*4*pi*domarea,1);
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
    spec=spec/domarea;
    specvar=specvar/domarea;
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

