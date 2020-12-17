function [invspecO,ls,rs]=invertCoreSpecMcLeod(lon,lat,rad,Br,rplanet,rsat,TH,rotcoord,meth,Lmax,Ltap,J,lrng,rstart)
  % [invspecO,ls,rs]=invertCoreSpecMcLeod(lon,lat,rad,Br,rplanet,rsat,TH,rotcoord,meth,Lmax,Ltap,J,lrng,rstart)
  %
  % Inverts for a power spectrum from raw data and fits it with a
  % McLeod spectrum to obtain source depth
  %
  % Note: When you use method 1, you could also use all three components
  % of the magnetic field and you could even simultaneosuly invert for an
  % external field. Would need to modify this code slightly for that.
  %
  % INPUT:
  %
  % data information:
  % lon          data locations longitude in degrees
  % lat          data locations latitude in degrees
  % rad          data locations radial positions
  % Br           radial component of your magnetic data
  % rplanet      radius of your planet / moon
  % rsat         average satellite altitude
  %
  % inversion parameters:
  % TH           semi-opening angle of your cap (degrees)
  % rotcoord     [lon, colat] of the center of your cap (degrees)
  % meth         LocTap (1) or TapLoc (2)?
  % Lmax         maximum spherical-harmonic degree for which you are inverting
  % Ltap         tapering bandwidth
  % J            number of altitude-cognizant Slepian functions to use
  % lrng         degrees to use for depth determination
  % rstart       Starting value for search for source radius. Shouldn't matter much
  %
  % OUTPUT:
  % invspecO    inverted power spectrum
  % ls          spehrical-harmonic degrees for the power spectrum
  % rs          estimated source radius from fitting McLeod spec
  %
  % Last modified by plattner-at-alumni.ethz.ch, 09/02/2020 
  

%%% Set the degrees
  ls = 0:Lmax;

  
%%%% Invert for local power spectrum
  if meth == 1
    % ... using Method LocTap
    coef = LocalIntField(Br(:),rad,(90-lat)*pi/180,lon*pi/180,TH,Lmax,J,rplanet,rsat);
    invspecO = localspectrum(coef2lmcosi(coef/sqrt(4*pi)),Ltap,TH,[],rotcoord,[],[],2,rplanet);% /spharea(TH,1);
    bias = false;   
  elseif meth == 2
    % ... using Method TapLoc
    invspecS_direct_dens =  localspectrumDataVarAlt2(Br, lon, lat, rad, rplanet, rsat, Lmax, Ltap, TH, J,[], rotcoord, [], 0)*spharea(TH,1)^2;
    invspecS_direct = invspecS_direct_dens.*(2*ls(:)+1).*(ls(:)+1)/rplanet^2*(rsat/rplanet)^2;
    invspecO = Simons2Olsen_spec(invspecS_direct);
    bias = true;
  else
    error('pick valid method')
  end
  

%%%% Estimate source depth
  rs = findDepthMinDiff_McLeod(invspecO,lrng,rplanet,rplanet,rstart,Ltap,Lmax,bias);
  

  
