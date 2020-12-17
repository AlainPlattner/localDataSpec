function coefs=loadEarthCoef(Lmax,datafile)
  % coefs=loadEarthCoef(Lmax,datafile)
  %
  % Loads coefficients from one of the CHAOS model files found on
  % https://www.space.dtu.dk/english/research/scientific_data_and_models/magnetic_field_models
  %
  defval('Lmax',20);  
  load(datafile)
  pp.dim = Lmax*(Lmax+2);
  cfs6 = reshape(pp.coefs,[],pp.pieces,pp.order);
  pp.coefs = reshape(cfs6(1:Lmax*(Lmax+2),:,:),[],pp.order);
  % Get coef by evaluating splines at time 2005
  % This is Jan 15 2005:
  t = 1841;
  coefs=[0;ppval(t,pp)];
