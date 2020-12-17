function varargout=LocalIntField_many(data,rad,cola,lon,dom,Lmax,J,rplanet,avgsat,rotcoord,loadJ,savename,niter,weights,lambda,coef0)
% [coef,condi,dataweights,fnpl]=LocalIntField_many(data,rad,cola,lon,dom,Lmax,J,rplanet,avgsat,rotcoord,loadJ,savename,niter,weights,lambda,coef0)
%
% Calculates potential crustal field at radius rplanet from local satellite
% data within chosen region dom.
%
% Works for radial conponent OR vector data (figures out depending on the 
% input data).
%
% This version solves for many data sets evaluated at the same locations.
%
% INPUT:
%
% Required:
%
% data      EITHER radial component data:
%              Matrix of data column vectors for Br for the points given as rad,cola,lon
%           OR vectorial data:
%              given as data{1} = rad component matrix of column vectors,
%              data{2} = colat component matrix of column vectors,
%              data{3} = lon component matrix of column vectors, 
%              or data is a matrix of size 3*length(rad) x number of data vectors, 
%              ordered as [rad;cola;lon]
% rad       radial position of satellite (planet radius + altitude)
% cola,lon  colatitude/longitude positions of the data values (both as
%           column values), 0<=cola<=pi; 0<=lon<=2pi
% dom       Slepian concentration target region either given as string 
%           (one of the region names)
%           OR integer (polar cap opening angle in degrees)
%           OR two polar cap opening angles for the ring in between them
%           OR [lon lat] an ordered list defining a closed curve [degrees]
%           OR several regions to add up/subtract (ONLY FOR VECTORIAL DATA): 
%           struct with
%           dom.name    for name of the combined region
%           dom.parts   for the cell array of names of the parts, 
%                       or cap opening angles
%                       or [cap,lon,colat] for rotated caps
%           dom.sign    for adding or subtracting 
%           Example: dom.parts{1}='namerica'; dom.parts{2}='samerica';
%                    dom.sign=[1,1]; dom.name='americas';
%                    dom.name='weirdRing'
%                    dom.parts{1}=30; dom.parts{2}=[5,5,10]; dom.sign=[1,-1]
%                    subtracts the ring of cTH=5, clon=5, ccola=10 from the
%                    larger polar cap
% Lmax      Maximum spherical harmonic degree
% J         How many Slepian functions should be used to calculate the
%           solution? More means more sensitive to noise but higher spatial
%           resolution
% rplanet   Planet radius to which the solution should be calculated
%
% Recommended:
%
% avgsat    Average satellite altitude for which the Slepian functions
%           should be calculated. The crustal magnetic field is of course 
%           calculated taking the true radial positions into account.
%
% Optional, depending on what you want to do:
%
% rotcoord  This option is needed if you work with spherical caps or rings
%           that are not centered at the north pole. Provide center as
%           [longitude colatitude] IN DEGREES [0<=360 0<=180].
% loadJ     For polar caps/rings that are not centered at north pole: You
%           might only want to calculate and save the first "loadJ"
%           evaluated Slepian functions. Default: loadJ=J
% savename  To make future calculations faster: Povide a name to save the 
%           evaluated Slepian functions. Default: no saving
% niter     Numer of iterations for iteratively reweighted residual
% weights   Weights you want to use in least squares, it is a vector, 
%           same size as either number of data locations or 3*number of 
%           data locations if you want to weight each component: 
%           [wBr wBth wBph]
% lambda    Damping factor [default = 0]
% coef0     reference model for damping: 
%           min( ||A*coef - dat||^2 + lambda*||coef-coef0||^2 )
%           [default = zeros]
% DAMPING CURRENTLY ONLY FOR RADIAL-ONLY INVERSION
%
% OUTPUT:
%
% coef          spherical-harmonic coefficients of the potential field
%               in ADDMON format
% condi         Conditioning number of the square matrix of evaluated 
%               Slepian functions
% dataweights   Final weights used in the iteratively reweighted least 
%               squares solution
% fnpl          path to saved evaluated matrix file
%
% Last modified by plattner-at-alumni.ethz.ch,  07/14/2020

defval('avgsat',[])
defval('rotcoord',[0 0])
defval('loadJ',J)
defval('savename',[])
defval('niter',10)
defval('weights',[])
defval('lambda',0)
defval('coef0',[])

if length(data)==3
    data=[data{1};data{2};data{3}];
end

if isempty(avgsat)
    avgsat=mean(rad);
end

% If we are working with radial derivative data only:
if size(data,1)==length(rad)
  if length(Lmax)>1
    warning('Bandpass for only radial data not yet implemented')
  end
  % Polar caps/rings?
    if (~ischar(dom)) & (rotcoord(1)~=0 || rotcoord(2)~=0)
        [G,~]=glmalphauptoJp(dom,Lmax,avgsat,rplanet,...
            rotcoord(1),rotcoord(2),0,loadJ);        
    else
        [G,~]=glmalphaup(dom,Lmax,avgsat,rplanet);
    end
    
    % Now, if available, load the previously calculated evaluated 
    % Slepian functions up to J=loadJ. For this make the dataname:
    filoc=fullfile(getenv('IFILES'),'SLEPEVAL');
    
    if length(Lmax)==1
    if ischar(dom)
        fnpl=sprintf('%s/%s-rad-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-rad-TH%i-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-rad-TH%i_%i-L%i-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom(1),dom(2),Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)>2
      if exist('octave_config_info')
	nam=builtin('hash','sha1',dom);
      else
	nam=hash(dom,'sha1');
      end
      fnpl=sprintf('%s/%s-rad-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
		   savename,nam,min(Lmax),avgsat,rplanet,loadJ);      
    end
    else
    if ischar(dom)
        fnpl=sprintf('%s/%s-rad-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-rad-TH%i-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-rad-TH%i_%i-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom(1),dom(2),min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)>2
       if exist('octave_config_info')
	nam=builtin('hash','sha1',dom);
      else
	nam=hash(dom,'sha1');
       end
       fnpl=sprintf('%s/%s-rad-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
		    savename,nam,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    end
    end
    
    
    % .. and load it
    if exist(fnpl,'file')==2 && ~isempty(savename)
        load(fnpl)
        fprintf('%s loaded by LocalInnerField\n',fnpl)
    else
        % We need to calculate it  
        MloadJ=rGscal(G(:,1:loadJ),cola,lon,rad,rplanet,1);
        % And save it if dataname is provided
        if ~isempty(savename) 
            % See if it's Octave or Matlab
            if (exist('octave_config_info'))
                % It's octave
                save(fnpl,'MloadJ')
            else
                % It's Matlab
                save(fnpl,'MloadJ','-v7.3')
            end
        end
    end
    % Now fit the data with the evaluate Slepian functions using
    % iteratively reweighted least squares.
       
    % First: Only select the Slepian functions we want:
    M=MloadJ(1:J,:);
    
    %%%% Loop over the different data given as column vectors
    for dati = 1:size(data,2)
    
    
        if ~isempty(weights)
            % Here I am multiplying the matrix with the weights
            M=M.*repmat(weights(:)',J,1);
            % and here the right hand side.
            data(:,dati)=weights(:).*data(:,dati);
        end    
        
        % Prepare damping if requested:
        lambda=abs(lambda);
        if isempty(coef0)
            slepcoef0 = zeros(J,1);
        else
            % transform the provided spherical-harmonic reference 
            % coefficients into Slepian reference coefficients
            slepcoef0 = (G(:,1:J))'*coef0;
        end

        % Start with a least squares step:
        MM=M*M' + lambda*eye(J);
        Md=M*data(:,dati) + lambda*slepcoef0;
        slepcoef=MM\Md;

        % And now iteratively reweighted residual calculation
        [slepcoef,dataweights]=itweighres(M,data(:,dati),slepcoef,niter,lambda,slepcoef0);
    
    
        % Turn the coefficients into spherical-harmonic coefficients: 
        coef{dati}=G(:,1:J)*slepcoef;       
   
    end
        
    % We are now done for the radial derivative case.

    % Now the full gradient case: This is almost the same:
elseif size(data,1)==3*length(rad)
    % Vectorial case  
    % First load or calculate the Slepian function coefficients
    if isstruct(dom)
        [H,~]=gradvecglmalphaup(dom,Lmax,avgsat,rplanet);
    elseif (~ischar(dom)) & (rotcoord(1)~=0 || rotcoord(2)~=0)        
        [H,~]=gradvecglmalphauptoJp(dom,Lmax,avgsat,rplanet,...
           rotcoord(1),rotcoord(2),0,loadJ);       
    else
        [H,~]=gradvecglmalphaup(dom,Lmax,avgsat,rplanet);
    end
    
    % Now, if available, load the previously calculated evaluated Slepian
    % functions up to J=loadJ
    filoc=fullfile(getenv('IFILES'),'SLEPEVAL');
    if length(Lmax)==1
    if ischar(dom)
      fnpl=sprintf('%s/%s-vec-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
                   savename,dom,Lmax,avgsat,rplanet,loadJ);
    elseif isstruct(dom)
      fnpl=sprintf('%s/%s-vec-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
                   savename,dom.name,Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && ~isstruct(dom)  && length(dom)==1
        fnpl=sprintf('%s/%s-vec-TH%i-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && ~isstruct(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-vec-TH%i_%i-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom(1),dom(2),Lmax,avgsat,rplanet,loadJ);    
    end
    else
    if ischar(dom)
        fnpl=sprintf('%s/%s-vec-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-vec-TH%i-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-vec-TH%i_%i-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom(1),dom(2),min(Lmax),max(Lmax),avgsat,rplanet,loadJ);    
    end
    end
    
    if exist(fnpl,'file')==2 && ~isempty(savename)
        load(fnpl)
        fprintf('%s loaded by LocalInnerField\n',fnpl)

    else
        % We need to calculate it  
        MloadJ=rGvec(H(:,1:loadJ),cola,lon,rad,rplanet,1);
        % And save it if dataname is provided
        if ~isempty(savename)
            % See if it's MAtlab or Octave
            if (exist('octave_config_info'))
                % It's Octave
                save(fnpl,'MloadJ')
            else
                % Matlab
                save(fnpl,'MloadJ','-v7.3')
            end
        end                
    end 
    
    % Now fit the data with the evaluate Slepian functions using
    % iteratively reweighted least squares.
       
    % First: Only select the Slepian functions we want:
    M=MloadJ(1:J,:);
        
    %%%% Loop over the different data given as column vectors
    for dati = 1:size(data,2)
    
        if ~isempty(weights)
            % if you want to weight each point, 
            % instead of each component of each point:
            %t1=tic;
            if length(weights)~=size(data,1)
                weights=[weights(:);weights(:);weights(:)];
            end 
        
            % Here I am reweighing the data points
            % Here the matrix
            M=M.*repmat(weights',J,1);                
            % Here the right hand side
            data(:,dati)=weights.*data(:,dati);
        end        
        
        % Start with a least squares step:
        MM=M*M';
        Md=M*data(:,dati);
        slepcoef=MM\Md;

        % And now iteratively reweighted residual calculation
        [slepcoef,dataweights]=itweighres(M,data(:,dati),slepcoef,niter);    
    
        % Turn the coefficients into spherical-harmonic coefficients: 
        coef{dati}=H(:,1:J)*slepcoef;   
        
    end
       
else
    error('Something is not right with the provided data')        
end

% Coefs are in ADDMOUT. Transform to ADDMON:
for i=1:size(data,2)
    coef{i}=out2on(coef{i},max(Lmax));

if nargout>2
    condi=cond(MM);
else
    condi=[];
end   


varns={coef,condi,dataweights,fnpl};
varargout=varns(1:nargout);

end

