%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: this script assumes that Rbar=1m, Bbar=1T to compute normalized 'Hat' quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(fullfile(pwd,'..','..','..','matlab'))% Add path to matlab directory to use functions from Matlab version of PERFECT
addpath(fullfile(pwd))

if exist('usePath')
  cd(usePath)
end

inputFileID = fopen(strcat(pwd,'/input.namelist'));
%inputArray = textscan(inputFileID, '%s %n', 'Delimiter', '=', 'CommentStyle', '!')
%celldisp(inputArray)
line = fgetl(inputFileID);
while ischar(line)
  split_line = strsplit(line,'=');
  if ~isscalar(split_line)
    % May be a parameter-value pair 
    if strcmp( strtrim(split_line{1}), 'psiDerivativeScheme' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      %value = strcat('int(',value,')')
      psiGridMode = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'psiMid' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      psiMid = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'psiDiameter' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      value = strrep(value,'d','e');
      psiDiameter = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'widthExtender' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      value = strrep(value,'d','e');
      widthExtender = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'leftBoundaryShift' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      value = strrep(value,'d','e');
      leftBoundaryShift = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'rightBoundaryShift' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      value = strrep(value,'d','e');
      rightBoundaryShift = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'NpsiPerDiameter' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      value = strrep(value,'d','e');
      NpsiPerDiameter = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'Ntheta' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      Ntheta = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'thetaDerivativeScheme' )
      value = strtrim(strsplit(split_line{2},'!'));
      value = value{1};
      thetaGridMode = str2num(value);
    end
    if strcmp( strtrim(split_line{1}), 'geometryFilename' )
      geometryFilename = strtrim(strsplit(split_line{2},'!'));
      geometryFilename = geometryFilename{1}(2:end-1); % (2:end-1) to get rid of " characters needed by Fortran
    end
  end
  line = fgetl(inputFileID);
end
fclose(inputFileID);

% Options for EFIT geometry reading
% .mat file must contain EFITFilename, topCropZ, bottomCropZ, innerCropR, outerCropR, polynomialFitDegreeForSmoothingEFITInPsi, numFourierModesInThetaToKeepInEFITGeometry and plotStuff
EFITOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate psi grid from options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Npsi = max(1,round(NpsiPerDiameter * (psiDiameter + 2*widthExtender - leftBoundaryShift + rightBoundaryShift))+1);
psiMin = psiMid - psiDiameter/2. - widthExtender + leftBoundaryShift;
psiMax = psiMid + psiDiameter/2. + widthExtender + rightBoundaryShift;

% Generate abscissae, quadrature weights, and derivative matrix for psi grid.
% Both abscissae and weights should be column vectors.
% Only really need the psi grid here
if Npsi>1
  switch psiGridMode
  case 1
    % Uniform grid, finite differences with a 3-point stencil
    scheme = 2;
    [psi, psiWeights, ddpsi, d2dpsi2] = differentiationMatricesForUniformGrid(Npsi, psiMin, psiMax, scheme);
  case 2
    % Uniform grid, finite differences with a 5-point stencil
    scheme = 12;
    [psi, psiWeights, ddpsi, d2dpsi2] = differentiationMatricesForUniformGrid(Npsi, psiMin, psiMax, scheme);
  otherwise
    error('Invalid psiGridMode!')
  end
else
  % If Npsi is 1, assume we are making the local approximation and so dI/dpsi is not needed, so ddpsi, d2dpsi2 are not needed
  psi = psiMid;
  psiWeights = 1.;
  ddpsi = 0.;
  d2dpsi2 = 0.;
end
psi = psi(:)';
psiWeights = psiWeights(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate theta grid from options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate abscissae, quadrature weights, and derivative matrix for theta grid.
% Both abscissae and weights should be row vectors.
switch thetaGridMode
case 0
  % Spectral uniform
  scheme = 20;
case 1
  % Uniform periodic 2nd order FD
  scheme = 0;
case 2
  % Uniform periodic 4th order FD
  scheme = 10;
otherwise
  error('Error! Invalid thetaGridMode')
end
[theta, thetaWeights, ddtheta, ~] = differentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load EFIT data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NPsi=1;
%psi = desiredPsi;
[thetaData, BData, BDotGradThetaData, IHat, qData, RData, as, R0, B0, psi0] = getGeometryFromEFITForSeveralFluxSurfaces(EFITFilename, psi, topCropZ, bottomCropZ, innerCropR, outerCropR, plotStuff);

%IHat = abs(IHat);
dIHatdpsi = (ddpsi * IHat')';

% Filter B(theta)
NThetaFilter = 100;
thetaForBFilter = linspace(0,2*pi,NThetaFilter+1);
thetaForBFilter(end)=[];
BModes = zeros(NThetaFilter, Npsi);
JModes = zeros(NThetaFilter, Npsi);
RModes = zeros(NThetaFilter, Npsi);
for psiIndex = 1:Npsi
  %temp = interp1(thetaData{psiIndex}, BData{psiIndex}/B0, thetaForBFilter, 'spline');
  temp = interp1(thetaData{psiIndex}, BData{psiIndex}, thetaForBFilter, 'spline');
  BModes(:,psiIndex) = fft(temp)/NThetaFilter;
  %temp = interp1(thetaData{psiIndex}, BDotGradThetaData{psiIndex}/B0, thetaForBFilter, 'spline');
  temp = interp1(thetaData{psiIndex}, BDotGradThetaData{psiIndex}, thetaForBFilter, 'spline');
  JModes(:,psiIndex) = fft(temp)/NThetaFilter;
  %temp = interp1(thetaData{psiIndex}, RData{psiIndex}/B0, thetaForBFilter, 'spline');
  temp = interp1(thetaData{psiIndex}, RData{psiIndex}, thetaForBFilter, 'spline');
  RModes(:,psiIndex) = fft(temp)/NThetaFilter;
end

numFourierModesInThetaToKeepInEFITGeometry = 5;

epsilon = as/R0;
Miller_A = 1./epsilon;
Miller_q = qData;
fprintf('Inverse aspect ratio derived from EFIT equilibrium: %g to %g\n',min(epsilon),max(epsilon))


BHat_beforeSmoothing = ones(Ntheta,1) * BModes(1,:);
JHat_beforeSmoothing = ones(Ntheta,1) * JModes(1,:);
RHat_beforeSmoothing = ones(Ntheta,1) * RModes(1,:);
keepUpDownAsymmetry = 1; % This variable should be either 1 or 0.
for m=1:numFourierModesInThetaToKeepInEFITGeometry
  BHat_beforeSmoothing = BHat_beforeSmoothing + 2*cos(m*theta)*real(BModes(m+1,:)) - keepUpDownAsymmetry*2*sin(m*theta)*imag(BModes(m+1,:));
  JHat_beforeSmoothing = JHat_beforeSmoothing + 2*cos(m*theta)*real(JModes(m+1,:)) - keepUpDownAsymmetry*2*sin(m*theta)*imag(JModes(m+1,:));
  RHat_beforeSmoothing = RHat_beforeSmoothing + 2*cos(m*theta)*real(RModes(m+1,:)) - keepUpDownAsymmetry*2*sin(m*theta)*imag(RModes(m+1,:));
end

% Smooth in the psi direction by fitting a polynomial:
BHat = zeros(Ntheta,Npsi);
JHat = zeros(Ntheta,Npsi);
RHat = zeros(Ntheta,Npsi);
dBHatDPsi = zeros(Ntheta,Npsi);
for itheta = 1:Ntheta
  if Npsi >= 5
    [p,S,mu] = polyfit(psi, BHat_beforeSmoothing(itheta,:), polynomialFitDegreeForSmoothingEFITInPsi);
    [y,~]= polyval(p,psi,S,mu);
    BHat(itheta,:) = y;
  elseif Npsi == 1
    BHat(itheta,:) = BHat_beforeSmoothing(itheta,:);
  else
    error('Invalid Npsi, must be 1 or greater than 5')
  end

  % Analytically differentiate the fitting polynomial:
  if Npsi >= 5
    [y,~]= polyval([0,polyder(p)],psi,S,mu);
    dBHatDPsi(itheta,:) = y / mu(2);
  elseif Npsi == 1
    dBHatDPsi(itheta,:) = zeros(Npsi);
  else
    error('Invalid Npsi, must be 1 or greater than 5');
  end

  if Npsi >= 5
    [p,S,mu] = polyfit(psi, JHat_beforeSmoothing(itheta,:), polynomialFitDegreeForSmoothingEFITInPsi);
    [y,~]= polyval(p,psi,S,mu);
    JHat(itheta,:) = y;
  elseif Npsi == 1
    JHat(itheta,:) = JHat_beforeSmoothing(itheta,:);
  else
    error('Invalid Npsi, must be 1 or greater than 5')
  end

  if Npsi >= 5
    [p,S,mu] = polyfit(psi, RHat_beforeSmoothing(itheta,:), polynomialFitDegreeForSmoothingEFITInPsi);
    [y,~]= polyval(p,psi,S,mu);
    RHat(itheta,:) = y;
  elseif Npsi == 1
    RHat(itheta,:) = RHat_beforeSmoothing(itheta,:);
  else
    error('Invalid Npsi, must be 1 or greater than 5')
  end
end

% Spectral uniform differentiation matrix:
scheme = 20;
[~, ~, ddthetaForBHat, ~] = differentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
dBHatdtheta = ddthetaForBHat * BHat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check potato orbit radius and print warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.17380166e-5.m^2.T^2=2*(Deuteron mass)*(1keV)/e^2 is rho^2*B^2.
fprintf('%s%f%s%f%s\n\n','The smallest minor radius is ',min(as),...
                   'm while the potato orbit radius for 1keV Deuterium is about ',...
                   (qData(1).^2*4.17380166e-5/B0.^2*R0).^(1./3.),'m.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save geometrical quantities into HDF5 file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geometryFilePath = strcat(pwd,'/',geometryFilename);
group = sprintf('%s%i%s%i%s','/Npsi',Npsi,'Ntheta',Ntheta,'/');
hdf5write(geometryFilePath, strcat(group,'psiMin'), psiMin)
hdf5write(geometryFilePath, strcat(group,'psiMax'), psiMax, 'WriteMode', 'append')
h5create(geometryFilePath,strcat(group,'BHat'),[Ntheta,Npsi])
h5write(geometryFilePath,strcat(group,'BHat'),BHat)
h5create(geometryFilePath,strcat(group,'dBHatdpsi'),[Ntheta,Npsi])
h5write(geometryFilePath,strcat(group,'dBHatdpsi'),dBHatDPsi)
h5create(geometryFilePath,strcat(group,'dBHatdtheta'),[Ntheta,Npsi])
h5write(geometryFilePath,strcat(group,'dBHatdtheta'),dBHatdtheta)
h5create(geometryFilePath,strcat(group,'JHat'),[Ntheta,Npsi])
h5write(geometryFilePath,strcat(group,'JHat'),JHat)
h5create(geometryFilePath,strcat(group,'IHat'),[Npsi])
h5write(geometryFilePath,strcat(group,'IHat'),IHat)
h5create(geometryFilePath,strcat(group,'dIHatdpsi'),[Npsi])
h5write(geometryFilePath,strcat(group,'dIHatdpsi'),dIHatdpsi)
% Write extra stuff into geometry files that is not used by PERFECT, but may be useful elsewhere
hdf5write(geometryFilePath, strcat(group,'R0'), R0, 'WriteMode', 'append')
h5create(geometryFilePath,strcat(group,'epsilon'),[Npsi])
h5write(geometryFilePath,strcat(group,'epsilon'),epsilon)
h5create(geometryFilePath,strcat(group,'RHat'),[Ntheta,Npsi])
h5write(geometryFilePath,strcat(group,'RHat'),RHat)
hdf5write(geometryFilePath,strcat(group,'psi0'),psi0, 'WriteMode', 'append')

quit
