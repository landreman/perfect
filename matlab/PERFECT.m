function PERFECT()
% PERFECT:
% Petestal and Edge Radially-global Fokker-Planck Evaluation of Collisional Transport
%
% Originally written in 2013 by Matt Landreman,
% MIT Plasma Science & Fusion Center

% **********************************************************************
% Program flow options
% **********************************************************************

runMode = 0;
% -1 = plot the input profiles, then stop. (For this option, the kinetic equation is not solved.)
%  0 = single run
%  1 = convergence scan
%  2 = convergence scan in local approximation
%  3 = scan dT/dpsi
%  4 = scan velocity-space structure of sources
%  5 = compare 3 different options for radial boundary conditions
%  6 = scan desiredFWHMInRhoTheta
%  7 = scan desiredU
%  8 = scan Nx
%  9 = scan psiDiameter
% 10 = scan exponent at fixed U
% 11 = scan U, and within that scan, scan exponent
% 12 = scan widthMultiplier
% 13 = scan Ntheta
% 14 = scan NpsiPerDiameter

%saveResults = true;
saveResults = false;
% If true, the matlab workspace is saved to a .mat file upon program completion.

filenameNote='test';
% This string is appended to the name of the workspace file which is
% generated when saveResults = true.

% **********************************************************************
% Geometry options
% **********************************************************************

geometry = 1;
% 0 = concentric circular flux surfaces
% 1 = Miller geometry
% 2 = Boozer poloidal angle, for comparison with SFINCS or other stellarator neoclassical codes
% 3 = up-down asymmetric
% 4 = EFIT. When using this option, the normalization length \bar{R} and
%     the normalization field \bar{B} are 1 meter and 1 Tesla respectively.


epsilon = 0.3;

% Beginning of parameters that only matter when geometry=1:
Miller_A = 3.17;
Miller_kappa = 1.66;
Miller_delta = 0.416;
Miller_s_delta = 1.37;
Miller_s_kappa = 0.70;
Miller_dRdr = -0.354;
Miller_q = 4;
% End of parameters that only matter when geometry=1:

% Beginning of parameters that only matter when geometry=3:
epsilonUDB = 0.0;
epsilonUDJacobian = 0;
% End of parameters that only matter when geometry=3:

% Beginning of parameters that matter only when geometry=4:
EFITFilename = 'C:\Users\landreman\Documents\MATLAB\DIII-D data\g145781.02141_583';  
topCropZ = 1;    
bottomCropZ = -1;

polynomialFitDegreeForSmoothingEFITInPsi = 4;
numFourierModesInThetaToKeepInEFITGeometry = 5;
% End of parameters that matter only when geometry=4:

% **********************************************************************
% Parameters that determine the input profiles 
% (n and T for each species, along with the electrostatic potential.)
% **********************************************************************

radialProfiles = 0;
% 0 = simple: no radial variation of I, B, or J
% 1 = more complicated
% 2 = profiles for comparison with analytic theory
% 3 = constant dPhi/dpsi
% 4 = Gaussian dPhi/dpsi, at constant original (wrong) analytic heat flux
% 5 = C-Mod like profiles
% 6 = Gaussian dPhi/dpsi, at constant corrected analytic heat flux
% 7 = Profiles for comparison with analytic theory, but keep analytic heat
%     flux constant.
% 8 = C-Mod data from Michael Churchill

CModDataFile = 'C:\Users\landreman\Documents\MATLAB\Churchill asymmetry measurements\hmode_1120803011_90_112.csv';

CModProfileStretchFactor = 4;

T_ped = 0.4;%0.8
T_wall = 0.4; %0.075;
T_Delta = 0.05;
T_linear = 2.8;

T_ped_2 = 0.4;
T_wall_2 = 0.5;
T_Delta_2 = 0.013;

n_ped = 0.9;
n_wall = 0.3;
n_Delta = 0.018;
%n_Delta = T_Delta;
n_linear = 1;

% Er offset and depth are in kV/m
Er_offset=0; %25
Er_depth = 50; %80
Er_Delta = n_Delta;
%Er_Delta = 0.01;

%forceElectrostaticConfinement = true;
forceElectrostaticConfinement = false;

linearCorner = 50;

psiMid = 0.97;
%psiMid = 0.85;
psiInnerExtender = 0.0;
%psiInnerExtender = 0.13;

setTPrimeToBalanceHeatFlux = true;
%setTPrimeToBalanceHeatFlux = false;

widthExtender = 0.0;
widthExtenders = linspace(1.8, 0.6, 7);

desiredU = 0.8;
desiredUs = linspace(0,1,16); %desiredUs(1)=[];
desiredFWHMInRhoTheta = 2; %2
%exponent = 0.994493;
exponent = 2;
%dTHatdpsiScalar = -0.7; %For U=0.7
%dTHatdpsiScalar = -0.4; %For U=1
dTHatdpsiScalar = -6; %0.9
detaHatdpsiScalar = 3;

% **********************************************************************
% Physics parameters
% **********************************************************************

Delta = 0.0006;
omega = 0.0014;
psiAHat = 0.03;

nu_r = 0.1;

%Delta = 0.001;
%Delta = 0.0006*(epsilon/0.01)^2/(desiredFWHMInRhoTheta/0.24)/30;
%omega=Delta;

% If BBar = 1T, RBar = 1m, mBar = m_Deuterium, and PhiBar = 1 kV:
%Delta = 0.00646;
%omega = 0.00323;
%nu_r = 0.00832;

%psiAHat = 0.1102; % for C-Mod.
%minorRadius = 0.2; % in meters. Used to put radial electric field in kV/m.
%minorRadius = 0.56; % in meters. Used to put radial electric field in kV/m.

% DIII-D parameters:
%Delta = 0.00193;
%omega = 0.000967;
%psiAHat = 0.0359;
minorRadius = 0.56; % in meters. Used to convert kv/m to psi.


%Delta = 0.0011;
%omega = 0.0014;
%Delta = 0.0003;
%Delta = 0.000001;
%Delta=0;
%omega = 0;

%psiAHat = 0.05;
%psiAHat = 3.333333333333333e-5;
%psiAHat = 3.333333333333333e-6;
%psiAHat = epsilon*epsilon/Miller_q;
%psiAHat = epsilon*epsilon*3;

%{
%nuPrimeMeaningful = 0.065;
%nuPrimeMeaningful =  0.0224719101123596d+0;
%nuPrimeMeaningful =  0.010;
nuPrimeMeaningful =  3e-5;
%nuPrimeMeaningful = 0.1;
%nuPrimeMeaningful = 1/178*4;
%nuPrimeMeaningful = 0.03 / (Miller_A^1.5);
nuPrime = nuPrimeMeaningful/Miller_q;
nuStar = nuPrime * (Miller_A^1.5);
%nuStar = 2;
%nuPrime = nuStar/(Miller_A^1.5);
%}

% Collisionality at the reference parameters:
% (equal to nuPrime / Miller_q)
%nu_r = 0.0139;
%nu_r = 0.0229666666666667;
%nu_r = 0.0216666666666667;
%nu_r = 0.0033333333333333333333;
%nu_r = 0.01 / Miller_q;
%nu_r = epsilon^(3/4) /  Miller_q
%nu_r = 1/30;

% The following parameter is used only to convert the flow to dimensional
% units:
vBar_kmPerSec = 310;


% **********************************************************************
% Species parameters
% **********************************************************************

% Here is an example for 1 species:
charges = 1;
masses = 1;
scalarNHats = 1;
scalarTHats = 1;

% Here is an example for 2 species:
%{
charges = [1, -1];
masses = [1, 2.72e-4];
scalarNHats = [1, 1];
scalarTHats = [1, 1];
%}


% **********************************************************************
% Options for terms in the kinetic equation
% **********************************************************************

makeLocalApproximation = true;
%makeLocalApproximation = false;
% If makeLocalApproximation is set to true, the conventional radially-local
% kinetic equation is solved instead of the radially global kinetic equation.
  
includeNewStreamingAndMirrorTerms = true;
%includeNewStreamingAndMirrorTerms = false;

includeddpsiTerm = true;
%includeddpsiTerm = false;

% **********************************************************************
% Options for the particle and energy sources
% **********************************************************************

sourcePoloidalVariation = 0;
% 0 = Souces are constant in theta
% 1 = Souces are \propto 1 + cos(theta)
% 2 = Souces are \propto 1 + source_a * cos(theta)
% 3 = Souces are \propto 1 + source_a * sin(theta)

source_a = 0.8;
% This parameter only matters when sourcePoloidalVariation = 2 or 3.

% **********************************************************************
% Numerical resolution parameters
% **********************************************************************

% This value, when multiplied by the range of normalized psi included in the simulation,
% gives the number of grid points in the radial direction.
% Memory and time requirements DO depend strongly on this parameter.
NpsiPerDiameterConverged = 100;

psiDiameterConverged = 0.25;

% Number of grid points in the poloidal direction.
% Memory and time requirements DO depend strongly on this parameter.
NthetaConverged = 11;

% Number of Legendre polynomials used to represent the distribution function.
% Memory and time requirements DO depend strongly on this parameter.
% The value of this parameter required for convergence depends strongly on
% the collisionality. At high collisionality, this parameter can be as low
% as ~ 5. At low collisionality, this parameter may need to be many 10s or
% even > 100 for convergence.
NxiConverged = 13; 

% Number of Legendre polynomials used to represent the Rosenbluth
% potentials. Except in exceptional circumstances, this number should be 4.
% Memory and time requirements do NOT depend strongly on this parameter. 
NLConverged = 4;

% Number of grid points in energy used to represent the distribution function.
% Memory and time requirements DO depend strongly on this parameter.
% This parameter almost always needs to be at least 5.
% Usually a value in the range 5-8 is plenty for convergence.
NxConverged = 6;

% Number of grid points in energy used to represent the Rosenbluth potentials.
% Memory and time requirements do NOT depend strongly on this parameter.
NxPotentialsPerVthConverged = 30;

% Maximum normalized speed for the Rosenbluth potential grid.
% Memory and time requirements do NOT depend strongly on this parameter.
% Typically a value of 5 is good.
xMax = 5;

% Tolerance used to define convergence of the Krylov solver.
% This parameter does not affect memory requirements but it does affect the
% time required for solution somewhat.
log10tolConverged = 5;

% Below are some parameters that need to be documented...

Nthetas = floor(linspace(15,40,6));
NLs = 2:2:4;
Nxis = floor(linspace(10,20,5));
Nxs=6:9;
NxPotentialsPerVths = linspace(5,50,9);

psiDiameters = psiDiameterConverged * logspace(0,0.5,2);
NpsiPerDiameters = [100,120,150];
NpsiIntervals = 1;


% **********************************************************************
% Other numerical options
% **********************************************************************

tryIterativeSolver = true;
%tryIterativeSolver = false;
% If this parameter is true, a Krylov solver is attempted.
% If this parameter is false, a sparse direct solver is used instead.
% You almost always want tryIterativeSolver = true,
% except perhaps for low-resolution runs in which the Krylov solver is not
% converging well.

tryDirectSolverIfIterativeSolversFail = false;

orderOfSolversToTry = [1, 3, 4, 5];
% 1 = GMRES
% 2 = BiCGStab
% 3 = BiCGStab(l)
% 4 = TFQMR
% 5 = CGS

% Below are some settings for the Krylov solvers.
% These settings are irrelevant if tryIterativeSolver = false.
maxIterations = 800;
restart = maxIterations; % Used only for GMRES.

psiGridMode = 2;
% This option determines the discretization scheme for the radial
% coordinate.
% 0 = Chebyshev grid and differentiation
% 1 = Uniform grid, finite difference derivatives with a 3-point stencil.
% 2 = Uniform grid, finite difference derivatives with a 5-point stencil.
% This parameter should almost always be 2.

% Note: the setting below only matters when preconditionerMethod_psi=1:
psiDerivativeForPreconditioner = 0;
% 0 = centered 2nd order finite difference
% 1 = 2-point upwinded difference

thetaGridMode = 2;
% 0 = uniform periodic spectral
% 1 = 2nd order uniform finite-difference
% 2 = 4th order uniform finite-difference

xGridScheme = 0;
% 0 = New polynomials
% 1 = Chebyshev on [0.1, 5]
% 2 = uniform, 5 point stencil
% 3 = uniform, 5 point stencil, set last 2 points of f to zero
% 4 = uniform with exp(x^2) weight, no boundary condition at xMax.

xMaxForDistribution = 4;

%imposeLocalSolutionWhereTrajectoriesAreParallel = true;
imposeLocalSolutionWhereTrajectoriesAreParallel = false;
tiny=1e-12;

% Adds dtheta/2 to all theta grid location
%shiftTheta = true;
shiftTheta = false;

upwindingSign = -1;

% **********************************************************************
% Preconditioner options
% **********************************************************************

preconditioner_species = 1;
% 0 = keep full species coupling
% 1 = drop all cross-species coupling

preconditionerMethod_x = 1;
% 0 = keep full x coupling
% 1 = keep only diagonal in x
% 2 = keep upper-triangular part in x
% 3 = keep lower-triangular part in x

preconditioner_x_min_L = 0;
% This is the lowest value of L at which the simplified x coupling is used
% in the preconditioner.
% Set it to 0 to use the simplified x coupling for all L in the
% preconditioner.

preconditionerMethod_psi = 0;
% 0 = keep full psi coupling
% 1 = use less accurate ddpsi
% 2 = drop ddpsi term, even at boundaries
% 3 = keep only diagonal of ddpsi term, which is only nonzero near the boundaries
% 4 = drop all global terms (those proportional to Delta or omega)

preconditionerMethod_theta = 0;
% 0 = keep full theta coupling
% 1 = use 3-point finite difference stencil

preconditionerMethod_xi = 0;
% 0 = keep full xi coupling
% 1 = drop terms that are +/- 2 from the diagonal in xi, so preconditioner
%     is tridiagonal in xi.

% **********************************************************************
% Plotting options:
% **********************************************************************

figureOffset=30;

drawIntroFigures = true;
%drawIntroFigures = false;

drawOutputFigures = true;
%drawOutputFigures = false;

%plotVelocitySpaceGrid = true;
plotVelocitySpaceGrid = false;

legendSize = 5;
% Font size for the legend in some plots.

% **********************************************************************
% **********************************************************************
%
% End of the main input parameters
%
% **********************************************************************
% **********************************************************************



% **********************************************************************
% Validate some of the input parameters:
% **********************************************************************

numSpecies = numel(masses);
if numSpecies ~= numel(charges)
    error('The number of species charges must equal the numer of species masses.')
end
if numSpecies ~= numel(scalarTHats)
    error('The number of species temperatures must equal the numer of species masses.')
end
if numSpecies ~= numel(scalarNHats)
    error('The number of species densities must equal the numer of species masses.')
end

if tiny<0
    error('tiny must be positive')
end

% **********************************************************************
% End of input parameter validation.
% **********************************************************************

dBdthetaResolutionMultiplier = 10;

% Compute a few quantities related to the Miller equilibrium:
Miller_x = asin(Miller_delta);
Miller_Q = Miller_kappa / (2*pi*Miller_A) * quadgk(@QIntegrand, 0, 2*pi);


q=0;
kq=0;
particleFlux=0;
density=0;
theta=0;
theta1D=linspace(0,2*pi,100);
iteration=0;
kParOutboard=0;
nOutboard=0;
psi = 0;
x_i=0;
ddpsi=0;
ddtheta=0;
ddx=0;
d2dx2=0;
soln=0;
BHat=0;
JHat=0;
dBHatdtheta=0;
dBHatDPsi=0;
dPhiHatdpsi=0;
THat=0;
dTHatdpsi=0;
IHat=0;
dIHatdpsi=0;
etaHat=0;
nHat=0;
deltaTMax = 0;
deltaTMin = 0;
deltanMax = 0;
deltanMin = 0;
nuPrimeLocal=0;
nuStarLocal=0;
residual=0;
didItConverge = true;
FSAKPar = 0;
heatFluxRelativeToLocalPlateauRegime =0;
FSAB2=0;
densityPerturbation=0;
potentialPerturbation=0;
flow=0;
pressurePerturbation=0;
temperaturePerturbation=0;
momentumFlux=0;
nHats=0;
THats=0;
etaHats=0;
detaHatdpsis=0;
dTHatdpsis=0;
PhiHat=0;
dPhiHatdpsi=0;
deltaTs=0;
deltans=0;
deltaetas=0;

if plotVelocitySpaceGrid
    figure(figureOffset+7)
    clf
end

switch runMode
    case {0,-1}
        % Single run
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('nu_r = %g\n',nu_r);
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        solveDKE();
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
    case 1
        % Convergence scan
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning convergence scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        scanStartTime = tic;
        
        baseParameters = [NpsiPerDiameterConverged; NthetaConverged; NxiConverged; NxConverged; NLConverged; NxPotentialsPerVthConverged; false; psiDiameterConverged; log10tolConverged];
        numRunsInScan = numel(baseParameters)+2;
        
        parameters = repmat(baseParameters, [1, numRunsInScan]);
        parameters(7, 1) = true;
        parameters(2, 3) = 2*NthetaConverged+1;
        parameters(3, 4) = 2*NxiConverged;
        parameters(4, 5) = NxConverged+1;
        parameters(4, 6) = floor(1.5*NxConverged);
        parameters(5, 7) = 2*NLConverged;
        parameters(6, 8) = 2*NxPotentialsPerVthConverged;
        parameters(1, 9) = 2*NpsiPerDiameterConverged;
        if geometry==4
            parameters(8, 10) = 0.8*psiDiameterConverged;
        else
            parameters(8, 10) = 0.8*psiDiameterConverged;
        end
        parameters(9, 11) = log10tolConverged + 1;
        legendText={'local approx','Base case','2x N\theta', '2x N\xi', 'Nx + 1','1.5x Nx', '2x NL', '2x NxPotentials','2x N\psi PerDiameter','2x \psi Diameter','tol/10'};
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        FSABFlowScan=cell(numRunsInScan,1);
        flowOutboardScan=cell(numRunsInScan,1);
        flowInboardScan=cell(numRunsInScan,1);
        heatFluxScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        particleSourceScan=cell(numRunsInScan,1);
        heatSourceScan=cell(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            fprintf('\n\nRun %d of %d: %s\n\n',runNum, numRunsInScan, legendText{runNum})
            NpsiPerDiameter = parameters(1, runNum);
            Ntheta = parameters(2, runNum);
            Nxi = parameters(3, runNum);
            Nx = parameters(4, runNum);
            NL = parameters(5, runNum);
            NxPotentialsPerVth = parameters(6, runNum);
            makeLocalApproximation = parameters(7, runNum);
            psiDiameter = parameters(8, runNum);
            tol = 10^(-parameters(9,runNum));
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            flowOutboardScan{runNum} = flowOutboard;
            flowInboardScan{runNum} = flowInboard;
            FSABFlowScan{runNum} = FSABFlow;
            heatFluxScan{runNum} = heatFlux;
            particleFluxScan{runNum} = particleFlux;
            particleSourceScan{runNum} = particleSource;
            heatSourceScan{runNum} = heatSource;
            didItConvergeScan(runNum) = didItConverge;
        end
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.9,0.6,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.9,0.6,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'--','-','-','-','-','-',':',':',':',':',':'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);

        for ispecies = 1:numSpecies
            figure(8+ispecies+figureOffset)
            clf
            numRows = 2;
            numCols = 4;
            plotNum = 1;
            
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, FSABFlowScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('<BV>')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, flowOutboardScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Outboard V_{||}')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, flowInboardScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Inboard V_{||}')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, heatFluxScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('heatFlux')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, particleFluxScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('particleFlux')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, particleSourceScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('particle source')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, heatSourceScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            legend(legendText)
            title('heat source')
            axis tight
        end
            
                
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 2
        % Convergence scan in the local approximation
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning convergence scans.\n')
%        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        makeLocalApproximation = true;
        psiDiameter = psiDiameterConverged;
        
        scanStartTime = tic;
        
        baseParameters = [NpsiPerDiameterConverged; NthetaConverged; NxiConverged; NxConverged; NLConverged; NxPotentialsPerVthConverged; log10tolConverged];
        numRunsInScan = numel(baseParameters)+1;
        
        parameters = repmat(baseParameters, [1, numRunsInScan]);
        parameters(2, 2) = 2*NthetaConverged+1;
        parameters(3, 3) = 2*NxiConverged;
        parameters(4, 4) = floor(1.5*NxConverged);
        parameters(5, 5) = 2*NLConverged;
        parameters(6, 6) = 2*NxPotentialsPerVthConverged;
        parameters(1, 7) = 2*NpsiPerDiameterConverged;
        parameters(7, 8) = log10tolConverged + 1;
        legendText={'Base case','2x N\theta', '2x N\xi', '1.5x Nx', '2x NL', '2x NxPotentials','2x N\psi PerDiameter','tol/10'};
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        FSABFlowScan=cell(numRunsInScan,1);
        flowOutboardScan=cell(numRunsInScan,1);
        flowInboardScan=cell(numRunsInScan,1);
        heatFluxScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        particleSourceScan=cell(numRunsInScan,1);
        heatSourceScan=cell(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            fprintf('\n\nRun %d of %d: %s\n\n',runNum, numRunsInScan, legendText{runNum})
            NpsiPerDiameter = parameters(1, runNum);
            Ntheta = parameters(2, runNum);
            Nxi = parameters(3, runNum);
            Nx = parameters(4, runNum);
            NL = parameters(5, runNum);
            NxPotentialsPerVth = parameters(6, runNum);
            tol = 10^(-parameters(7,runNum));
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            flowOutboardScan{runNum} = flowOutboard;
            flowInboardScan{runNum} = flowInboard;
            FSABFlowScan{runNum} = FSABFlow;
            heatFluxScan{runNum} = heatFlux;
            particleFluxScan{runNum} = particleFlux;
            particleSourceScan{runNum} = particleSource;
            heatSourceScan{runNum} = heatSource;
            didItConvergeScan(runNum) = didItConverge;
        end
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.9,0.6,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.9,0.6,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'--','-','-','-','-','-',':',':',':',':',':'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);

        for ispecies = 1:numSpecies
            figure(8+ispecies+figureOffset)
            clf
            numRows = 2;
            numCols = 4;
            plotNum = 1;
            
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, FSABFlowScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('<BV_{||}>')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, flowOutboardScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Outboard V_{||}')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, flowInboardScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Inboard V_{||}')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, heatFluxScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('heatFlux')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, particleFluxScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('particleFlux')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, particleSourceScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('particle source')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, heatSourceScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            legend(legendText)
            title('heat source')
        end
        
                
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 3
        % Scan dT/dpsi
        
        dTHatdpsis = -logspace(-1,0.5,5);
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning convergence scans.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInScan = numel(dTHatdpsis);
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        kParOutboardScan=cell(numRunsInScan,1);
        kParInboardScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        qScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        residualScan=zeros(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            fprintf('\n\nRun %d of %d: dTHatdpsi = %g\n\n',runNum, numRunsInScan, dTHatdpsi)
            dTHatdpsiScalar = dTHatdpsis(runNum);
            legendText{runNum} = ['dTHat/d\psi = ',num2str(dTHatdpsiScalar)];
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            kParOutboardScan{runNum} = kParOutboard;
            kParInboardScan{runNum} = kParInboard;
            nOutboardScan{runNum} = nOutboard;
            qScan{runNum} = q;
            particleFluxScan{runNum} = particleFlux;
            residualScan(runNum) = residual;
            didItConvergeScan(runNum) = didItConverge;
        end
        
        
        figure(8+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParInboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=\pi')
        
        subplot(numRows, numCols, 3)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, 4)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, qScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('q')
        
        subplot(numRows, numCols, 5)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleFluxScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 6
        % Scan desiredFWHMInRhoTheta
        
        desiredFWHMInRhoThetas = linspace(1,3,5);
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning convergence scans.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInScan = numel(desiredFWHMInRhoThetas);
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        kParOutboardScan=cell(numRunsInScan,1);
        kParInboardScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        qScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        residualScan=zeros(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            fprintf('\n\nRun %d of %d: desiredFWHMInRhoThetas = %g\n\n',runNum, numRunsInScan, desiredFWHMInRhoTheta)
            desiredFWHMInRhoTheta = desiredFWHMInRhoThetas(runNum);
            Delta = 0.0006*(epsilon/0.01)/(desiredFWHMInRhoTheta/0.24)/30;
            legendText{runNum} = ['desiredFWHMInRhoTheta = ',num2str(desiredFWHMInRhoTheta)];
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            kParOutboardScan{runNum} = kParOutboard;
            kParInboardScan{runNum} = kParInboard;
            nOutboardScan{runNum} = nOutboard;
            qScan{runNum} = q;
            particleFluxScan{runNum} = particleFlux;
            residualScan(runNum) = residual;
            didItConvergeScan(runNum) = didItConverge;
        end
        
        
        figure(8+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParInboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=\pi')
        
        subplot(numRows, numCols, 3)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, 4)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, qScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('q')
        
        subplot(numRows, numCols, 5)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleFluxScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 7
        % Scan desiredU
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of U.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInScan = numel(desiredUs);
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        kParOutboardScan=cell(numRunsInScan,1);
        kParInboardScan=cell(numRunsInScan,1);
        FSAKParScan=cell(numRunsInScan,1);
        LHSOfKEquationScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        qScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        residualScan=zeros(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        LHSOfKEquationAtPsiMidScan = zeros(numRunsInScan,1);
        heatFluxRelativeToLocalPlateauRegimeScan = cell(numRunsInScan,1);
        qNormalizedToPlateauInUScan = zeros(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            desiredU = desiredUs(runNum);
            fprintf('\n\nRun %d of %d: desiredU = %g\n\n',runNum, numRunsInScan, desiredU)
            legendText{runNum} = ['U = ',num2str(desiredU)];
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            kParOutboardScan{runNum} = kParOutboard;
            kParInboardScan{runNum} = kParInboard;
            nOutboardScan{runNum} = nOutboard;
            qScan{runNum} = q;
            particleFluxScan{runNum} = particleFlux;
            residualScan(runNum) = residual;
            didItConvergeScan(runNum) = didItConverge;
            FSAKParScan{runNum} = FSAKPar;
            LHSOfKEquationScan{runNum} = LHSOfKEquation;
            LHSOfKEquationAtPsiMidScan(runNum) = LHSOfKEquationAtPsiMid;
            heatFluxRelativeToLocalPlateauRegimeScan{runNum} = heatFluxRelativeToLocalPlateauRegime;
            qNormalizedToPlateauInUScan(runNum) = interp1(psi,heatFluxRelativeToLocalPlateauRegime,psiMid,'cubic');
        end
        
        
        figure(8+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParInboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=\pi')
        
        subplot(numRows, numCols, 3)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, 4)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, qScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('q')
        
        subplot(numRows, numCols, 5)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleFluxScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case {8,13,14}
        % 8 = Scan Nx
        % 13 = Scan Ntheta
        % 14 = Scan NpsiPerDiameter
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of Nx.\n')
%        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        switch runMode
            case 8
                Ntheta=NthetaConverged;
                NpsiPerDiameter = NpsiPerDiameterConverged;
                numRunsInScan = numel(Nxs);
            case 13
                Nx=NxConverged;
                NpsiPerDiameter = NpsiPerDiameterConverged;
                numRunsInScan = numel(Nthetas);
            case 14
                Nx=NxConverged;
                Ntheta=NthetaConverged;
                numRunsInScan = numel(NpsiPerDiameters);
        end
        scanStartTime = tic;
        
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        FSABFlowScan=cell(numRunsInScan,1);
        flowOutboardScan=cell(numRunsInScan,1);
        flowInboardScan=cell(numRunsInScan,1);
        heatFluxScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        particleSourceScan=cell(numRunsInScan,1);
        heatSourceScan=cell(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            switch runMode
                case 8
                    Nx = Nxs(runNum);
                    legendText{runNum} = ['Nx = ',num2str(Nx)];
                    fprintf('\n\nRun %d of %d: Nx = %d\n\n',runNum, numRunsInScan, Nx)
                case 13
                    Ntheta = Nthetas(runNum);
                    legendText{runNum} = ['N\theta = ',num2str(Ntheta)];
                    fprintf('\n\nRun %d of %d: Ntheta = %d\n\n',runNum, numRunsInScan, Ntheta)
                case 14
                    NpsiPerDiameter = NpsiPerDiameters(runNum);
                    legendText{runNum} = ['NpsiPerDiameter = ',num2str(NpsiPerDiameter)];
                    fprintf('\n\nRun %d of %d: NpsiPerDiameter = %d\n\n',runNum, numRunsInScan, NpsiPerDiameter)
            end
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            flowOutboardScan{runNum} = flowOutboard;
            flowInboardScan{runNum} = flowInboard;
            FSABFlowScan{runNum} = FSABFlow;
            heatFluxScan{runNum} = heatFlux;
            particleFluxScan{runNum} = particleFlux;
            particleSourceScan{runNum} = particleSource;
            heatSourceScan{runNum} = heatSource;
            didItConvergeScan(runNum) = didItConverge;
        end
        
        
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.9,0.6,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.9,0.6,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'--','-','-','-','-','-',':',':',':',':',':'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);

        for ispecies = 1:numSpecies
            figure(8+ispecies+figureOffset)
            clf
            numRows = 2;
            numCols = 4;
            plotNum = 1;
            
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, FSABFlowScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('<BV>')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, flowOutboardScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Outboard V_{||}')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, flowInboardScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Inboard V_{||}')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, heatFluxScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('heatFlux')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, particleFluxScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('particleFlux')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, particleSourceScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('particle source')
            axis tight
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum+1;
            for runNum=1:numRunsInScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum}, heatSourceScan{runNum}(ispecies,:), linespecs{index},'Color',colors(index,:))
                hold on
            end
            legend(legendText)
            title('heat source')
            axis tight
        end
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 9
        % Scan psiDiameter
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of psiDiameter.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInScan = numel(psiDiameters);
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        kParOutboardScan=cell(numRunsInScan,1);
        kParInboardScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        qScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        residualScan=zeros(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            psiDiameter = psiDiameters(runNum);
            legendText{runNum} = ['psiDiameter = ',num2str(psiDiameter)];
            fprintf('\n\nRun %d of %d: psiDiameter = %g\n\n',runNum, numRunsInScan, psiDiameter)
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            kParOutboardScan{runNum} = kParOutboard;
            kParInboardScan{runNum} = kParInboard;
            nOutboardScan{runNum} = nOutboard;
            qScan{runNum} = q;
            particleFluxScan{runNum} = particleFlux;
            residualScan(runNum) = residual;
            didItConvergeScan(runNum) = didItConverge;
        end
        
        
        figure(8+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParInboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=\pi')
        
        subplot(numRows, numCols, 3)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, 4)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, qScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('q')
        
        subplot(numRows, numCols, 5)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleFluxScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 10
        % Scan exponent
        
        minExponent = 0.85;
        maxExponent = 1.15;
        exponents = linspace(minExponent, maxExponent,4);
        heatSourcesAtPsiMid = zeros(size(exponents));
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of exponent.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInScan = numel(exponents)+1;
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        kParOutboardScan=cell(numRunsInScan,1);
        kParInboardScan=cell(numRunsInScan,1);
        FSAKParScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        qScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        heatSources=cell(numRunsInScan,1);
        particleSources=cell(numRunsInScan,1);
        residualScan=zeros(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            if runNum < numRunsInScan
                exponent = exponents(runNum);
            else
                p = polyfit(exponents, heatSourcesAtPsiMid, 2);
                myRoots=roots(p);
                
                figure(200)
                clf
                plot(exponents, heatSourcesAtPsiMid,'.-')
                xlabel('exponent')
                ylabel('heat source at \psi Mid')
                numValidRoots = 0;
                for i=1:numel(myRoots)
                    if myRoots(i)>= 0 && myRoots(i) < 2
                        numValidRoots = numValidRoots + 1;
                        bestExponent = myRoots(i);
                    end
                end
                if numValidRoots < 1 || numValidRoots > 1
                    beep
                    fprintf('\n***\nWarning: No good exponent was found between 0 and 2.\n***\n')
                    exponent = 1;
                else
                    exponent = bestExponent;
                end
            end
            fprintf('\n\nRun %d of %d: exponent = %g\n\n',runNum, numRunsInScan, exponent)
            legendText{runNum} = ['exponent = ',num2str(exponent)];
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            kParOutboardScan{runNum} = kParOutboard;
            kParInboardScan{runNum} = kParInboard;
            FSAKParScan{runNum} = FSAKPar;
            nOutboardScan{runNum} = nOutboard;
            qScan{runNum} = q;
            particleFluxScan{runNum} = particleFlux;
            heatSources{runNum} = yEnergySourceProfile;
            particleSources{runNum} = yParticleSourceProfile;
            residualScan(runNum) = residual;
            didItConvergeScan(runNum) = didItConverge;
            heatSourcesAtPsiMid(runNum) = interp1(psi,yEnergySourceProfile, psiMid,'cubic');
        end
        
        
        figure(8+figureOffset)
        clf
        numRows = 3;
        numCols = 3;
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);
        
        plotNum = 1;
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParInboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('k_{||} at \theta=\pi')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, FSAKParScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('<k_{||}>')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, heatSources{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('Heat sources')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleSources{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('Particle sources')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, qScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        %legend(legendText)
        title('q')
        
        subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleFluxScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 11
        % Scan U, and within that scan, scan exponent
        
        minExponent = 0.85;
        maxExponent = 1.15;
        exponents = linspace(minExponent, maxExponent,4);
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of exponent.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInExponentScan = numel(exponents)+1;
        numRunsInUScan = numel(desiredUs);
        kParInUScan = zeros(numRunsInUScan,1);
        qNormalizedToPlateauInUScan = zeros(numRunsInUScan,1);
        heatSourcesAtPsiMid = zeros(numRunsInExponentScan, numRunsInUScan);
        
        legendText=cell(numRunsInExponentScan,1);
        
        psiScan=cell(numRunsInExponentScan,numRunsInUScan);
        thetaScan=cell(numRunsInExponentScan,numRunsInUScan);
        kParOutboardScan=cell(numRunsInExponentScan,numRunsInUScan);
        kParInboardScan=cell(numRunsInExponentScan,numRunsInUScan);
        FSAKParScan=cell(numRunsInExponentScan,numRunsInUScan);
        nOutboardScan=cell(numRunsInExponentScan,numRunsInUScan);
        qScan=cell(numRunsInExponentScan,numRunsInUScan);
        particleFluxScan=cell(numRunsInExponentScan,numRunsInUScan);
        heatSources=cell(numRunsInExponentScan,numRunsInUScan);
        particleSources=cell(numRunsInExponentScan,numRunsInUScan);
        residualScan=zeros(numRunsInExponentScan,numRunsInUScan);
        didItConvergeScan=true(numRunsInExponentScan,numRunsInUScan);
        
        overallRunNum=0;
        for URunNum = 1:numRunsInUScan
            desiredU = desiredUs(URunNum);
            for runNum = 1:numRunsInExponentScan
                if runNum < numRunsInExponentScan
                    exponent = exponents(runNum);
                else
                    p = polyfit(exponents(:), heatSourcesAtPsiMid(1:end-1,URunNum), 2);
                    myRoots=roots(p);
                    
                    figure(200)
                    clf
                    plot(exponents(:), heatSourcesAtPsiMid(1:end-1,URunNum),'.-')
                    xlabel('exponent')
                    ylabel('heat source at \psi Mid')
                    numValidRoots = 0;
                    for i=1:numel(myRoots)
                        if myRoots(i)>= 0 && myRoots(i) < 2
                            numValidRoots = numValidRoots + 1;
                            bestExponent = myRoots(i);
                        end
                    end
                    if numValidRoots < 1 || numValidRoots > 1
                        beep
                        fprintf('\n***\nWarning: No good exponent was found between 0 and 2.\n***\n')
                        exponent = 1;
                    else
                        exponent = bestExponent;
                    end
                end
                overallRunNum = overallRunNum+1;
                fprintf('\n\nRun %d of %d: U = %g, exponent = %g\n\n',overallRunNum, numRunsInExponentScan*numRunsInUScan, desiredU, exponent)
                legendText{runNum} = ['exponent = ',num2str(exponent)];
                
                solveDKE()
                
                psiScan{runNum,URunNum} = psi;
                thetaScan{runNum,URunNum} = theta;
                kParOutboardScan{runNum,URunNum} = kParOutboard;
                kParInboardScan{runNum,URunNum} = kParInboard;
                FSAKParScan{runNum,URunNum} = FSAKPar;
                nOutboardScan{runNum,URunNum} = nOutboard;
                qScan{runNum,URunNum} = q;
                qRelativeToPlateauScan{runNum,URunNum} = heatFluxRelativeToLocalPlateauRegime;
                particleFluxScan{runNum,URunNum} = particleFlux;
                heatSources{runNum,URunNum} = yEnergySourceProfile;
                particleSources{runNum,URunNum} = yParticleSourceProfile;
                residualScan(runNum,URunNum) = residual;
                didItConvergeScan(runNum,URunNum) = didItConverge;
                heatSourcesAtPsiMid(runNum,URunNum) = interp1(psi,yEnergySourceProfile, psiMid,'cubic');
                if runNum == numRunsInExponentScan
                    kParInUScan(URunNum) = interp1(psi,FSAKPar, psiMid,'cubic');
                    qNormalizedToPlateauInUScan(URunNum) = interp1(psi,heatFluxRelativeToLocalPlateauRegime, psiMid,'cubic');
                end
            end
            
            
            figure(8+figureOffset)
            clf
            numRows = 3;
            numCols = 4;
            
            colors = [0,0,0;
                1,0,0;
                0,0.7,0;
                0,0,1;
                0.6,0.6,0.6;
                1,0,0;
                0,0.7,0;
                0,0,1;
                0.6,0.6,0.6;
                0.7,0.4,0;
                0,0,0;
                0,0.7,1;
                1,0,0];
            linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
            maxPlotStyles=min([numel(linespecs),size(colors,1)]);
            
            plotNum = 1;
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, kParOutboardScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('k_{||} at \theta=0')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, kParInboardScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('k_{||} at \theta=\pi')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, FSAKParScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('<k_{||}>')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, nOutboardScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('n_1 at \theta=0')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, nOutboardScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('n_1 at \theta=0')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, heatSources{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Heat sources')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, particleSources{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('Particle sources')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, qScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('q')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, qRelativeToPlateauScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            %legend(legendText)
            title('q / q_{analytic local plateau}')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            for runNum=1:numRunsInExponentScan
                index = 1+mod(runNum-1,maxPlotStyles);
                plot(psiScan{runNum,URunNum}, particleFluxScan{runNum,URunNum}, linespecs{index},'Color',colors(index,:))
                hold on
            end
            legend(legendText)
            title('particle flux')
        end
        
        figure(300)
        clf
        numRows=2;
        numCols = 1;
        plotNum=1;
        
        U = linspace(0,1);
        J = (1+4*U.^2 + 6*U.^4 + 12*U.^6) ./ (1+2*U.^2 + 2*U.^4);
        k=-0.5*J;
        U2 = U.^2;
        P = exp(-U2) .* (1 + 4*U.^2 + 8*U.^4 + (4*(4*U.^6 + U.^8)/3)) ./ (1+2*U.^2 + 2*U.^4);
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        plot(desiredUs, kParInUScan,'or')
        hold on
        plot(U, k)
        legend('numerical','analytic')
        xlabel('U')
        ylabel('k_{||}')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        plot(desiredUs, qNormalizedToPlateauInUScan,'or')
        hold on
        plot(U, P)
        legend('numerical','analytic')
        xlabel('U')
        ylabel('P (heat flux relative to local plateau regime)')
        
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 4
        % scan of source types
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of source types.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        scanStartTime = tic;
        
        numRunsInScan = numSourceOptions*(numSourceOptions-1)/2;
        legendText=cell(numRunsInScan,1);
        
        Npsi = NpsiConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        kParOutboardScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        
        runNum=1;
        for source1 = 1:(numSourceOptions-1)
            for source2 = (source1+1):numSourceOptions
                fprintf('\n\nRun %d of %d: source1 = %d, source2 = %d\n\n',runNum, numRunsInScan, source1, source2)
                legendText{runNum} = ['s1=',num2str(source1),' s2=',num2str(source2)];
                solveDKE()
                kParOutboardScan{runNum} = kParOutboard;
                nOutboardScan{runNum} = nOutboard;
                runNum = runNum+1;
            end
        end
        
        figure(8+figureOffset)
        clf
        numRows = 1;
        numCols = 2;
        
        colors = [1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1];
        linespecs = {'-','-','-','-',':',':',':'};
        numPlotStyles = min([numel(linespecs), size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            plot(psi, kParOutboardScan{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            plot(psi, nOutboardScan{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 5
        % scan of options for radial boundary conditions
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of source types.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        scanStartTime = tic;
        
        numRunsInScan = 3;
        legendText={'Dirichlet BCs everywhere','Partial Dirichlet, parallel=t','Partial Dirichlet, parallel=f'};
        DirichletEverywheres = [true,false,false];
        imposeParallels = [true, true, false];
        
        Npsi = NpsiConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        kParOutboardScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        particleSources=cell(numRunsInScan,1);
        heatSources=cell(numRunsInScan,1);
        particleFluxes=cell(numRunsInScan,1);
        
        for runNum=1:numRunsInScan
            imposeDirichletConditionsEverywhereAtRadialBoundaries = DirichletEverywheres(runNum);
            imposeLocalSolutionWhereTrajectoriesAreParallel = imposeParallels(runNum);
            fprintf('\n\nRun %d of %d\n\n',runNum, numRunsInScan)
            solveDKE()
            kParOutboardScan{runNum} = kParOutboard;
            nOutboardScan{runNum} = nOutboard;
            particleSources{runNum} = yParticleSourceProfile;
            heatSources{runNum} = yEnergySourceProfile;
            particleFluxes{runNum} = particleFlux;
        end
        
        figure(8+figureOffset)
        clf
        numRows = 2;
        numCols = 2;
        
        colors = [1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1];
        linespecs = {'-','-','-','-',':',':',':'};
        numPlotStyles = min([numel(linespecs), size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            plot(psi, kParOutboardScan{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        %{
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            plot(psi, nOutboardScan{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        %}
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            plot(psi, particleFluxes{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        subplot(numRows, numCols, 3)
        for runNum=1:numRunsInScan
            plot(psi, particleSources{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('particle source')
        
        subplot(numRows, numCols, 4)
        for runNum=1:numRunsInScan
            plot(psi, heatSources{runNum}, linespecs{1+mod(runNum-1,numPlotStyles)},'Color',colors(1+mod(runNum-1,numPlotStyles),:))
            hold on
        end
        legend(legendText)
        title('heat source')
        
        
        %{
stringForTop=sprintf('A = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, Legendre modal, Landreman polynomial colocation in x, solutionMethod = %d',Miller_A, nuPrime,nuStar, thetaGridMode, solutionMethod);
    annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',12,'LineStyle','none','String',stringForTop);
        %}
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
    case 12
        % Scan widthMultiplier
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning scan of widthMultiplier.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        psiDiameter = psiDiameterConverged;
        NpsiPerDiameter = NpsiPerDiameterConverged;
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        scanStartTime = tic;
        
        numRunsInScan = numel(widthMultipliers);
        
        legendText=cell(numRunsInScan,1);
        
        psiScan=cell(numRunsInScan,1);
        thetaScan=cell(numRunsInScan,1);
        kParOutboardScan=cell(numRunsInScan,1);
        kParInboardScan=cell(numRunsInScan,1);
        FSAKParScan=cell(numRunsInScan,1);
        LHSOfKEquationScan=cell(numRunsInScan,1);
        nOutboardScan=cell(numRunsInScan,1);
        qScan=cell(numRunsInScan,1);
        UScan=cell(numRunsInScan,1);
        particleFluxScan=cell(numRunsInScan,1);
        residualScan=zeros(numRunsInScan,1);
        didItConvergeScan=true(numRunsInScan,1);
        LHSOfKEquationAtPsiMidScan = zeros(numRunsInScan,1);
        heatFluxRelativeToLocalPlateauRegimeScan = cell(numRunsInScan,1);
        qNormalizedToPlateauInUScan = zeros(numRunsInScan,1);
        
        for runNum = 1:numRunsInScan
            widthExtender = widthExtenders(runNum);
            fprintf('\n\nRun %d of %d: widthExtender = %g\n\n',runNum, numRunsInScan, widthExtender)
            legendText{runNum} = ['widthExtender = ',num2str(widthExtender)];
            
            solveDKE()
            
            psiScan{runNum} = psi;
            thetaScan{runNum} = theta;
            kParOutboardScan{runNum} = kParOutboard;
            kParInboardScan{runNum} = kParInboard;
            nOutboardScan{runNum} = nOutboard;
            qScan{runNum} = q;
            UScan{runNum} = U;
            particleFluxScan{runNum} = particleFlux;
            residualScan(runNum) = residual;
            didItConvergeScan(runNum) = didItConverge;
            FSAKParScan{runNum} = FSAKPar;
            LHSOfKEquationScan{runNum} = LHSOfKEquation;
            LHSOfKEquationAtPsiMidScan(runNum) = LHSOfKEquationAtPsiMid;
            heatFluxRelativeToLocalPlateauRegimeScan{runNum} = heatFluxRelativeToLocalPlateauRegime;
            qNormalizedToPlateauInUScan(runNum) = interp1(psi,heatFluxRelativeToLocalPlateauRegime,psiMid,'cubic');
        end
        
        
        figure(8+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        
        colors = [0,0,0;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            1,0,0;
            0,0.7,0;
            0,0,1;
            0.6,0.6,0.6;
            0.7,0.4,0;
            0,0,0;
            0,0.7,1;
            1,0,0];
        linespecs = {'.-','.-','.-','.-','.:','.:','.:','.:'};
        maxPlotStyles=min([numel(linespecs),size(colors,1)]);
        
        subplot(numRows, numCols, 1)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=0')
        
        subplot(numRows, numCols, 2)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, kParInboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('k_{||} at \theta=\pi')
        
        subplot(numRows, numCols, 3)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, nOutboardScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('n_1 at \theta=0')
        
        subplot(numRows, numCols, 4)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, qScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('q')
        
        subplot(numRows, numCols, 5)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            plot(psiScan{runNum}, particleFluxScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('particle flux')
        
        subplot(numRows, numCols, 6)
        for runNum=1:numRunsInScan
            index = 1+mod(runNum-1,maxPlotStyles);
            alpha = -2*UScan{runNum}.^2;
            plot(psiScan{runNum}, LHSOfKEquationScan{runNum} + alpha .* FSAKParScan{runNum}, linespecs{index},'Color',colors(index,:))
            hold on
        end
        legend(legendText)
        title('LHS of k equation')
        
        if saveResults
            temp=dbstack;
            nameOfThisProgram=sprintf('%s',temp.file);
            filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
            filename=[filenameBase,'.mat'];
            save(filename)
        end
        
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',toc(scanStartTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    otherwise
        error('Invalid runMode')
end

    function solveDKE()
        
        startTimeForThisRun=tic;
        
        sqrtpi=sqrt(pi);
        iteration = iteration+1;
        
        
        % Order of the rows of the matrix and of the RHS:
        % --------------------------------
        % for j = 1:Npsi
        %   for s = 1:numSpecies
        %     for x = dx to xMax-dx  (Nx-2 points)
        %       for L = 0:(NL-1)
        %         for itheta = 1:Ntheta
        %           Enforce the drift-kinetic equation, or set f_1 to the local solution.
        % for j = 1:Npsi
        %   for s = 1:numSpecies
        %     Force <n_1> = 0
        %     Force <p_1> = 0
        
        
        % Order of the vector of unknowns & of columns in the matrix:
        % --------------------------------
        % for j = 1:Npsi
        %   for s = 1:numSpecies
        %     for x = dx to xMax-dx  (Nx-2 points)
        %       for L = 0:(NL-1)
        %         for itheta = 1:Ntheta
        %           f_1
        % for j = 1:Npsi
        %   for s = 1:numSpecies
        %     Particle source
        %     Heat source
        
        if mod(Ntheta,2)==0
            % Force Ntheta to be odd.
            Ntheta = Ntheta+1;
        end
        
        psiDiameterFinal = psiDiameter + widthExtender*2;
        
        psiMax = psiMid + psiDiameterFinal/2;
        psiMin = psiMid - psiDiameterFinal/2 - psiInnerExtender;
        Npsi = max([5, round(NpsiPerDiameter * (psiMax-psiMin)) + 1]);
        
        fprintf('*************************************************\n')
        fprintf('NpsiPerDiameter = %g, Npsi = %d, Ntheta = %d,  Nxi = %d,  Nx = %d, NL = %d,  NxPtentialsPerVth = %g, tol = %g\n',NpsiPerDiameter,Npsi,Ntheta,Nxi,Nx,NL,NxPotentialsPerVth,tol)
        fprintf('psiDiameterFinal = %g, psiMin = %g, psiMax = %g, Delta = %g, omega = %g\n',psiDiameterFinal, psiMin, psiMax, Delta, omega)
        
        tic
        
        % Generate abscissae, quadrature weights, and derivative matrix for theta grid.
        % Both abscissae and weights should be column vectors.
        switch psiGridMode
            case 0
                % Chebyshev
                [psi, psiWeights, ddpsi] = multiChebyshevWeightsAndDifferentiation(Npsi, psiMin, psiMax, NpsiIntervals);
                d2dpsi2 = ddpsi*ddpsi;
            case 1
                % 2nd order finite difference
                scheme = 2;
                [psi, psiWeights, ddpsi, d2dpsi2] = differentiationMatricesForUniformGrid(Npsi, psiMin, psiMax, scheme);
            case 2
                % 4th order finite difference
                scheme = 12;
                [psi, psiWeights, ddpsi, d2dpsi2] = differentiationMatricesForUniformGrid(Npsi, psiMin, psiMax, scheme);
            otherwise
                error('Invalid psiGridMode!')
        end
        psi = psi(:)';
        psiWeights = psiWeights(:)';
        
        switch psiDerivativeForPreconditioner
            case 0
                %error('Need to write code here')
                if psiGridMode == 0
                    error('This psiDerivativeForPreconditioner setting does not work for psiGridMode=0')
                    return
                end
                scheme = 2;
                [psiPreconditioner, ~, ddpsiPreconditioner, d2dpsi2Preconditioner] = differentiationMatricesForUniformGrid(Npsi, psiMin, psiMax, scheme);
                if any(psiPreconditioner(:) ~= psi(:))
                    error('Something went wrong here...')
                end
            case 1
                d = psi(:) - circshift(psi(:), [1,0]);
                ddpsiPreconditioner = diag(1./d) - diag(1./d(2:end), -1);
                ddpsiPreconditioner(1,1)=0;
            otherwise
                error('Invalid psiDerivativeForPreconditioner')
        end

        
        
        
        
        
        
        switch radialProfiles
            case 0
                THat = 1 - (psi-1);
                dTHatdpsi = -1 * ones(size(psi));
                IHat = ones(size(psi));
                dIHatdpsi = zeros(size(psi));
                etaHat = ones(size(psi));
            case 1
                THat = exp(1-psi);
                dTHatdpsi = -exp(1-psi);
                IHat =  1+0.1*sin(psi);
                dIHatdpsi = 0.1*cos(psi);
                etaHat = exp((psi-1)*1.7);
        end
        
        switch radialProfiles
            case {0,1}
                middle = (psiMin+psiMax)/2;
                nHat = 1 - 0.5*erf((psi-middle)/0.03);
                
                PhiHat = Delta/(2*omega)*THat .* log(etaHat./nHat);
                if psiGridMode > 0
                    scheme = 12;
                    [psi, psiWeights, ddpsiForPhi, ~] = differentiationMatricesForUniformGrid(Npsi, psiMin, psiMax, scheme);
                else
                    ddpsiForPhi = ddpsi;
                end
                dPhiHatdpsi = (ddpsiForPhi * PhiHat(:))';
                if omega==0
                    dPhiHatdpsi = zeros(size(psi));
                end
            case {2,3}
                IHat = ones(size(psi));
                dIHatdpsi = zeros(size(psi));
                ss=150;
                r0=0.07 + widthExtender;
                %r0 = 0.07;
                rampAmplitude = desiredU * psiAHat / omega;
                if setTPrimeToBalanceHeatFlux
                    % Nonlinearly solve for T(psi) and n(psi) such that
                    % n^exponent * dT/dpsi is constant, for the given
                    % eta(psi) and PhiHat(psi):
                    
                    numIterations = 10;
                    
                    NpsiFine = 200;
                    NpsiIntervals=1;
                    %[psiFine, ~, ddpsiForT] = multiChebyshevWeightsAndDifferentiation(NpsiFine, psiMin, psiMax, NpsiIntervals);
                    
                    scheme = 12;
                    [psiFine, ~, ddpsiForT, ~] = differentiationMatricesForUniformGrid(NpsiFine, psiMin-(1e-10), psiMax+(1e-10), scheme);
                    psiFine = psiFine(:)';
                    
                    PhiHat = compute_Phi(psiFine-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    THat = ones(size(psiFine));
                    etaHat = 1 + (psiFine-psiMid)*detaHatdpsiScalar;
                    ddpsiForT(1,:)=0;
                    ddpsiForT(1,1)=1;
                    
                    figure(100)
                    clf
                    
                    for i=1:numIterations
                        plot(psiFine,THat)
                        hold on
                        oldTHat = THat;
                        nHat = etaHat .* exp(-2*omega/Delta*PhiHat./THat);
                        dTHatdpsi = dTHatdpsiScalar ./ (nHat.^exponent);
                        rhsForT = dTHatdpsi(:);
                        rhsForT(1) = 0;
                        THat = ddpsiForT \ rhsForT;
                        if mod(NpsiFine,2)==1
                            THat = THat - THat((NpsiFine+1)/2);
                        else
                            THat = THat - 0.5*(THat(NpsiFine/2) + THat(NpsiFine/2 + 1));
                        end
                        THat = (1+THat)';
                        
                        %THat = (THat + oldTHat)/2;
                    end
                    %etaHat = nHat .* exp(2*omega/Delta*PhiHat./THat);
                    
                    %THat = chebint2(THat,psiMin, psiMax, psi)';
                    %dTHatdpsi = chebint2(dTHatdpsi,psiMin, psiMax, psi)';
                    %nHat = chebint2(nHat,psiMin, psiMax, psi)';
                    
                    coarsifyMatrix = m20121127_02_makeHighOrderInterpolationMatrix(psiFine,psi,0,'f');
                    THat = (coarsifyMatrix * THat')';
                    dTHatdpsi = (coarsifyMatrix * dTHatdpsi')';
                    nHat = (coarsifyMatrix * nHat')';
                    
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    PhiHat = compute_Phi(psi-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psi-psiMid, ss, r0, rampAmplitude);
                else
                    THat = 1 + dTHatdpsiScalar * (psi-psiMid);
                    dTHatdpsi = dTHatdpsiScalar * ones(size(psi));
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    
                    PhiHat = compute_Phi(psi-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psi-psiMid, ss, r0, rampAmplitude);
                    
                    nHat = etaHat .* exp(-2*omega/Delta*PhiHat./THat);
                end
            case {4}
                IHat = ones(size(psi));
                dIHatdpsi = zeros(size(psi));
                ss = 25;
                r0=0;
                rampAmplitude = desiredU * psiAHat / omega;
                FSAB2 = 1.00005000374406;
                if setTPrimeToBalanceHeatFlux
                    % Nonlinearly solve for T(psi) and n(psi) such that
                    % n^exponent * dT/dpsi is constant, for the given
                    % eta(psi) and PhiHat(psi):
                    
                    numIterations = 10;
                    
                    NpsiFine = 200;
                    NpsiIntervals=1;
                    %[psiFine, ~, ddpsiForT] = multiChebyshevWeightsAndDifferentiation(NpsiFine, psiMin, psiMax, NpsiIntervals);
                    
                    scheme = 12;
                    [psiFine, ~, ddpsiForT, ~] = differentiationMatricesForUniformGrid(NpsiFine, psiMin-(1e-10), psiMax+(1e-10), scheme);
                    psiFine = psiFine(:)';
                    
                    matrixForT = ddpsiForT;
                    matrixForT(1,:)=0;
                    matrixForT(1,1)=1;
                    
                    PhiHat = compute_Phi(psiFine-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psiFine-psiMid, ss, r0, rampAmplitude);
                    THat = ones(size(psiFine));
                    etaHat = 1 + (psiFine-psiMid)*detaHatdpsiScalar;
                    IHatFine = ones(size(psiFine));
                    rhsForTPrime = dTHatdpsiScalar*ones(NpsiFine,1);
                    
                    THats = THat;
                    for i=1:numIterations
                        nHat = etaHat .* exp(-2*omega/Delta*PhiHat./THat);
                        U = omega * IHatFine .* dPhiHatdpsi ./ (psiAHat * sqrt(FSAB2 .* THat)); %A more accurate U using the real FSAB2 will be calculated later.
                        
                        % Here is my original formula for the analytic heat
                        % flux, which I now believe to be wrong:
                        %matrixForTPrime = diag((1/3)*(4*U.^8+16*U.^6+24*U.^4+12*U.^2+3)./(2*U.^4+2*U.^2+1).*exp(-U.*U).*nHat.*(THat.^1.5)) ...
                        %    * (eye(NpsiFine) - diag(Delta/psiAHat*U.*sqrt(THat/FSAB2).*IHatFine)*ddpsiForT);
                        
                        % Here is my corrected formula for the analytic
                        % heat flux:
                        %matrixForTPrime = diag((1/3)*(4*U.^8+16*U.^6+24*U.^4+12*U.^2+3)./(2*U.^4+2*U.^2+1).*exp(-U.*U).*nHat.*(THat.^1.5)) * ddpsiForT;
                        
                        %dTHatdpsi = dTHatdpsiScalar ./ (nHat.^exponent);
                        %dTHatdpsi = dTHatdpsiScalar ./ (nHat.^2);
                        dTHatdpsi = matrixForTPrime \ rhsForTPrime;
                        
                        %{
                        assignin('base','pm',PhiHat)
                        assignin('base','nm',nHat)
                        assignin('base','Um',U)
                        assignin('base','mm',matrixForTPrime)
                        assignin('base','rm',rhsForTPrime)
                        return
                        %}
                        
                        rhsForT = dTHatdpsi;
                        rhsForT(1)=0;
                        THat = matrixForT \ rhsForT;
                        
                        THatAtPsiMin = interp1(psiFine,THat,psiMid,'cubic');
                        THat = (THat - THatAtPsiMin + 1)';
                        THats = [THats;THat];
                    end
                    
                    figure(2)
                    clf
                    plot(psiFine, THats,'.-')
                    
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    %etaHat = nHat .* exp(2*omega/Delta*PhiHat./THat);
                    dTHatdpsi = dTHatdpsi';
                    
                    %{
                    THat = chebint2(THat,psiMin, psiMax, psi)';
                    dTHatdpsi = chebint2(dTHatdpsi,psiMin, psiMax, psi)';
                    nHat = chebint2(nHat,psiMin, psiMax, psi)';
                    %}
                    coarsifyMatrix = m20121127_02_makeHighOrderInterpolationMatrix(psiFine,psi,0,'f');
                    THat = (coarsifyMatrix * THat')';
                    dTHatdpsi = (coarsifyMatrix * dTHatdpsi')';
                    nHat = (coarsifyMatrix * nHat')';
                    
                    assignin('base','tm',THat)
                    assignin('base','dtm',dTHatdpsi)
                    assignin('base','nm',nHat)
                    
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    PhiHat = compute_Phi(psi-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psi-psiMid, ss, r0, rampAmplitude);
                else
                    THat = 1 + dTHatdpsiScalar * (psi-psiMid);
                    dTHatdpsi = dTHatdpsiScalar * ones(size(psi));
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    
                    PhiHat = compute_Phi(psi-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psi-psiMid, ss, r0, rampAmplitude);
                    
                    nHat = etaHat .* exp(-2*omega/Delta*PhiHat./THat);
                end
            case {6,7}
                IHat = ones(size(psi));
                dIHatdpsi = zeros(size(psi));
                if radialProfiles==6
                    ss = 25;
                    r0=0;
                else
                    ss=150;
                    r0=0.07 + widthExtender;
                end
                rampAmplitude = desiredU * psiAHat / omega;
                FSAB2 = 1.00005000374406;
                if setTPrimeToBalanceHeatFlux
                    % Nonlinearly solve for T(psi) and n(psi) such that
                    % n^exponent * dT/dpsi is constant, for the given
                    % eta(psi) and PhiHat(psi):
                    
                    numIterations = 10;
                    
                    NpsiFine = 200;
                    NpsiIntervals=1;
                    %[psiFine, ~, ddpsiForT] = multiChebyshevWeightsAndDifferentiation(NpsiFine, psiMin, psiMax, NpsiIntervals);
                    
                    scheme = 12;
                    [psiFine, ~, ddpsiForT, ~] = differentiationMatricesForUniformGrid(NpsiFine, psiMin-(1e-10), psiMax+(1e-10), scheme);
                    psiFine = psiFine(:)';
                    
                    ddpsiForT(1,:)=0;
                    ddpsiForT(1,1)=1;
                    
                    PhiHat = compute_Phi(psiFine-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psiFine-psiMid, ss, r0, rampAmplitude);
                    THat = ones(size(psiFine));
                    etaHat = 1 + (psiFine-psiMid)*detaHatdpsiScalar;
                    IHatFine = ones(size(psiFine));
                    
                    THats = THat;
                    analyticHeatFluxAtPsiMid = (1/3)*(4*desiredU.^8+16*desiredU.^6+24*desiredU.^4+12*desiredU.^2+3)...
                        ./(2*desiredU.^4+2*desiredU.^2+1).*exp(-desiredU.*desiredU);
                    for i=1:numIterations
                        nHat = etaHat .* exp(-2*omega/Delta*PhiHat./THat);
                        U = omega * IHatFine .* dPhiHatdpsi ./ (psiAHat * sqrt(FSAB2 .* THat)); %A more accurate U using the real FSAB2 will be calculated later.
                        
                        % Here is my corrected formula for the analytic
                        % heat flux:
                        analyticHeatFluxPrefactor = (1/3)*(4*U.^8+16*U.^6+24*U.^4+12*U.^2+3)./(2*U.^4+2*U.^2+1).*exp(-U.*U).*nHat.*(THat.^1.5);
                        
                        dTHatdpsi = dTHatdpsiScalar * analyticHeatFluxAtPsiMid ./ analyticHeatFluxPrefactor;
                        rhsForT = dTHatdpsi(:);
                        rhsForT(1) = 0;
                        THat = ddpsiForT \ rhsForT;
                        if mod(NpsiFine,2)==1
                            THat = THat - THat((NpsiFine+1)/2);
                        else
                            THat = THat - 0.5*(THat(NpsiFine/2) + THat(NpsiFine/2 + 1));
                        end
                        THat = (1+THat)';

                        %THatAtPsiMin = interp1(psiFine,THat,psiMid,'cubic');
                        %THat = (THat - THatAtPsiMin + 1)';
                        THats = [THats;THat];
                    end
                    
                    figure(2)
                    clf
                    plot(psiFine, THats,'.-')
                    
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    %etaHat = nHat .* exp(2*omega/Delta*PhiHat./THat);
                    dTHatdpsi = dTHatdpsi';
                    
                    %{
                    THat = chebint2(THat,psiMin, psiMax, psi)';
                    dTHatdpsi = chebint2(dTHatdpsi,psiMin, psiMax, psi)';
                    nHat = chebint2(nHat,psiMin, psiMax, psi)';
                    %}
                    coarsifyMatrix = m20121127_02_makeHighOrderInterpolationMatrix(psiFine,psi,0,'f');
                    THat = (coarsifyMatrix * THat')';
                    dTHatdpsi = (coarsifyMatrix * dTHatdpsi)';
                    nHat = (coarsifyMatrix * nHat')';
                    
                    assignin('base','tm',THat)
                    assignin('base','dtm',dTHatdpsi)
                    assignin('base','nm',nHat)
                    
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    PhiHat = compute_Phi(psi-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psi-psiMid, ss, r0, rampAmplitude);
                else
                    THat = 1 + dTHatdpsiScalar * (psi-psiMid);
                    dTHatdpsi = dTHatdpsiScalar * ones(size(psi));
                    etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar;
                    
                    PhiHat = compute_Phi(psi-psiMid, ss, r0, rampAmplitude) - compute_Phi(0, ss, r0, rampAmplitude);
                    dPhiHatdpsi = compute_dPhidr(psi-psiMid, ss, r0, rampAmplitude);
                    
                    nHat = etaHat .* exp(-2*omega/Delta*PhiHat./THat);
                end
            case 5
                %THat = T_wall + (T_ped-T_wall)*(tanh((psiMid-psi)/T_Delta)+1)/2 + (1-psi)*T_linear;
                %dTHatdpsi = (T_wall-T_ped)./(2*T_Delta*cosh((psiMid-psi)/T_Delta).^2) - T_linear;
                THat = T_wall + (T_ped-T_wall)*(tanh((psiMid-psi)/T_Delta)+1)/2 + T_linear*log(exp(linearCorner*(1-psi))+1)/linearCorner;
                dTHatdpsi = (T_wall-T_ped)./(2*T_Delta*cosh((psiMid-psi)/T_Delta).^2) - T_linear./(exp(linearCorner*(psi - 1)) + 1);
                dIHatdpsi = zeros(size(psi));
                
                IHat = ones(size(psi));
                etaHat = psi*0;
                
                PhiHat = minorRadius/2*(-Er_offset*(psi-psiMid) + (Er_depth+Er_offset)*Er_Delta*tanh((psi-psiMid)/Er_Delta));
                dPhiHatdpsi = minorRadius/2*(-Er_offset + (Er_depth+Er_offset)./(cosh((psi-psiMid)/Er_Delta).^2));
                
                nHat = n_wall + (n_ped-n_wall)*(tanh((psiMid-psi)/n_Delta)+1)/2 + n_linear*log(exp(linearCorner*(1-psi))+1)/linearCorner;
                
            case 8
                % C-mod profiles from Michael Churchill
                
                data = importdata(CModDataFile);
                data_psi = data.data(:,1);
                data_n = data.data(:,2)/(1e20);
                data_T = data.data(:,3)/1000;

                TShiftInPsi = 0.032; %0.02; 
                nShiftInPsi = 0.042; %0.03; 
                THat = interp1(data_psi, data_T, (psi-psiMid)/CModProfileStretchFactor + psiMid+ TShiftInPsi,'cubic');
                nHat = interp1(data_psi, data_n, (psi-psiMid)/CModProfileStretchFactor + psiMid + nShiftInPsi, 'cubic');
                
                % Don't let THat get too small:
                T0 = 0.05;
                w = 0.1;
                THat = 0.0 + THat - (w*log(tanh((THat - T0)/w) + 1))/2;
                dTHatdpsi = (ddpsi * THat')';
                
                PhiHat = minorRadius/2*(-Er_offset*(psi-psiMid) + (Er_depth+Er_offset)*Er_Delta*tanh((psi-psiMid)/Er_Delta));
                dPhiHatdpsi = minorRadius/2*(-Er_offset + (Er_depth+Er_offset)./(cosh((psi-psiMid)/Er_Delta).^2));

                % These next few lines are irrelevant if EFIT geometry is
                % used:
                IHat = ones(1,Npsi);
                dIHatdpsi = zeros(1,Npsi);
            otherwise
                error('Invalid setting for radialProfiles')
        end
        
        if forceElectrostaticConfinement
            PhiHat = -Delta/(2*omega)*THat.*log(nHat);
            dPhiHatdpsi = (ddpsi * PhiHat')';
        end
        
        % For now, for simplicity, force the density and temperature
        % profiles to be the same for all species, up to a constant:
        if false %radialProfiles==5 && numSpecies==2
            THats = [T_wall + (T_ped-T_wall)*(tanh((psiMid-psi)/T_Delta)+1)/2;
                T_wall_2 + (T_ped_2-T_wall_2)*(tanh((psiMid-psi)/T_Delta_2)+1)/2];
            dTHatdpsis = [(T_wall-T_ped)./(2*T_Delta*cosh((psiMid-psi)/T_Delta).^2);
                (T_wall_2-T_ped_2)./(2*T_Delta_2*cosh((psiMid-psi)/T_Delta_2).^2)];
                
        else
            THats = kron(scalarTHats(:), THat);
            dTHatdpsis = kron(scalarTHats(:), dTHatdpsi);
        end
        nHats = kron(scalarNHats(:), nHat);
        dnHatdpsis = (ddpsi * nHats')';
        etaHats = zeros(size(THats));
        for ispecies = 1:numSpecies
            etaHats(ispecies,:) =  nHats(ispecies,:) .* exp(2*charges(ispecies)*omega/Delta*PhiHat./THats(ispecies,:));
            detaHatdpsis = (ddpsi * etaHats')';
        end
        
        %{
        etaHats = kron(scalarNHats(:), etaHat);
        detaHatdpsis = (ddpsi * etaHats')';
        nHats = zeros(size(THats));
        for ispecies = 1:numSpecies
            nHats(ispecies,:) =  etaHats(ispecies,:) .* exp(-2*charges(ispecies)*omega/Delta*PhiHat./THats(ispecies,:));
            dnHatdpsis = (ddpsi * nHats')';
        end
        %}
        
        %{
        PhiHat
        THats
        dTHatdpsis
        nHats
        etaHats
        dPhiHatdpsi
        dnHatdpsis
        %}
        
        
        
        %{
        c1 = 6.8;
        c2 = 0.7;
        %phi = -0.9 + c2 * sin(c1*psi);
        %dphidpsi = c1 * c2 * cos(c1*psi);
        %dphidpsi = 1./cosh((psi-middle)/0.03);
        dphidpsi = (0.3*(1-cos(2*pi*(psi-psiMax)/(psiMax-psiMin)))).^2;
        phi = zeros(size(nHat));
        %}
        
        
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
        
        if shiftTheta
            dtheta = theta(2)-theta(1);
            theta = theta + dtheta/2;
        end
        %theta = theta';
        %thetaWeights = thetaWeights';
        
        
        if geometry == 4
            % Load EFIT data
            
            plotStuff = true;
            %plotStuff = false;
            
            %NPsi=1;
            %psi = desiredPsi;
            [thetaData, BData, BDotGradThetaData, IHat qData, as, R0, B0] = getGeometryFromEFITForSeveralFluxSurfaces(EFITFilename, psi, topCropZ, bottomCropZ, plotStuff);
            
            %IHat = abs(IHat);
            dIHatdpsi = (ddpsi * IHat')';
            
            % Filter B(theta)
            NThetaFilter = 100;
            thetaForBFilter = linspace(0,2*pi,NThetaFilter+1);
            thetaForBFilter(end)=[];
            BModes = zeros(NThetaFilter, Npsi);
            JModes = zeros(NThetaFilter, Npsi);
            for psiIndex = 1:Npsi
                %temp = interp1(thetaData{psiIndex}, BData{psiIndex}/B0, thetaForBFilter, 'spline');
                temp = interp1(thetaData{psiIndex}, BData{psiIndex}, thetaForBFilter, 'spline');
                BModes(:,psiIndex) = fft(temp)/NThetaFilter;
                %temp = interp1(thetaData{psiIndex}, BDotGradThetaData{psiIndex}/B0, thetaForBFilter, 'spline');
                temp = interp1(thetaData{psiIndex}, BDotGradThetaData{psiIndex}, thetaForBFilter, 'spline');
                JModes(:,psiIndex) = fft(temp)/NThetaFilter;
            end
            
            numFourierModesInThetaToKeepInEFITGeometry = 5;
            
            epsilon = as/R0;
            Miller_A = 1./epsilon;
            Miller_q = qData;
            fprintf('Inverse aspect ratio derived from EFIT equilibrium: %g to %g\n',min(epsilon),max(epsilon))
            
            
            BHat_beforeSmoothing = ones(Ntheta,1) * BModes(1,:);
            JHat_beforeSmoothing = ones(Ntheta,1) * JModes(1,:);
            keepUpDownAsymmetry = 1; % This variable should be either 1 or 0.
            for m=1:numFourierModesInThetaToKeepInEFITGeometry
                BHat_beforeSmoothing = BHat_beforeSmoothing + 2*cos(m*theta)*real(BModes(m+1,:)) - keepUpDownAsymmetry*2*sin(m*theta)*imag(BModes(m+1,:));
                JHat_beforeSmoothing = JHat_beforeSmoothing + 2*cos(m*theta)*real(JModes(m+1,:)) - keepUpDownAsymmetry*2*sin(m*theta)*imag(JModes(m+1,:));
            end

            % Smooth in the psi direction by fitting a polynomial:
            BHat = zeros(Ntheta,Npsi);
            JHat = zeros(Ntheta,Npsi);
            dBHatDPsi = zeros(Ntheta,Npsi);
            for itheta = 1:Ntheta
                [p,S,mu] = polyfit(psi, BHat_beforeSmoothing(itheta,:), polynomialFitDegreeForSmoothingEFITInPsi);
                [y,~]= polyval(p,psi,S,mu);
                BHat(itheta,:) = y;
                
                % Analytically differentiate the fitting polynomial:
                [y,~]= polyval([0,polyder(p)],psi,S,mu);
                dBHatDPsi(itheta,:) = y / mu(2);
                
                [p,S,mu] = polyfit(psi, JHat_beforeSmoothing(itheta,:), polynomialFitDegreeForSmoothingEFITInPsi);
                [y,~]= polyval(p,psi,S,mu);
                JHat(itheta,:) = y;
            end
            
            % Spectral uniform differentiation matrix:
            scheme = 20;
            [~, ~, ddthetaForBHat, ~] = differentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
            dBHatdtheta = ddthetaForBHat * BHat;
            
        else
            % We're not using an EFIT reconstruction.
            
            BHat1D = b(theta);
            
            switch radialProfiles
                case {0,2,3,4,5,6,7,8}
                    BHat = BHat1D * ones(1,Npsi);
                    JHat = (BHat1D ./ oneOverqRbDotGradTheta(theta)) * ones(1,Npsi);
                    dBHatDPsi = zeros(size(BHat));
                    factorFordBHatdtheta = ones(1,Npsi);
                case 1
                    BHat = BHat1D * sqrt(psi(:)');
                    %JHat = BHat ./(( 1./ oneOverqRbDotGradTheta(theta)) * ones(1,Npsi));
                    JHat = BHat ./(oneOverqRbDotGradTheta(theta) * ones(1,Npsi));
                    dBHatDPsi = BHat1D * (0.5./sqrt(psi(:)'));
                    factorFordBHatdtheta = sqrt(psi(:)');
                otherwise
                    error('Invalid setting for radialProfiles')
            end
            
            switch radialProfiles
                case {2,3,4,5,6,7,8}
                    JHat = JHat / Miller_q;
                case {0,1}
                otherwise
                    error('Invalid setting for radialProfiles')
            end
            
            if geometry==1
                % Miller: too hard to analytically differentiate b(theta) so do
                % it numerically:
                scheme = 20;
                [thetaFine, ~, ddthetaFine, ~] = differentiationMatricesForUniformGrid(Ntheta*dBdthetaResolutionMultiplier, 0, 2*pi, scheme);
                dbdthetaFine = ddthetaFine * b(thetaFine);
                dBHatdtheta = dbdthetaFine(1:dBdthetaResolutionMultiplier:end);
            else
                dBHatdtheta = dbdtheta(theta);
            end
            dBHatdtheta = dBHatdtheta * factorFordBHatdtheta;
        end
        
        if preconditionerMethod_theta==1
            % Uniform periodic 2nd order FD
            scheme = 0;
            [~, ~, ddthetaForPreconditioner, ~] = differentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
        else
            ddthetaForPreconditioner=ddtheta;
        end
        
        %{
        BHatMax = max(BHat,[],1);
        BHatMin = min(BHat,[],1);
        deltaTMax = -Delta * IHat ./ (psiAHat * BHatMin .* sqrt(THat)) .* dTHatdpsi;
        deltaTMin = -Delta * IHat ./ (psiAHat * BHatMax .* sqrt(THat)) .* dTHatdpsi;
        dnHatdpsi = (ddpsi*(nHat'))';
        deltanMax = -Delta * IHat ./ (psiAHat * BHatMin .* nHat) .* sqrt(THat) .* dnHatdpsi;
        deltanMin = -Delta * IHat ./ (psiAHat * BHatMax .* nHat) .* sqrt(THat) .* dnHatdpsi;
        %}
        
        VPrimeB0OverR0 = ((1./JHat')*thetaWeights)';
        FSAB2 = ((BHat.*BHat./JHat)'*thetaWeights)' ./ VPrimeB0OverR0;
        U = omega * IHat .* dPhiHatdpsi ./ (psiAHat * sqrt(FSAB2 .* THat));
        
        typicalB = sqrt(FSAB2);
        
        %nuPrimeLocal = nuPrimeMeaningful * nHat ./ (THat.*THat);
        %nuStarLocal = nuPrimeLocal * (Miller_A^1.5);
        nuPrimeLocal = nu_r * Miller_q .* nHat ./ (THat.*THat);
        nuStarLocal = nuPrimeLocal .* (Miller_A' .^1.5);
        
        deltaTs = zeros(numSpecies, Npsi);
        deltaetas = zeros(numSpecies, Npsi);
        deltans = zeros(numSpecies, Npsi);
        for ispecies = 1:numSpecies
            prefactor = -Delta*sqrt(THats(ispecies,:)*masses(ispecies)).*IHat./(psiAHat*charges(ispecies)*typicalB);
            deltaTs(ispecies,:) = prefactor./THats(ispecies,:).*dTHatdpsis(ispecies,:);
            deltaetas(ispecies,:) = prefactor./etaHats(ispecies,:).*detaHatdpsis(ispecies,:);
            deltans(ispecies,:) = prefactor./nHats(ispecies,:).*dnHatdpsis(ispecies,:);
        end
        
        drdpsi = psiAHat/Delta*sqrt(FSAB2./THat)./IHat;
        drdpsi = drdpsi(:);
        ddpsiForr = ddpsi;
        ddpsiForr(1,:)=0;
        ddpsiForr(1,1)=1;
        r = ddpsiForr \ drdpsi;
        %r = r - max(r)/2;
        % shift r so r=0 at psiMid:
        if (mod(Npsi,2)==1)
            r = r - r((Npsi+1)/2);
        else
            r = r - (r(Npsi/2) + r(Npsi/2+1))/2;
        end
        
        % Generate abscissae, quadrature weights, and derivative matrix for x grid.
        % Both abscissae and weights should be row vectors.
        switch xGridScheme
            case 0
                powerOfX = 0;
                scale = 1;
                pointAtZero = false;
                [x_i, ddx, d2dx2, xWeights_i] = spectralNodesWeightsAndDifferentiationMatricesForV(Nx, powerOfX, scale, pointAtZero);
            case 1
                % Chebyshev
                xMin=0.1;
                xMaxPolynomial = 5;
                NxIntervals = 1;
                [x_i, xWeights_i, ddx] = multiChebyshevWeightsAndDifferentiation(Nx, xMin, xMaxPolynomial, NxIntervals);
                d2dx2 = ddx*ddx;
            case 2
                xMin = xMaxForDistribution / (Nx+1) /5;
                scheme = 12;
                [x_i, xWeights_i, ddx, d2dx2] = differentiationMatricesForUniformGrid(Nx+1, xMin, xMaxForDistribution, scheme);
                x_i(end)=[];
                xWeights_i(end)=[];
                ddx = ddx(1:(end-1), 1:(end-1));
                d2dx2 = d2dx2(1:(end-1), 1:(end-1));
            case 3
                % Same as 2, but set the last 2 points of f to zero.
                xMin = xMaxForDistribution / (Nx+1) /5;
                scheme = 12;
                [x_i, xWeights_i, ddx, d2dx2] = differentiationMatricesForUniformGrid(Nx+2, xMin, xMaxForDistribution, scheme);
                x_i(end-1:end)=[];
                xWeights_i(end-1:end)=[];
                ddx = ddx(1:(end-2), 1:(end-2));
                d2dx2 = d2dx2(1:(end-2), 1:(end-2));
            case 4
                xMin = xMaxForDistribution / (Nx+1) /5;
                scheme = 12;
                [x_i, xWeights_i, dgdx, d2gdx2] = differentiationMatricesForUniformGrid(Nx, xMin, xMaxForDistribution, scheme);
                exp_px2 = exp(x_i.^2);
                exp_mx2 = exp(-x_i.^2);
                ddx = diag(exp_mx2)*dgdx*diag(exp_px2) - 2*diag(x_i);
                d2dx2 = diag(exp_mx2)*d2gdx2*diag(exp_px2) - 4*diag(x_i.*exp_mx2)*dgdx*diag(exp_px2) + diag(4*x_i.^2-2);
            otherwise
                error('Invalid xGridScheme')
        end
        
        function y=weight(x)
            x2temp=x.*x;
            y=exp(-x2temp);
        end
        
        
        x_i = x_i(:)';
        xWeights_i = xWeights_i(:)';
        
        xMaxPolynomial = max(x_i);
        xMaxPotentials = max([xMax, xMaxPolynomial]);
        NxPotentials = ceil(xMax * NxPotentialsPerVth);
        % Uniform, higher order FD
        scheme = 12;
        xMin=0;
        [xPotentials, xWeightsPotentials, ddxPotentials, d2dx2Potentials] = differentiationMatricesForUniformGrid(NxPotentials, xMin, xMaxPotentials, scheme);
        xPotentials = xPotentials(:)';
        regridPolynomialToUniform = polynomialInterpolationMatrix(x_i,xPotentials,weight(x_i),weight(xPotentials));
        
        
        if plotVelocitySpaceGrid
            figure(figureOffset+7)
            plot(xPotentials,zeros(size(xPotentials))+iteration,'.r')
            hold on
            plot(x_i, zeros(size(x_i))+iteration,'o')
            title('Velocity-space grid')
            xlabel('x')
            ylabel('iteration')
        end
        
        if drawIntroFigures
            colors = [1,0,0;
                0,0.7,0;
                0,0,1;
                0,0,0];
                
            figure(6+figureOffset)
            clf
            numRows = 4;
            numCols = 5;
            numContours = 10;
            
            figureNum=1;
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            contourf(psi, theta, BHat, numContours, 'EdgeColor','none')
            colorbar
            xlabel('\psi_N')
            ylabel('\theta')
            title('BHat')
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            contourf(psi, theta, dBHatdtheta, numContours, 'EdgeColor','none')
            colorbar
            xlabel('\psi_N')
            ylabel('\theta')
            title('d BHat / d \theta')
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            contourf(psi, theta, dBHatDPsi, numContours, 'EdgeColor','none')
            colorbar
            xlabel('\psi_N')
            ylabel('\theta')
            title('d BHat / d \psi_N')
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            contourf(psi, theta, JHat, numContours, 'EdgeColor','none')
            colorbar
            xlabel('\psi_N')
            ylabel('\theta')
            title('JHat')
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,THats(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('\psi_N')
            title('THat')
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            axis tight
            ylim([0,Inf])
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,dTHatdpsis(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
                %plot(psi,(ddpsi*THats(ispecies,:)')','x-','Color',colors(ispecies,:),'DisplayName','')
            end
            xlabel('\psi_N')
            title('d THat / d \psi_N')
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            axis tight
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,IHat,'.-')
            xlabel('\psi_N')
            title('IHat')
            axis tight
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,dIHatdpsi,'.-')
            xlabel('\psi_N')
            title('d IHat / d \psi_N')
            axis tight
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,nHats(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('\psi_N')
            title('nHat')
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            axis tight
            ylim([0,Inf])
            
            %{
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(r,nHats(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('r')
            title('nHat')
            axis tight
            ylim([0,Inf])
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            %}
            
            %{
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            if radialProfiles==4
                analyticHeatFlux = (1/3)*(4*U.^8+16*U.^6+24*U.^4+12*U.^2+3)./(2*U.^4+2*U.^2+1).*exp(-U.*U).*nHat.*(THat.^1.5)...
                    .*(dTHatdpsi - Delta/psiAHat*sqrt(THat./FSAB2).*IHat.*U.*(ddpsi*dTHatdpsi')');
                
                plot(psi,analyticHeatFlux,'.-')
                title('Analytic radial heat flux')
            else
                plot(psi,nHat.^exponent .* dTHatdpsi,'.-')
                title(['nHat^{',num2str(exponent),'} dTHat/d\psi'])
            end
            xlabel('\psi_N')
            %}
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,etaHats(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('\psi_N')
            title('\eta Hat')
            axis tight
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,-dPhiHatdpsi*2/minorRadius,'.-')
            xlabel('\psi_N')
            title('E_r (kV/m)')
            hold on
            axis tight
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,dPhiHatdpsi,'.-')
            xlabel('\psi_N')
            title('d \Phi Hat / d \psi_N')
            hold on
            plot(psi, (ddpsi*PhiHat')', 'x-r')
            axis tight
            
            if radialProfiles > 1
                subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
                plot(psi,PhiHat,'.-')
                xlabel('\psi_N')
                title('\Phi Hat')
                axis tight
            end
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,nuPrimeLocal,'.-')
            xlabel('\psi_N')
            hold on
            plot(psi,nuStarLocal,'.-r')
            legend('\nu''', '\nu_*')
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            semilogy(psi,nuPrimeLocal,'.-')
            xlabel('\psi_N')
            hold on
            plot(psi,nuStarLocal,'.-r')
            legend('\nu''', '\nu_*')
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,deltaTs(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('\psi_N')
            title('\delta_T = \rho_\theta/r_T  (Should be << 1)')
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,deltans(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('\psi_N')
            title('\delta_n = \rho_\theta/r_n')
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            for ispecies = 1:numSpecies
                plot(psi,deltaetas(ispecies,:),'.-','Color',colors(ispecies,:),'DisplayName',['species ',num2str(ispecies)])
                hold on
            end
            xlabel('\psi_N')
            title('\delta_\eta = \rho_\theta/r_\eta')
            if numSpecies>1
                legend show
                set(legend(),'FontSize',legendSize)
            end
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,r,'.-')
            xlabel('\psi')
            ylabel('r')
            axis tight
            
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,U,'.-')
            xlabel('\psi')
            ylabel('U')

            %{
            subplot(numRows,numCols,figureNum); figureNum = figureNum+1;
            plot(psi,FSAB2,'.-')
            xlabel('\psi')
            ylabel('<B^2>')
            %}
        end
        
        figure(80)
        clf
        set(gcf,'Color','w','Units','in','Position',[1,1,7,7])
        numPlots = 8;
    
        % Create textbox
        annotation(gcf,'textbox',...
            [0.227832512315271 0.95 0.105911330049261 0.05],...
            'String',{'Inputs'},...
            'FontSize',15,...
            'FitBoxToText','off',...
            'EdgeColor','none');
        
        % Create textbox
        annotation(gcf,'textbox',...
            [0.665320197044335 0.95 0.105911330049261 0.05],...
            'String',{'Outputs'},...
            'FontSize',15,...
            'FitBoxToText','off',...
            'EdgeColor','none');
        
        leftMargin = 0.09;
        topMargin = 0.06;
        bottomMargin = 0.08;
        interPlotMargin = 0.01;
        
        plotWidth = 0.38;
        plotHeight = (1-topMargin-bottomMargin+interPlotMargin)/numPlots-interPlotMargin;
        plotNum = 0;
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,1000*THats(1,:),'.-')
        axis tight
        set(gca,'XTickLabel',[])
        ylabel('T_i (eV)')
        ylim([0,Inf])
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,nHats(1,:),'.-')
        axis tight
        set(gca,'XTickLabel',[])
        ylabel('n (10^{20} m^{-3})')
        ylim([0,Inf])
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,-dPhiHatdpsi*2/minorRadius,'.-')
        ylabel('E_r (kV/m)')
        axis tight
        set(gca,'XTickLabel',[])
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,etaHats(1,:),'.-')
        axis tight
        set(gca,'XTickLabel',[])
        ylabel('\eta (10^{20} m^{-3})')
        ylim([0,Inf])
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,deltaTs(1,:),'.-')
        axis tight
        set(gca,'XTickLabel',[])
        ylabel('\delta_{Ti} = \rho_\theta/r_{Ti}')
        xlabel('\Psi_N')
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,deltans(1,:),'.-')
        axis tight
        set(gca,'XTickLabel',[])
        ylabel('\delta_n = \rho_\theta/r_n')
        xlabel('\Psi_N')
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        plot(psi,deltaetas(1,:),'.-')
        set(gca,'XTickLabel',[])
        axis tight
        ylabel('\delta_\eta = \rho_\theta/r_\eta')
        
        plotNum = plotNum+1;
        subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
        semilogy(psi,nuStarLocal,'.-r')
        hold on
        semilogy(psi,nuPrimeLocal,'.-')
        axis tight
        ylim([0.1,100])
        set(gca,'YMinorTick','on','YTick',[0.1,1,10,100])
        h=legend('\nu_*', '\nu''');
        set(h,'FontSize',8,'box','off')
        xlabel('\psi_N')

        
        localMatrixSize = numSpecies * Ntheta * Nx * Nxi;
        globalMatrixSize = localMatrixSize * Npsi + numSpecies*Npsi*2;
        %estimated_nnz = floor(1.1*(4*Nx*Nx*Nxi*Ntheta + 5*Nxi*nnz(ddtheta)*Nx+5*Nxi*Nx*Ntheta));
        estimated_nnz = floor(1.1*(4*Nx*Nx*Nxi*Ntheta*numSpecies + 5*Nxi*nnz(ddtheta)*Nx*numSpecies + 5*Nxi*Nx*Ntheta*numSpecies + numSpecies*numSpecies*Nx*Nx*Nxi*Ntheta));
        fprintf('localMatrixSize: %d.   globalMatrixSize: %d.\n',localMatrixSize,globalMatrixSize)
        
        if iteration==1 && drawIntroFigures && geometry==1
            figure(figureOffset+4);
            clf
            numRows=2;
            numCols=3;
            
            subplot(numRows,numCols,1)
            plot(RHat(theta1D), ZHat(theta1D),'r')
            xlabel('R/R_0')
            ylabel('Z/R_0')
            title('Flux surface shape')
            axis equal
            
            subplot(numRows,numCols,2)
            plot(theta1D, BPoloidal(theta1D))
            xlabel('\theta')
            ylabel('B_p / B_0')
            xlim([0, 2*pi])
            
            subplot(numRows,numCols,3)
            plot(theta1D, b(theta1D))
            xlabel('\theta')
            ylabel('b')
            xlim([0, 2*pi])
            
            subplot(numRows,numCols,4)
            plot(theta, dBHatdtheta,'.')
            xlabel('\theta')
            ylabel('db/d\theta')
            xlim([0, 2*pi])
            hold on
            bb=b(theta1D);
            dtheta=theta1D(2)-theta1D(1);
            plot(theta1D, (circshift(bb,[0,-1]) - circshift(bb,[0,1]))/(2*dtheta))
            
            subplot(numRows,numCols,5)
            plot(theta1D, BDotGradTheta(theta1D))
            xlabel('\theta')
            ylabel('R_0 / B_0 * B dot grad \theta')
            xlim([0, 2*pi])
        end
        
        if runMode==-1
            return
        end
        
        % Begin timer for matrix construction:
        tic
        constructionStartTime=tic;
        
        
        % *************************************************************
        % *************************************************************
        %
        % Build the right-hand side:
        %
        % *************************************************************
        % *************************************************************
        
        x2=x_i.*x_i;
        expx2=exp(-x2);
        rhs=zeros(localMatrixSize,Npsi);
        
        for ipsi = 1:Npsi
            for ispecies = 1:numSpecies
                for itheta = 1:Ntheta
                    stuffToAdd = masses(ispecies)^2 / (2*pi*sqrtpi*charges(ispecies)*sqrt(THats(ispecies, ipsi))) ...
                        * nHats(ispecies, ipsi) * IHat(ipsi) * JHat(itheta, ipsi) / (BHat(itheta,ipsi)^3 * psiAHat) ...
                        * dBHatdtheta(itheta, ipsi) * x2 .* expx2 .* ...
                        (1/nHats(ispecies, ipsi) * dnHatdpsis(ispecies, ipsi) + 2*charges(ispecies)*omega/(THats(ispecies, ipsi)*Delta) * dPhiHatdpsi(ipsi) ...
                        + (x2-3/2)/THats(ispecies, ipsi)*dTHatdpsis(ispecies, ipsi));
                    
                    L = 0;
                    indices = (ispecies-1)*Ntheta*Nxi*Nx + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    rhs(indices, ipsi) = stuffToAdd * (4/3);
                    
                    L = 2;
                    indices = (ispecies-1)*Ntheta*Nxi*Nx + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    rhs(indices, ipsi) = stuffToAdd * (2/3);
                end
            end
        end
        
        assignin('base','rhsm',rhs)
                
        % *************************************************************
        % *************************************************************
        %
        % Done with the right-hand side.
        % Next, build the local matrices:
        %
        % *************************************************************
        % *************************************************************
        
        sparseCreatorIndex=1;
        sparseCreator_i=0;
        sparseCreator_j=0;
        sparseCreator_s=0;
        resetSparseCreator()
        
        localMatrices = cell(Npsi,1);
        localMatricesForPreconditioner = cell(Npsi,1);
        
        fprintf('Building local matrices: ')
        
        if tryIterativeSolver
            matricesToMake = [0,1];
        else
            matricesToMake = 0;
        end
        % Options for whichMatrix:
        % 0 = full matrix
        % 1 = preconditioner
        
        for whichMatrix = matricesToMake
            if whichMatrix==0
                ipsiMin = 0;
                ipsiMax = Npsi+1;
                ddthetaToUse = ddtheta;
            else
                ipsiMin = 1;
                ipsiMax = Npsi;
                ddthetaToUse = ddthetaForPreconditioner;
            end
            for ipsiIncludingBoundaries = ipsiMin:ipsiMax
                ipsi = ipsiIncludingBoundaries;
                if ipsi==0
                    ipsi=1;
                end
                if ipsi==Npsi+1
                    ipsi=Npsi;
                end
                fprintf('%d ',ipsi)
                if mod(ipsi,30)==0
                    fprintf('\n')
                end
                
                localDelta = Delta;
                localOmega = omega;
                if ipsiIncludingBoundaries == 0 || ipsiIncludingBoundaries==Npsi+1
                    localDelta = 0;
                    localOmega = 0;
                end
                if whichMatrix==1 && preconditionerMethod_psi==4
                    localDelta=0;
                    localOmega=0;
                end
                if makeLocalApproximation
                    localDelta=0;
                    localOmega=0;
                end
                
                for ispecies = 1:numSpecies
                    speciesFactor = sqrt(masses(ispecies))/charges(ispecies);
                    thetaPartOfLocalStreamingTerm = sqrt(THats(ispecies,ipsi))*diag(JHat(:,ipsi)./BHat(:,ipsi))*ddthetaToUse;
                    thetaPartOfLocalMirrorTerm = -0.5*sqrt(THats(ispecies,ipsi))*JHat(:,ipsi) .* dBHatdtheta(:,ipsi)./(BHat(:,ipsi).^2);
                    thetaPartOfGlobalDDXiTerm = 1/(2*psiAHat) * JHat(:,ipsi) .* dBHatdtheta(:,ipsi) ./ (BHat(:,ipsi).^3);
                    thetaPartOfOffDiagonalGlobalDDThetaTerm = speciesFactor*localDelta*THats(ispecies,ipsi)/psiAHat*diag(JHat(:,ipsi)./(BHat(:,ipsi).^2) ...
                        .*(0.5*IHat(ipsi)*dBHatDPsi(:,ipsi)./BHat(:,ipsi) - dIHatdpsi(ipsi)))*ddthetaToUse;
                    poloidalExBTerm = sqrt(masses(ispecies))*localOmega*IHat(ipsi)/psiAHat*dPhiHatdpsi(ipsi)*diag(JHat(:,ipsi)./(BHat(:,ipsi).^2))*ddthetaToUse;
                    
                    for ix=1:Nx
                        
                        x=x_i(ix);
                        
                        thetaAndXPartOfGlobalDDXiTerm = sqrt(masses(ispecies))*(localOmega*IHat(ipsi)*dPhiHatdpsi(ipsi) + localDelta/charges(ispecies)*x*x*THats(ispecies,ipsi)*dIHatdpsi(ipsi)) * thetaPartOfGlobalDDXiTerm;
                        
                        for L=0:(Nxi-1)
                            rowIndices = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + (1:Ntheta);
                            
                            % super-diagonals:
                            if L<(Nxi-1) % && whichMatrixToMake ~= 3
                                colIndices = rowIndices + Ntheta;
                                % Streaming term
                                addSparseBlock(rowIndices, colIndices, x*thetaPartOfLocalStreamingTerm*(L+1)/(2*L+3));
                                % Mirror term
                                addToSparse(rowIndices, colIndices, x*thetaPartOfLocalMirrorTerm*(L+1)*(L+2)/(2*L+3));
                            end
                            
                            % Sub-diagonals:
                            if L>0
                                colIndices = rowIndices - Ntheta;
                                % Streaming term
                                addSparseBlock(rowIndices, colIndices, x*thetaPartOfLocalStreamingTerm*L/(2*L-1));
                                % Mirror term
                                addToSparse(rowIndices, colIndices, -x*thetaPartOfLocalMirrorTerm*(L-1)*L/(2*L-1));
                            end
                            
                            if includeNewStreamingAndMirrorTerms
                                % Terms that are diagonal in L:
                                colIndices = rowIndices;
                                addToSparse(rowIndices, colIndices, thetaAndXPartOfGlobalDDXiTerm*L*(L+1)/((2*L+3)*(2*L - 1)));
                                stuffToAdd = poloidalExBTerm + speciesFactor * localDelta*x*x*THats(ispecies,ipsi)/(psiAHat*(2*L+3)*(2*L-1)) ...
                                    * diag(JHat(:,ipsi)./(BHat(:,ipsi).*BHat(:,ipsi)).*( (3*L*L+3*L-2)*IHat(ipsi)*dBHatDPsi(:,ipsi)./BHat(:,ipsi) - (2*L*L+2*L-1)*dIHatdpsi(ipsi) ))*ddthetaToUse;
                                addSparseBlock(rowIndices, colIndices, stuffToAdd)
                                
                                if (whichMatrix==0 || preconditionerMethod_xi==0)
                                    % super-super-diagonals:
                                    if L<(Nxi-2)
                                        % ell = L+2:
                                        colIndices = rowIndices + Ntheta*2;
                                        addSparseBlock(rowIndices, colIndices, (L+2)*(L+1)/((2*L+5)*(2*L+3))*x*x*thetaPartOfOffDiagonalGlobalDDThetaTerm)
                                        addToSparse(rowIndices, colIndices, thetaAndXPartOfGlobalDDXiTerm*(L+1)*(L+2)*(L+3)/((2*L+5)*(2*L+3)));
                                    end
                                    
                                    % Sub-sub-diagonals:
                                    if L>1
                                        % ell = L-2:
                                        colIndices = rowIndices - Ntheta*2;
                                        addSparseBlock(rowIndices, colIndices, (L-1)*L/((2*L-3)*(2*L-1))*x*x*thetaPartOfOffDiagonalGlobalDDThetaTerm)
                                        addToSparse(rowIndices, colIndices, thetaAndXPartOfGlobalDDXiTerm*(-L)*(L-1)*(L-2)/((2*L-3)*(2*L-1)));
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Now add the collisionless xDot term:
                for ispecies = 1:numSpecies
                    xPartOfDDXTerm_original = sqrt(masses(ispecies))*diag(localDelta/charges(ispecies)*dTHatdpsis(ispecies,ipsi)*0.5*(x_i.^3) + localOmega*dPhiHatdpsi(ipsi)*x_i) * ddx;
                    thetaPartOfDDXTerm = IHat(ipsi)/(2*psiAHat)*JHat(:,ipsi).*dBHatdtheta(:,ipsi)./(BHat(:,ipsi).^3);
                    for L=0:(Nxi-1)
                        xPartOfDDXTerm = xPartOfDDXTerm_original;
                        if whichMatrix==1  && L >= preconditioner_x_min_L
                            switch preconditionerMethod_x
                                case 0
                                    % Nothing to do here.
                                case 1
                                    xPartOfDDXTerm = diag(diag(xPartOfDDXTerm));
                                case 2
                                    xPartOfDDXTerm = triu(xPartOfDDXTerm);
                                case 3
                                    xPartOfDDXTerm = tril(xPartOfDDXTerm);
                                otherwise
                                    error('Invalid preconditionerMethod_x')
                            end
                        end
                        
                        for itheta=1:Ntheta
                            rowIndices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                            
                            % Terms that are diagonal in L:
                            colIndices = rowIndices;
                            addSparseBlock(rowIndices, colIndices, 2*(3*L*L+3*L-2)/((2*L+3)*(2*L-1))*thetaPartOfDDXTerm(itheta)*xPartOfDDXTerm);
                            
                            if (whichMatrix==0 || preconditionerMethod_xi==0)
                                % super-super-diagonals:
                                if L<(Nxi-2)
                                    % ell = L+2:
                                    colIndices = rowIndices + Ntheta*2;
                                    addSparseBlock(rowIndices, colIndices, (L+2)*(L+1)/((2*L+5)*(2*L+3))*thetaPartOfDDXTerm(itheta)*xPartOfDDXTerm);
                                end
                                
                                % Sub-sub-diagonals:
                                if L>1
                                    % ell = L-2:
                                    colIndices = rowIndices - Ntheta*2;
                                    addSparseBlock(rowIndices, colIndices, L*(L-1)/((2*L-3)*(2*L-1))*thetaPartOfDDXTerm(itheta)*xPartOfDDXTerm);
                                end
                            end
                        end
                    end
                end

                
                % ******************************************************
                % Build the collision operator:
                % ******************************************************
                
                xWith0s = [0, xPotentials(2:(end-1)), 0];
                M21 = 4*pi*diag(xWith0s.^2) * regridPolynomialToUniform;
                xWith0s = [0, xPotentials(2:(end-1)), 0];
                M32 = -2*diag(xWith0s.^2);
                LaplacianTimesX2WithoutL = diag(xPotentials.^2)*d2dx2Potentials + 2*diag(xPotentials)*ddxPotentials;
                                
                erfs=erf(x_i);
                x2 = x_i.*x_i;
                x3 = x2.*x_i;
                expx2 = exp(-x_i.*x_i);
                
                
                CE = zeros(Nx, Nx, numSpecies);
                nuD = zeros(Nx, numSpecies);
                regridSpecies = zeros(Nx, Nx, numSpecies, numSpecies);
                M12IncludingX0 = zeros(Nx, NxPotentials, numSpecies, numSpecies, NL);
                M13IncludingX0 = zeros(Nx, NxPotentials, numSpecies, numSpecies, NL);
                for speciesA = 1:numSpecies
                    for speciesB = 1:numSpecies
                        speciesFactorTest = 3*sqrtpi/4*nHats(speciesB,ipsi) * charges(speciesA)*charges(speciesA)*charges(speciesB)*charges(speciesB)/(THats(speciesA,ipsi)^(3/2));
                        xb = x_i * sqrt(THats(speciesA,ipsi)*masses(speciesB)/(THats(speciesB,ipsi)*masses(speciesA)));
                        erfs = erf(xb);
                        xb2  = xb.*xb;
                        xb3 = xb2.*xb;
                        expxb2 = exp(-xb2);
                        Psi = (erfs - 2/sqrtpi*xb .* expxb2) ./ (2*xb2);
                        nuD(:,speciesA) = nuD(:,speciesA) + (speciesFactorTest * (erfs - Psi) ./ (x_i.^3))';
                        coefficientOfd2dx2 = Psi./x_i;
                        coefficientOfddx = -2*THats(speciesA,ipsi)*masses(speciesB)/(THats(speciesB,ipsi)*masses(speciesA))*Psi*(1-masses(speciesA)/masses(speciesB)) ...
                            + (erfs - Psi)./(x_i.*x_i);
                        diagonalPartOfCE = 4/sqrtpi*THats(speciesA,ipsi)/THats(speciesB,ipsi)*sqrt(THats(speciesA,ipsi)*masses(speciesB)/(THats(speciesB,ipsi)*masses(speciesA))) .* expxb2;
                        CE(:,:,speciesA) = CE(:,:,speciesA) + speciesFactorTest*(diag(coefficientOfd2dx2)*d2dx2 + diag(coefficientOfddx)*ddx + diag(diagonalPartOfCE));
                        
                        if speciesA==speciesB
                            regridSpecies(:,:,speciesA,speciesB) = eye(Nx);
                        else
                            regridSpecies(:,:,speciesA,speciesB) = polynomialInterpolationMatrix(x_i,xb,weight(x_i),weight(xb));
                        end
                        
                        speciesFactorField = nHats(speciesA,ipsi) * charges(speciesA)*charges(speciesA)*charges(speciesB)*charges(speciesB)...
                            * masses(speciesA) * THats(speciesB,ipsi)/(THats(speciesA,ipsi)^(5/2) * masses(speciesB));
                        for L=0:(NL-1)
                            regridUniformToPolynomial = makeHighOrderUniformRegriddingMatrix(xPotentials,xb,L,'H');
                            M12IncludingX0(:,:,speciesA, speciesB, L+1) = nu_r * 3/(2*pi)*speciesFactorField*diag(expx2)* regridUniformToPolynomial...
                                * (diag(xPotentials*(1-masses(speciesA)/masses(speciesB)))*ddxPotentials + eye(NxPotentials)) ;
                            regridUniformToPolynomial = makeHighOrderUniformRegriddingMatrix(xPotentials,xb,L,'G');
                            M13IncludingX0(:,:,speciesA, speciesB, L+1) = -nu_r * 3/(2*pi) * speciesFactorField * diag(x2.*expx2) * regridUniformToPolynomial* d2dx2Potentials;
                        end
                    end
                end
                                
                for L=0:(Nxi-1)
                    if L <= (NL-1)
                        % Add Rosenbluth potential stuff
                        
                        M22 = LaplacianTimesX2WithoutL-L*(L+1)*eye(NxPotentials);
                        % Add Dirichlet or Neumann boundary condition for
                        % potentials at x=0:
                        if L==0
                            M22(1,:)=ddxPotentials(1,:);
                        else
                            M22(1,:) = 0;
                            M22(1,1) = 1;
                            %M12(:,1) = 0;
                            %M13(:,1) = 0;
                        end
                        M33 = M22;
                        
                        % Add Robin boundary condition for potentials at x=xMax:
                        M22(NxPotentials,:) = xMaxPotentials*ddxPotentials(NxPotentials,:);
                        M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1;
                        
                        % My original b.c.:
                        M33(NxPotentials,:) = xMaxPotentials*xMaxPotentials*d2dx2Potentials(NxPotentials,:) + (2*L+1)*xMaxPotentials*ddxPotentials(NxPotentials,:);
                        M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1);
                        
                        if L~=0
                            M22(NxPotentials,1)=0;
                            M33(NxPotentials,1)=0;
                        end
                        
                        M22BackslashM21 = M22 \ M21;
                        M33BackslashM32 = M33 \ M32;
                        
                    end
                    
                    for speciesA = 1:numSpecies
                        if whichMatrix == 0
                            % We're making the main matrix.
                            speciesBToUse = 1:numSpecies;
                        else
                            % We're making the preconditioner.
                            switch preconditioner_species
                                case 0
                                    % Full inter-species coupling
                                    speciesBToUse = 1:numSpecies;
                                case 1
                                    % No inter-species coupling
                                    speciesBToUse = speciesA;
                                otherwise
                                    error('Invalid preconditioner_species')
                            end
                        end
                        for speciesB = speciesBToUse
                            % Add CD
                            CD = 3*nHats(speciesA,ipsi)*charges(speciesA)*charges(speciesA)*charges(speciesB)*charges(speciesB)...
                                * masses(speciesA)/(masses(speciesB)*THats(speciesA,ipsi)*sqrt(THats(speciesA,ipsi))) ...
                                * diag(expx2) * regridSpecies(:,:,speciesA, speciesB);
                            
                            if speciesA == speciesB
                                M11 = -nu_r * (-0.5*diag(nuD(:,speciesA))*L*(L+1) + CE(:,:,speciesA) + CD);
                            else
                                M11 = -nu_r * CD;
                            end
                            
                            %if false
                            if L <= (NL-1)
                                % Add Rosenbluth potential stuff
                                
                                M13 = M13IncludingX0(:,:,speciesA, speciesB, L+1);
                                M12 = M12IncludingX0(:,:,speciesA, speciesB, L+1);
                                
                                % Add Dirichlet or Neumann boundary condition for
                                % potentials at x=0:
                                if L~=0
                                    M12(:,1) = 0;
                                    M13(:,1) = 0;
                                end
                                
                                KWithoutThetaPart = M11 -  (M12 - M13 * M33BackslashM32) * M22BackslashM21;
                            else
                                KWithoutThetaPart = M11;
                            end
                            
                            % The lines below are invoked to make the local preconditioner.
                            if whichMatrix == 1 && L >= preconditioner_x_min_L
                                switch preconditionerMethod_x
                                    case 0
                                        % Nothing to do here.
                                    case 1
                                        KWithoutThetaPart = diag(diag(KWithoutThetaPart));
                                    case 2
                                        KWithoutThetaPart = triu(KWithoutThetaPart);
                                    case 3
                                        KWithoutThetaPart = tril(KWithoutThetaPart);
                                    otherwise
                                        error('Invalid preconditionerMethod_x')
                                end
                                
                            end
                            
                            for itheta=1:Ntheta
                                rowIndices = (speciesA-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                                colIndices = (speciesB-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                                addSparseBlock(rowIndices, colIndices, KWithoutThetaPart)
                            end
                            
                            
                        end
                    end
                    
                end
            
                
                switch whichMatrix
                    case 0
                        if ipsiIncludingBoundaries==0
                            leftLocalMatrix = createSparse(localMatrixSize, localMatrixSize);
                        elseif ipsiIncludingBoundaries==Npsi+1
                            rightLocalMatrix = createSparse(localMatrixSize, localMatrixSize);
                        else
                            localMatrices{ipsi} = createSparse(localMatrixSize, localMatrixSize);
                        end
                    case 1
                        localMatricesForPreconditioner{ipsi} = createSparse(localMatrixSize, localMatrixSize);
                    otherwise
                        error('Invalid whichMatrix')
                end
            end
            
            % *****************************************************
            % Done building the local matrices.
            % Now build the big global matrix:
            % *****************************************************
            
            estimated_nnz = nnz(localMatrices{2})*Npsi + 3*nnz(ddpsi)*Nx*Nxi*Ntheta*numSpecies ...
                + 2*numSpecies*Nx*Ntheta + 2*numSpecies*Nx*Ntheta;
            %                + Npsi + max([nnz(ddpsi),nnz(d2dpsi2)])*2*Nx*Ntheta + Npsi*Nx*Ntheta;
            resetSparseCreator()
            
            % For the interior points, add everything but the d/dpsi term:
            for ipsi = 1:Npsi
                switch whichMatrix
                    case 0
                        [rowIndices, colIndices, values] = find(localMatrices{ipsi});
                    case 1
                        [rowIndices, colIndices, values] = find(localMatricesForPreconditioner{ipsi});
                    otherwise
                        error('Invalid whichMatrix')
                end
                rowIndices = rowIndices + (ipsi-1)*localMatrixSize;
                colIndices = colIndices + (ipsi-1)*localMatrixSize;
                addToSparse(rowIndices, colIndices, values)
            end
            
            % Finally, add the d/dpsi term:
            includeddpsiTermThisTime = includeddpsiTerm;
            if makeLocalApproximation
                includeddpsiTermThisTime = false;
            end
            if whichMatrix==1
                switch preconditionerMethod_psi
                    case 0
                        ddpsiToUse = ddpsi;
                        d2dpsi2ToUse = d2dpsi2;
                    case 1
                        ddpsiToUse = ddpsiPreconditioner;
                        %d2dpsi2ToUse = d2dpsi2Preconditioner;
                        d2dpsi2ToUse = d2dpsi2;
                    case {2,4}
                        includeddpsiTermThisTime = false;
                    case 3
                        ddpsiToUse = diag(diag(ddpsi));
                        d2dpsi2ToUse = diag(diag(d2dpsi2));
                    otherwise
                end
            else
                ddpsiToUse = ddpsi;
                d2dpsi2ToUse = d2dpsi2;
            end
            if includeddpsiTermThisTime
                for ispecies = 1:numSpecies
                    for itheta = 1:Ntheta
                        if psiGridMode==3 || psiGridMode == 4
                            % psiDot below is not the full time derivative of
                            % psi, but should have the same sign.
                            psiDot = IHat(1) * JHat(itheta,1) * dBHatdtheta(itheta,1) / psiAHat;
                            if psiDot * upwindingSign > 0
                                ddpsiToUse = ddpsiLeft;
                            else
                                ddpsiToUse = ddpsiRight;
                            end
                        end
                        for ix = 1:Nx
                            for L=0:(Nxi-1)
                                rowIndices = ((1:Npsi)-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta;
                                
                                stuffToAddWithoutTheXiPart = -sqrt(masses(ispecies))/charges(ispecies)*0.5*Delta/psiAHat*x2(ix) ...
                                    * diag(JHat(itheta,:) .* THats(ispecies,:) .* IHat .* dBHatdtheta(itheta,:) ./ (BHat(itheta,:).^3)) * ddpsiToUse;
                                stuffToAddWithoutTheXiPart = stuffToAddWithoutTheXiPart(1:Npsi, :);
                                
                                % Stuff that is diagonal in L:
                                ell = L;
                                colIndices = ((1:Npsi)-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + ell*Ntheta + itheta;
                                stuffToAdd = 2*(3*L*L+3*L-2)/((2*L+3)*(2*L-1)) * stuffToAddWithoutTheXiPart;
                                addSparseBlock(rowIndices, colIndices, stuffToAdd)
                                
                                if (whichMatrix==0 || preconditionerMethod_xi==0)
                                    % Super-super diagonal in L:
                                    if L<(Nxi-2)
                                        ell = L+2;
                                        colIndices = ((1:Npsi)-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + ell*Ntheta + itheta;
                                        addSparseBlock(rowIndices, colIndices, (L+2)*(L+1)/((2*L+5)*(2*L+3)) * stuffToAddWithoutTheXiPart)
                                    end
                                    
                                    % Sub-sub diagonal in L:
                                    if L>1
                                        ell = L-2;
                                        colIndices = ((1:Npsi)-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + ell*Ntheta + itheta;
                                        addSparseBlock(rowIndices, colIndices, (L-1)*L/((2*L-3)*(2*L-1)) * stuffToAddWithoutTheXiPart)
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            
            %**********************************************
            % Require < density > = 0
            %**********************************************
            
            xPart = (x_i.^2) .* xWeights_i;
            xiPart = zeros(1,Nxi);
            xiPart(1) = 1;
            
            for ispecies = 1:numSpecies
                for ipsi = 1:Npsi
                    thetaPart = thetaWeights(:) ./ (JHat(:,ipsi));
                    stuffToAdd = kron(xPart, kron(xiPart, thetaPart'));
                    rowIndex = localMatrixSize*Npsi + (ipsi-1)*2*numSpecies + (ispecies-1)*2 + 1;
                    colIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (1:(Nx*Nxi*Ntheta));
                    addSparseBlock(rowIndex, colIndices, stuffToAdd)
                end
            end
            
            %**********************************************
            % Require < pressure > = 0
            %**********************************************
            
            xPart = (x_i.^4) .* xWeights_i;
            xiPart = zeros(1,Nxi);
            xiPart(1) = 1;
            
            for ispecies = 1:numSpecies
                for ipsi = 1:Npsi
                    thetaPart = thetaWeights(:) ./ (JHat(:,ipsi));
                    stuffToAdd = kron(xPart, kron(xiPart, thetaPart'));
                    %rowIndex = localMatrixSize*Npsi + Npsi + ipsi;
                    %colIndices = (ipsi-1)*localMatrixSize + (1:localMatrixSize);
                    rowIndex = localMatrixSize*Npsi + (ipsi-1)*2*numSpecies + (ispecies-1)*2 + 2;
                    colIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (1:(Nx*Nxi*Ntheta));
                    addSparseBlock(rowIndex, colIndices, stuffToAdd)
                end
            end
            
            
            %**********************************************
            % Add particle source:
            %**********************************************
            
            switch sourcePoloidalVariation
                case 0
                    sourceThetaProfile = ones(size(theta));
                case 1
                    sourceThetaProfile = 1+cos(theta);
                case 2
                    sourceThetaProfile = 1 + source_a * cos(theta);
                case 3
                    sourceThetaProfile = 1 + source_a * sin(theta);
                otherwise
                    error('Invalid sourcePoloidalVariation')
            end
            
            xPartOfSource = (x2-5/2) .* exp(-x2);
            L=0;
            
            for ispecies = 1:numSpecies
                for ipsi=1:Npsi
                    for itheta=1:Ntheta
                        rowIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                        colIndex = localMatrixSize*Npsi + (ipsi-1)*2*numSpecies + (ispecies-1)*2 + 1;
                        addSparseBlock(rowIndices, colIndex, sourceThetaProfile(itheta) * xPartOfSource(:))
                    end
                end
            end
            
            %**********************************************
            % Add heat source:
            %**********************************************
            
            xPartOfSource = (x2-3/2) .* exp(-x2);
            L=0;
            
            for ispecies = 1:numSpecies
                for ipsi=1:Npsi
                    for itheta=1:Ntheta
                        rowIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                        colIndex = localMatrixSize*Npsi + (ipsi-1)*2*numSpecies + (ispecies-1)*2 + 2;
                        addSparseBlock(rowIndices, colIndex, sourceThetaProfile(itheta) * xPartOfSource(:))
                    end
                end
            end
            
            switch whichMatrix
                case 0
                    matrix = createSparse(globalMatrixSize, globalMatrixSize);
                case 1
                    preconditionerMatrix = createSparse(globalMatrixSize, globalMatrixSize);
            end
            
        end
        
        fprintf('Done. Time to build matrices: %g sec.\n',toc(constructionStartTime))
        
        fprintf('Predicted nnz for global matrix: %d.    Actual nnz: %d\n',estimated_nnz,nnz(matrix))

        % *****************************************************
        % Zero out rows of the global matrix in which trajectories
        % enter the domain:
        % *****************************************************
        
        for ipsi=2:(Npsi-1)
            indices = (ipsi-1)*localMatrixSize + (1:localMatrixSize);
            addToSparse(indices, indices, ones(size(indices)))
        end
        indices = Npsi*localMatrixSize + (1:(2*Npsi*numSpecies));
        addToSparse(indices, indices, ones(size(indices)))
        
        for ispecies = 1:numSpecies
            ipsi=1;
            indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (1:(Nx*Nxi*Ntheta));
            signOfPsiDot = -IHat(ipsi)*JHat(:,ipsi).*dBHatdtheta(:,ipsi)/(psiAHat*charges(ispecies));
            if imposeLocalSolutionWhereTrajectoriesAreParallel
                trajectoriesEnteringDomain = signOfPsiDot > -tiny;
            else
                trajectoriesEnteringDomain = signOfPsiDot > tiny;
            end
            trajectoriesEnteringDomainBig = logical(kron(true(Nxi*Nx,1),trajectoriesEnteringDomain));
            indices = indices(~trajectoriesEnteringDomainBig);
            addToSparse(indices, indices, ones(size(indices)))
            
            ipsi=Npsi;
            indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta + (1:(Nx*Nxi*Ntheta));
            signOfPsiDot = -IHat(ipsi)*JHat(:,ipsi).*dBHatdtheta(:,ipsi)/(psiAHat*charges(ispecies));
            if imposeLocalSolutionWhereTrajectoriesAreParallel
                trajectoriesEnteringDomain = signOfPsiDot < tiny;
            else
                trajectoriesEnteringDomain = signOfPsiDot < -tiny;
            end
            trajectoriesEnteringDomainBig = logical(kron(true(Nxi*Nx,1),trajectoriesEnteringDomain));
            indices = indices(~trajectoriesEnteringDomainBig);
            addToSparse(indices, indices, ones(size(indices)))
        end
        
        rowsInWhichDKEIsEnforced = createSparse(globalMatrixSize, globalMatrixSize);
        
        onesToAddToDiagonal = speye(globalMatrixSize) - rowsInWhichDKEIsEnforced;
        
        matrix = rowsInWhichDKEIsEnforced*matrix + onesToAddToDiagonal;
        if tryIterativeSolver
            preconditionerMatrix = rowsInWhichDKEIsEnforced * preconditionerMatrix + onesToAddToDiagonal;
        end
        
        % *****************************************************
        % This is the end of building the matrices.
        % Next, do some preparation for the global solve:
        % *****************************************************
        
        fprintf('Solving for distribution function at endpoints: left...')
        solutionLeft = leftLocalMatrix \ rhs(:,1);
        fprintf('  right...')
        solutionRight = rightLocalMatrix \ rhs(:,Npsi);
        fprintf('  Done.\n')
        
        % Subtract any density or pressure in the boundary solutions.
        % This is kind of a hack - it would be preferable to just put the
        % corresponding constraints in the local matrix.
        sqrtPiOver4 = xWeights_i*((x_i.^2.*expx2)');
        sqrtPiTimes3Over8 = xWeights_i*((x_i.^4.*expx2)');
        L=0;
        for iteration = 1:10
            for ispecies = 1:numSpecies
                ipsi = 1;
                localDensity = zeros(Ntheta,1);
                localPressure = zeros(Ntheta,1);
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    localDensity(itheta) = x_i.*x_i.*xWeights_i * solutionLeft(indices);
                    localPressure(itheta) = x_i.*x_i.*x_i.*x_i.*xWeights_i * solutionLeft(indices);
                end
                FSADensity = thetaWeights'*(localDensity./JHat(:,ipsi))/VPrimeB0OverR0(ipsi);
                FSAPressure = thetaWeights'*(localPressure./JHat(:,ipsi))/VPrimeB0OverR0(ipsi);
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    solutionLeft(indices) = solutionLeft(indices) - FSADensity*expx2(:)/sqrtPiOver4 - FSAPressure*expx2(:).*(x_i(:).^2-3/2)/sqrtPiTimes3Over8;
                end
                
                ipsi = Npsi;
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    localDensity(itheta) = x_i.*x_i.*xWeights_i * solutionRight(indices);
                    localPressure(itheta) = x_i.*x_i.*x_i.*x_i.*xWeights_i * solutionRight(indices);
                end
                FSADensity = thetaWeights'*(localDensity./JHat(:,ipsi))/VPrimeB0OverR0(ipsi);
                FSAPressure = thetaWeights'*(localPressure./JHat(:,ipsi))/VPrimeB0OverR0(ipsi);
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    solutionRight(indices) = solutionRight(indices) - FSADensity*expx2(:)/sqrtPiOver4 - FSAPressure*expx2(:).*(x_i(:).^2-3/2)/sqrtPiTimes3Over8;
                end
                
            end
        end
        
        
        %solutionLeft = solutionLeft*0;
        %solutionRight = solutionRight*0;
        
        assignin('base','rhslm',rhs(:,1))
        assignin('base','slm',solutionLeft)
        assignin('base','srm',solutionRight)
        assignin('base','llmm',leftLocalMatrix)
        assignin('base','rlmm',rightLocalMatrix)
        
        
        % Build the rhs for the global solve.
        % It's the same as the rhs we already built, except for boundary
        % conditions:
        globalRhs = rhs;
        
        for ispecies = 1:numSpecies
            speciesPart = false(numSpecies,1);
            speciesPart(ispecies) = true;
            
            ipsi=1;
            signOfPsiDot = -IHat(ipsi)*JHat(:,ipsi).*dBHatdtheta(:,ipsi)/(psiAHat*charges(ispecies));
            if imposeLocalSolutionWhereTrajectoriesAreParallel
                trajectoriesEnteringDomain = signOfPsiDot > -tiny;
            else
                trajectoriesEnteringDomain = signOfPsiDot > tiny;
            end
            trajectoriesEnteringDomainBig = logical(kron(speciesPart, kron(true(Nxi*Nx,1),trajectoriesEnteringDomain)));
            globalRhs(trajectoriesEnteringDomainBig,ipsi) = solutionLeft(trajectoriesEnteringDomainBig);
            
            ipsi=Npsi;
            signOfPsiDot = -IHat(ipsi)*JHat(:,ipsi).*dBHatdtheta(:,ipsi)/(psiAHat*charges(ispecies));
            if imposeLocalSolutionWhereTrajectoriesAreParallel
                trajectoriesEnteringDomain = signOfPsiDot < tiny;
            else
                trajectoriesEnteringDomain = signOfPsiDot < -tiny;
            end
            trajectoriesEnteringDomainBig = logical(kron(speciesPart, kron(true(Nxi*Nx,1),trajectoriesEnteringDomain)));
            globalRhs(trajectoriesEnteringDomainBig,ipsi) = solutionRight(trajectoriesEnteringDomainBig);
        end
        
        globalRhs = reshape(globalRhs,[localMatrixSize*Npsi,1]);
        globalRhs = [globalRhs; zeros(2*Npsi*numSpecies,1)];
                
        if tryIterativeSolver
            fprintf('LU-factorizing preconditioner...')
            tic
            [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditionerMatrix);
            fprintf('done.  Took %g seconds.\n',toc)
        end
        
        
        function solnVector=preconditioner(rhsVector)
            solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
        end
        
        % ***********************************************************
        % ***********************************************************
        %
        % Begin the main solve!
        %
        % ***********************************************************
        % ***********************************************************
        
        fprintf('Beginning global solve!  This might take a while...\n')
        
        assignin('base','rhsm',globalRhs)
        assignin('base','mm',matrix)
        if tryIterativeSolver
            assignin('base','pm',preconditionerMatrix)
        end
        
        [soln, didItConverge, residual] ...
            = solveLinearSystem(matrix, globalRhs, @preconditioner, ...
            tryIterativeSolver, orderOfSolversToTry, tol, maxIterations, restart, ...
            figureOffset+3, tryDirectSolverIfIterativeSolversFail);
        
        assignin('base','sm',soln)
        
        % **************************************************************
        % **************************************************************
        %
        % Calculate moments of the distribution function
        % and other output quantities:
        %
        % **************************************************************
        % **************************************************************
        
        solnMatrix = reshape(soln(1:(localMatrixSize*Npsi)),[localMatrixSize,Npsi]);
        
        bigCharges = kron(charges(:), ones(1,Npsi));
        bigMasses = kron(masses(:), ones(1,Npsi));
        bigIHat = kron(ones(numSpecies, 1), IHat);
        
        densityPerturbation = zeros(Ntheta, Npsi, numSpecies);
        flow = zeros(Ntheta, Npsi, numSpecies);
        pressurePerturbation = zeros(Ntheta, Npsi, numSpecies);
        particleFluxBeforeThetaIntegral = zeros(Ntheta, Npsi, numSpecies);
        momentumFluxBeforeThetaIntegral = zeros(Ntheta, Npsi, numSpecies);
        heatFluxBeforeThetaIntegral = zeros(Ntheta, Npsi, numSpecies);
        
        densityIntegralWeight = (x_i.^2);
        densityFactors = 4*pi*Delta*THats.^(3/2) ./ (nHats .* (bigMasses.^(3/2)));
        
        flowIntegralWeight = (x_i.^3);
        flowFactors = 4*pi/3*(THats.^2)./(nHats .* bigMasses .* bigMasses);
        %flowFactors = 4*pi/3*psiAHat * FSAB2 .* THat .* THat ./ (nHat .* IHat .* dTHatdpsi);
        
        pressureIntegralWeight = (x_i.^4);
        pressureFactors = 8*pi*Delta*THats.^(3/2) ./ (3 * nHats .* (bigMasses.^(3/2)));
        
        particleFluxIntegralWeight = (x_i.^4);
        particleFluxFactors = -bigMasses .* bigIHat ./ bigCharges .* ((THats./bigMasses).^(5/2));
        
        momentumFluxIntegralWeight = (x_i.^5);
        momentumFluxFactors = -bigMasses .* (bigIHat.^2) ./ bigCharges .* ((THats./bigMasses).^3);
        
        heatFluxIntegralWeight = (x_i.^6);
        heatFluxFactors = -bigMasses .* THats .* bigIHat ./ bigCharges .* ((THats./bigMasses).^(5/2));
        
        L=0;
        for ispecies = 1:numSpecies
            for ipsi = 1:Npsi
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    fSlice = solnMatrix(indices, ipsi);
                    densityPerturbation(itheta, ipsi, ispecies) = densityFactors(ispecies, ipsi)  * xWeights_i*(densityIntegralWeight' .* fSlice);
                    pressurePerturbation(itheta, ipsi, ispecies) = pressureFactors(ispecies, ipsi)  * xWeights_i*(pressureIntegralWeight' .* fSlice);
                    particleFluxBeforeThetaIntegral(itheta, ipsi, ispecies) = particleFluxFactors(ispecies, ipsi) * (8/3)*xWeights_i*(particleFluxIntegralWeight' .* fSlice);
                    heatFluxBeforeThetaIntegral(itheta, ipsi, ispecies) = heatFluxFactors(ispecies, ipsi) * (8/3)*xWeights_i*(heatFluxIntegralWeight' .* fSlice);
                end
            end
        end
        
        L=1;
        for ispecies = 1:numSpecies
            for ipsi = 1:Npsi
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    fSlice = solnMatrix(indices, ipsi);
                    flow(itheta, ipsi, ispecies) = flowFactors(ispecies, ipsi)  * xWeights_i*(flowIntegralWeight' .* fSlice);
                    momentumFluxBeforeThetaIntegral(itheta, ipsi, ispecies) = (16/15)*momentumFluxFactors(ispecies, ipsi)  * xWeights_i*(momentumFluxIntegralWeight' .* fSlice);
                end
            end
        end
        
        L=2;
        for ispecies = 1:numSpecies
            for ipsi = 1:Npsi
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    fSlice = solnMatrix(indices, ipsi);
                    particleFluxBeforeThetaIntegral(itheta, ipsi, ispecies) = particleFluxBeforeThetaIntegral(itheta, ipsi, ispecies) ...
                        + particleFluxFactors(ispecies, ipsi) * (4/15)*xWeights_i*(particleFluxIntegralWeight' .* fSlice);
                    heatFluxBeforeThetaIntegral(itheta, ipsi, ispecies) = heatFluxBeforeThetaIntegral(itheta, ipsi, ispecies) ...
                        + heatFluxFactors(ispecies, ipsi) * (4/15)*xWeights_i*(heatFluxIntegralWeight' .* fSlice);
                end
            end
        end
        
        L=3;
        for ispecies = 1:numSpecies
            for ipsi = 1:Npsi
                for itheta = 1:Ntheta
                    indices = (ispecies-1)*Nx*Nxi*Ntheta + ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    fSlice = solnMatrix(indices, ipsi);
                    momentumFluxBeforeThetaIntegral(itheta, ipsi, ispecies) = momentumFluxBeforeThetaIntegral(itheta, ipsi, ispecies) ...
                        + (4/35)*momentumFluxFactors(ispecies, ipsi) * xWeights_i*(momentumFluxIntegralWeight' .* fSlice);
                end
            end
        end
        
        %flowBeforeThetaIntegral = zeros(size(particleFluxBeforeThetaIntegral));
        FSADensityPerturbation = zeros(numSpecies,Npsi);
        FSABFlow = zeros(numSpecies,Npsi);
        FSAPressurePerturbation = zeros(numSpecies,Npsi);
        particleFlux = zeros(numSpecies,Npsi);
        momentumFlux = zeros(numSpecies,Npsi);
        heatFlux = zeros(numSpecies,Npsi);
        particleSource = zeros(numSpecies,Npsi);
        heatSource = zeros(numSpecies,Npsi);
        Machs = zeros(Ntheta,Npsi,numSpecies);
        for ispecies = 1:numSpecies
            particleFluxBeforeThetaIntegral(:,:,ispecies) = particleFluxBeforeThetaIntegral(:,:,ispecies) .* dBHatdtheta ./ (BHat.^3);
            momentumFluxBeforeThetaIntegral(:,:,ispecies) = momentumFluxBeforeThetaIntegral(:,:,ispecies) .* dBHatdtheta ./ (BHat.^4);
            heatFluxBeforeThetaIntegral(:,:,ispecies) = heatFluxBeforeThetaIntegral(:,:,ispecies) .* dBHatdtheta ./ (BHat.^3);
            %flowBeforeThetaIntegral(:,:,ispecies) = BHat ./ JHat .* flow(:,:,ispecies);
            
            temp = densityPerturbation(:,:,ispecies) ./ JHat;
            FSADensityPerturbation(ispecies,:) = (temp'*thetaWeights)' ./ VPrimeB0OverR0;
            
            temp = BHat .* flow(:,:,ispecies) ./ JHat;
            FSABFlow(ispecies,:) = (temp'*thetaWeights)' ./ VPrimeB0OverR0;
            
            temp = pressurePerturbation(:,:,ispecies) ./ JHat;
            FSAPressurePerturbation(ispecies,:) = (temp'*thetaWeights)' ./ VPrimeB0OverR0;
            
            particleFlux(ispecies,:) = (particleFluxBeforeThetaIntegral(:,:,ispecies)'*thetaWeights)';
            momentumFlux(ispecies,:) = (momentumFluxBeforeThetaIntegral(:,:,ispecies)'*thetaWeights)';
            heatFlux(ispecies,:) = (heatFluxBeforeThetaIntegral(:,:,ispecies)'*thetaWeights)';
            
            indices = Npsi*localMatrixSize + ((1:Npsi)-1)*numSpecies*2 + (ispecies-1)*2 + 1;
            particleSource(ispecies,:) = soln(indices);
            
            indices = Npsi*localMatrixSize + ((1:Npsi)-1)*numSpecies*2 + (ispecies-1)*2 + 2;
            heatSource(ispecies,:) = soln(indices);
            
            Machs(:,:,ispecies) = Delta*flow(:,:,ispecies).*sqrt(masses(ispecies)./(ones(Ntheta,1)*THats(ispecies,:)));
        end
   
        flowOutboard = squeeze(flow(1,:,:))';
        if mod(Ntheta,2)==0
            flowInboard = squeeze(flow(Ntheta/2,:,:))';
        else
            index = (Ntheta+1)/2;
            flowInboard = 0.5*squeeze(flow(index,:,:) + flow(index+1,:,:))';
        end
        
        if numSpecies==1
            flowOutboard = flowOutboard';
            flowInboard = flowInboard';
        end
        
        if numSpecies==1
            potentialPerturbation = Delta/(4*omega)*(ones(Ntheta,1)*THats(1,:)) .* densityPerturbation;
        end
        
        temperaturePerturbation = pressurePerturbation - densityPerturbation;
        
        %{
        %ySourceProfile = soln((globalMatrixSize-Npsi+1):globalMatrixSize);
        yParticleSourceProfile = soln(Npsi*localMatrixSize + (1:Npsi));
        %yMomentumSourceProfile = soln(Npsi*localMatrixSize + (1:Npsi) + Npsi);
        yEnergySourceProfile = soln(Npsi*localMatrixSize + (1:Npsi) + Npsi);

        kParOutboard = kpar(1,:);
        if mod(Ntheta,2)==0
            kParInboard = kpar(Ntheta/2+1, :);
        else
            index = (Ntheta+1)/2;
            kParInboard = (kpar(index,:) + kpar(index+1,:))/2;
        end
        
        FSAKPar = ((kpar./JHat)'*thetaWeights)' ./ VPrimeB0OverR0;
        LHSOfKEquation = FSAKPar - omega*Delta/(psiAHat^2)*IHat.*IHat.*dPhiHatdpsi./FSAB2 .* (ddpsi*FSAKPar')';
        LHSOfKEquationAtPsiMid = interp1(psi, LHSOfKEquation, psiMid,'cubic');
        
        nOutboard = densityPerturbation(1,:);
        particleFlux = thetaWeights(:)' * particleFluxBeforeThetaAverage;
        q = thetaWeights(:)' * qBeforeThetaAverage;
        heatFluxRelativeToLocalPlateauRegime = -4*Miller_q*sqrtpi*psiAHat*FSAB2.*q./(3*epsilon*epsilon*VPrimeB0OverR0.*nHat.*(THat.^1.5).*(IHat.^2).*dTHatdpsi);
        kq = q;
        %}
        
        if drawOutputFigures
            for ispecies = 1:numSpecies
                figure(ispecies+6+figureOffset)
                clf
                numContours = 15;
                numRows = 3;
                numCols = 5;
                
                plotNum = 1;
                
                %{
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            contourf(psi, theta, kpar, numContours,'EdgeColor','none')
            xlabel('\psi_N')
            ylabel('\theta')
            title('k_{||}')
            colorbar
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            plot(psi, kParOutboard,'.-r')
            hold on
            plot(psi, kParInboard, '.-b')
            plot(psi, FSAKPar, '.-','Color',[0,0.7,0])
            xlabel('\psi_N')
            ylabel('k_{||}')
            legend('Outboard','Inboard','FSA')
            
            subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
            plot(psi, LHSOfKEquation,'.-')
            xlabel('\psi_N')
            title('LHS of k equation')
                %}
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, densityPerturbation(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('\Delta density')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, flow(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('flow')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, Machs(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('Mach number')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, pressurePerturbation(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('\Delta pressure')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, temperaturePerturbation(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('\Delta temperature')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, particleFluxBeforeThetaIntegral(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('particle flux before \theta average')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, momentumFluxBeforeThetaIntegral(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('momentumFluxBeforeThetaIntegral')
                colorbar
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                contourf(psi, theta, heatFluxBeforeThetaIntegral(:,:,ispecies), numContours,'EdgeColor','none')
                xlabel('\psi_N')
                ylabel('\theta')
                title('heatFluxBeforeThetaIntegral')
                colorbar
                
                if numSpecies==1
                    subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                    contourf(psi, theta, potentialPerturbation(:,:,ispecies), numContours,'EdgeColor','none')
                    xlabel('\psi_N')
                    ylabel('\theta')
                    title('potentialPerturbation')
                    colorbar
                else
                    subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                    plot(psi,FSADensityPerturbation(ispecies,:), '.-')
                    xlabel('\psi_N')
                    title('FSADensityPerturbation')
                end
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                plot(psi,FSABFlow(ispecies,:), '.-')
                xlabel('\psi_N')
                title('FSABFlow')
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                plot(psi,FSAPressurePerturbation(ispecies,:), '.-')
                xlabel('\psi_N')
                title('FSAPressurePerturbation')
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                plot(psi,particleFlux(ispecies,:), '.-')
                xlabel('\psi_N')
                title('<particle flux>')
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                plot(psi,momentumFlux(ispecies,:), '.-')
                xlabel('\psi_N')
                title('<momentum flux>')
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                plot(psi,heatFlux(ispecies,:), '.-')
                xlabel('\psi_N')
                title('heatFlux')
                
                %subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                %plot(psi,heatFluxRelativeToLocalPlateauRegime, '.-')
                %xlabel('\psi_N')
                %title('q / q_{local plateau regime}')
                
                subplot(numRows, numCols, plotNum); plotNum = plotNum + 1;
                plot(psi,particleSource(ispecies,:),'.-')
                hold on
                plot(psi,heatSource(ispecies,:),'.-r')
                legend('particles','energy')
                xlabel('\psi_N')
                ylabel('Source y')
                
            end

            figure(80)
            numPlots = 3;
            numContours = 15;
            ispecies=1;
        
            
            leftMargin = 0.55;
            topMargin = 0.09;
            bottomMargin = 0.08;
            interPlotMargin = 0.06;
            
            plotWidth = 0.42;
            plotHeight = (1-topMargin-bottomMargin+interPlotMargin)/numPlots-interPlotMargin;
            plotNum = 0;
            
            plotNum = plotNum+1;
            subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
            contourf(psi, theta, densityPerturbation(:,:,ispecies), numContours,'EdgeColor','none')
            ylabel('\theta')
            set(gca,'XTickLabel',[])
            title('\Delta density')
            colorbar

            plotNum = plotNum+1;
            subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
            contourf(psi, theta, flow(:,:,ispecies)*Delta*vBar_kmPerSec, numContours,'EdgeColor','none')
            ylabel('\theta')
            set(gca,'XTickLabel',[])
            title('Parallel flow (km/s)')
            colorbar

            plotNum = plotNum+1;
            subplot('Position',[leftMargin, 1-topMargin+interPlotMargin - plotNum*(plotHeight+interPlotMargin), plotWidth, plotHeight])
            contourf(psi, theta, pressurePerturbation(:,:,ispecies)-densityPerturbation(:,:,ispecies), numContours,'EdgeColor','none')
            ylabel('\theta')
            xlabel('\psi_N')
            %set(gcf,'XTickLabel',[])
            title('\Delta temperature')
            colorbar


        end
        
        % *********************************************************
        % *********************************************************
        %
        % Below are utilities related to creating sparse matrices:
        %
        % *********************************************************
        % *********************************************************
        
        function resetSparseCreator()
            sparseCreatorIndex=1;
            sparseCreator_i=zeros(estimated_nnz,1);
            sparseCreator_j=zeros(estimated_nnz,1);
            sparseCreator_s=zeros(estimated_nnz,1);
        end
        
        function addToSparse(i,j,s)
            n=numel(i);
            if n ~= numel(j)
                error('Error A');
            end
            if n ~= numel(s)
                error('Error B');
            end
            if any(i<1)
                error('Error Q: i<1');
            end
            if any(j<1)
                error('Error Q: j<1');
            end
            sparseCreator_i(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = i;
            sparseCreator_j(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = j;
            sparseCreator_s(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = s;
            sparseCreatorIndex = sparseCreatorIndex+n;
            if sparseCreatorIndex > estimated_nnz
                fprintf('Warning! estimated_nnz is too small.\n')
            end
        end
        
        function addSparseBlock(rowIndices, colIndices, block)
            s=size(block);
            if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
                s
                size(rowIndices)
                size(colIndices)
                error('Error in addSparseBlock!')
            end
            [rows, cols, values] = find(block);
            addToSparse(rowIndices(rows),colIndices(cols),values)
        end
        
        function sparseMatrix = createSparse(sizeForThisMatrix1, sizeForThisMatrix2)
            sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), sizeForThisMatrix1, sizeForThisMatrix2);
            resetSparseCreator()
        end
    end

% *********************************************************
% *********************************************************
%
% Below are functions related to the magnetic geometry
%
% *********************************************************
% *********************************************************

    function zz = RHat(theta)
        zz = 1 + (1/Miller_A)*cos(theta + Miller_x*sin(theta));
    end

    function zz = ZHat(theta)
        zz = (Miller_kappa/Miller_A)*sin(theta);
    end

    function zz = QIntegrand(theta)
        zz = ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) .* (1+Miller_x*cos(theta)) .* sin(theta) ...
            + cos(theta) .* (Miller_dRdr + cos(theta + Miller_x *sin(theta)) - Miller_s_delta*sin(theta + Miller_x*sin(theta)) .* sin(theta))) ./ RHat(theta);
    end

    function zz = BPoloidal(theta)
        zz = Miller_Q/(Miller_kappa*Miller_q)*sqrt((sin(theta+Miller_x*sin(theta)) .* (1+Miller_x*cos(theta))).^2 + (Miller_kappa*cos(theta)).^2) ...
            ./ (RHat(theta) .* ( cos(Miller_x*sin(theta)) + Miller_dRdr*cos(theta) + (Miller_s_kappa-Miller_s_delta*cos(theta) + (1+Miller_s_kappa)*Miller_x*cos(theta)) .* sin(theta) .* sin(theta + Miller_x*sin(theta))));
    end

    function zz = BDotGradTheta(theta)
        zz = - Miller_A*Miller_Q./(Miller_kappa*Miller_q*RHat(theta) .* ...
            ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) .* (1+Miller_x*cos(theta)) .* sin(theta) ...
            + cos(theta) .* (Miller_dRdr + cos(theta + Miller_x *sin(theta)) - Miller_s_delta*sin(theta + Miller_x*sin(theta)) .* sin(theta))));
    end


    function zz = b(theta)
        switch geometry
            case 0
                zz = 1./(1 + cos(theta)/Miller_A);
            case 1
                zz = sqrt(BPoloidal(theta).^2 + 1./((RHat(theta)).^2));
            case 2
                zz = 1 + cos(theta)/Miller_A;
            case 3
                zz = 1 - cos(theta)/Miller_A + epsilonUDB*sin(theta);
            otherwise
                error('Invalid geometry')
        end
    end

    function zz = dbdtheta(theta)
        switch geometry
            case 0
                zz = sin(theta) ./(Miller_A * (1 + cos(theta)/Miller_A).^2);
            case 1
                error('Program should not get here!')
            case 2
                zz = - sin(theta)/Miller_A;
            case 3
                zz = sin(theta)/Miller_A + epsilonUDB * cos(theta);
            otherwise
                error('Invalid geometry')
        end
    end

    function zz = oneOverqRbDotGradTheta(theta)
        switch geometry
            case 0
                zz = ones(size(theta));
            case 1
                zz = b(theta) ./ (Miller_q*BDotGradTheta(theta));
            case 2
                zz = 1./b(theta);
            case 3
                zz = 1 + epsilonUDJacobian*sin(theta);
            otherwise
                error('Invalid geometry')
        end
    end

    function pp=compute_Phi(r, ss, r0, rampAmplitude)
        switch radialProfiles
            case {2,7}
                e=exp(1);
                pp = rampAmplitude * ( ...
                    -e^(r0*ss)*log(-(sqrt(e^(2*r0*ss) - 4)*e^(r0*ss) - 2*e.^(-r*ss + r0*ss) - e^(2*r0*ss))./(sqrt(e^(2*r0*ss) - 4)*e^(r0*ss) + 2*e.^(-r*ss + r0*ss) + e^(2*r0*ss)))/(sqrt(e^(2*r0*ss) - 4)*ss));
            case 3
                pp = rampAmplitude * r;
            case {4,6}
                pp = rampAmplitude * erf(r*ss) * sqrt(pi)/(2*ss);
        end
    end
    function pp=compute_dPhidr(r, ss, r0, rampAmplitude)
        switch radialProfiles
            case {2,7}
                pp = rampAmplitude * ( ...
                    1./(1 + exp(ss*(r-r0)) + exp(-ss*(r+r0))));
            case 3
                pp = rampAmplitude * ones(size(r));
            case {4,6}
                pp = rampAmplitude * exp(-(r*ss).^2);
        end
    end

end