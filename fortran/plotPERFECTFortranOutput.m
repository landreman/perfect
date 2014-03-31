function plotPERFECTFortranOutput()

filename='PERFECTOutput.h5';

colors = [1,0,0;
    0.8,0.6,0;
    0,0.7,0;
    0,0.8,0.9;
    0,0,1;
    1,0,1;
    0.6,0.6,0.6;
    0,0,0];


figureOffset = 10;

info = h5info(filename);
fprintf('Fields saved in this HDF5 file:\n')
for i=1:numel(info.Groups(1).Datasets)
    fprintf('  %s\n',info.Groups(1).Datasets(i).Name)
end

numRuns = numel(info.Groups);

programMode = h5read(filename,'/programMode');

switch programMode
    case 1
        legendText={};
    case {2,3,0}
        
        legendText = cell(numRuns,1);
        location  = getLocationString(1);
        
        baseNpsiPerDiameter = h5read(filename,[location,'NpsiPerDiameter']);
        basePsiDiameter = h5read(filename,[location,'psiDiameter']);
        baseWidthExtender = h5read(filename,[location,'widthExtender']);
        baseNtheta = h5read(filename,[location,'Ntheta']);
        baseNxi = h5read(filename,[location,'Nxi']);
        baseNL = h5read(filename,[location,'NL']);
        baseNx = h5read(filename,[location,'Nx']);
        baseNxPotentialsPerVth = h5read(filename,[location,'NxPotentialsPerVth']);
        baseXMax = h5read(filename,[location,'xMax']);
        baseSolverTolerance = h5read(filename,[location,'solverTolerance']);
        makeLocalApproximation = h5read(filename,[location,'makeLocalApproximation']);
        j=1;
        if programMode==3
            text = 'local result';
            location  = getLocationString(1);
            didItConverge = h5read(filename,[location,'didItConverge']);
            if didItConverge<0
                text = [text,' (not converged)'];
            end
            legendText{1} = text;
            j=2;
        end
        text = 'base case';
        location  = getLocationString(j);
        didItConverge = h5read(filename,[location,'didItConverge']);
        if didItConverge<0
            text = [text,' (not converged)'];
        end
        legendText{j} = text;
        
        for i=(j+1):numRuns
            text='';
            location  = getLocationString(i);
            NpsiPerDiameter = h5read(filename,[location,'NpsiPerDiameter']);
            psiDiameter = h5read(filename,[location,'psiDiameter']);
            widthExtender = h5read(filename,[location,'widthExtender']);
            Ntheta = h5read(filename,[location,'Ntheta']);
            Nxi = h5read(filename,[location,'Nxi']);
            NL = h5read(filename,[location,'NL']);
            Nx = h5read(filename,[location,'Nx']);
            NxPotentialsPerVth = h5read(filename,[location,'NxPotentialsPerVth']);
            xMax = h5read(filename,[location,'xMax']);
            solverTolerance = h5read(filename,[location,'solverTolerance']);
            didItConverge = h5read(filename,[location,'didItConverge']);
            makeLocalApproximation = h5read(filename,[location,'makeLocalApproximation']);
            if NpsiPerDiameter ~= baseNpsiPerDiameter
                text = [text,'NpsiPerDiameter=',num2str(NpsiPerDiameter),'  '];
            end
            if psiDiameter ~= basePsiDiameter
                text = [text,'psiDiameter=',num2str(psiDiameter),'  '];
            end
            if widthExtender ~= baseWidthExtender
                text = [text,'widthExtender=',num2str(widthExtender),'  '];
            end
            if Ntheta ~= baseNtheta
                text = [text,'Ntheta=',num2str(Ntheta),'  '];
            end
            if Nxi ~= baseNxi
                text = [text,'Nxi=',num2str(Nxi),'  '];
            end
            if NL ~= baseNL
                text = [text,'NL=',num2str(NL),'  '];
            end
            if Nx ~= baseNx
                text = [text,'Nx=',num2str(Nx),'  '];
            end
            if NxPotentialsPerVth ~= baseNxPotentialsPerVth
                text = [text,'NxPotentialsPerVth=',num2str(NxPotentialsPerVth),'  '];
            end
            if xMax ~= baseXMax
                text = [text,'xMax=',num2str(xMax),'  '];
            end
            if solverTolerance ~= baseSolverTolerance
                text = [text,'solverTolerance=',num2str(solverTolerance),'  '];
            end
            if didItConverge<0
                text = [text,'(not converged)'];
            end
            legendText{i} = text;
        end
    case 4
        legendText = {'local approx','particle source, not upwinding','momentum source, not upwinding','particle source, upwinding','momentum source, upwinding'};
    case 5
        % U scan
        legendText = cell(numRuns,1);
        for runNum=1:numRuns
            desiredU = h5read(filename,[getLocationString(runNum),'/desiredU']);
            legendText{runNum} = ['U=',num2str(desiredU)];
        end
    otherwise
        error('Unrecognized programMode')
end

%{
if numRuns==1
    programMode = 1;
    % Single run
else
    programMode = 2;
    % Convergence scan
end
%}

for i=1:numRuns
    location  = getLocationString(i);
    didItConverge = h5read(filename,[location,'didItConverge']);
    if didItConverge<0
        beep
        fprintf('Warning: run %d did not converge.\n',i)
    end
end

% Find the run with highest Npsi and display details of that run.
runWithMaxNpsi = 1;
location  = getLocationString(1);
maxNpsi = h5read(filename,[location,'Npsi']);
for i=2:numRuns
    location  = getLocationString(i);
    Npsi = h5read(filename,[location,'Npsi']);
    if Npsi >= maxNpsi
        runWithMaxNpsi = i;
        maxNpsi = Npsi;
    end
end

location  = getLocationString(runWithMaxNpsi);
masses = h5read(filename,[location,'masses']);
numSpecies = numel(masses);
%{
if numSpecies>1
    error('This m-file is not designed for runs with >1 species')
end
%}

plotInputs(runWithMaxNpsi)
%plotInputs(2)

% ************************************************************************
% Quantities that depend on both theta and psi, only shown for the run with
% highest Npsi.
% ************************************************************************


psi = h5read(filename,[location,'psi']);
theta = h5read(filename,[location,'theta']);
Delta = h5read(filename,[location,'Delta']);
particleSourceProfile = h5read(filename,[location,'particleSourceProfile']);
heatSourceProfile = h5read(filename,[location,'heatSourceProfile']);
kPar = h5read(filename,[location,'kPar']);
flow = h5read(filename,[location,'flow']);
THats = h5read(filename,[location,'THat']);
densityPerturbation = h5read(filename,[location,'densityPerturbation']);
pressurePerturbation = h5read(filename,[location,'pressurePerturbation']);
kParOutboard = h5read(filename,[location,'kParOutboard']);
kParInboard = h5read(filename,[location,'kParInboard']);
FSAKPar = h5read(filename,[location,'FSAKPar']);
particleFlux = h5read(filename,[location,'particleFlux']);
particleFluxBeforeThetaIntegral = h5read(filename,[location,'particleFluxBeforeThetaIntegral']);
momentumFlux = h5read(filename,[location,'momentumFlux']);
momentumFluxBeforeThetaIntegral = h5read(filename,[location,'momentumFluxBeforeThetaIntegral']);
heatFlux = h5read(filename,[location,'heatFlux']);
heatFluxBeforeThetaIntegral = h5read(filename,[location,'heatFluxBeforeThetaIntegral']);
%{
kThetaWith3PointStencil = h5read(filename,[location,'kThetaWith3PointStencil']);
kThetaWith5PointStencil = h5read(filename,[location,'kThetaWith5PointStencil']);
kThetaOutboardWith3PointStencil = h5read(filename,[location,'kThetaOutboardWith3PointStencil']);
kThetaOutboardWith5PointStencil = h5read(filename,[location,'kThetaOutboardWith5PointStencil']);
kThetaInboardWith3PointStencil = h5read(filename,[location,'kThetaInboardWith3PointStencil']);
kThetaInboardWith5PointStencil = h5read(filename,[location,'kThetaInboardWith5PointStencil']);
PhiTermInKTheta = h5read(filename,[location,'PhiTermInKTheta']);
pPerpTermInKThetaBeforePsiDerivative = h5read(filename,[location,'pPerpTermInKThetaBeforePsiDerivative']);
pPerpTermInKThetaWith3PointStencil = h5read(filename,[location,'pPerpTermInKThetaWith3PointStencil']);
pPerpTermInKThetaWith5PointStencil = h5read(filename,[location,'pPerpTermInKThetaWith5PointStencil']);
%}

for ispecies = 1:numSpecies
    figure(numSpecies+figureOffset+ispecies)
    
    clf
    numRows = 3;
    numCols = 4;
    plotNum=1;
    
    Ntheta = numel(theta);
    bigTHat = ones(Ntheta,1)*THats(ispecies,:);
    Mach = Delta * squeeze(flow(ispecies,:,:)) ./ sqrt(bigTHat);
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,squeeze(kPar(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('k_{||}')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,squeeze(flow(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('flow')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,Mach,numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('Mach number')
    
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    %contours = linspace(-1e-4,1e-4,20);
    %contourf(psi,theta,densityPerturbation,contours,'EdgeColor','none')
    contourf(psi,theta,squeeze(densityPerturbation(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('density perturbation')
    %set(gca,'CLim',[-2e-4,2e-4])
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,squeeze(pressurePerturbation(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('pressure perturbation')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    plot(psi,kParOutboard(ispecies,:),'.-')
    hold on
    plot(psi,kParInboard(ispecies,:),'.-r')
    plot(psi,FSAKPar(ispecies,:),'.-','Color',[0.8,0.7,0])
    xlabel('\psi')
    ylabel('k_{||}')
    legend('k_{||} outboard','k_{||} inboard','<k_{||}>')
    axis tight
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,squeeze(particleFluxBeforeThetaIntegral(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('particle flux before \theta integral')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,squeeze(momentumFluxBeforeThetaIntegral(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('momentum flux before \theta integral')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    contourf(psi,theta,squeeze(heatFluxBeforeThetaIntegral(ispecies,:,:)),numContours,'EdgeColor','none')
    colorbar
    xlabel('\psi')
    ylabel('\theta')
    title('heat flux before \theta integral')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    plot(psi,particleSourceProfile(ispecies,:),'.-')
    hold on
    plot(psi,heatSourceProfile(ispecies,:),'.-r')
    xlabel('\psi')
    ylabel('source magnitude')
    legend('particles','heat')
    axis tight
    
    titleString = sprintf('Outputs for species %d: mass=%g, charge = %g',ispecies,masses(ispecies),charges(ispecies));
    
    annotation(gcf,'textbox',...
        [0 0.95 1 0.05],'LineStyle','none','HorizontalAlignment','center',...
        'String',titleString,...
        'FitBoxToText','off');
    
    
end

% ************************************************************************
% Profiles, comparing runs if there is more than 1 run.
% ************************************************************************

for ispecies = 1:numSpecies
    figure(2*numSpecies+figureOffset+ispecies)
    clf
    numRows = 3;
    numCols = 3;
    
    for runNum = 1:numRuns
        location = getLocationString(runNum);
        psi = h5read(filename,[location,'psi']);
        theta = h5read(filename,[location,'theta']);
        particleSourceProfile = h5read(filename,[location,'particleSourceProfile']);
        heatSourceProfile = h5read(filename,[location,'heatSourceProfile']);
        kParOutboard = h5read(filename,[location,'kParOutboard']);
        kParInboard = h5read(filename,[location,'kParInboard']);
        FSAKPar = h5read(filename,[location,'FSAKPar']);
        particleFlux = h5read(filename,[location,'particleFlux']);
        momentumFlux = h5read(filename,[location,'momentumFlux']);
        heatFlux = h5read(filename,[location,'heatFlux']);
        densityPerturbation = h5read(filename,[location,'densityPerturbation']);
        Ntheta = numel(theta);
        densityPerturbationAtPsiMid = zeros(Ntheta,1);
        
        psiMid = 0.9;
        for itheta = 1:Ntheta
            densityPerturbationAtPsiMid(itheta) = interp1(psi,squeeze(densityPerturbation(ispecies,itheta,:)), psiMid);
        end
        
        %{
    kThetaOutboardWith3PointStencil = h5read(filename,[location,'kThetaOutboardWith3PointStencil']);
    kThetaOutboardWith5PointStencil = h5read(filename,[location,'kThetaOutboardWith5PointStencil']);
    kThetaInboardWith3PointStencil = h5read(filename,[location,'kThetaInboardWith3PointStencil']);
    kThetaInboardWith5PointStencil = h5read(filename,[location,'kThetaInboardWith5PointStencil']);
    LHSOfKParEquation = h5read(filename,[location,'LHSOfKParEquation']);
        %}
        %{
        psiMid = (min(psi)+max(psi))/2;
        U = h5read(filename,[location,'U']);
        if max(abs(psi))>0
            UMid = interp1(psi(:),U(:),psiMid,'cubic');
        end
        alpha = -4*UMid.*UMid;
        %    LHSOfKParEquation = LHSOfKParEquation + alpha * FSAKPar;
        %}
        
        colorIndex = 1 + mod(runNum-1, size(colors,1));
        
        plotNum = 1;
        
        myLinespec = '.-';
        myColor = colors(colorIndex,:);
        if runNum > size(colors,1)
            myLinespec = '.:';
        end
        if runNum > 2*size(colors,1)
            myLinespec = '.--';
        end
        if programMode==3 && runNum==1
            myLinespec = ':k';
            myColor = [0,0,0];
        end
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,kParOutboard(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,kParInboard(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,FSAKPar(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,particleSourceProfile(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,heatSourceProfile(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,particleFlux(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,momentumFlux(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,heatFlux(ispecies,:),myLinespec,'Color',myColor)
        hold on
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(theta,densityPerturbationAtPsiMid,myLinespec,'Color',myColor)
        hold on
        
        fftd = fft(densityPerturbationAtPsiMid);
        sinAmplitude = -2*imag(fftd(2))/Ntheta;
        cosAmplitude = 2*real(fftd(2))/Ntheta;
        thetaFine = linspace(0,2*pi);
        plot(thetaFine,sinAmplitude*sin(thetaFine)+cosAmplitude*cos(thetaFine),'-','Color',myColor)
        
    end
    
    plotNum = 1;
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('k_{||} at outboard side')
    axis tight
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('k_{||} at inboard side')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('<k_{||}>')
    
    legend(legendText)
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('particle source magnitude')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('heat source magnitude')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('particle flux')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('momentum flux')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\psi')
    ylabel('heat flux')
    axis tight
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
    xlabel('\theta')
    ylabel('density perturbation at  \psi_{mid}')
    axis tight
   
    titleString = sprintf('Outputs for species %d: mass=%g, charge = %g',ispecies,masses(ispecies),charges(ispecies));
    
    annotation(gcf,'textbox',...
        [0 0.95 1 0.05],'LineStyle','none','HorizontalAlignment','center',...
        'String',titleString,...
        'FitBoxToText','off');

    
end

% --------------------------------------------------------

    function l = getLocationString(runNum)
        l = sprintf('/run%3d/',runNum);
    end

    function plotInputs(runNum)
        location  = getLocationString(runNum);
        
        masses = h5read(filename,[location,'masses']);
        charges = h5read(filename,[location,'charges']);
        psi = h5read(filename,[location,'psi']);
        r = h5read(filename,[location,'r']);
        exponent = h5read(filename,[location,'exponent']);
        theta = h5read(filename,[location,'theta']);
        THat = h5read(filename,[location,'THat']);
        dTHatdpsi = h5read(filename,[location,'d(THat)d(psi)']);
        IHat = h5read(filename,[location,'IHat']);
        dIHatdpsi = h5read(filename,[location,'d(IHat)d(psi)']);
        nHat = h5read(filename,[location,'nHat']);
        etaHat = h5read(filename,[location,'etaHat']);
        PhiHat = h5read(filename,[location,'PhiHat']);
        dPhiHatdpsi = h5read(filename,[location,'d(PhiHat)d(psi)']);
        BHat = h5read(filename,[location,'BHat']);
        dBHatdtheta = h5read(filename,[location,'d(BHat)d(theta)']);
        dBHatdpsi = h5read(filename,[location,'d(BHat)d(psi)']);
        JHat = h5read(filename,[location,'JHat']);
        nuPrimeProfile = h5read(filename,[location,'nuPrimeProfile']);
        nuStarProfile = h5read(filename,[location,'nuStarProfile']);
        deltaN = h5read(filename,[location,'deltaN']);
        deltaT = h5read(filename,[location,'deltaT']);
        deltaEta = h5read(filename,[location,'deltaEta']);
        %{
        deltaNAtBMax = h5read(filename,[location,'deltaNAtBMax']);
        deltaNAtBMin = h5read(filename,[location,'deltaNAtBMin']);
        deltaTAtBMax = h5read(filename,[location,'deltaTAtBMax']);
        deltaTAtBMin = h5read(filename,[location,'deltaTAtBMin']);
        deltaEtaAtBMax = h5read(filename,[location,'deltaEtaAtBMax']);
        deltaEtaAtBMin = h5read(filename,[location,'deltaEtaAtBMin']);
        %}
        U = h5read(filename,[location,'U']);
        
        for ispecies = 1:numSpecies
            figure(figureOffset+ispecies)
            clf
            numRows = 4;
            numCols = 4;
            numContours = 12;
            plotNum = 1;
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            contourf(psi,theta,BHat,numContours,'EdgeColor','none')
            colorbar
            xlabel('\psi')
            ylabel('\theta')
            title('BHat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            contourf(psi,theta,dBHatdtheta,numContours,'EdgeColor','none')
            colorbar
            xlabel('\psi')
            ylabel('\theta')
            title('dBHat/d\theta')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            contourf(psi,theta,dBHatdpsi,numContours,'EdgeColor','none')
            colorbar
            xlabel('\psi')
            ylabel('\theta')
            title('dBHat/d\psi')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            contourf(psi,theta,JHat,numContours,'EdgeColor','none')
            colorbar
            xlabel('\psi')
            ylabel('\theta')
            title('JHat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,THat(ispecies,:),'.-')
            xlabel('\psi')
            ylabel('THat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,dTHatdpsi(ispecies,:),'.-')
            xlabel('\psi')
            ylabel('d(THat)/d\psi')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,IHat,'.-')
            xlabel('\psi')
            ylabel('IHat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,(nHat(ispecies,:).^exponent).*dTHatdpsi(ispecies,:),'.-')
            xlabel('\psi')
            ylabel(['n^',num2str(exponent),'*(dT/d\psi)'])
            
            %{
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,dIHatdpsi,'.-')
        xlabel('\psi')
        ylabel('d(IHat)/d\psi')
            %}
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,nHat(ispecies,:),'.-')
            xlabel('\psi')
            ylabel('nHat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,etaHat(ispecies,:),'.-')
            xlabel('\psi')
            ylabel('\eta Hat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,PhiHat,'.-')
            xlabel('\psi')
            ylabel('\Phi Hat')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,dPhiHatdpsi,'.-')
            xlabel('\psi')
            ylabel('d(\Phi Hat)/d\psi')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,U(ispecies,:),'.-')
            xlabel('\psi')
            ylabel('U')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(psi,nuStarProfile(ispecies,:),'.-')
            hold on
            plot(psi,nuPrimeProfile(ispecies,:),'.-r')
            legend('\nu_*','\nu''')
            xlabel('\psi')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            semilogy(psi,nuStarProfile(ispecies,:),'.-')
            hold on
            semilogy(psi,nuPrimeProfile(ispecies,:),'.-r')
            legend('\nu_*','\nu''')
            xlabel('\psi')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
            plot(r(ispecies,:),deltaN(ispecies,:),'.-')
            hold on
            plot(r(ispecies,:),deltaT(ispecies,:),'.-r')
            plot(r(ispecies,:),deltaEta(ispecies,:),'.-','Color',[0,0.7,0])
            xlabel('r')
            legend('\delta_n','\delta_T','\delta_\eta')
            
            %{
        subplot(numRows,numCols,plotNum); plotNum = plotNum+1;
        plot(psi,deltaNAtBMin,'.:')
        hold on
        plot(psi,deltaTAtBMin,'.:r')
        plot(psi,deltaEtaAtBMin,'.:','Color',[0,0.7,0])
        plot(psi,deltaNAtBMax,'.-')
        plot(psi,deltaTAtBMax,'.-r')
        plot(psi,deltaEtaAtBMax,'.-','Color',[0,0.7,0])
        xlabel('\psi')
        legend('\delta_n (Bmin)','\delta_T (BMin)','\delta_\eta (BMin)','\delta_n (Bmax)','\delta_T (BMax)','\delta_\eta (BMax)')
            %}
            
            titleString = sprintf('Inputs for species %d: mass=%g, charge = %g',ispecies,masses(ispecies),charges(ispecies));
            
            annotation(gcf,'textbox',...
                [0 0.95 1 0.05],'LineStyle','none','HorizontalAlignment','center',...
                'String',titleString,...
                'FitBoxToText','off');

        end
        
    end

end