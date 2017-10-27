function [thetaSurfs, BPSurfs, BDotGradThetaSurfs, ISurfs, qSurfs, RSurfs, ZSurfs, as, R0, Z0, B0, psi0] = getGeometryFromEFITForSeveralFluxSurfaces(filename, desiredPsiNs, topCropZ, bottomCropZ, innerCropR, outerCropR, plotStuff, saveSeparatrix, extrapolateBeyondPsiN, extrapolatePsiNInterval, rmaxExpansionFactor)

% 'BPSurfs' and 'B0' should have units of Tesla.
% Both 'as' and 'R0' should have units of meters.
% 'BDotGradThetaSurfs' should have units of Tesla / meters.

% Good values for topCropZ and bottomCropZ for Alcator C-Mod are 0.4 and -0.4.

if (any(desiredPsiNs>1) && extrapolateBeyondPsiN<0) || any(desiredPsiNs<0)
    error('All desired psi_N values must be between 0 and 1.')
end

efit = read_eqdsk(filename);

psiN = (efit.psi_grid-efit.psiaxis)/(efit.psiedge-efit.psiaxis);
psiN2D = (efit.psi-efit.psiaxis)/(efit.psiedge-efit.psiaxis);
psi0 = (efit.psiedge-efit.psiaxis); % Needed to set psiAHat for PERFECT

R_grid = efit.R_grid;
Z_grid = efit.Z_grid;

if extrapolateBeyondPsiN>0
  % extrapolate psiN along lines of constant poloidal angle theta
  % in order to create closed flux surfaces for a buffer zone

  % create grids
  [Rg,Zg] = meshgrid(R_grid-efit.Raxis,Z_grid-efit.Zaxis);
  [theta_RZ,r_RZ] = cart2pol(Rg,Zg);
  rmax = max(r_RZ(:))*rmaxExpansionFactor;
  % make theta grid be uniformly spaced between -pi and pi, including the points at -pi and pi and some extra points to make cubic interpolation OK around +-pi
  %[theta,r] = meshgrid(linspace(-pi*(1-1/nfine),pi*(1-1/nfine),nfine),linspace(0,rmax,nfine));
  %[theta,r] = meshgrid(linspace(-pi-2*2*pi/(nfine+1),pi+2*2*pi/(nfine+1),nfine+5),linspace(0,rmax,nfine));
  [theta,r] = meshgrid(linspace(-pi-2*2*pi/(nfinetheta+1),pi+2*2*pi/(nfinetheta+1),nfinetheta+5),linspace(0,rmax,nfiner));
  [R_rtheta,Z_rtheta] = pol2cart(theta,r);

  % interpolate psiN onto r-theta grid
  psiN2D_rtheta = interp2(Rg,Zg,psiN2D,R_rtheta,Z_rtheta,'cubic');

  % extrapolate along constant theta lines
  for itheta=1:nfinetheta+5
    iextrap = find(psiN2D_rtheta(:,itheta)>extrapolateBeyondPsiN, 1)-1;
    %igrad = find(psiN2D_rtheta(:,itheta)>(extrapolateBeyondPsiN-extrapolatePsiNInterval), 1)-1;
    %if igrad<1
    %  igrad=1;
    %end
    %psiN2D_rtheta((iextrap+1):end,itheta) = psiN2D_rtheta(iextrap,itheta) + (psiN2D_rtheta(iextrap,itheta)-psiN2D_rtheta(igrad,itheta))*(1:(nfiner-iextrap))/(iextrap-igrad);
    thisinds = 1 : find(psiN2D_rtheta(:,itheta)>1,1)-1;
    thispsiN = psiN2D_rtheta(thisinds,itheta);
    thisr = r(thisinds,itheta);
    r1 = interp1(thispsiN,thisr,extrapolateBeyondPsiN,'pchip');
    % extrapolatePsiNInterval should be small for a continuous
    % dpsi/dr
    % TODO: perhaps respect derivative scheme to be used.
    % Might not matter after smoothing?      
    r2 = interp1(thispsiN,thisr,extrapolateBeyondPsiN-extrapolatePsiNInterval,'pchip');
    gradient = extrapolatePsiNInterval/(r1-r2);
    psiN2D_rtheta((iextrap+1):end,itheta) = extrapolateBeyondPsiN + (r(iextrap+1:end,itheta)-r1)*gradient;
  end

  % Make new R-Z grid to make sure values interpolated from r-theta do not fall off the edge of the grid
  rmax = max(r(:));
  R_grid = linspace(efit.Raxis-rmax,efit.Raxis+rmax,nfiner);
  Z_grid = linspace(efit.Zaxis-rmax,efit.Zaxis+rmax,nfiner);
  [Rg,Zg] = meshgrid(R_grid-efit.Raxis,Z_grid-efit.Zaxis);
  % interpolate back to R-Z grid
  [theta_RZ,r_RZ] = cart2pol(Rg,Zg);
  psiN2D = interp2(theta,r,psiN2D_rtheta,theta_RZ,r_RZ,'cubic');

  % recalculate unnormalized psi, used for computing magnetic field components
  efit.psi = (efit.psiedge-efit.psiaxis)*psiN2D+efit.psiaxis;
  
  % set I to constant in extrapolated region
  constantI = false;
  constantDIHatDPsi = true;
  iextrap2 = find(psiN>extrapolateBeyondPsiN,1)-1;
  if constantI
      efit.T(iextrap2:end) = efit.T(iextrap2);
  elseif constantDIHatDPsi
      % TODO: perhaps respect derivative scheme to be used.
      % Might not matter after smoothing?
      DI = efit.T(iextrap2) - efit.T(iextrap2-1);
      Dpsi = psiN(iextrap2) - psiN(iextrap2-1);
      DIDpsi = DI/Dpsi;
      efit.T(iextrap2:end) = efit.T(iextrap2) + DIDpsi * (psiN(iextrap2:end) - psiN(iextrap2));
  else
      disp('Invalid I extrapolation!')
  end
else
  valueForCropping = max(max(psiN2D));
  psiN2D(Z_grid > topCropZ, :) = valueForCropping;
  psiN2D(Z_grid < bottomCropZ, :) = valueForCropping;
  psiN2D(:, R_grid < innerCropR) = valueForCropping;
  psiN2D(:, R_grid > outerCropR) = valueForCropping;
end

qSurfs = interp1(psiN, efit.q, desiredPsiNs, 'pchip');
ISurfs = interp1(psiN, efit.T, desiredPsiNs, 'pchip');

R0 = efit.Raxis;
Z0 = efit.Zaxis;
B0 = efit.B0EXP;

N = numel(desiredPsiNs);

if plotStuff
    fig0 = figure('Visible','off');
    %contour(R_grid, Z_grid, efit.psi, efit.psiaxis+efit.psiedge*linspace(.01,1.1,110))
    %contour(R_grid, Z_grid, (efit.psi-efit.psiaxis)/(efit.psiedge-efit.psiaxis), linspace(.01,1.1,1.1*N))
    if extrapolateBeyondPsiN>0
      contour(R_grid, Z_grid, psiN2D, linspace(.01,max(desiredPsiNs),1.1*N))
    else
      contour(R_grid, Z_grid, psiN2D, linspace(.01,1.1,1.1*N))
    end
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    if extrapolateBeyondPsiN<=0
      % add box showing crop values
      line([innerCropR,innerCropR,outerCropR,outerCropR,innerCropR],[bottomCropZ,topCropZ,topCropZ,bottomCropZ,bottomCropZ])
    else
      % add circle showing the interpolation grid boundary
      plot(efit.Raxis+rmax*cos(theta),efit.Zaxis+rmax*sin(theta))
    end
    axis equal
    xlabel('R (m)')
    ylabel('Z (m)')
    title('\psi_N')
    print(fig0,'EFITgeometry_testFig0','-dpdf')
end
disp('done initial plotStuff')

scheme = 12;
[R, ~, ddR, d2dR2] = differentiationMatricesForUniformGrid(efit.nrbox, efit.rboxleft, efit.rboxleft+efit.rboxlength, scheme);
[Z, ~, ddZ, d2dZ2] = differentiationMatricesForUniformGrid(efit.nzbox, -efit.zboxlength/2, efit.zboxlength/2, scheme);

[R2D, Z2D] = meshgrid(R, Z);
theta = atan2(Z2D - efit.Zaxis, R2D-efit.Raxis);

psiNoNAN = efit.psi;
psiNoNAN(isnan(psiNoNAN)) = 0.;
dpsidZ = ddZ * psiNoNAN;
dpsidR = (ddR * (psiNoNAN'))';

BZ = dpsidR./R2D;
BR = -dpsidZ./R2D;

BPol = sqrt(dpsidZ.^2 + dpsidR.^2)./ R2D;
I2D = interp1(efit.psi_grid, efit.T, psiNoNAN,'pchip',efit.T(end));
BTor = I2D ./ R2D;

B = sqrt(BPol.^2 + BTor.^2);

hypotenuse = (R2D-R0).^2 + (Z2D-Z0).^2;
dthetadZ = (R2D-R0) ./ hypotenuse;
dthetadR = -(Z2D-Z0) ./ hypotenuse;

thetaSurfs = cell(N,1);
BPSurfs = cell(N,1);
RSurfs = cell(N,1);
ZSurfs = cell(N,1);
BDotGradThetaSurfs = cell(N,1);
Rss = cell(N,1);
Zss = cell(N,1);
as = zeros(N,1);

BDotGradTheta = BR .* dthetadR + BZ .* dthetadZ;
%bDotGradTheta = BDotGradTheta ./ B;
    
for i=1:N
    desiredPsiN = desiredPsiNs(i);
    c = contourc(R_grid, Z_grid, psiN2D, [desiredPsiN, desiredPsiN]);
    
    if size(c,2) < 10
        error(['There do not appear to be enough points in the contour for psi_N = ',num2str(desiredPsiN)])
    end
    if size(c,2) ~= c(2,1)+1
        error('There is more than 1 separate curve in the contour.')
    end
    
    % Skip point 2 since it's repeated at the end:
    Rs = c(1,3:end);
    Zs = c(2,3:end);
    Rss{i} = Rs;
    Zss{i} = Zs;
    
    %thetaSurf = interp2(R, Z, theta, Rs, Zs,'spline');
    % Recalculate thetas from the Rs and Zs on the contour instead of interpolating.
    % Gives better results when smoothing/interpolating onto poloidal angle grid later
    thetaSurf = atan2(Zs-efit.Zaxis, Rs-efit.Raxis);

    BPSurf = interp2(R, Z, BPol, Rs, Zs,'cubic');
    BDotGradThetaSurf = interp2(R, Z, BDotGradTheta, Rs, Zs,'cubic');
    
    % Sometimes, a point can fall in the tiny sliver of the R-Z plane where theta decreases
    % from 2pi back to 0 as you move clockwise, instead of the main region
    % where theta increases as you move clockwise. We must remove these
    % points.
    thetaIncreasing = (thetaSurf - circshift(thetaSurf, [0,1])) > 0;
    % Written for theta (mostly) decreasing around the contour, if theta generally increases, then need to take the complement of theteIncreasing
    if sum(thetaIncreasing)>numel(thetaIncreasing)/2.
      thetaIncreasing = ~thetaIncreasing;
    end
    % One such point will always be present, and if there is only 1 point it is not a problem,
    % but if there are more than 1 such point, remove all but the last:
    thetaIncreasing(find(thetaIncreasing,1,'last')) = false;
    thetaSurf(thetaIncreasing) = [];
    BPSurf(thetaIncreasing) = [];
    BDotGradThetaSurf(thetaIncreasing) = [];
    Rs(thetaIncreasing) = [];
    Zs(thetaIncreasing) = [];
    
    % Make 3 copies:
    thetaSurf = [thetaSurf-2*pi, thetaSurf, thetaSurf+2*pi];
    BPSurf = [BPSurf, BPSurf, BPSurf];
    BDotGradThetaSurf = [BDotGradThetaSurf, BDotGradThetaSurf, BDotGradThetaSurf];
    RSurf = [Rs, Rs, Rs];
    ZSurf = [Zs, Zs, Zs];
    
    % Sort:
    [thetaSurf, permutation] = sort(thetaSurf');
    BPSurf = BPSurf(permutation)';
    BDotGradThetaSurf = BDotGradThetaSurf(permutation)';
    RSurf = RSurf(permutation)';
    ZSurf = ZSurf(permutation)';
    
    %bSurf = BSurf / abs(efit.B0EXP);
    
    thetaSurfs{i} = thetaSurf;
    BPSurfs{i} = BPSurf;
    RSurfs{i} = RSurf;
    ZSurfs{i} = ZSurf;
    BDotGradThetaSurfs{i} = BDotGradThetaSurf;
    as(i) = (max(Rs)-min(Rs))/2;
end

if plotStuff
    %fig1 = figure(1)
    fig1 = figure('Visible','off');
    clf
    
    numRows=3;
    numCols=4;
    
    subplot(numRows,numCols,[1, (numCols+1)])
    numContours=20;
    %contourf(R_grid, Z_grid, psiN2D, numContours)
    contourf(R_grid, Z_grid, efit.psi, numContours)
    hold on
    %contour(R_grid, Z_grid, psiN2D, [1, 1],'Color','r')
    for i=1:N
        plot(Rss{i},Zss{i},':c')
    end
    plot(efit.Raxis, efit.Zaxis,'xw')
    plot(efit.R_LCFS,efit.Z_LCFS,'r')
    axis equal
    colorbar
    xlabel('R (m)')
    ylabel('Z (m)')
    title('\psi (Tm^2)')
    
    subplot(numRows,numCols,2)
    plot(psiN,efit.q,'.-')
    hold on
    plot(desiredPsiNs, qSurfs,'or')
    xlabel('\psi_N')
    ylabel('q')
    ylim([0, 7])
    
    subplot(numRows,numCols,3)
    plot(psiN,efit.T)
    xlabel('\psi_N')
    ylabel('I (Tm)')
    
    subplot(numRows,numCols,4)
    plot(psiN,efit.p)
    xlabel('\psi_N')
    ylabel('p (Pa)')
    
    subplot(numRows,numCols,6)
    contourf(R, Z, BPol, 20)
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    title('B_{pol} (T)')
    
    subplot(numRows,numCols,7)
    contourf(R, Z, BTor, 20)
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    title('B_{tor} (T)')
    
    subplot(numRows,numCols,8)
    contourf(R, Z, B, 20)
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    title('B (T)')
    
    subplot(numRows,numCols,9)
    contourf(R, Z, BR, 20)
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    title('B_R (T)')
    
    subplot(numRows,numCols,10)
    contourf(R, Z, BZ, 20)
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    title('B_Z (T)')
    
    subplot(numRows,numCols,11)
    contourf(R, Z, BDotGradTheta, 20)
    colorbar
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    title('B.\nabla\theta (T/m)')
    
    %fig3 = figure(3)
    fig3 = figure('visible','off');
    clf
    
    numRows=2;
    numCols=1;
    
    colors=[...
        0,0,0;...
        1,0.5,0;...
        0,0.8,0;...
        0,0.5,1;...
        0.8, 0.8, 0;...
        0.5,0.5,0.5;...
        0,1,1;
        0.4,0,0.8;
        1,0,0;...
        1, 0, 1;...
        0,0,1    ];
    numColors = size(colors,1);
    
    subplot(numRows,numCols,1)
    for i=1:N
        colorNum = mod(i-1,numColors)+1;
        plot(thetaSurfs{i}, BPSurfs{i},'.-','Color',colors(colorNum,:))
        hold on
    end
    xlabel('\theta')
    ylabel('B_poloidal (T)')
    
    subplot(numRows,numCols,2)
    for i=1:N
        colorNum = mod(i-1,numColors)+1;
        plot(thetaSurfs{i}, BDotGradThetaSurfs{i},'.-','Color',colors(colorNum,:))
        hold on
    end
    xlabel('\theta')
    ylabel('B.\nabla\theta (T/m)')

    %waitfor(fig1);
    %waitfor(fig3);
    set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[40,30])
    print(fig1,'EFITgeometry_testFig1','-dpdf','-bestfit')
    print(fig3,'EFITgeometry_testFig3','-dpdf')
    
end

if saveSeparatrix

  separatrixFileName = 'EFITseparatrix.dat';
  fileID = fopen(separatrixFileName,'w');
  fprintf(fileID,'#Separatrix contour read from EFIT file %s\n',filename);
  fprintf(fileID,'#R / m\t\tZ / m\n');
  fclose(fileID);
  dlmwrite(separatrixFileName,[efit.R_LCFS,efit.Z_LCFS],'delimiter','\t','precision',16,'-append');

end
