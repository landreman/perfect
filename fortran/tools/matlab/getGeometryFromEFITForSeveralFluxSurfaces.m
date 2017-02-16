function [thetaSurfs, BPSurfs, BDotGradThetaSurfs, ISurfs, qSurfs, RSurfs, ZSurfs, as, R0, Z0, B0, psi0] = getGeometryFromEFITForSeveralFluxSurfaces(filename, desiredPsiNs, topCropZ, bottomCropZ, innerCropR, outerCropR, plotStuff, saveSeparatrix)

% 'BPSurfs' and 'B0' should have units of Tesla.
% Both 'as' and 'R0' should have units of meters.
% 'BDotGradThetaSurfs' should have units of Tesla / meters.

% Good values for topCropZ and bottomCropZ for Alcator C-Mod are 0.4 and -0.4.

if any(desiredPsiNs>1) || any(desiredPsiNs<0)
    error('All desired psi_N values must be between 0 and 1.')
end

efit = read_eqdsk(filename);

psiN = (efit.psi_grid-efit.psiaxis)/(efit.psiedge-efit.psiaxis);
psiN2D = (efit.psi-efit.psiaxis)/(efit.psiedge-efit.psiaxis);
psi0 = (efit.psiedge-efit.psiaxis); % Needed to set psiAHat for PERFECT

qSurfs = interp1(psiN, efit.q, desiredPsiNs, 'spline');
ISurfs = interp1(psiN, efit.T, desiredPsiNs, 'spline');

valueForCropping = max(max(psiN2D));
psiN2D(efit.Z_grid > topCropZ, :) = valueForCropping;
psiN2D(efit.Z_grid < bottomCropZ, :) = valueForCropping;
psiN2D(:, efit.R_grid < innerCropR) = valueForCropping;
psiN2D(:, efit.R_grid > outerCropR) = valueForCropping;

R0 = efit.Raxis;
Z0 = efit.Zaxis;
B0 = efit.B0EXP;

scheme = 12;
[R, ~, ddR, d2dR2] = differentiationMatricesForUniformGrid(efit.nrbox, efit.rboxleft, efit.rboxleft+efit.rboxlength, scheme);
[Z, ~, ddZ, d2dZ2] = differentiationMatricesForUniformGrid(efit.nzbox, -efit.zboxlength/2, efit.zboxlength/2, scheme);

[R2D, Z2D] = meshgrid(R, Z);
theta = atan2(Z2D - efit.Zaxis, R2D-efit.Raxis);

dpsidZ = ddZ * efit.psi;
dpsidR = (ddR * (efit.psi'))';

BZ = dpsidR./R2D;
BR = -dpsidZ./R2D;

BPol = sqrt(dpsidZ.^2 + dpsidR.^2)./ R2D;
I2D = interp1(efit.psi_grid, efit.T, efit.psi,'spline',efit.T(end));
BTor = I2D ./ R2D;

B = sqrt(BPol.^2 + BTor.^2);

hypotenuse = (R2D-R0).^2 + (Z2D-Z0).^2;
dthetadZ = (R2D-R0) ./ hypotenuse;
dthetadR = -(Z2D-Z0) ./ hypotenuse;

N = numel(desiredPsiNs);
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
    c = contourc(efit.R_grid, efit.Z_grid, psiN2D, [desiredPsiN, desiredPsiN]);
    
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

    BPSurf = interp2(R, Z, BPol, Rs, Zs,'spline');
    BDotGradThetaSurf = interp2(R, Z, BDotGradTheta, Rs, Zs,'spline');
    
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
    fig0 = figure('Visible','off');
    %contour(efit.R_grid, efit.Z_grid, efit.psi, efit.psiaxis+efit.psiedge*linspace(.01,1.1,110))
    contour(efit.R_grid, efit.Z_grid, (efit.psi-efit.psiaxis)/efit.psiedge, linspace(.01,1.1,110))
    hold on
    plot(efit.Raxis,efit.Zaxis,'xk')
    plot(efit.R_LCFS,efit.Z_LCFS,'k')
    axis equal
    xlabel('R (m)')
    ylabel('Z (m)')
    title('\psi_N')
    print(fig0,'EFITgeometry_testFig0','-dpdf')

    %fig1 = figure(1)
    fig1 = figure('Visible','off');
    clf
    
    numRows=3;
    numCols=4;
    
    subplot(numRows,numCols,[1, (numCols+1)])
    numContours=20;
    contourf(efit.R_grid, efit.Z_grid, psiN2D, numContours)
    hold on
    %contour(efit.R_grid, efit.Z_grid, psiN2D, [1, 1],'Color','r')
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
    title('B_{pol} (T)')
    
    subplot(numRows,numCols,7)
    contourf(R, Z, BTor, 20)
    colorbar
    title('B_{tor} (T)')
    
    subplot(numRows,numCols,8)
    contourf(R, Z, B, 20)
    colorbar
    title('B (T)')
    
    subplot(numRows,numCols,9)
    contourf(R, Z, BR, 20)
    colorbar
    title('B_R (T)')
    
    subplot(numRows,numCols,10)
    contourf(R, Z, BZ, 20)
    colorbar
    title('B_Z (T)')
    
    subplot(numRows,numCols,11)
    contourf(R, Z, BDotGradTheta, 20)
    colorbar
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
    print(fig1,'EFITgeometry_testFig1','-dpdf')
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
