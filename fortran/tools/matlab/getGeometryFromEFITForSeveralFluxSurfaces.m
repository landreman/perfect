function [thetaSurfs, BSurfs, BDotGradThetaSurfs, ISurfs, qSurfs, RSurfs, as, R0, B0] = getGeometryFromEFITForSeveralFluxSurfaces(filename, desiredPsiNs, topCropZ, bottomCropZ, innerCropR, outerCropR, plotStuff)

% 'BSurfs' and 'B0' should have units of Tesla.
% Both 'as' and 'R0' should have units of meters.
% 'BDotGradThetaSurfs' should have units of Tesla / meters.

% Good values for topCropZ and bottomCropZ for Alcator C-Mod are 0.4 and -0.4.

if any(desiredPsiNs>1) || any(desiredPsiNs<0)
    error('All desired psi_N values must be between 0 and 1.')
end

efit = read_eqdsk(filename);

psiN = (efit.psi_grid-efit.psiaxis)/(efit.psiedge-efit.psiaxis);
psiN2D = (efit.psi-efit.psiaxis)/(efit.psiedge-efit.psiaxis);

qSurfs = interp1(psiN, efit.q, desiredPsiNs, 'spline');
ISurfs = interp1(psiN, efit.T, desiredPsiNs, 'spline');

valueForCropping = max(max(psiN2D));
psiN2D(efit.Z_grid > topCropZ, :) = valueForCropping;
psiN2D(efit.Z_grid < bottomCropZ, :) = valueForCropping;
psiN2D(efit.R_grid < innerCropR, :) = valueForCropping;
psiN2D(efit.R_grid > outerCropR, :) = valueForCropping;

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
BSurfs = cell(N,1);
RSurfs = cell(N,1);
BDotGradThetaSurfs = cell(N,1);
Rss = cell(N,1);
Zss = cell(N,1);
as = zeros(N,1);

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
    
    BDotGradTheta = BR .* dthetadR + BZ .* dthetadZ;
    %bDotGradTheta = BDotGradTheta ./ B;
    
    thetaSurf = interp2(R, Z, theta, Rs, Zs);
    BSurf = interp2(R, Z, B, Rs, Zs);
    BDotGradThetaSurf = interp2(R, Z, BDotGradTheta, Rs, Zs);
    
    % Sometimes, a point can fall in the tiny sliver of the R-Z plane where theta decreases
    % from 2pi back to 0 as you move clockwise, instead of the main region
    % where theta increases as you move clockwise. We must remove these
    % points.
    thetaIncreasing = (thetaSurf - circshift(thetaSurf, [0,1])) > 0;
    % One such point will always be present, and if there is only 1 point it is not a problem,
    % but if there are more than 1 such point, remove all but the last:
    thetaIncreasing(find(thetaIncreasing,1,'last')) = false;
    thetaSurf(thetaIncreasing) = [];
    BSurf(thetaIncreasing) = [];
    BDotGradThetaSurf(thetaIncreasing) = [];
    
    % Make 3 copies:
    thetaSurf = [thetaSurf-2*pi, thetaSurf, thetaSurf+2*pi];
    BSurf = [BSurf, BSurf, BSurf];
    BDotGradThetaSurf = [BDotGradThetaSurf, BDotGradThetaSurf, BDotGradThetaSurf];
    RSurf = [Rs, Rs, Rs];
    
    % Sort:
    [thetaSurf, permutation] = sort(thetaSurf');
    BSurf = BSurf(permutation)';
    BDotGradThetaSurf = BDotGradThetaSurf(permutation)';
    RSurf = RSurf(permutation)';
    
    %bSurf = BSurf / abs(efit.B0EXP);
    
    thetaSurfs{i} = thetaSurf;
    BSurfs{i} = BSurf;
    RSurfs{i} = RSurf
    BDotGradThetaSurfs{i} = BDotGradThetaSurf;
    as(i) = (max(Rs)-min(Rs))/2;
end

if plotStuff
    fig1 = figure(1)
    clf
    
    numRows=3;
    numCols=4;
    
    subplot(numRows,numCols,[1, (numCols+1)])
    numContours=20;
    contourf(efit.R_grid, efit.Z_grid, psiN2D, numContours)
    hold on
    contour(efit.R_grid, efit.Z_grid, psiN2D, [1, 1],'Color','r')
    for i=1:N
        plot(Rss{i},Zss{i},':c')
    end
    plot(efit.Raxis, efit.Zaxis,'xw')
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
    
    fig3 = figure(3)
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
        plot(thetaSurfs{i}, BSurfs{i},'.-','Color',colors(colorNum,:))
        hold on
    end
    xlabel('\theta')
    ylabel('B (T)')
    
    subplot(numRows,numCols,2)
    for i=1:N
        colorNum = mod(i-1,numColors)+1;
        plot(thetaSurfs{i}, BDotGradThetaSurfs{i},'.-','Color',colors(colorNum,:))
        hold on
    end
    xlabel('\theta')
    ylabel('B.\nabla\theta (T/m)')

    waitfor(fig1);
    waitfor(fig3);
    
end

