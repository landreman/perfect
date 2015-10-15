module moments

  use globalVariables
  use grids

#include <finclude/petsckspdef.h>

  implicit none

contains

  !**************************************************************************
  !**************************************************************************
  ! 
  !  Calculate moments of the solution:
  !
  !**************************************************************************
  !**************************************************************************
  subroutine calculateMoments(soln)

    Vec, intent(in) :: soln

    PetscErrorCode :: ierr
    PetscScalar, dimension(:), allocatable :: densityFactors, densityIntegralWeights
    PetscScalar, dimension(:), allocatable :: flowFactors, flowIntegralWeights
    PetscScalar, dimension(:), allocatable :: pressureFactors, pressureIntegralWeights
    PetscScalar, dimension(:), allocatable :: particleFluxFactors, particleFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: momentumFluxFactors, momentumFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: heatFluxFactors, heatFluxIntegralWeights
    Vec :: solnOnProc0
    PetscScalar :: speciesFactor
    VecScatter :: VecScatterContext
    PetscScalar, pointer :: solnArray(:)
    integer, dimension(:), allocatable :: indices
    integer :: ix, itheta, ipsi, L, index
    integer :: ispecies
    integer :: ixi

    ! First, send the entire solution vector to the master process:
    call VecScatterCreateToZero(soln, VecScatterContext, solnOnProc0, ierr)
    CHKERRQ(ierr)
    call VecScatterBegin(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(VecScatterContext, ierr)
    CHKERRQ(ierr)


    if (masterProcInSubComm) then
       ! All computation of moments of the distribution function is then done on the master process:

       allocate(particleSourceProfile(numSpecies,Npsi))
       allocate(heatSourceProfile(numSpecies,Npsi))

       allocate(densityPerturbation(numSpecies,Ntheta,Npsi))
       allocate(flow(numSpecies,Ntheta,Npsi))
       allocate(kPar(numSpecies,Ntheta,Npsi))
       allocate(pressurePerturbation(numSpecies,Ntheta,Npsi))
       allocate(particleFluxBeforeThetaIntegral(numSpecies,Ntheta,Npsi))
       allocate(momentumFluxBeforeThetaIntegral(numSpecies,Ntheta,Npsi))
       allocate(heatFluxBeforeThetaIntegral(numSpecies,Ntheta,Npsi))

       allocate(FSADensityPerturbation(numSpecies,Npsi))
       allocate(kParOutboard(numSpecies,Npsi))
       allocate(kParInboard(numSpecies,Npsi))
       allocate(FSAKPar(numSpecies,Npsi))
       allocate(flowOutboard(numSpecies,Npsi))
       allocate(flowInboard(numSpecies,Npsi))
       allocate(FSABFlow(numSpecies,Npsi))
       allocate(FSAPressurePerturbation(numSpecies,Npsi))
       allocate(particleFlux(numSpecies,Npsi))
       allocate(momentumFlux(numSpecies,Npsi))
       allocate(heatFlux(numSpecies,Npsi))

       allocate(deltaFOutboard(numSpecies,Npsi,NxUniform,NxiUniform))
       allocate(fullFOutboard(numSpecies,Npsi,NxUniform,NxiUniform))

       !       allocate(LHSOfKParEquation(Npsi))
       !       allocate(kThetaWith3PointStencil(Ntheta,Npsi))
       !       allocate(kThetaWith5PointStencil(Ntheta,Npsi))
       !       allocate(kThetaOutboardWith3PointStencil(Npsi))
       !       allocate(kThetaInboardWith3PointStencil(Npsi))
       !       allocate(kThetaOutboardWith5PointStencil(Npsi))
       !       allocate(kThetaInboardWith5PointStencil(Npsi))
       !       allocate(PhiTermInKTheta(Ntheta,Npsi))
       !       allocate(pPerpTermInKThetaWith3PointStencil(Ntheta,Npsi))
       !       allocate(pPerpTermInKThetaWith5PointStencil(Ntheta,Npsi))
       !       allocate(pPerpTermInKThetaBeforePsiDerivative(Ntheta,Npsi))
       !       allocate(ddpsiForKTheta(Npsi,Npsi))

       allocate(flowFactors(Npsi))
       allocate(densityFactors(Npsi))
       allocate(pressureFactors(Npsi))
       allocate(particleFluxFactors(Npsi))
       allocate(momentumFluxFactors(Npsi))
       allocate(heatFluxFactors(Npsi))
       !       allocate(pPerpTermInKThetaFactors(Npsi))

       allocate(flowIntegralWeights(Nx))
       allocate(densityIntegralWeights(Nx))
       allocate(pressureIntegralWeights(Nx))
       allocate(particleFluxIntegralWeights(Nx))
       allocate(momentumFluxIntegralWeights(Nx))
       allocate(heatFluxIntegralWeights(Nx))

       densityIntegralWeights = x*x
       flowIntegralWeights = x*x*x
       pressureIntegralWeights = x*x*x*x
       particleFluxIntegralWeights = x*x*x*x
       momentumFluxIntegralWeights = x*x*x*x*x
       heatFluxIntegralWeights = x*x*x*x*x*x

       allocate(indices(Nx))

       ! Convert the PETSc vector into a normal Fortran array:
       call VecGetArrayF90(solnOnProc0, solnArray, ierr)
       CHKERRQ(ierr)

       do ispecies = 1,numSpecies
          densityFactors = Delta*4*pi*THats(ispecies,:)*sqrtTHats(ispecies,:) &
               / (nHats(ispecies,:)*masses(ispecies)*sqrt(masses(ispecies)))
          flowFactors = 4*pi/(three*nHats(ispecies,:)) * ((THats(ispecies,:)/masses(ispecies)) ** 2)
          pressureFactors = Delta*8*pi/(three*nHats(ispecies,:)) * ((THats(ispecies,:)/masses(ispecies)) ** (1.5d+0))
          particleFluxFactors = -masses(ispecies) / charges(ispecies) * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
          momentumFluxFactors = -masses(ispecies) / charges(ispecies) * IHat*IHat * ((THats(ispecies,:)/masses(ispecies)) ** 3)
          heatFluxFactors = -masses(ispecies) / charges(ispecies) * THats(ispecies,:) &
               * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
          !       pPerpTermInKThetaFactors = THat ** (5/two)

          ! The final elements of the solution vector correspond to the source profiles:
          do ipsi=1,Npsi
             particleSourceProfile(ispecies,ipsi) = solnArray(localMatrixSize*Npsi + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 1)
             heatSourceProfile(ispecies,ipsi) = solnArray(localMatrixSize*Npsi + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 2)
          end do

          L = 0
          do ipsi=1,Npsi
             do itheta=1,Ntheta
                indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                     + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

                densityPerturbation(ispecies,itheta,ipsi) = dot_product(xWeights, densityIntegralWeights * solnArray(indices)) &
                     * densityFactors(ipsi)

                pressurePerturbation(ispecies,itheta,ipsi) = dot_product(xWeights, pressureIntegralWeights * solnArray(indices)) &
                     * pressureFactors(ipsi)

                particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = (8/three) * particleFluxFactors(ipsi) &
                     * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices))

                heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = (8/three) * heatFluxFactors(ipsi) &
                     * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices))

                !             pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) = &
                !                  (4/three) * pPerpTermInKThetaFactors(ipsi) &
                !                  * dot_product(xWeights, pressureIntegralWeights * solnArray(indices))

             end do
          end do

          L = 1
          do ipsi=1,Npsi
             do itheta=1,Ntheta
                indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                     + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

                flow(ispecies,itheta,ipsi) = dot_product(xWeights, flowIntegralWeights * solnArray(indices)) &
                     * flowFactors(ipsi)

                momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = ((16d+0)/15) * momentumFluxFactors(ipsi) &
                     * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices))

             end do
          end do

          L = 2
          do ipsi=1,Npsi
             do itheta=1,Ntheta
                indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                     + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

                particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                     + (four/15) * particleFluxFactors(ipsi) &
                     * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices))

                heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                     + (four/15) * heatFluxFactors(ipsi) &
                     * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices))

                !             pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) = &
                !                  pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) &
                !                  - ((4d+0)/15) * pPerpTermInKThetaFactors(ipsi) &
                !                  * dot_product(xWeights, pressureIntegralWeights * solnArray(indices))

             end do
          end do

          L = 3
          do ipsi=1,Npsi
             do itheta=1,Ntheta
                indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                     + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + itheta

                momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                     + (four/35) * momentumFluxFactors(ipsi) &
                     * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices))

             end do
          end do

          do itheta=1,Ntheta
             kpar(ispecies,itheta,:) = FSABHat2/(BHat(itheta,:)*BHat(itheta,:)*dTHatdpsis(ispecies,:)) &
                  * (2*charges(ispecies)*psiAHat*BHat(itheta,:)/IHat*flow(ispecies,itheta,:) + dTHatdpsis(ispecies,:) &
                  + THats(ispecies,:)/nHats(ispecies,:)*dnHatdpsis(ispecies,:) + 2*charges(ispecies)*omega/Delta*dPhiHatdpsi)
          end do

          particleFluxBeforeThetaIntegral(ispecies,:,:) = particleFluxBeforeThetaIntegral(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat)

          momentumFluxBeforeThetaIntegral(ispecies,:,:) = momentumFluxBeforeThetaIntegral(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat * BHat)

          heatFluxBeforeThetaIntegral(ispecies,:,:) = heatFluxBeforeThetaIntegral(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat)

   !!$         if (psiDerivativeScheme == 0) then
   !!$            ddpsiForKTheta = ddpsiLeft
   !!$         else
   !!$            ! centered finite differences, no upwinding, 3-point stencil
   !!$            scheme = 2
   !!$            call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForKTheta, d2dpsi2)
   !!$         end if
   !!$         !       pPerpTermInKThetaWith3PointStencil = transpose(matmul(ddpsiForKTheta, &
   !!$         !            transpose(pPerpTermInKThetaBeforePsiDerivative)))
   !!$
   !!$         if (psiDerivativeScheme == 0) then
   !!$            ddpsiForKTheta = ddpsiLeft
   !!$         else
   !!$            ! centered finite differences, no upwinding, 5-point stencil
   !!$            scheme = 12
   !!$            call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForKTheta, d2dpsi2)
   !!$         end if
   !!$         !       pPerpTermInKThetaWith5PointStencil = transpose(matmul(ddpsiForKTheta, &
   !!$         !            transpose(pPerpTermInKThetaBeforePsiDerivative)))

          do ipsi = 1,Npsi
             particleFlux(ispecies,ipsi) = dot_product(thetaWeights, particleFluxBeforeThetaIntegral(ispecies,:,ipsi))
             momentumFlux(ispecies,ipsi) = dot_product(thetaWeights, momentumFluxBeforeThetaIntegral(ispecies,:,ipsi))
             heatFlux(ispecies,ipsi) = dot_product(thetaWeights, heatFluxBeforeThetaIntegral(ispecies,:,ipsi))

             !          pPerpTermInKThetaWith3PointStencil(:,ipsi) = pPerpTermInKThetaWith3PointStencil(:,ipsi) &
             !               * 2 * pi * delta * FSABHat2(ipsi) / (nHat(ipsi) * BHat(:,ipsi) * BHat(:,ipsi) * dTHatdpsi(ipsi))
             !
             !          pPerpTermInKThetaWith5PointStencil(:,ipsi) = pPerpTermInKThetaWith5PointStencil(:,ipsi) &
             !               * 2 * pi * delta * FSABHat2(ipsi) / (nHat(ipsi) * BHat(:,ipsi) * BHat(:,ipsi) * dTHatdpsi(ipsi))
             !
             !          PhiTermInKTheta(:,ipsi) = densityPerturbation(:,ipsi) * 2 * omega * dphidpsi(ipsi) * FSABHat2(ipsi) &
             !               / (delta * BHat(:,ipsi) * BHat(:,ipsi) * dTHatdpsi(ipsi))

             FSADensityPerturbation(ispecies,ipsi) = dot_product(thetaWeights, &
                  densityPerturbation(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

             FSABFlow(ispecies,ipsi) = dot_product(thetaWeights, flow(ispecies,:,ipsi)*BHat(:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

             FSAkPar(ispecies,ipsi) = dot_product(thetaWeights, kPar(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

             FSAPressurePerturbation(ispecies,ipsi) = dot_product(thetaWeights, &
                  pressurePerturbation(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)
          end do

          !       LHSOfKParEquation = FSAKPar - omega*delta/(psiAHat*psiAHat) * IHat*IHat * dphidpsi / FSABHat2 &
          !            * matmul(ddpsiLeft,FSAKPar)

          !       kThetaWith3PointStencil = kpar + PhiTermInKTheta + pPerpTermInKThetaWith3PointStencil
          !       kThetaWith5PointStencil = kpar + PhiTermInKTheta + pPerpTermInKThetaWith5PointStencil

          kParOutboard(ispecies,:) = kPar(ispecies,1,:)
          flowOutboard(ispecies,:) = flow(ispecies,1,:)
          !       kThetaOutboardWith3PointStencil = kThetaWith3PointStencil(1,:)
          !       kThetaOutboardWith5PointStencil = kThetaWith5PointStencil(1,:)
          if (mod(Ntheta,2)==0) then
             kParInboard(ispecies,:) = kPar(ispecies,Ntheta/2+1,:)
             flowInboard(ispecies,:) = flow(ispecies,Ntheta/2+1,:)
             !          kThetaInboardWith3PointStencil = kThetaWith3PointStencil(Ntheta/2+1,:)
             !          kThetaInboardWith5PointStencil = kThetaWith5PointStencil(Ntheta/2+1,:)
          else
             index = (Ntheta+1)/2
             kParInboard(ispecies,:) = (kPar(ispecies,index,:) + kPar(ispecies,index+1,:))*oneHalf
             flowInboard(ispecies,:) = (flow(ispecies,index,:) + flow(ispecies,index+1,:))*oneHalf
             !          kThetaInboardWith3PointStencil = (kThetaWith3PointStencil(index,:) + kThetaWith3PointStencil(index+1,:))*oneHalf
             !          kThetaInboardWith5PointStencil = (kThetaWith5PointStencil(index,:) + kThetaWith5PointStencil(index+1,:))*oneHalf
          end if

       end do

       LegendresOnXiUniform_m1 = 1
       deltaFOutboard = 0
       do L = 0,(Nxi-1)
          ! Recursively evaluate Legendre polynomials on a uniform grid in xi.
          ! The results will be used to map the distribution function from the modal discretization
          ! to a uniform grid.
          if (L == 0) then
             LegendresOnXiUniform = 1
          elseif (L == 1) then
             LegendresOnXiUniform = xiUniform
          else
             LegendresOnXiUniform_m2 = LegendresOnXiUniform_m1
             LegendresOnXiUniform_m1 = LegendresOnXiUniform
             LegendresOnXiUniform = ((2*L-1)*xiUniform * LegendresOnXiUniform_m1 - (L-1)*LegendresOnXiUniform_m2)/L
         end if
         
         do ispecies=1,numSpecies
            do ipsi = 1,Npsi
               indices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + [(ix-1, ix=1,Nx)]*Nxi*Ntheta + L*Ntheta + thetaIndexForOutboard

               do ixi = 1,NxiUniform
                  deltaFOutboard(ispecies,ipsi,:,ixi) = deltaFOutboard(ispecies,ipsi,:,ixi) + &
                       LegendresOnXiUniform(ixi) * matmul(regridPolynomialToUniformForDiagnostics, solnArray(indices))
               end do
            end do
         end do
      end do

      ! Now build the full f from delta f:
      do ispecies = 1,numSpecies
         do ipsi = 1,Npsi
            speciesFactor = nHats(ispecies,ipsi)/(pi*sqrtpi*((THats(ispecies,ipsi)/masses(ispecies)) ** (1.5d+0)))
            do ixi = 1,NxiUniform
               do ix = 1,NxUniform
                  fullFOutboard(ispecies,ipsi,ix,ixi) = deltaFOutboard(ispecies,ipsi,ix,ixi) * Delta &
                       + speciesFactor * exp(-xUniform(ix)*xUniform(ix)) ! This line adds the Maxwellian F_M
               end do
            end do
         end do
      end do

      deallocate(indices)

      call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
      CHKERRQ(ierr)

      call VecDestroy(solnOnProc0, ierr) !dubious
   end if

  end subroutine calculateMoments

end module moments
