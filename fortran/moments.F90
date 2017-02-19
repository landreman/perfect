module moments

  use globalVariables
  use grids
  use writeHDF5Output ! used to write output for debugging
  use geometry ! used to get poloidal and toroidal flows (only Miller supported ATM)

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

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
    PetscScalar, dimension(:), allocatable :: temp, pPerpTermInVpFactors
    PetscScalar, dimension(:), allocatable :: particleFluxPotentialFactors, expx2
    
    
    Vec :: solnOnProc0
    PetscScalar :: speciesFactor
    VecScatter :: VecScatterContext
    PetscScalar, pointer :: solnArray(:)
    PetscScalar, dimension(:), allocatable :: solnAtL
    integer, dimension(:), allocatable :: indices
    integer :: ix, itheta, ipsi, L, index
    integer :: ispecies, isources, iextraSources, ispeciesIndepSources
    integer :: ixi
    integer :: ispecies2
    integer, dimension(:), allocatable :: indices2
    PetscScalar :: temp1

    PetscScalar :: oneOverAbsNablaTheta,oneOverAbsNablaPhi,signBP
    

    ! First, send the entire solution vector to the master process:
    call VecScatterCreateToZero(soln, VecScatterContext, solnOnProc0, ierr)
    CHKERRQ(ierr)
    call VecScatterBegin(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(VecScatterContext, ierr)
    CHKERRQ(ierr)


    if (masterProcInSubComm) then
       ! All computation of moments of the distribution function is then done on the master process:

       allocate(sourceProfile(Nsources,Nspecies,Npsi-NpsiSourcelessLeft-NpsiSourcelessRight))
       allocate(densityPerturbation(Nspecies,Ntheta,Npsi))
       allocate(flow(Nspecies,Ntheta,Npsi))
       allocate(pPerpTermInVp(Nspecies,Ntheta,Npsi))
       allocate(pPerpTermInVpBeforePsiDerivative(Nspecies,Ntheta,Npsi))
       allocate(toroidalFlow(Nspecies,Ntheta,Npsi))
       allocate(poloidalFlow(Nspecies,Ntheta,Npsi))
       allocate(kPar(Nspecies,Ntheta,Npsi))
       allocate(pressurePerturbation(Nspecies,Ntheta,Npsi))
       allocate(particleFluxBeforeThetaIntegral(Nspecies,Ntheta,Npsi))
       allocate(momentumFluxBeforeThetaIntegral(Nspecies,Ntheta,Npsi))
       allocate(heatFluxBeforeThetaIntegral(Nspecies,Ntheta,Npsi))

       allocate(FSADensityPerturbation(Nspecies,Npsi))
       allocate(kParOutboard(Nspecies,Npsi))
       allocate(kParInboard(Nspecies,Npsi))
       allocate(FSAKPar(Nspecies,Npsi))
       allocate(flowOutboard(Nspecies,Npsi))
       allocate(flowInboard(Nspecies,Npsi))
       allocate(FSAToroidalFlow(Nspecies,Npsi))
       allocate(FSAPoloidalFlow(Nspecies,Npsi))
       allocate(FSAFlow(Nspecies,Npsi))
       allocate(FSABFlow(Nspecies,Npsi))
       allocate(FSAPressurePerturbation(Nspecies,Npsi))
       allocate(particleFlux(Nspecies,Npsi))
       allocate(momentumFlux(Nspecies,Npsi))
       allocate(heatFlux(Nspecies,Npsi))
       
       allocate(deltaFOutboard(Nspecies,Npsi,NxUniform,NxiUniform))
       allocate(fullFOutboard(Nspecies,Npsi,NxUniform,NxiUniform))

       allocate(particleFluxPotentialBeforeSums(Nspecies,Ntheta,Npsi))
       allocate(particleFluxPotential(Nspecies,Npsi))

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
       allocate(pPerpTermInVpFactors(Npsi))
       allocate(particleFluxPotentialFactors(Npsi))

       allocate(flowIntegralWeights(Nx))
       allocate(densityIntegralWeights(Nx))
       allocate(pressureIntegralWeights(Nx))
       allocate(particleFluxIntegralWeights(Nx))
       allocate(momentumFluxIntegralWeights(Nx))
       allocate(heatFluxIntegralWeights(Nx))

       

       ! to calculate poloidal and toroidal flow and flow differentials
       allocate(temp(Npsi))
       ! and sums needed for the potential contribution to particle fluxes
       allocate(expx2(Nx))
       expx2 = exp(-x*x)
       
       allocate(extraSourceProfile(NextraSources,Nspecies,Npsi-NpsiSourcelessLeft-NpsiSourcelessRight))
       allocate(speciesIndepSourceProfile(NspeciesIndepSources,Nspecies,Npsi))
       
       densityIntegralWeights = x*x
       flowIntegralWeights = x*x*x
       pressureIntegralWeights = x*x*x*x
       particleFluxIntegralWeights = x*x*x*x
       momentumFluxIntegralWeights = x*x*x*x*x
       heatFluxIntegralWeights = x*x*x*x*x*x

       
       
       allocate(indices(Nx))

       allocate(indices2(Nx))


       ! Convert the PETSc vector into a normal Fortran array:
       call VecGetArrayF90(solnOnProc0, solnArray, ierr)
       CHKERRQ(ierr)
       
       do ispecies = 1,Nspecies
          densityFactors = Delta*4*pi*THats(ispecies,:)*sqrtTHats(ispecies,:) &
               / (nHats(ispecies,:)*masses(ispecies)*sqrt(masses(ispecies)))
          flowFactors = 4*pi/(three*nHats(ispecies,:)) * ((THats(ispecies,:)/masses(ispecies)) ** 2)
          pressureFactors = Delta*8*pi/(three*nHats(ispecies,:)) * ((THats(ispecies,:)/masses(ispecies)) ** (1.5d+0))
          particleFluxFactors = -masses(ispecies) / charges(ispecies) * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
          momentumFluxFactors = -masses(ispecies) / charges(ispecies) * IHat*IHat * ((THats(ispecies,:)/masses(ispecies)) ** 3)
          heatFluxFactors = -masses(ispecies) / charges(ispecies) * THats(ispecies,:) &
               * IHat * ((THats(ispecies,:)/masses(ispecies)) ** (5/two))
          pPerpTermInVpFactors = (8*pi/3)*(THats(ispecies,:)/masses(ispecies))**(5/two)
          !       pPerpTermInKThetaFactors = THat ** (5/two)
          particleFluxPotentialFactors = -particleFluxFactors * charges(ispecies) * nHats(ispecies,:) * masses(ispecies)**(3/two) &
               /(THats(ispecies,:)**(5/two) * 2**(3/two) * pi*sqrtpi)

          ! The final elements of the solution vector correspond to the source profiles:
          do ipsi=lowestEnforcedIpsi,highestEnforcedIpsi
             do isources = 1,Nsources
                sourceProfile(isources,ispecies,ipsi - lowestEnforcedIpsi + 1) = &
                     solnArray(getIndexSources(isources,ispecies,ipsi) + 1)
             end do
          end do

          do iextraSources = 1,NextraSources
             do ipsi=lowestEnforcedIpsi,highestEnforcedIpsi
                extraSourceProfile(iextraSources,ispecies,ipsi - lowestEnforcedIpsi + 1) = &
                     extraSourceSpeciesPart(iextraSources,ispecies)*solnArray(getIndexExtraSources(iextraSources,ipsi) + 1)
             end do
          end do

          do ispeciesIndepSources = 1,NspeciesIndepSources
             do ipsi=1,Npsi
                speciesIndepSourceProfile(ispeciesIndepSources,ispecies,ipsi - lowestEnforcedIpsi + 1) = &
                     speciesIndepSourceSpeciesPart(ispeciesIndepSources,ispecies) &
                     * solnArray(getIndexSpeciesIndepSources(ispeciesIndepSources,ipsi) + 1)
             end do
          end do
          
          L = 0
          do ipsi=1,Npsi
             do itheta=1,Ntheta
	        indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

                densityPerturbation(ispecies,itheta,ipsi) = dot_product(xWeights, densityIntegralWeights * solnArray(indices+1)) &
                     * densityFactors(ipsi)

                pressurePerturbation(ispecies,itheta,ipsi) = dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1)) &
                     * pressureFactors(ipsi)

                particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = (8/three) * particleFluxFactors(ipsi) &
                     * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices+1))

                particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) = zero
                do ispecies2 = 1, Nspecies
                   indices2 = [(getIndex(ispecies2,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]
                   particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) = particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) &
                        + (8/three) * charges(ispecies2) &
                        * dot_product(xWeights, expx2 * particleFluxIntegralWeights * solnArray(indices2+1))
                end do

                heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = (8/three) * heatFluxFactors(ipsi) &
                     * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices+1))

                pPerpTermInVpBeforePsiDerivative(ispecies,itheta,ipsi) = pPerpTermInVpFactors(ipsi) &
                     * dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1)) 

                
                !             pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) = &
                !                  (4/three) * pPerpTermInKThetaFactors(ipsi) &
                !                  * dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1))

             end do
          end do

          L = 1
          do ipsi=1,Npsi
             do itheta=1,Ntheta
	        indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

                flow(ispecies,itheta,ipsi) = dot_product(xWeights, flowIntegralWeights * solnArray(indices+1)) &
                     * flowFactors(ipsi)

                momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = ((16d+0)/15) * momentumFluxFactors(ipsi) &
                     * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices+1))

             end do
          end do

          L = 2
          do ipsi=1,Npsi
             do itheta=1,Ntheta
	     	indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

                particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = particleFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                     + (four/15) * particleFluxFactors(ipsi) &
                     * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices+1))

                temp1 = zero
                do ispecies2 = 1, Nspecies
                   indices2 = [(getIndex(ispecies2,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]
                   particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) = particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) &
                        + (four/15) * charges(ispecies2) &
                        * dot_product(xWeights, expx2 * particleFluxIntegralWeights * solnArray(indices2+1))

                   temp1 = temp1 + charges(ispecies2)**2 * nHats(ispecies2,ipsi)/THats(ispecies2,ipsi)
                end do
                particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) = particleFluxPotentialBeforeSums(ispecies,itheta,ipsi) &
                     * particleFluxPotentialFactors(ipsi)/temp1

                heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = heatFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                     + (four/15) * heatFluxFactors(ipsi) &
                     * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices+1))

                pPerpTermInVpBeforePsiDerivative(ispecies,itheta,ipsi) = pPerpTermInVpBeforePsiDerivative(ispecies,itheta,ipsi) &
                     - (one/five) * pPerpTermInVpFactors(ipsi) &
                     * dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1))

                !             pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) = &
                !                  pPerpTermInKThetaBeforePsiDerivative(itheta,ipsi) &
                !                  - ((4d+0)/15) * pPerpTermInKThetaFactors(ipsi) &
                !                  * dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1))

             end do
          end do

          L = 3
          do ipsi=1,Npsi
             do itheta=1,Ntheta
	        indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

                momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) = momentumFluxBeforeThetaIntegral(ispecies,itheta,ipsi) &
                     + (four/35) * momentumFluxFactors(ipsi) &
                     * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices+1))

             end do
          end do

          do itheta=1,Ntheta
             kpar(ispecies,itheta,:) = FSABHat2/(BHat(itheta,:)*BHat(itheta,:)*dTHatdpsis(ispecies,:)) &
                  * (2*charges(ispecies)*psiAHatArray(:)*BHat(itheta,:)/IHat*flow(ispecies,itheta,:) + dTHatdpsis(ispecies,:) &
                  + THats(ispecies,:)/nHats(ispecies,:)*dnHatdpsis(ispecies,:) + 2*charges(ispecies)*omega/Delta*dPhiHatdpsi)
          end do

          particleFluxBeforeThetaIntegral(ispecies,:,:) = particleFluxBeforeThetaIntegral(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat)

          momentumFluxBeforeThetaIntegral(ispecies,:,:) = momentumFluxBeforeThetaIntegral(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat * BHat)

          heatFluxBeforeThetaIntegral(ispecies,:,:) = heatFluxBeforeThetaIntegral(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat)

          particleFluxPotentialBeforeSums(ispecies,:,:) = particleFluxPotentialBeforeSums(ispecies,:,:) &
               * dBHatdtheta / (BHat * BHat * BHat)

          do itheta=1,Ntheta
             pPerpTermInVp(ispecies,itheta,:) = matmul(ddpsiLeft, pPerpTermInVpBeforePsiDerivative(ispecies,itheta,:))
             toroidalFlow(ispecies,itheta,:) = pPerpTermInVp(ispecies,itheta,:) 
             poloidalFlow(ispecies,itheta,:) = pPerpTermInVp(ispecies,itheta,:)
             
             toroidalFlow(ispecies,itheta,:)=(Delta/(2*psiAHat)) &
                  *(masses(ispecies))/(charges(ispecies)*BHat(itheta,:)**2*nHats(ispecies,:)) &
                  * (-BPHat(itheta,:)**2) * RHat(itheta,:) * toroidalFlow(ispecies,itheta,:)
             poloidalFlow(ispecies,itheta,:)=(Delta/(2*psiAHat)) &
                  *(masses(ispecies))/(charges(ispecies)*BHat(itheta,:)**2*nHats(ispecies,:))&
                  *BPHat(itheta,:)*IHat(:)* poloidalFlow(ispecies,itheta,:) 

             toroidalFlow(ispecies,itheta,:) = toroidalFlow(ispecies,itheta,:) &
                  + (BTHat(itheta,:)/BHat(itheta,:))*flow(ispecies,itheta,:)
             poloidalFlow(ispecies,itheta,:) = poloidalFlow(ispecies,itheta,:) &
                  + (BPHat(itheta,:)/BHat(itheta,:))*flow(ispecies,itheta,:)

             toroidalFlow(ispecies,itheta,:) = toroidalFlow(ispecies,itheta,:) -omega/(Delta*psiAHat)*dPhiHatdpsi &
                  *densityPerturbation(ispecies,itheta,:)*(BPHat(itheta,:)**2*RHat(itheta,:))/(BHat(itheta,:)**2)
             poloidalFlow(ispecies,itheta,:) = poloidalFlow(ispecies,itheta,:) +omega/(Delta*psiAHat)*dPhiHatdpsi &
                  *densityPerturbation(ispecies,itheta,:)*(BPHat(itheta,:)*IHat)/(BHat(itheta,:)**2)

             !since we are not using that variable for anything else
             temp = dnHatdpsis(ispecies,:)/nHats(ispecies,:) +dTHatdpsis(ispecies,:)/THats(ispecies,:) &
                  + (2*omega*charges(ispecies)/(Delta*THats(ispecies,:)))*dPhiHatdpsi
             toroidalFlow(ispecies,itheta,:) = toroidalFlow(ispecies,itheta,:) &
                  -THats(ispecies,:)/(2*psiAHat*charges(ispecies)*BHat(itheta,:)**2)*temp*BPHat(itheta,:)**2*RHat(itheta,:)
             poloidalFlow(ispecies,itheta,:) = poloidalFlow(ispecies,itheta,:) &
                  +THats(ispecies,:)/(2*psiAHat*charges(ispecies)*BHat(itheta,:)**2)*temp*BPHat(itheta,:)*IHat(:)
          end do
          
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
             particleFluxPotential(ispecies,ipsi) = dot_product(thetaWeights, particleFluxPotentialBeforeSums(ispecies,:,ipsi))

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

             FSAFlow(ispecies,ipsi) = dot_product(thetaWeights, flow(ispecies,:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

             FSAToroidalFlow(ispecies,ipsi) = dot_product(thetaWeights, toroidalFlow(ispecies,:,ipsi)/JHat(:,ipsi)) /VPrimeHat(ipsi)

             FSAPoloidalFlow(ispecies,ipsi) = dot_product(thetaWeights, poloidalFlow(ispecies,:,ipsi)/JHat(:,ipsi)) /VPrimeHat(ipsi)

             
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
       allocate(solnAtL(Nx))
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
         
         do ispecies=1,Nspecies
            do ipsi = 1,Npsi
	       indices = [(getIndex(ispecies,ix,L,thetaIndexForOutboard,ipsi), ix=min_x_for_L(L),Nx)]
               solnAtL(1:min_x_for_L(L)-1) = 0d0
               solnAtL(min_x_for_L(L):Nx) = solnArray(indices+1)

               do ixi = 1,NxiUniform
                  deltaFOutboard(ispecies,ipsi,:,ixi) = deltaFOutboard(ispecies,ipsi,:,ixi) + &
                       LegendresOnXiUniform(ixi) * matmul(regridPolynomialToUniformForDiagnostics, solnAtL)
               end do
            end do
         end do
      end do

      ! Now build the full f from delta f:
      do ispecies = 1,Nspecies
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

      deallocate(solnAtL)
      deallocate(indices)

      call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
      CHKERRQ(ierr)

      call VecDestroy(solnOnProc0, ierr) !dubious
    end if

  end subroutine calculateMoments

  !**************************************************************************
  !**************************************************************************
  ! 
  !  Calculate moments of a local solution, e.g. for a boundary:
  !
  !**************************************************************************
  !**************************************************************************
  subroutine calculateLocalMoments(soln,ipsi,filename)
    Vec, intent(in) :: soln
    integer, intent(in) :: ipsi
    character(len=*), intent(in) :: filename

    PetscErrorCode :: ierr
    PetscScalar :: flowFactors, densityFactors, pressureFactors,&
         particleFluxFactors, momentumFluxFactors, heatFluxFactors
    PetscScalar, dimension(:), allocatable :: densityIntegralWeights
    PetscScalar, dimension(:), allocatable :: flowIntegralWeights
    PetscScalar, dimension(:), allocatable :: pressureIntegralWeights
    PetscScalar, dimension(:), allocatable :: particleFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: momentumFluxIntegralWeights
    PetscScalar, dimension(:), allocatable :: heatFluxIntegralWeights
    PetscScalar, dimension(:,:), allocatable :: this_densityPerturbation, this_flow, this_kPar,&
                                                this_pPerpTermInVpBeforePsiDerivative,this_pPerpTermInVp,&
                                                this_poloidalFlow,this_toroidalFlow,&
                                                this_pressurePerturbation, this_particleFluxBeforeThetaIntegral,&
                                                this_momentumFluxBeforeThetaIntegral, this_heatFluxBeforeThetaIntegral
    PetscScalar, dimension(:), allocatable :: this_FSADensityPerturbation, this_kParOutboard,&
                                              this_kParInboard, this_FSAKPar, this_flowOutboard,&
                                              this_flowInboard, this_FSAFlow, this_FSABFlow, this_FSAPressurePerturbation,&
                                              this_particleFlux, this_momentumFlux, this_heatFlux
    PetscScalar, dimension(:,:,:), allocatable :: this_deltaFOutboard, this_fullFOutboard
    PetscScalar :: speciesFactor
    !! Vec :: solnOnProc0
    !! VecScatter :: VecScatterContext
    PetscScalar, pointer :: solnArray(:)
    PetscScalar, dimension(:), allocatable :: solnAtL
    integer, dimension(:), allocatable :: indices
    integer :: ix, itheta, L, index
    integer :: ispecies
    integer :: ixi
    
    PetscScalar :: this_temp
    ! should only be called by one processor
    !! ! First, send the entire solution vector to the master process:
    !! call VecScatterCreateToZero(soln, VecScatterContext, solnOnProc0, ierr)
    !! CHKERRQ(ierr)
    !! call VecScatterBegin(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !! call VecScatterEnd(VecScatterContext, soln, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !! call VecScatterDestroy(VecScatterContext, ierr)
    !! CHKERRQ(ierr)


    !! if (masterProcInSubComm) then
       ! All computation of moments of the distribution function is then done on the master process:

       allocate(this_densityPerturbation(Nspecies,Ntheta))
       allocate(this_flow(Nspecies,Ntheta))
       allocate(this_pPerpTermInVpBeforePsiDerivative(Nspecies,Ntheta))
       allocate(this_pPerpTermInVp(Nspecies,Ntheta))
       allocate(this_poloidalFlow(Nspecies,Ntheta))
       allocate(this_toroidalFlow(Nspecies,Ntheta))
       allocate(this_kPar(Nspecies,Ntheta))
       allocate(this_pressurePerturbation(Nspecies,Ntheta))
       allocate(this_particleFluxBeforeThetaIntegral(Nspecies,Ntheta))
       allocate(this_momentumFluxBeforeThetaIntegral(Nspecies,Ntheta))
       allocate(this_heatFluxBeforeThetaIntegral(Nspecies,Ntheta))

       allocate(this_FSADensityPerturbation(Nspecies))
       allocate(this_kParOutboard(Nspecies))
       allocate(this_kParInboard(Nspecies))
       allocate(this_FSAKPar(Nspecies))
       allocate(this_flowOutboard(Nspecies))
       allocate(this_flowInboard(Nspecies))
       allocate(this_FSAFlow(Nspecies))
       allocate(this_FSABFlow(Nspecies))
       allocate(this_FSAPressurePerturbation(Nspecies))
       allocate(this_particleFlux(Nspecies))
       allocate(this_momentumFlux(Nspecies))
       allocate(this_heatFlux(Nspecies))

       allocate(this_deltaFOutboard(Nspecies,NxUniform,NxiUniform))
       allocate(this_fullFOutboard(Nspecies,NxUniform,NxiUniform))

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
       !! call VecGetArrayF90(solnOnProc0, solnArray, ierr)
       call VecGetArrayF90(soln, solnArray, ierr)
       CHKERRQ(ierr)
       do ispecies = 1,Nspecies
          densityFactors = Delta*4*pi*THats(ispecies,ipsi)*sqrtTHats(ispecies,ipsi) &
               / (nHats(ispecies,ipsi)*masses(ispecies)*sqrt(masses(ispecies)))
          flowFactors = 4*pi/(three*nHats(ispecies,ipsi)) * ((THats(ispecies,ipsi)/masses(ispecies)) ** 2)
          pressureFactors = Delta*8*pi/(three*nHats(ispecies,ipsi)) * ((THats(ispecies,ipsi)/masses(ispecies)) ** (1.5d+0))
          particleFluxFactors = -masses(ispecies) / charges(ispecies) * IHat(ipsi)&
                                * ((THats(ispecies,ipsi)/masses(ispecies)) ** (5/two))
          momentumFluxFactors = -masses(ispecies) / charges(ispecies) * IHat(ipsi)*IHat(ipsi)&
                                * ((THats(ispecies,ipsi)/masses(ispecies)) ** 3)
          heatFluxFactors = -masses(ispecies) / charges(ispecies) * THats(ispecies,ipsi) &
               * IHat(ipsi) * ((THats(ispecies,ipsi)/masses(ispecies)) ** (5/two))
          
          L = 0
          do itheta=1,Ntheta
	     indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

             this_densityPerturbation(ispecies,itheta) = dot_product(xWeights, densityIntegralWeights * solnArray(indices+1)) &
                  * densityFactors

             this_pressurePerturbation(ispecies,itheta) = dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1)) &
                  * pressureFactors

             this_particleFluxBeforeThetaIntegral(ispecies,itheta) = (8/three) * particleFluxFactors &
                  * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices+1))

             this_heatFluxBeforeThetaIntegral(ispecies,itheta) = (8/three) * heatFluxFactors &
                  * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices+1))

             this_pPerpTermInVpBeforePsiDerivative(ispecies,itheta) &
                 = dot_product(xWeights, pressureIntegralWeights * solnArray(indices+1)) &
                  * 2*pi*(4/three)*THats(ispecies,ipsi)**(5.0/2.0)/(masses(ispecies)**(5.0/2.0))
          end do

          L = 1
          do itheta=1,Ntheta
	     indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

             this_flow(ispecies,itheta) = dot_product(xWeights, flowIntegralWeights * solnArray(indices+1)) &
                  * flowFactors

             this_momentumFluxBeforeThetaIntegral(ispecies,itheta) = ((16d+0)/15) * momentumFluxFactors &
                  * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices+1))
          end do

          L = 2
          do itheta=1,Ntheta
	     indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

             this_particleFluxBeforeThetaIntegral(ispecies,itheta) = this_particleFluxBeforeThetaIntegral(ispecies,itheta) &
                  + (four/15) * particleFluxFactors &
                  * dot_product(xWeights, particleFluxIntegralWeights * solnArray(indices+1))

             this_heatFluxBeforeThetaIntegral(ispecies,itheta) = this_heatFluxBeforeThetaIntegral(ispecies,itheta) &
                  + (four/15) * heatFluxFactors &
                  * dot_product(xWeights, heatFluxIntegralWeights * solnArray(indices+1))

             this_pPerpTermInVpBeforePsiDerivative(ispecies,itheta) &
                  = this_pPerpTermInVpBeforePsiDerivative(ispecies,itheta) &
                  - dot_product(xWeights, pressureIntegralWeights*solnArray(indices+1)) &
                  * 2*pi*(four/15)*THats(ispecies,ipsi)**(5.0/2.0)/(masses(ispecies)**(5.0/2.0))
          end do

          L = 3
          do itheta=1,Ntheta
	     indices = [(getIndex(ispecies,ix,L,itheta,ipsi), ix=min_x_for_L(L),Nx)]

             this_momentumFluxBeforeThetaIntegral(ispecies,itheta) = this_momentumFluxBeforeThetaIntegral(ispecies,itheta) &
                  + (four/35) * momentumFluxFactors &
                  * dot_product(xWeights, momentumFluxIntegralWeights * solnArray(indices+1))
          end do

          do itheta=1,Ntheta
             this_kPar(ispecies,itheta) = FSABHat2(ipsi)/(BHat(itheta,ipsi)*BHat(itheta,ipsi)*dTHatdpsis(ispecies,ipsi)) &
                  * (2*charges(ispecies)*psiAHatArray(ipsi)*BHat(itheta,ipsi)/IHat(ipsi)*this_flow(ispecies,itheta)&
                  + dTHatdpsis(ispecies,ipsi) &
                  + THats(ispecies,ipsi)/nHats(ispecies,ipsi)*dnHatdpsis(ispecies,ipsi)&
                  + 2*charges(ispecies)*omega/Delta*dPhiHatdpsi(ipsi))
          end do

          do itheta=1,Ntheta
             this_pPerpTermInVp(ispecies,itheta) = 0
             this_toroidalFlow(ispecies,itheta) = this_pPerpTermInVp(ispecies,itheta) 
             this_poloidalFlow(ispecies,itheta) = this_pPerpTermInVp(ispecies,itheta)
             
             this_toroidalFlow(ispecies,itheta) &
                  = (Delta/(2*psiAHat))*(masses(ispecies))/(charges(ispecies)*BHat(itheta,ipsi)**2*nHats(ispecies,ipsi)) &
                  *BPHat(itheta,ipsi)**2 * RHat(itheta,ipsi) * this_toroidalFlow(ispecies,itheta)
             this_poloidalFlow(ispecies,itheta) &
                  = (Delta/(2*psiAHat))*(masses(ispecies))/(charges(ispecies)*BHat(itheta,ipsi)**2*nHats(ispecies,ipsi))&
                  *(-1)*BPHat(itheta,ipsi)*IHat(ipsi)* this_poloidalFlow(ispecies,itheta) 

             this_toroidalFlow(ispecies,itheta) = this_toroidalFlow(ispecies,itheta) &
                  + (BTHat(itheta,ipsi)/BHat(itheta,ipsi))*this_flow(ispecies,itheta)
             this_poloidalFlow(ispecies,itheta) = this_poloidalFlow(ispecies,itheta) &
                  + (BPHat(itheta,ipsi)/BHat(itheta,ipsi))*this_flow(ispecies,itheta)

             this_toroidalFlow(ispecies,itheta) = this_toroidalFlow(ispecies,itheta) -omega/(Delta*psiAHat)*dPhiHatdpsi(ipsi) &
                  *this_densityPerturbation(ispecies,itheta)*(BPHat(itheta,ipsi)**2*RHat(itheta,ipsi))/(BHat(itheta,ipsi)**2)
             this_poloidalFlow(ispecies,itheta) = this_poloidalFlow(ispecies,itheta) +omega/(Delta*psiAHat)*dPhiHatdpsi(ipsi) &
                  *this_densityPerturbation(ispecies,itheta)*(BPHat(itheta,ipsi)*IHat(ipsi))/(BHat(itheta,ipsi)**2)

             !since we are not using that variable for anything else
             this_temp = dnHatdpsis(ispecies,ipsi)/nHats(ispecies,ipsi) +dTHatdpsis(ispecies,ipsi)/THats(ispecies,ipsi) &
                  + (2*omega*charges(ispecies)/(Delta*THats(ispecies,ipsi)))*dPhiHatdpsi(ipsi)
             this_toroidalFlow(ispecies,itheta) = this_toroidalFlow(ispecies,itheta) &
                  -THats(ispecies,ipsi)/(2*psiAHat*charges(ispecies)*BHat(itheta,ipsi)**2)*this_temp*BPHat(itheta,ipsi)**2 &
                  *RHat(itheta,ipsi)
             this_poloidalFlow(ispecies,itheta) = this_poloidalFlow(ispecies,itheta) +THats(ispecies,ipsi)&
                  /(2*psiAHat*charges(ispecies)*BHat(itheta,ipsi)**2)*this_temp*BPHat(itheta,ipsi)*IHat(ipsi)
          end do

          
          this_particleFluxBeforeThetaIntegral(ispecies,:) = this_particleFluxBeforeThetaIntegral(ispecies,:) &
               * dBHatdtheta(:,ipsi) / (BHat(:,ipsi) * BHat(:,ipsi) * BHat(:,ipsi))

          this_momentumFluxBeforeThetaIntegral(ispecies,:) = this_momentumFluxBeforeThetaIntegral(ispecies,:) &
               * dBHatdtheta(:,ipsi) / (BHat(:,ipsi) * BHat(:,ipsi) * BHat(:,ipsi))

          this_heatFluxBeforeThetaIntegral(ispecies,:) = this_heatFluxBeforeThetaIntegral(ispecies,:) &
               * dBHatdtheta(:,ipsi) / (BHat(:,ipsi) * BHat(:,ipsi) * BHat(:,ipsi))

          this_particleFlux(ispecies) = dot_product(thetaWeights, this_particleFluxBeforeThetaIntegral(ispecies,:))
          this_momentumFlux(ispecies) = dot_product(thetaWeights, this_momentumFluxBeforeThetaIntegral(ispecies,:))
          this_heatFlux(ispecies) = dot_product(thetaWeights, this_heatFluxBeforeThetaIntegral(ispecies,:))

          this_FSADensityPerturbation(ispecies) = dot_product(thetaWeights, &
               this_densityPerturbation(ispecies,:)/JHat(:,ipsi)) / VPrimeHat(ipsi)

          this_FSAFlow(ispecies) = dot_product(thetaWeights, this_flow(ispecies,:)/JHat(:,ipsi)) / VPrimeHat(ipsi)

          
          this_FSABFlow(ispecies) = dot_product(thetaWeights, this_flow(ispecies,:)*BHat(:,ipsi)/JHat(:,ipsi)) / VPrimeHat(ipsi)

          this_FSAkPar(ispecies) = dot_product(thetaWeights, this_kPar(ispecies,:)/JHat(:,ipsi)) / VPrimeHat(ipsi)

          this_FSAPressurePerturbation(ispecies) = dot_product(thetaWeights, &
               this_pressurePerturbation(ispecies,:)/JHat(:,ipsi)) / VPrimeHat(ipsi)

          this_kParOutboard(ispecies) = this_kPar(ispecies,1)
          this_flowOutboard(ispecies) = this_flow(ispecies,1)
          if (mod(Ntheta,2)==0) then
             this_kParInboard(ispecies) = this_kPar(ispecies,Ntheta/2+1)
             this_flowInboard(ispecies) = this_flow(ispecies,Ntheta/2+1)
          else
             index = (Ntheta+1)/2
             this_kParInboard(ispecies) = (this_kPar(ispecies,index) + this_kPar(ispecies,index+1))*oneHalf
             this_flowInboard(ispecies) = (this_flow(ispecies,index) + this_flow(ispecies,index+1))*oneHalf
          end if

       end do

       LegendresOnXiUniform_m1 = 1
       this_deltaFOutboard = 0
       allocate(solnAtL(Nx))
       do L = 0,(Nxi-1)
          ! Recursively evaluate Legendre polynomials on a uniform grid in xi.
          ! The results will be used to map the distribution function from the modal discretization
          ! to a uniform grid.
          if (L == 0) then
             LegendresOnXiUniform = 1
          else if (L == 1) then
             LegendresOnXiUniform = xiUniform
          else
             LegendresOnXiUniform_m2 = LegendresOnXiUniform_m1
             LegendresOnXiUniform_m1 = LegendresOnXiUniform
             LegendresOnXiUniform = ((2*L-1)*xiUniform * LegendresOnXiUniform_m1 - (L-1)*LegendresOnXiUniform_m2)/L
         end if
         
         do ispecies=1,Nspecies
	    indices = [(getIndex(ispecies,ix,L,thetaIndexForOutboard,ipsi), ix=min_x_for_L(L),Nx)]
            solnAtL(1:min_x_for_L(L)-1) = 0d0
            solnAtL(min_x_for_L(L):Nx) = solnArray(indices+1)

            do ixi = 1,NxiUniform
               this_deltaFOutboard(ispecies,:,ixi) = this_deltaFOutboard(ispecies,:,ixi) + &
                    LegendresOnXiUniform(ixi) * matmul(regridPolynomialToUniformForDiagnostics, solnAtL)
            end do
         end do
      end do

      ! Now build the full f from delta f:
      do ispecies = 1,Nspecies
         speciesFactor = nHats(ispecies,ipsi)/(pi*sqrtpi*((THats(ispecies,ipsi)/masses(ispecies)) ** (1.5d+0)))
         do ixi = 1,NxiUniform
            do ix = 1,NxUniform
               this_fullFOutboard(ispecies,ix,ixi) = this_deltaFOutboard(ispecies,ix,ixi) * Delta &
                    + speciesFactor * exp(-xUniform(ix)*xUniform(ix)) ! This line adds the Maxwellian F_M
            end do
         end do
      end do

      deallocate(solnAtL)
      deallocate(indices)

      !! call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)
      call VecRestoreArrayF90(soln, solnArray, ierr)
      CHKERRQ(ierr)

    !! end if

    call openDebugOutputFile(filename)

    call writeDebugArray(this_densityPerturbation,"densityPerturbation")
    call writeDebugArray(this_flow,"flow")
    call writeDebugArray(this_pPerpTermInVpBeforePsiDerivative,"pPerpTermInVpBeforePsiDerivative")
    call writeDebugArray(this_kPar,"kPar")
    call writeDebugArray(this_pressurePerturbation,"pressurePerturbation")
    call writeDebugArray(this_particleFluxBeforeThetaIntegral,"particleFluxBeforeThetaIntegral")
    call writeDebugArray(this_momentumFluxBeforeThetaIntegral,"momentumFluxBeforeThetaIntegral")
    call writeDebugArray(this_heatFluxBeforeThetaIntegral,"heatFluxBeforeThetaIntegral")
    call writeDebugArray(this_FSADensityPerturbation,"FSADensityPerturbation")
    call writeDebugArray(this_kParOutboard,"kParOutboard")
    call writeDebugArray(this_kParInboard,"kParInboard")
    call writeDebugArray(this_FSAKPar,"FSAKPar")
    call writeDebugArray(this_flowOutboard,"flowOutboard")
    call writeDebugArray(this_flowInboard,"flowInboard")
    call writeDebugArray(this_FSABFlow,"FSABFlow")
    call writeDebugArray(this_FSAFlow,"FSAFlow")
    call writeDebugArray(this_FSAPressurePerturbation,"FSAPressurePerturbation")
    call writeDebugArray(this_particleFlux,"particleFlux")
    call writeDebugArray(this_momentumFlux,"momentumFlux")
    call writeDebugArray(this_heatFlux,"heatFlux")
    if (outputScheme == 2) then
      call writeDebugArray(this_deltaFOutboard,"deltaFOutboard")
      call writeDebugArray(this_fullFOutboard,"fullFOutboard")
    end if

    call closeDebugOutputFile()

  end subroutine calculateLocalMoments

  !function Miller_QQ()
!
 !   implicit none
  !  integer :: i
   ! integer, parameter :: NThetaIntegral = 100
!
 !   PetscScalar :: theta, Miller_QQ
!
 !   Miller_QQ = 0
  !           do i=1,NThetaIntegral
   !             Miller_QQ = Miller_QQ + QQIntegrand(2*pi*i/NThetaIntegral)
    !         end do
     !        Miller_QQ = Miller_kappa / (2*pi*Miller_A) * (Miller_QQ * 2*pi/NThetaIntegral)

  !end function Miller_QQ

end module moments
