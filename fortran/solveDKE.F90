! For compilers that do not include the error function erf(x), the line
! below should be un-commented, and you will need to link to GSL:
!#define USE_GSL_ERF

#include <finclude/petsckspdef.h>
#include <finclude/petscdmdadef.h>

#include "PETScVersions.F90"
  
subroutine solveDKE()

  use geometry
  use globalVariables
  use grids
  use petscdmda
  use profiles
  use sparsify

  implicit none

  PetscErrorCode :: ierr
  Vec :: rhs, rhsLeft, rhsRight, soln, solnLeft, solnRight, solnOnProc0
  Mat :: matrix, preconditionerMatrix
  Mat :: leftMatrix, leftPreconditionerMatrix
  Mat :: rightMatrix, rightPreconditionerMatrix
  PetscViewer MatlabOutput
  PetscScalar, dimension(:), allocatable :: diagonalOfXDot, sourceThetaPart
  PetscScalar, dimension(:), allocatable :: rSingleSpecies
  PetscScalar, dimension(:,:), allocatable :: thetaPartOfStreamingTerm, xPartOfXDot
  integer :: i, j, ix, itheta, ipsi, L, index
  integer :: ispecies, iSpeciesA, iSpeciesB
  integer :: ithetaRow, ithetaCol, scheme, ipsiMinInterior, ipsiMaxInterior, ell
  PetscScalar, dimension(:,:), allocatable :: fToFInterpolationMatrix
  PetscScalar, dimension(:,:), allocatable :: thetaPartMatrix
  PetscScalar, dimension(:), allocatable :: thetaPartOfMirrorTerm
  PetscScalar, dimension(:,:), allocatable :: potentialsToFInterpolationMatrix
  PetscScalar :: VPrime, dtheta, T32
  integer, dimension(:), allocatable :: indices, rowIndices, colIndices
  PetscScalar, dimension(:,:), allocatable :: M11, M21, M32, LaplacianTimesX2WithoutL, nuDHat
  PetscScalar, dimension(:,:,:,:), allocatable :: CECD
  PetscScalar, dimension(:), allocatable :: erfs, xb, expxb2, Psi_Chandra
  PetscScalar, dimension(:,:), allocatable :: KWithoutThetaPart, M22, M33, M12, M13
  PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32, ddpsiToUse
  PetscScalar, dimension(:,:,:), allocatable :: M22BackslashM21s, M33BackslashM32s
  integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
  integer :: LAPACKInfo
  PetscScalar :: xDotFactor, LFactor, temp, temp1, temp2, speciesFactor, speciesFactor2
  PetscScalar, dimension(:), allocatable :: thetaPartOfPsiDot
  PetscScalar, dimension(:,:), allocatable :: localddpsiToUse, everythingButLInPsiDot, ddpsiForKTheta
  PetscScalar :: signOfPsiDot, xPartOfSource
  integer :: sourceToUse, rowIndex, colIndex
  logical :: upwinding, includeddpsiTermThisTime
  integer :: constraintToAdd
  PetscScalar, dimension(:), allocatable :: dnHatdpsi, detaHatdpsi
  PetscScalar, dimension(:), allocatable :: xAndThetaPartOfConstraint
  PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermDiagonal1
  PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermDiagonal2
  PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermDiagonal3
  PetscScalar, dimension(:,:), allocatable :: spatialPartOfStreamingTermOffDiagonal
  PetscScalar, dimension(:), allocatable :: densityFactors, densityIntegralWeights
  PetscScalar, dimension(:), allocatable :: flowFactors, flowIntegralWeights
  PetscScalar, dimension(:), allocatable :: pressureFactors, pressureIntegralWeights
  PetscScalar, dimension(:), allocatable :: particleFluxFactors, particleFluxIntegralWeights
  PetscScalar, dimension(:), allocatable :: momentumFluxFactors, momentumFluxIntegralWeights
  PetscScalar, dimension(:), allocatable :: heatFluxFactors, heatFluxIntegralWeights
  !    PetscScalar, dimension(:), allocatable :: pPerpTermInKThetaFactors
  PetscLogDouble :: time1, time2, startTime
  KSP :: KSPInstance, KSPBoundary
  PC :: preconditionerContext, PCBoundary
  KSPConvergedReason :: reason
  PetscScalar, pointer :: solnArray(:)
  integer :: ipsiMinForThisTheta, ipsiMaxForThisTheta, NpsiToUse
  VecScatter :: VecScatterContext
  integer :: ixi, newIndex, oldIndex
  PetscScalar, dimension(:,:), allocatable :: tempMatrix, tempMatrix2
  logical :: makeLocalApproximationOriginal
  integer :: whichMatrix, whichMatrixMin, rowIndexArray(1), tempInt1, tempInt2, keepXCoupling
  PetscScalar :: singleValueArray(1), stuffToAdd, sqrtMass
  Mat :: permutationMatrix, tempMat
  Vec :: tempVec
  double precision :: myMatInfo(MAT_INFO_SIZE)
  integer :: NNZMain, NNZPreconditioner
  !    PetscScalar :: scaleFactor, maxXForDistribution
  PetscScalar, dimension(:,:), allocatable :: sqrtTHats

  ! *******************************************************************************
  ! *******************************************************************************
  !
  ! Start solveDKE()
  !
  ! *******************************************************************************
  ! *******************************************************************************

  call PetscTime(time1, ierr)
  startTime = time1

  if ((.not. isAParallelDirectSolverInstalled) .and. (numProcsInSubComm > 1)) then
     if (masterProcInSubComm) then
        print *,"Error! Neither mumps nor superlu_dist appears to be installed,"
        print *," yet you have asked for a matrix to be distributed across processsors."
     end if
     stop
  end if

  upwinding = .false.
  call create_grids(upwinding)

   ! *******************************************************************************
   ! *******************************************************************************
   !
   ! Evaluate input physical quantities on the (psi, theta) grids:
   !
   ! *******************************************************************************
   ! *******************************************************************************

   ! Initialize the magnetic geometry
   allocate(BHat(Ntheta,Npsi))
   allocate(dBHatdpsi(Ntheta,Npsi))
   allocate(dBHatdtheta(Ntheta,Npsi))
   allocate(JHat(Ntheta,Npsi))
   allocate(IHat(Npsi))
   allocate(dIHatdpsi(Npsi))

   ! This subroutine to initialize the magnetic geometry is in geometry.F90.
   ! It must fill the following arrays:
   ! BHat(Npsi,Ntheta)
   ! dBHatdpsi(Npsi,Ntheta)
   ! dBHatdtheta(Npsi,Ntheta)
   ! JHat(Npsi,Ntheta)
   ! IHat(Npsi)
   ! dIHatdpsi(Npsi)
   call computeMagneticQuantitiesOnGrids()

   ! Initialize some arrays that can be calculated from the magnetic geometry.
   allocate(VPrimeHat(Npsi))
   allocate(FSABHat2(Npsi))
   allocate(typicalB(Npsi))

   do i=1,Npsi
      VPrimeHat(i) = dot_product(thetaWeights, 1/JHat(:,i))
      FSABHat2(i) = dot_product(thetaWeights, BHat(:,i) * BHat(:,i) / JHat(:,i)) / VPrimeHat(i)
      typicalB(i) = sqrt(FSABHat2(i))
   end do

   ! Initialize the radial physics profiles.
   allocate(PhiHat(Npsi))
   allocate(dPhiHatdpsi(Npsi))
   allocate(THats(numSpecies,Npsi))
   allocate(dTHatdpsis(numSpecies,Npsi))
   allocate(etaHats(numSpecies,Npsi))
   allocate(detaHatdpsis(numSpecies,Npsi))
   allocate(nHats(numSpecies,Npsi))
   allocate(dnHatdpsis(numSpecies,Npsi))

   ! This subroutine to initialize the radial profiles is in profiles.F90.
   ! It must fill the following arrays:
   ! PhiHat(Npsi)
   ! dPhiHatdpsi(Npsi)
   ! THats(numSpecies,Npsi)
   ! dTHatdpsis(numSpecies,Npsi)
   ! nHats(numSpecies,Npsi)
   ! dNHatdpsis(numSpecies,NPsi)
   ! etaHats(numSpecies,Npsi)
   ! detaHatdpsis(numSpecies,Npsi)
   call initializeProfiles()

   ! Initialize some arrays that can be computed from the radial physics profiles.
   allocate(nuPrimeProfile(numSpecies,Npsi))
   allocate(nuStarProfile(numSpecies,Npsi))
   allocate(deltaN(numSpecies,Npsi))
   allocate(deltaT(numSpecies,Npsi))
   allocate(deltaEta(numSpecies,Npsi))
   allocate(U(numSpecies,Npsi))
   allocate(r(numSpecies,Npsi))
   allocate(rSingleSpecies(Npsi))

   do ispecies = 1,numSpecies

      do i=1,Npsi
         deltaT(ispecies,i) = abs(delta*sqrt(masses(ispecies))*IHat(i)/(psiAHat*typicalB(i) &
              *charges(ispecies)*sqrt(THats(ispecies,i)))*dTHatdpsis(ispecies,i))
         deltaN(ispecies,i) = abs(delta*sqrt(masses(ispecies)*THats(ispecies,i))*IHat(i) &
              / (psiAHat*charges(ispecies)*typicalB(i)*nHats(ispecies,i)) * dnHatdpsis(ispecies,i))
         deltaEta(ispecies,i) = abs(delta*sqrt(masses(ispecies)*THats(ispecies,i))*IHat(i) &
              / (psiAHat*charges(ispecies)*typicalB(i)*etaHats(ispecies,i)) * detaHatdpsis(ispecies,i))
      end do

      do i=1,Npsi
         nuPrimeProfile(ispecies,i) = nu_r * Miller_q * nHats(ispecies,i) / (THats(ispecies,i)*THats(ispecies,i))
         nuStarProfile(ispecies,i) = nuPrimeProfile(ispecies,i) / (epsil*sqrt(epsil))
      end do

      U(ispecies,:) = omega*IHat*dPhiHatdpsi/psiAHat*sqrt(masses(ispecies)/(FSABHat2*THats(ispecies,:)))

      ! Next, compute r.
      ! Store dr/dpsi in the variable r, since LAPACK will over-write dr/dpsi with r in a few lines:
      rSingleSpecies = psiAHat / delta * sqrt(FSABHat2 / THats(ispecies,:)) / IHat

      ! Re-create ddpsi_accurate, since it is over-written in the loop.
      ! centered finite differences, no upwinding, 5-point stencil
      scheme = 12
      call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsi_accurate, d2dpsi2)
      ! Change first row of ddpsi matrix so matrix is nonsingular:
      ddpsi_accurate(1,:) = 0
      ddpsi_accurate(1,1) = 1

      allocate(IPIV(Npsi))
      ! The command below overwrites both ddpsi_accurate and r:
#if defined(PETSC_USE_REAL_SINGLE)
      call SGESV(Npsi, 1, ddpsi_accurate, Npsi, IPIV, rSingleSpecies, Npsi, LAPACKInfo)
#else
      call DGESV(Npsi, 1, ddpsi_accurate, Npsi, IPIV, rSingleSpecies, Npsi, LAPACKInfo)
#endif
      if (LAPACKInfo /= 0) then
         print *,"LAPACK error 2!!  Info = ",LAPACKInfo
         stop
       end if
      deallocate(IPIV)
      ! Finally, shift r so its value is 0 at psiMid:
      if (mod(Npsi,2)==1) then
         rSingleSpecies = rSingleSpecies - rSingleSpecies((Npsi+1)/2)
      else
         rSingleSpecies = rSingleSpecies - (rSingleSpecies(Npsi/2) + rSingleSpecies(Npsi/2+1))/2
      end if

      r(ispecies,:) = rSingleSpecies
   end do

   deallocate(rSingleSpecies)
   allocate(sqrtTHats(numSpecies,Npsi))
   sqrtTHats = sqrt(THats)


   ! *********************************************************
   ! *********************************************************
   !
   ! Now build the main matrix, as well as local matrices for 
   ! the left and right boundaries.
   !
   ! *********************************************************
   ! *********************************************************

   ! *********************************************************
   ! Allocate matrices:
   ! *********************************************************

   allocate(xb(Nx))
   allocate(expxb2(Nx))
   allocate(erfs(Nx))
   allocate(Psi_Chandra(Nx))
   allocate(nuDHat(numSpecies, Nx))
   allocate(fToFInterpolationMatrix(Nx,Nx))
   allocate(potentialsToFInterpolationMatrix(Nx, NxPotentials))
   allocate(CECD(numSpecies, numSpecies, Nx, Nx))

   allocate(M21(NxPotentials, Nx))
   allocate(M32(NxPotentials, NxPotentials))
   allocate(M22BackslashM21(NxPotentials, Nx))
   allocate(M33BackslashM32(NxPotentials, NxPotentials))
   allocate(M22BackslashM21s(NL,NxPotentials, Nx))
   allocate(M33BackslashM32s(NL,NxPotentials, NxPotentials))
   allocate(LaplacianTimesX2WithoutL(NxPotentials, NxPotentials))

   allocate(M12(Nx,NxPotentials))
   allocate(M13(Nx,NxPotentials))
   allocate(M22(NxPotentials,NxPotentials))
   allocate(M33(NxPotentials,NxPotentials))

   allocate(M11(Nx,Nx))
   allocate(KWithoutThetaPart(Nx,Nx))
   allocate(IPIV(NxPotentials))

   !predictedNNZForEachRowOfTotalMatrix = 5*(Nx + 5*3 + 5*5 + 3*Nx + 5 +2)
   !predictedNNZForEachRowOfTotalMatrixBoundary = (2*Ntheta + 2 + Nx)

   makeLocalApproximationOriginal = makeLocalApproximation

   if (useIterativeSolver) then
      whichMatrixMin = 0
   else
      whichMatrixMin = 1
   end if
   do whichMatrix = whichMatrixMin,1
      ! When whichMatrix = 0, build the preconditioner.
      ! When whichMatrix = 1, build the final matrix.

      makeLocalApproximation = makeLocalApproximationOriginal
      if (whichMatrix==0 .and. preconditioner_psi==4) then
         makeLocalApproximation = .true.
      end if

      ! In  preallocateMatrix, the last parameter is 0 for the global matrix, or 1 for the local matrices.
      call preallocateMatrix(matrix, 0)
      if (procThatHandlesLeftBoundary) then
         call preallocateMatrix(leftMatrix, 1)
      end if
      if (procThatHandlesRightBoundary) then
         call preallocateMatrix(rightMatrix, 1)
      end if

      ! Sometimes PETSc complains if any of the diagonal elements are not set.
      ! Therefore, set the entire diagonal to 0 to be safe.
      if (masterProcInSubComm) then
         do i=1,matrixSize
            call MatSetValue(matrix, i-1, i-1, zero, ADD_VALUES, ierr)
         end do
      end if
      if (procThatHandlesLeftBoundary) then
         do i=1,localMatrixSize
            call MatSetValue(leftMatrix, i-1, i-1, zero, ADD_VALUES, ierr)
         end do
      end if
      if (procThatHandlesRightBoundary) then
         do i=1,localMatrixSize
            call MatSetValue(rightMatrix, i-1, i-1, zero, ADD_VALUES, ierr)
         end do
      end if


      ! *********************************************************
      ! Add the streaming and mirror terms which persist even 
      ! in the (delta,omega)=(0,0) limit:
      ! *********************************************************

      if (whichMatrix==1 .or. preconditioner_theta==0) then
         ddthetaToUse = ddtheta
      else
         ddthetaToUse = ddtheta_preconditioner
      end if

      allocate(thetaPartOfMirrorTerm(Ntheta))
      allocate(thetaPartMatrix(Ntheta,Ntheta))
      allocate(thetaPartOfStreamingTerm(Ntheta,Ntheta))

      allocate(rowIndices(Ntheta))
      allocate(colIndices(Ntheta))
      do ispecies = 1,numSpecies
         do ipsi = ipsiMin, ipsiMax
            thetaPartOfMirrorTerm = -oneHalf * sqrtTHats(ispecies,ipsi) &
                 * JHat(:,ipsi) * dBHatdtheta(:,ipsi) &
                 / (BHat(:,ipsi) * BHat(:,ipsi))
            do itheta=1,Ntheta
               thetaPartOfStreamingTerm(itheta,:) = JHat(itheta,ipsi)*sqrtTHats(ispecies,ipsi) &
                    / BHat(itheta,ipsi)*ddthetaToUse(itheta,:)
            end do
            do ix=1,Nx
               do L=0,Nxi-1
                  rowIndices = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + [(i, i=1,Ntheta)] - 1
                  if (L < Nxi-1) then
                     ! Super-diagonal in L:
                     colIndices = rowIndices + Ntheta

                     ! Add streaming term:
                     thetaPartMatrix = x(ix)*(L+1)/(two*L+3)*thetaPartOfStreamingTerm

                     ! Add mirror term:
                     ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                     do itheta=1,Ntheta
                        thetaPartMatrix(itheta,itheta) = x(ix)*(L+1)*(L+2)/(two*L+3)*thetaPartOfMirrorTerm(itheta)
                     end do

                     ! Petsc uses a transposed format relative to Fortran:
                     thetaPartMatrix = transpose(thetaPartMatrix)
                     ! Put values in matrix, 
                     if (ipsi==1) then
                        call MatSetValuesSparse(leftMatrix, Ntheta, rowIndices, Ntheta, colIndices, &
                             thetaPartMatrix, ADD_VALUES, ierr)
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                / (psiAHat*charges(ispecies))
                           if (signOfPsiDot < -thresh) then
                              call MatSetValuesSparse(matrix, 1, rowIndices(itheta), Ntheta, colIndices, &
                                   thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     elseif (ipsi==Npsi) then
                        call MatSetValuesSparse(rightMatrix, Ntheta, rowIndices, Ntheta, colIndices, &
                             thetaPartMatrix, ADD_VALUES, ierr)
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                / (psiAHat*charges(ispecies))
                           if (signOfPsiDot > thresh) then
                              rowIndexArray = (ipsi-1)*localMatrixSize+rowIndices(itheta)
                              call MatSetValuesSparse(matrix, 1, rowIndexArray, &
                                   Ntheta, (ipsi-1)*localMatrixSize+colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     else
                        call MatSetValuesSparse(matrix, Ntheta, (ipsi-1)*localMatrixSize+rowIndices, Ntheta, &
                             (ipsi-1)*localMatrixSize+colIndices, &
                             thetaPartMatrix, ADD_VALUES, ierr)
                     end if
                  end if
                  if (L>0) then
                     ! Sub-diagonal in L:
                     colIndices = rowIndices - Ntheta

                     ! Streaming term
                     thetaPartMatrix = x(ix)*L/(two*L-1)*thetaPartOfStreamingTerm

                     ! Mirror term
                     ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                     do itheta=1,Ntheta
                        thetaPartMatrix(itheta,itheta) = -x(ix)*L*(L-1)/(two*L-1)*thetaPartOfMirrorTerm(itheta)
                     end do

                     thetaPartMatrix = transpose(thetaPartMatrix)

                     ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
                     if (ipsi==1) then
                        call MatSetValuesSparse(leftMatrix, Ntheta, rowIndices, Ntheta, colIndices, &
                             thetaPartMatrix, ADD_VALUES, ierr)
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi)&
                                / (psiAHat*charges(ispecies))
                           if (signOfPsiDot < -thresh) then
                              call MatSetValuesSparse(matrix, 1, rowIndices(itheta), Ntheta, colIndices, &
                                   thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     elseif (ipsi==Npsi) then
                        call MatSetValuesSparse(rightMatrix, Ntheta, rowIndices, Ntheta, colIndices, &
                             thetaPartMatrix, ADD_VALUES, ierr)
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                / (psiAHat*charges(ispecies))
                           if (signOfPsiDot > thresh) then
                              rowIndexArray = (ipsi-1)*localMatrixSize+rowIndices(itheta)
                              call MatSetValuesSparse(matrix, 1, rowIndexArray, &
                                   Ntheta, (ipsi-1)*localMatrixSize+colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     else
                        call MatSetValuesSparse(matrix, Ntheta,  (ipsi-1)*localMatrixSize+rowIndices, &
                             Ntheta, (ipsi-1)*localMatrixSize+colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                     end if
                  end if
               end do
            end do
         end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(thetaPartOfStreamingTerm)

      ! *********************************************************
      ! Add the streaming and mirror terms which vanish
      ! in the (delta,omega)=(0,0) limit:
      ! *********************************************************
      if (.not. makeLocalApproximation) then

         allocate(spatialPartOfStreamingTermDiagonal1(Ntheta,Ntheta))
         allocate(spatialPartOfStreamingTermDiagonal2(Ntheta,Ntheta))
         allocate(spatialPartOfStreamingTermDiagonal3(Ntheta,Ntheta))
         allocate(spatialPartOfStreamingTermOffDiagonal(Ntheta,Ntheta))

         allocate(rowIndices(Ntheta))
         allocate(colIndices(Ntheta))

         do ispecies = 1,numSpecies
            sqrtMass = sqrt(masses(ispecies))
            do ipsi = ipsiMin, ipsiMax
               do itheta=1,Ntheta
                  spatialPartOfStreamingTermDiagonal1(itheta,:) = &
                       sqrtMass * omega*JHat(itheta,ipsi)*IHat(ipsi)*dPhiHatdpsi(ipsi) &
                       / (psiAHat*BHat(itheta,ipsi)*BHat(itheta,ipsi)) &
                       * ddthetaToUse(itheta,:)

                  spatialPartOfStreamingTermDiagonal2(itheta,:) = &
                       sqrtMass/charges(ispecies)*delta*THats(ispecies,ipsi)*JHat(itheta,ipsi) &
                       /(psiAHat*BHat(itheta,ipsi)*BHat(itheta,ipsi)) &
                       * IHat(ipsi)*dBHatdpsi(itheta,ipsi)/BHat(itheta,ipsi) &
                       * ddthetaToUse(itheta,:)

                  spatialPartOfStreamingTermDiagonal3(itheta,:) = &
                       sqrtMass/charges(ispecies)*delta*THats(ispecies,ipsi)*JHat(itheta,ipsi) &
                       /(psiAHat*BHat(itheta,ipsi)*BHat(itheta,ipsi)) &
                       * dIHatdpsi(ipsi) &
                       * ddthetaToUse(itheta,:)

                  spatialPartOfStreamingTermOffDiagonal(itheta,:) = &
                       sqrtMass/charges(ispecies)*delta*THats(ispecies,ipsi)*JHat(itheta,ipsi) &
                       / (BHat(itheta,ipsi)*BHat(itheta,ipsi)*psiAHat) &
                       * (IHat(ipsi)/(two*BHat(itheta,ipsi))*dBHatdpsi(itheta,ipsi) - dIHatdpsi(ipsi)) &
                       * ddthetaToUse(itheta,:)
               end do
               do ix=1,Nx
                  thetaPartOfMirrorTerm = sqrt(masses(ispecies))*(omega*dPhiHatdpsi(ipsi)*IHat(ipsi) &
                       + delta*x2(ix)*THats(ispecies,ipsi)/charges(ispecies)*dIHatdpsi(ipsi)) &
                       * JHat(:,ipsi)*dBHatdtheta(:,ipsi) / (two*psiAHat*(BHat(:,ipsi) ** 3))

                  do L=0,Nxi-1
                     rowIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                          + (ix-1)*Nxi*Ntheta  + L*Ntheta + [(i, i=1,Ntheta)] - 1

                     ! Term that is diagonal in L:
                     colIndices = rowIndices

                     ! Add streaming term:
                     thetaPartMatrix = spatialPartOfStreamingTermDiagonal1 &
                          + x2(ix) * ((3*L*L+3*L-2)*spatialPartOfStreamingTermDiagonal2 &
                          - (2*L*L+2*L-1)*spatialPartOfStreamingTermDiagonal3) / ((two*L+3)*(two*L-1))

                     ! Add mirror term:
                     ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                     do itheta=1,Ntheta
                        thetaPartMatrix(itheta,itheta) = L*(L+1)/((two*L-1)*(two*L+3))&
                             * thetaPartOfMirrorTerm(itheta)
                     end do

                     thetaPartMatrix = transpose(thetaPartMatrix)

                     ! Put values in matrix
                     if (ipsi>1 .and. ipsi<Npsi) then
                        call MatSetValuesSparse(matrix, Ntheta, rowIndices, &
                             Ntheta, colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                     else
                        do itheta=1,Ntheta
                           signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                / (psiAHat*charges(ispecies))
                           if ((ipsi==1 .and. signOfPsiDot<-thresh) .or. (ipsi==Npsi .and. signOfPsiDot>thresh)) then
                              call MatSetValuesSparse(matrix, 1, rowIndices(itheta), &
                                   Ntheta, colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                           end if
                        end do
                     end if

                     ! End of term that is diagonal in L.

                     if ((L < Nxi-2) .and. (whichMatrix==1 .or. preconditioner_xi==0)) then
                        ! Super-super-diagonal in L:
                        colIndices = rowIndices + Ntheta*2

                        ! Add streaming term:
                        thetaPartMatrix = x2(ix)*(L+2)*(L+1)/((two*L+5)*(two*L+3)) &
                             * spatialPartOfStreamingTermOffDiagonal

                        ! Add mirror term:
                        ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                        do itheta=1,Ntheta
                           thetaPartMatrix(itheta,itheta) = (L+3)*(L+2)*(L+1)/((two*L+5)*(two*L+3)) &
                                * thetaPartOfMirrorTerm(itheta)
                        end do

                        thetaPartMatrix = transpose(thetaPartMatrix)

                        ! Put values in matrix
                        if (ipsi>1 .and. ipsi<Npsi) then
                           call MatSetValuesSparse(matrix, Ntheta, rowIndices, &
                                Ntheta, colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                        else
                           do itheta=1,Ntheta
                              signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                   / (psiAHat*charges(ispecies))
                              if ((ipsi==1 .and. signOfPsiDot<-thresh) .or. (ipsi==Npsi .and. signOfPsiDot>thresh)) then
                                 call MatSetValuesSparse(matrix, 1, rowIndices(itheta), &
                                      Ntheta, colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                              end if
                           end do
                        end if
 !!$                   ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
 !!$                   call MatSetValuesSparse(preconditionerMatrix, Ntheta, rowIndices, &
 !!$                        Ntheta, colIndices, transpose(thetaPartMatrix), ADD_VALUES, ierr)
                     end if

                     if ((L>1) .and. (whichMatrix==1 .or. preconditioner_xi==0)) then
                        ! Sub-sub-diagonal in L:
                        colIndices = rowIndices - Ntheta*2

                        ! Streaming term
                        thetaPartMatrix = x2(ix)*L*(L-1)/((two*L-3)*(two*L-1)) &
                             * spatialPartOfStreamingTermOffDiagonal

                        ! Mirror term
                        ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                        do itheta=1,Ntheta
                           thetaPartMatrix(itheta,itheta) = -L*(L-1)*(L-2)/((two*L-3)*(two*L-1)) &
                                * thetaPartOfMirrorTerm(itheta)
                        end do

                        thetaPartMatrix = transpose(thetaPartMatrix)

                        ! Put values in matrix
                        if (ipsi>1 .and. ipsi<Npsi) then
                           call MatSetValuesSparse(matrix, Ntheta, rowIndices, &
                                Ntheta, colIndices, thetaPartMatrix, ADD_VALUES, ierr)
                        else
                           do itheta=1,Ntheta
                              signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                                   / (psiAHat*charges(ispecies))
                              if ((ipsi==1 .and. signOfPsiDot<-thresh) .or. (ipsi==Npsi .and. signOfPsiDot>thresh)) then
                                 call MatSetValuesSparse(matrix, 1, rowIndices(itheta), &
                                      Ntheta, colIndices, thetaPartMatrix(:,itheta), ADD_VALUES, ierr)
                              end if
                           end do
                        end if
 !!$                   ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
 !!$                   call MatSetValuesSparse(preconditionerMatrix, Ntheta, rowIndices, &
 !!$                        Ntheta, colIndices, transpose(thetaPartMatrix), ADD_VALUES, ierr)
                     end if
                  end do
               end do
            end do
         end do
         deallocate(rowIndices)
         deallocate(colIndices)
         deallocate(spatialPartOfStreamingTermDiagonal1)
         deallocate(spatialPartOfStreamingTermDiagonal2)
         deallocate(spatialPartOfStreamingTermDiagonal3)
         deallocate(spatialPartOfStreamingTermOffDiagonal)
      end if
      deallocate(thetaPartOfMirrorTerm)
      deallocate(thetaPartMatrix)

      ! *********************************************************
      ! Add the collisionless d/dx term:
      ! *********************************************************

      keepXCoupling = 1
      if (whichMatrix==0 .and. preconditioner_x==1) then
         keepXCoupling = 0
      end if
      ! When keepXCoupling==1, keep the full x coupling.
      ! When keepXCoupling==0, drop everything off-diagonal in x.

      if (.not. makeLocalApproximation) then

         allocate(rowIndices(Nx))
         allocate(colIndices(Nx))
         allocate(xPartOfXDot(Nx,Nx))
         allocate(diagonalOfXDot(Nx))
         do ipsi=ipsiMin,ipsiMax
            do L=0,(Nxi-1)
               if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                  ddxToUse = ddxPreconditioner
               else
                  ddxToUse = ddx
               end if
               do ispecies = 1,numSpecies
                  do ix=1,Nx
                     xPartOfXDot(ix,:) = x(ix)*(delta*dTHatdpsis(ispecies,ipsi)*x2(ix)/(two*charges(ispecies))&
                          + omega*dPhiHatdpsi(ipsi)) * sqrt(masses(ispecies)) * ddxToUse(ix,:)
                  end do
                  xPartOfXDot = transpose(xPartOfXDot)  ! PETSc uses the opposite convention of Fortran
                  do itheta=1,Ntheta
                     signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                          / (psiAHat*charges(ispecies))
                     if (((ipsi > 1) .and. (ipsi < Npsi)) .or. ((ipsi == 1) .and. (signOfPsiDot < -thresh)) &
                          .or. ((ipsi==Npsi) .and. (signOfPsiDot > thresh))) then
                        ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                        ! so impose the kinetic equation here.

                        xDotFactor = JHat(itheta,ipsi) * IHat(ipsi) * dBHatdtheta(itheta,ipsi) &
                             / (two*psiAHat*(BHat(itheta,ipsi) ** 3))
                        rowIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                             + [(ix-1,ix=1,Nx)]*Ntheta*Nxi + L*Ntheta + itheta - 1

                        ! Term that is diagonal in L:
                        colIndices = rowIndices
                        LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))*xDotFactor
                        call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                             LFactor*xPartOfXDot, ADD_VALUES, ierr)

                        if (whichMatrix==1 .or. preconditioner_xi==0) then
                           ! Term that is super-super-diagonal in L:
                           if (L<(Nxi-2)) then
                              colIndices = rowIndices + 2*Ntheta
                              LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3))*xDotFactor
                              call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                                   LFactor*xPartOfXDot, ADD_VALUES, ierr)
                           end if

                           ! Term that is sub-sub-diagonal in L:
                           if (L>1) then
                              colIndices = rowIndices - 2*Ntheta
                              LFactor = L*(L-1)/((two*L-3)*(2*L-1))*xDotFactor
                              call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                                   LFactor*xPartOfXDot, ADD_VALUES, ierr)
                           end if
                        end if

                     end if
                  end do
               end do
            end do
         end do
         deallocate(rowIndices)
         deallocate(colIndices)
         deallocate(xPartOfXDot)
         deallocate(diagonalOfXDot)
      end if

      ! *********************************************************
      ! Add the collisionless d/dpsi term:
      ! *********************************************************

      includeddpsiTermThisTime = includeddpsiTerm
      if (makeLocalApproximation) then
         includeddpsiTermThisTime = .false.
      end if
      if (whichMatrix==0) then
         if (preconditioner_psi==2 .or. preconditioner_psi == 4) then
            includeddpsiTermThisTime = .false.
         end if
      end if
      if (includeddpsiTermThisTime) then
         allocate(colIndices(Npsi))
         allocate(thetaPartOfPsiDot(Npsi))
         allocate(ddpsiToUse(Npsi,Npsi))

         if (whichMatrix==1) then
            ddpsiToUse = ddpsiLeft
         else
            select case (preconditioner_psi)
            case (0)
               ddpsiToUse = ddpsiLeft
            case (1)
               ddpsiToUse = ddpsiForPreconditioner
            case (2,4)
               print *,"Error! Program should not get here."
               stop
            case (3)
               ! Keep only the diagonal of ddpsi, which is typically nonzero only near the boundaries.
               ddpsiToUse = 0
               do i=1,Npsi
                  ddpsiToUse(i,i) = ddpsiLeft(i,i)
               end do
            case default
               print *,"Error! Invalid setting for preconditioner_psi."
               stop
            end select
         end if

         !allocate(localddpsiToUse(localNpsiInterior,Npsi))
         do ispecies = 1,numSpecies
            do itheta = 1, Ntheta
               thetaPartOfPsiDot = -oneHalf * sqrt(masses(ispecies)) * delta * JHat(itheta,:) &
                    * IHat(:) * THats(ispecies,:) * dBHatdtheta(itheta,:) &
                    / (charges(ispecies) * psiAHat * (BHat(itheta,:) ** 3))
               if (upwinding .and. (maxval(thetaPartOfPsiDot)>0) .and. (minval(thetaPartOfPsiDot)<0)) then
                  print *,"Warning: psiDot at itheta =",itheta,&
                       " changes sign with psi, so upwinding is not well-defined."
               end if
               ipsiMinForThisTheta = ipsiMin
               ipsiMaxForThisTheta = ipsiMax
               if (procThatHandlesLeftBoundary .and. (thetaPartOfPsiDot(1) .ge. 0)) then
                  ipsiMinForThisTheta = 2
               end if
               if (procThatHandlesRightBoundary .and. (thetaPartOfPsiDot(Npsi) .le. 0)) then
                  ipsiMaxForThisTheta = Npsi-1
               end if

               NpsiToUse = ipsiMaxForThisTheta - ipsiMinForThisTheta + 1
               allocate(localddpsiToUse(NpsiToUse,Npsi))
               allocate(rowIndices(NpsiToUse))
               allocate(everythingButLInPsiDot(Npsi,NpsiToUse))
               ! The line below must be modified if we want to use upwinded derivatives:
               localddpsiToUse = ddpsiToUse(ipsiMinForThisTheta:ipsiMaxForThisTheta,:)

               !signOfPsiDot = thetaPartOfPsiDot(1)
               !if (signOfPsiDot > 0) then
               !   localddpsiToUse = localddpsiLeftInterior
               !else
               !   localddpsiToUse = localddpsiRightInterior
               !end if
               do ix = 1, Nx
                  do ipsi=1,NpsiToUse
                     ! Build everythingButLInPsiDot as the transpose of what it would normally be
                     ! because of the difference between Fortran and Petsc convention.
                     everythingButLInPsiDot(:,ipsi) = x2(ix) * thetaPartOfPsiDot(ipsi-1+ipsiMinForThisTheta) &
                          * localddpsiToUse(ipsi,:)
                  end do
                  do L=0,(Nxi-1)
                     !rowIndices = [(ipsi-1, ipsi=ipsiMinInterior,ipsiMaxInterior)]*localMatrixSize &
                     rowIndices = [(ipsi-1, ipsi=ipsiMinForThisTheta,ipsiMaxForThisTheta)]*localMatrixSize &
                          + (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Ntheta*Nxi + L*Ntheta + itheta - 1

                     ! Term that is diagonal in L:
                     ell = L
                     colIndices = [(ipsi-1, ipsi=1,Npsi)]*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                          + (ix-1)*Ntheta*Nxi + ell*Ntheta + itheta - 1
                     LFactor = two*(3*L*L+3*L-2)/((two*L+3)*(2*L-1))
                     call MatSetValuesSparse(matrix, NpsiToUse, rowIndices, Npsi, colIndices, &
                          LFactor*everythingButLInPsiDot, ADD_VALUES, ierr)

                     if (whichMatrix==1 .or. preconditioner_xi==0) then
                        ! Term that is super-super-diagonal in L:
                        if (L<(Nxi-2)) then
                           ell = L+2
                           colIndices = [(ipsi-1, ipsi=1,Npsi)]*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                                + (ix-1)*Ntheta*Nxi + ell*Ntheta + itheta - 1
                           LFactor = (L+1)*(L+2)/((two*L+5)*(2*L+3))
                           call MatSetValuesSparse(matrix, NpsiToUse, rowIndices, Npsi, colIndices, &
                                LFactor*everythingButLInPsiDot, ADD_VALUES, ierr)
                        end if

                        ! Term that is sub-sub-diagonal in L:
                        if (L>1) then
                           ell = L-2
                           colIndices = [(ipsi-1, ipsi=1,Npsi)]*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                                + (ix-1)*Ntheta*Nxi + ell*Ntheta + itheta - 1
                           LFactor = L*(L-1)/((two*L-3)*(2*L-1))
                           call MatSetValuesSparse(matrix, NpsiToUse, rowIndices, Npsi, colIndices, &
                                LFactor*everythingButLInPsiDot, ADD_VALUES, ierr)
                        end if
                     end if

                  end do
               end do
               deallocate(localddpsiToUse)
               deallocate(rowIndices)
               deallocate(everythingButLInPsiDot)
            end do
         end do
         !deallocate(rowIndices)
         deallocate(colIndices)
         deallocate(thetaPartOfPsiDot)
         deallocate(ddpsiToUse)
      end if

      ! *********************************************************
      ! In preparation for adding the collision operator,
      ! create several matrices which will be needed.
      ! *********************************************************

      allocate(rowIndices(Nx))
      allocate(colIndices(Nx))
      allocate(tempMatrix(Nx, NxPotentials))
      allocate(tempMatrix2(NxPotentials, NxPotentials))


      if (whichMatrix==0 .and. preconditioner_x==5) then
         ! For other preconditioner_x values, the simplified x derivative
         ! is implemented later instead of here.
         ddxToUse = ddxPreconditioner
         d2dx2ToUse = d2dx2Preconditioner
         ! Since I set ddxToUse here, the same ddxToUse is used for all values of L.
         ! This means that preconditioner_x_min_L does not behave as you might expect when
         ! preconditioner_x=5.
      else
         ddxToUse = ddx
         d2dx2ToUse = d2dx2
      end if

      ! First assemble rows 2 and 3 of the block linear system, since they
      ! are independent of psi and independent of species.

      M32 = zero
      M21 = 4*pi*regridPolynomialToUniform
      do i=2,NxPotentials-1
         M21(i,:) = M21(i,:)*xPotentials(i)*xPotentials(i)
         M32(i,i) = -2*xPotentials(i)*xPotentials(i)
      end do
      M21(1,:)=zero
      M21(NxPotentials,:)=zero
      M32(1,:)=zero
      M32(NxPotentials,:)=zero
      do i=1,NxPotentials
         LaplacianTimesX2WithoutL(i,:) = xPotentials(i)*xPotentials(i)*d2dx2Potentials(i,:) &
              + 2 * xPotentials(i) * ddxPotentials(i,:)
      end do

      do L=0,(NL-1)
         M22 = LaplacianTimesX2WithoutL
         do i=1,NxPotentials
            M22(i,i) = M22(i,i) - L*(L+1)
         end do

         ! Add Dirichlet or Neumann boundary condition for potentials at x=0:
         if (L==0) then
            M22(1,:)=ddxPotentials(1,:)
         else
            M22(1,:) = 0
            M22(1,1) = 1
         end if
         M33 = M22;

         ! Add Robin boundary condition for potentials at x=xMax:
         M22(NxPotentials,:) = xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
         M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1

         ! Boundary condition for G:
         M33(NxPotentials,:) = xMaxNotTooSmall*xMaxNotTooSmall*d2dx2Potentials(NxPotentials,:) &
              + (2*L+1)*xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
         M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1)

         if (L /= 0) then
            M22(NxPotentials,1)=0
            M33(NxPotentials,1)=0
         end if

         ! Call LAPACK subroutine DGESV to solve a linear system
         ! Note: this subroutine changes M22 and M33!
         M22BackslashM21 = M21  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
         call SGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#else
         call DGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#endif
         if (LAPACKInfo /= 0) then
            print *, "Error in LAPACK call: info = ", LAPACKInfo
            stop
         end if
         M33BackslashM32 = M32  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
         call SGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#else
         call DGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#endif
         if (LAPACKInfo /= 0) then
            print *, "Error in LAPACK call: info = ", LAPACKInfo
            stop
         end if

         M33BackslashM32s(L+1,:,:) = M33BackslashM32
         M22BackslashM21s(L+1,:,:) = M22BackslashM21
      end do

      do ipsi = ipsiMin, ipsiMax

         nuDHat = zero
         CECD = zero
         ! Before adding the collision operator, we must loop over both species
         ! to build several terms in the operator.
         ! row is species a, column is species b
         do iSpeciesA = 1,numSpecies
            do iSpeciesB = 1,numSpecies
               speciesFactor = sqrt(THats(iSpeciesA,ipsi)*masses(iSpeciesB) &
                    / (THats(iSpeciesB,ipsi) * masses(iSpeciesA)))
               xb =  x * speciesFactor
               expxb2 = exp(-xb*xb)
               do ix=1,Nx
                  ! erf is vectorized in gfortran but not pathscale
                  temp1 = xb(ix)
#ifdef USE_GSL_ERF
                  call erf(temp1, temp2)
#else
                  temp2 = erf(temp1)
#endif
                  erfs(ix) = temp2
               end do
               Psi_Chandra = (erfs - 2/sqrtpi * xb * expxb2) / (2*xb*xb)

               T32 = THats(iSpeciesA,ipsi) * sqrt(THats(iSpeciesA,ipsi))

               ! Build the pitch-angle scattering frequency:
               nuDHat(iSpeciesA, :) =  nuDHat(iSpeciesA, :) &
                    + (three*sqrtpi/four) / T32 &
                    * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) &
                    * nHats(iSpeciesB,ipsi)*(erfs - Psi_Chandra)/(x*x*x)

               ! Given a vector of function values on the species-B grid, multiply the vector
               ! by this regridding matrix to obtain its values on the species-A grid:
               if (iSpeciesA /= iSpeciesB) then
                  call polynomialInterpolationMatrix(Nx, Nx, x, xb, expx2, &
                       expxb2, fToFInterpolationMatrix)
               else
                  fToFInterpolationMatrix = zero
                  do i=1,Nx
                     fToFInterpolationMatrix(i, i) = one
                  end do
               end if

               ! Using the resulting interpolation matrix,
               ! add CD (the part of the field term independent of Rosenbluth potentials.
               ! CD is dense in the species indices.

               speciesFactor = 3 * nHats(iSpeciesA,ipsi)  * masses(iSpeciesA)/masses(iSpeciesB) &
                    * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) / T32

               do ix=1,Nx
                  CECD(iSpeciesA, iSpeciesB, ix, :) = CECD(iSpeciesA, iSpeciesB, ix, :) &
                       + speciesFactor * expx2(ix) * fToFInterpolationMatrix(ix, :)
               end do

               ! Done adding CD. Now add energy scattering (CE).
               ! Unlike CD, CE is diagonal in the species index.

               speciesFactor = 3*sqrtpi/four * nHats(iSpeciesB,ipsi)  &
                    * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) / T32

               do ix=1,Nx
                  !Now add the d2dx2 and ddx terms in CE:
                  !CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                  CECD(iSpeciesA, iSpeciesA, ix, :) = CECD(iSpeciesA, iSpeciesA, ix, :) &
                       + speciesFactor * (Psi_Chandra(ix)/x(ix)*d2dx2ToUse(ix,:) &
                       + (-2*THats(iSpeciesA,ipsi)*masses(iSpeciesB)/(THats(iSpeciesB,ipsi)*masses(iSpeciesA)) &
                       * Psi_Chandra(ix)*(1-masses(iSpeciesA)/masses(iSpeciesB)) &
                       + (erfs(ix)-Psi_Chandra(ix))/x2(ix)) * ddxToUse(ix,:))

                  ! Lastly, add the part of CE for which f is not differentiated:
                  ! CE is diagonal in the species indices, so use iSpeciesA for both indices in CECD:
                  CECD(iSpeciesA, iSpeciesA, ix, ix) = CECD(iSpeciesA, iSpeciesA, ix, ix) &
                       + speciesFactor *4/sqrtpi*THats(iSpeciesA,ipsi)/THats(iSpeciesB,ipsi) &
                       *sqrt(THats(iSpeciesA,ipsi)*masses(iSpeciesB)/(THats(iSpeciesB,ipsi)*masses(iSpeciesA))) &
                       * expxb2(ix)

               end do

            end do
         end do


         ! *****************************************************************
         ! Now we are ready to add the collision operator to the main matrix.
         ! *****************************************************************

         do L=0, Nxi-1
            do iSpeciesB = 1,numSpecies
               do iSpeciesA = 1,numSpecies
                  if (iSpeciesA==iSpeciesB .or. whichMatrix==1) then

                     ! Build M11
                     ! Eventually un-remark the next line:
                     M11 = -nu_r * CECD(iSpeciesA, iSpeciesB,:,:)
                     if (iSpeciesA == iSpeciesB) then
                        do i=1,Nx
                           M11(i,i) = M11(i,i) - nu_r * (-oneHalf*nuDHat(iSpeciesA,i)*L*(L+1))
                        end do
                     end if

                     if (L < NL) then
                     !   if (.false.) then
                        ! Add Rosenbluth potential terms.

                        speciesFactor2 = sqrt(THats(iSpeciesA,ipsi)*masses(iSpeciesB) &
                             / (THats(iSpeciesB,ipsi) * masses(iSpeciesA)))

                        ! Build M13:
                        scheme = 2
                        call interpolationMatrix(NxPotentials, Nx, xPotentials, x*speciesFactor2, &
                             potentialsToFInterpolationMatrix, scheme, L)

                        speciesFactor = -nu_r * 3/(2*pi)*nHats(iSpeciesA,ipsi) &
                             * charges(iSpeciesA)*charges(iSpeciesA)*charges(iSpeciesB)*charges(iSpeciesB) &
                             / (THats(iSpeciesA,ipsi) * sqrt(THats(iSpeciesA,ipsi))) &
                             * THats(iSpeciesB,ipsi)*masses(iSpeciesA)/(THats(iSpeciesA,ipsi)*masses(iSpeciesB))

                        tempMatrix = matmul(potentialsToFInterpolationMatrix, d2dx2Potentials)
                        do i=1,Nx
                           M13(i, :) = speciesFactor*expx2(i)*x2(i)*tempMatrix(i,:)
                        end do

                        ! Build M12:
                        scheme = 1
                        call interpolationMatrix(NxPotentials, Nx, xPotentials, x*speciesFactor2, &
                             potentialsToFInterpolationMatrix, scheme, L)

                        temp = 1-masses(iSpeciesA)/masses(iSpeciesB)
                        do i=1,NxPotentials
                           tempMatrix2(i,:) = temp*xPotentials(i)*ddxPotentials(i,:)
                           tempMatrix2(i,i) = tempMatrix2(i,i) + one
                        end do
                        tempMatrix = matmul(potentialsToFInterpolationMatrix, tempMatrix2)
                        do i=1,Nx
                           M12(i,:) = -speciesFactor*expx2(i)*tempMatrix(i,:)
                        end do

                        ! Possibly add Dirichlet boundary condition for potentials at x=0:
                        if (L /= 0) then
                           M12(:,1) = 0
                           M13(:,1) = 0
                        end if

                        !KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                        KWithoutThetaPart = M11 - matmul(M12 - matmul(M13, M33BackslashM32s(L+1,:,:)),&
                             M22BackslashM21s(L+1,:,:))

                     else
                        KWithoutThetaPart = M11;
                     end if

                     if (whichMatrix==0 .and. L >= preconditioner_x_min_L) then
                        ! We're making the preconditioner, so simplify the x part of the matrix if desired.
                        select case (preconditioner_x)
                        case (0)
                           ! Do nothing.
                        case (1)
                           ! Keep only diagonal in x:
                           do i=1,Nx
                              do j=1,Nx
                                 if (i /= j) then
                                    KWithoutThetaPart(i,j) = zero
                                 end if
                              end do
                           end do
                        case (2)
                           ! Keep only upper-triangular part:
                           do i=2,Nx
                              do j=1,(i-1)
                                 KWithoutThetaPart(i,j) = zero
                              end do
                           end do
                        case (3,5)
                           ! Keep only tridiagonal part:
                           do i=1,Nx
                              do j=1,Nx
                                 if (abs(i-j)>1) then
                                    KWithoutThetaPart(i,j) = zero
                                 end if
                              end do
                           end do
                        case (4)
                           ! Keep only the diagonal and super-diagonal:
                           do i=1,Nx
                              do j=1,Nx
                                 if (i /= j .and. j /= (i+1)) then
                                    KWithoutThetaPart(i,j) = zero
                                 end if
                              end do
                           end do
                        case default
                           print *,"Error! Invalid preconditioner_x"
                           stop
                        end select

                     end if

                     ! PETSc and Fortran use row-major vs column-major:
                     KWithoutThetaPart = transpose(KWithoutThetaPart)

                     do itheta=1,Ntheta
                        signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                             / (psiAHat*charges(iSpeciesA))
                        if ((ipsi > 1 .and. ipsi < Npsi) .or. (ipsi == 1 .and. signOfPsiDot < -thresh) &
                             .or. (ipsi==Npsi .and. signOfPsiDot > thresh)) then
                           ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                           ! so impose the kinetic equation here.

                           rowIndices = (ipsi-1)*localMatrixSize + (iSpeciesA-1)*Nx*Nxi*Ntheta &
                                + [(i, i=0,Nx-1)]*Nxi*Ntheta + L*Ntheta + itheta - 1
                           colIndices = (ipsi-1)*localMatrixSize + (iSpeciesB-1)*Nx*Nxi*Ntheta &
                                + [(i, i=0,Nx-1)]*Nxi*Ntheta + L*Ntheta + itheta - 1
                           call MatSetValuesSparse(matrix, Nx, rowIndices, Nx, colIndices, &
                                KWithoutThetaPart, ADD_VALUES, ierr)
                        end if
                     end do

                     if (procThatHandlesLeftBoundary .and. ipsi==1) then
                        do itheta=1,Ntheta
                           rowIndices = (iSpeciesA-1)*Nx*Nxi*Ntheta + [(i, i=0,Nx-1)]*Nxi*Ntheta &
                                + L*Ntheta + itheta - 1
                           colIndices = (iSpeciesB-1)*Nx*Nxi*Ntheta + [(i, i=0,Nx-1)]*Nxi*Ntheta &
                                + L*Ntheta + itheta - 1
                           call MatSetValuesSparse(leftMatrix, Nx, rowIndices, Nx, colIndices, &
                                KWithoutThetaPart, ADD_VALUES, ierr)
                        end do
                     end if

                     if (procThatHandlesRightBoundary .and. ipsi==Npsi) then
                        do itheta=1,Ntheta
                           rowIndices = (iSpeciesA-1)*Nx*Nxi*Ntheta + [(i, i=0,Nx-1)]*Nxi*Ntheta &
                                + L*Ntheta + itheta - 1
                           colIndices = (iSpeciesB-1)*Nx*Nxi*Ntheta + [(i, i=0,Nx-1)]*Nxi*Ntheta &
                                + L*Ntheta + itheta - 1
                           call MatSetValuesSparse(rightMatrix, Nx, rowIndices, Nx, colIndices, &
                                KWithoutThetaPart, ADD_VALUES, ierr)
                        end do
                     end if

                  end if
               end do
            end do
         end do
      end do

      deallocate(rowIndices)
      deallocate(colIndices)
      deallocate(tempMatrix)
      deallocate(tempMatrix2)


      ! *******************************************************************************
      ! *******************************************************************************
      !
      ! Done adding collision operator.
      !
      ! *******************************************************************************
      ! *******************************************************************************

      ! *******************************************************************************
      ! Put a 1 on the matrix diagonal where appropriate to enforce the radial boundary condition
      ! *******************************************************************************

      if (procThatHandlesLeftBoundary) then
         ipsi=1
         do ispecies = 1,numSpecies
            do itheta=1,Ntheta
               signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                    /(psiAHat*charges(ispecies))
               if (signOfPsiDot > -thresh) then
                  do ix=1,Nx
                     do L=0,(Nxi-1)
                        index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta+L*Ntheta+itheta - 1
                        call MatSetValueSparse(matrix, index, index, one, ADD_VALUES, ierr)
                     end do
                  end do
               end if
            end do
         end do
      end if
      if (procThatHandlesRightBoundary) then
         ipsi = Npsi
         do ispecies = 1,numSpecies
            do itheta=1,Ntheta
               signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                    / (psiAHat*charges(ispecies))
               if (signOfPsiDot < thresh) then
                  do ix=1,Nx
                     do L=0,(Nxi-1)
                        index = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                             + (ix-1)*Nxi*Ntheta+L*Ntheta+itheta - 1
                        call MatSetValueSparse(matrix, index, index, one, ADD_VALUES, ierr)
                     end do
                  end do
               end if
            end do
         end do

      end if

      ! *******************************************************************************
      ! Add sources:
      ! *******************************************************************************

      allocate(sourceThetaPart(Ntheta))
      select case (sourcePoloidalVariation)
      case (0)
         sourceThetaPart = 1
      case (1)
         do i=1,Ntheta
            sourceThetaPart(i) = 1 + cos(theta(i))
         end do
      case default
         print *,"Error! Invalid sourcePoloidalVariation."
         stop
      end select

      ! Add particle source:
      ! S = f_M * (x^2 - 5/2)  (Provides particles but no heat or momentum)
      L = 0
      do ix=1,Nx
         xPartOfSource = (x2(ix)-5/two)*exp(-x2(ix))
         do ipsi=ipsiMin,ipsiMax
            do ispecies = 1,numSpecies
               do itheta=1,Ntheta
                  signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                       / (psiAHat*charges(ispecies))
                  if ((ipsi > 1 .and. ipsi < Npsi) .or. (ipsi == 1 .and. signOfPsiDot < -thresh) &
                       .or. (ipsi==Npsi .and. signOfPsiDot > thresh)) then
                     ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                     ! so impose the kinetic equation here.

                     rowIndex = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                          + (ix-1)*Ntheta*Nxi + L*Ntheta + itheta - 1
                     colIndex = Npsi*localMatrixSize + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 1 - 1
                     call MatSetValueSparse(matrix, rowIndex, colIndex, &
                          sourceThetaPart(itheta) * xPartOfSource, ADD_VALUES, ierr)
                  end if
               end do
            end do
         end do
      end do

      ! Add heat source:
      ! S = f_M * (x^2 - 3/2)  (Provides heat but no particles or momentum)
      L = 0
      do ix=1,Nx
         xPartOfSource = (x2(ix)-3/two)*exp(-x2(ix))
         do ipsi=ipsiMin,ipsiMax
            do ispecies = 1,numSpecies
               do itheta=1,Ntheta
                  signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                       / (psiAHat*charges(ispecies))
                  if ((ipsi > 1 .and. ipsi < Npsi) .or. (ipsi == 1 .and. signOfPsiDot < -thresh) &
                       .or. (ipsi==Npsi .and. signOfPsiDot > thresh)) then
                     ! We're either in the interior, or on a boundary point at which trajectories leave the domain,
                     ! so impose the kinetic equation here.

                     rowIndex = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                          + (ix-1)*Ntheta*Nxi + L*Ntheta + itheta - 1
                     colIndex = Npsi*localMatrixSize + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 2 - 1
                     call MatSetValueSparse(matrix, rowIndex, colIndex, &
                          sourceThetaPart(itheta) * xPartOfSource, ADD_VALUES, ierr)
                  end if
               end do
            end do
         end do
      end do

      deallocate(sourceThetaPart)

      ! *******************************************************************************
      ! Add constraints:
      ! *******************************************************************************

      if (procThatHandlesRightBoundary) then
         ! The processor that owns the right-most psi point handles the constraints.

         L=0
         allocate(colIndices(Ntheta))
         allocate(xAndThetaPartOfConstraint(Ntheta))

         ! Enforce <n_1> = 0 at each psi
         do ispecies = 1,numSpecies
            do ipsi = 1, Npsi
               rowIndexArray = Npsi*localMatrixSize + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 1 - 1
               do ix = 1, Nx
                  xAndThetaPartOfConstraint = xWeights(ix)*x2(ix) * thetaWeights / JHat(:,ipsi)
                  colIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                       + (ix-1)*Nxi*Ntheta + L*Ntheta + [(itheta,itheta=1,Ntheta)]-1
                  call MatSetValuesSparse(matrix, 1, rowIndexArray, Ntheta, colIndices, xAndThetaPartOfConstraint,&
                       ADD_VALUES, ierr)
               end do
            end do
         end do

         ! Enforce <p_1> = 0 at each psi
         do ispecies = 1,numSpecies
            do ipsi = 1, Npsi
               rowIndexArray = Npsi*localMatrixSize + (ipsi-1)*numSpecies*2 + (ispecies-1)*2 + 2 - 1
               do ix = 1, Nx
                  xAndThetaPartOfConstraint = xWeights(ix)*x2(ix)*x2(ix) * thetaWeights / JHat(:,ipsi)
                  colIndices = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                       + (ix-1)*Nxi*Ntheta + L*Ntheta + [(itheta,itheta=1,Ntheta)]-1
                  call MatSetValuesSparse(matrix, 1, rowIndexArray, Ntheta, colIndices, xAndThetaPartOfConstraint,&
                       ADD_VALUES, ierr)
               end do
            end do
         end do

         deallocate(colIndices)
         deallocate(xAndThetaPartOfConstraint)
      end if

      ! *******************************************************************************
      ! Done inserting values into the matrices.
      ! Now finalize the matrices:
      ! *******************************************************************************

      call PetscTime(time2, ierr)
      if (masterProcInSubComm) then
         if (whichMatrix==0) then
            print *,"[",myCommunicatorIndex,"] Time to pre-assemble preconditioner matrix: ", time2-time1, " seconds."
         else
            print *,"[",myCommunicatorIndex,"] Time to pre-assemble matrix: ", time2-time1, " seconds."
         end if
      end if
      call PetscTime(time1, ierr)

      call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      if (procThatHandlesLeftBoundary) then
         call MatAssemblyBegin(leftMatrix, MAT_FINAL_ASSEMBLY, ierr)
      end if
      if (procThatHandlesRightBoundary) then
         call MatAssemblyBegin(rightMatrix, MAT_FINAL_ASSEMBLY, ierr)
      end if
      CHKERRQ(ierr)
      call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
      if (procThatHandlesLeftBoundary) then
         call MatAssemblyEnd(leftMatrix, MAT_FINAL_ASSEMBLY, ierr)
      end if
      if (procThatHandlesRightBoundary) then
         call MatAssemblyEnd(rightMatrix, MAT_FINAL_ASSEMBLY, ierr)
      end if
      CHKERRQ(ierr)

      if (whichMatrix==0) then
         preconditionerMatrix = matrix
         if (procThatHandlesLeftBoundary) then
            leftPreconditionerMatrix = leftMatrix
         end if
         if (procThatHandlesRightBoundary) then
            rightPreconditionerMatrix = rightMatrix
         end if
      end if

      call PetscTime(time2, ierr)
      if (masterProcInSubComm) then
         if (whichMatrix==0) then
            print *,"[",myCommunicatorIndex,"] Time to assemble preconditioner matrices: ", time2-time1, " seconds."
         else
            print *,"[",myCommunicatorIndex,"] Time to assemble matrices: ", time2-time1, " seconds."
         end if
      end if
      call PetscTime(time1, ierr)

   end do

   makeLocalApproximation = makeLocalApproximationOriginal

   ! *******************************************************************************
   ! *******************************************************************************
   ! Delete this next section eventually
   ! *******************************************************************************
   ! *******************************************************************************
 !!$        call PetscViewerASCIIOpen(MPIComm, &
 !!$             & "boundaries.m",&
 !!$             & MatlabOutput, ierr)
 !!$        CHKERRQ(ierr)
 !!$        call PetscViewerSetFormat(MatlabOutput, PETSC_VIEWER_ASCII_MATLAB, ierr)
 !!$        CHKERRQ(ierr)
 !!$
 !!$        !call PetscObjectSetName(leftPreconditionerMatrix, "leftPreconditionerMatrix", ierr)
 !!$        !call MatView(leftPreconditionerMatrix, MatlabOutput, ierr)
 !!$        call PetscObjectSetName(leftMatrix, "leftMatrix", ierr)
 !!$        call MatView(leftMatrix, MatlabOutput, ierr)
 !!$
 !!$        !call PetscObjectSetName(rightPreconditionerMatrix, "rightPreconditionerMatrix", ierr)
 !!$        !call MatView(rightPreconditionerMatrix, MatlabOutput, ierr)
 !!$        call PetscObjectSetName(rightMatrix, "rightMatrix", ierr)
 !!$        call MatView(rightMatrix, MatlabOutput, ierr)
 !!$
 !!$        call PetscViewerDestroy(MatlabOutput, ierr)

   ! *******************************************************************************
   ! *******************************************************************************
   !
   ! Create the right-hand side vector
   !
   ! *******************************************************************************
   ! *******************************************************************************

   call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr)
   CHKERRQ(ierr)
   if (procThatHandlesLeftBoundary) then
      ! This process handles the left boundary, so solve for the local solution there.
      call VecCreateSeq(MPI_COMM_SELF, localMatrixSize, rhsLeft, ierr)
   end if
   if (procThatHandlesRightBoundary) then
      ! This process handles the right boundary, so solve for the local solution there.
      call VecCreateSeq(MPI_COMM_SELF, localMatrixSize, rhsRight, ierr)
   end if
   CHKERRQ(ierr)

   do ispecies = 1, numSpecies
      do ipsi = ipsiMin, ipsiMax
         do itheta = 1, Ntheta
            do ix = 1, Nx

               stuffToAdd = masses(ispecies)*masses(ispecies) * nHats(ispecies,ipsi) * IHat(ipsi) &
                    * JHat(itheta,ipsi) * dBHatdtheta(itheta,ipsi) * x2(ix) * expx2(ix) &
                    /(2*pi*sqrtpi*charges(ispecies)*sqrtTHats(ispecies,ipsi)*psiAHat &
                    * (BHat(itheta,ipsi) ** 3)) &
                    * (dnHatdpsis(ispecies,ipsi)/nHats(ispecies, ipsi) &
                    + 2*charges(ispecies)/THats(ispecies,ipsi)*omega/Delta*dPhiHatdpsi(ipsi) &
                    + (x2(ix)-3/two)/THats(ispecies,ipsi)*dTHatdpsis(ispecies,ipsi))

               ! It's okay to assign the rhs even at ipsi=1 and ipsi=Npsi because we will
               ! over-write these values later.

               L = 0
               LFactor = 4/three
               index = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta - 1
               !call VecSetValues(rhs, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               call VecSetValue(rhs, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta &
                    + L*Ntheta + itheta - 1
               if (ipsi==1) then
                  ! This is the left boundary
                  !call VecSetValues(rhsLeft, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                  call VecSetValue(rhsLeft, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               elseif (ipsi==Npsi) then
                  ! This is the right boundary
                  !call VecSetValues(rhsRight, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                  call VecSetValue(rhsRight, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               end if

               L = 2
               LFactor = 2/three
               index = (ipsi-1)*localMatrixSize + (ispecies-1)*Nx*Nxi*Ntheta &
                    + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta - 1
               !call VecSetValues(rhs, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               call VecSetValue(rhs, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta &
                    + L*Ntheta + itheta - 1
               if (ipsi==1) then
                  ! This is the left boundary
                  !call VecSetValues(rhsLeft, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                  call VecSetValue(rhsLeft, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               elseif (ipsi==Npsi) then
                  ! This is the right boundary
                  !call VecSetValues(rhsRight, 1, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
                  call VecSetValue(rhsRight, index, LFactor*stuffToAdd, INSERT_VALUES, ierr)
               end if

            end do
         end do
      end do
   end do

   CHKERRQ(ierr)
   call deallocate_initialization_grid_arrays

   ! *********************************************************************************************
   ! If this process handles the left or right boundary, solve the local kinetic equation there:
   ! *********************************************************************************************
   if (procThatHandlesLeftBoundary) then
      ! This process handles the left boundary, so solve the local kinetic equation there.

      call VecAssemblyBegin(rhsLeft, ierr)
      call VecAssemblyEnd(rhsLeft, ierr)

      call KSPCreate(MPI_COMM_SELF, KSPBoundary, ierr)
      if (useIterativeSolver) then
         ! Use iterative solver
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4:
         call KSPSetOperators(KSPBoundary, leftMatrix, leftPreconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
         ! Syntax for PETSc version 3.5 and later
         call KSPSetOperators(KSPBoundary, leftMatrix, leftPreconditionerMatrix, ierr)
#endif
         call KSPGetPC(KSPBoundary, PCBoundary, ierr)
         call PCSetType(PCBoundary, PCLU, ierr)
         call KSPSetType(KSPBoundary, KSPBCGSL, ierr)
         call KSPSetTolerances(KSPBoundary, solverTolerance, PETSC_DEFAULT_REAL, &
              PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
         call KSPSetFromOptions(KSPBoundary, ierr)
         call KSPMonitorSet(KSPBoundary, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
      else
         ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4:
         call KSPSetOperators(KSPBoundary, leftMatrix, leftMatrix, SAME_PRECONDITIONER, ierr)
#else
         ! Syntax for PETSc version 3.5 and later
         call KSPSetOperators(KSPBoundary, leftMatrix, leftMatrix, ierr)
#endif
         call KSPGetPC(KSPBoundary, PCBoundary, ierr)
         call PCSetType(PCBoundary, PCLU, ierr)
         call KSPSetType(KSPBoundary, KSPPREONLY, ierr)
         call KSPSetFromOptions(KSPBoundary, ierr)
      end if

      call VecDuplicate(rhsLeft, solnLeft, ierr)
      CHKERRQ(ierr)
      select case (leftBoundaryScheme)
      case (0)
         print *,"[",myCommunicatorIndex,"] Setting f_1 at left boundary to 0."
         call VecSet(solnLeft, 0, ierr)
      case (1)
         print *,"[",myCommunicatorIndex,"] Proc",myRank, &
              " is solving local kinetic equation at left boundary ..."
         if (solveSystem) then
            call KSPSolve(KSPBoundary, rhsLeft, solnLeft, ierr)
         end if
         CHKERRQ(ierr)

         call PetscTime(time2, ierr)
         print *,"[",myCommunicatorIndex,"] Done solving for left boundary.  Time to solve: ", time2-time1, " seconds."
         call PetscTime(time1, ierr)

         if (useIterativeSolver) then
            call KSPGetConvergedReason(KSPBoundary, reason, ierr)
            if (reason>0) then
               print *,"[",myCommunicatorIndex,"] Solution at left boundary converged!  KSPConvergedReason = ", reason
               !didItConverge = integerToRepresentTrue
            else
               print *,"[",myCommunicatorIndex,"] Solution at left boundary did not converge :(   KSPConvergedReason = ", reason
               didItConverge = integerToRepresentFalse
               if (solveSystem) then
                  stop
               end if
            end if
         else
            !didItConverge = integerToRepresentTrue
         end if
      case default
         print *,"Error! Invalid setting for leftBoundaryScheme"
         stop
      end select

      ! Where trajectories enter the domain, copy solnLeft to the global rhs:
      call VecGetArrayF90(solnLeft, solnArray, ierr)
      ipsi = 1
      do ispecies=1,numSpecies
         do itheta=1,Ntheta
            signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                 / (psiAHat*charges(ispecies))
            if (signOfPsiDot > -thresh) then
               do L=0,(Nxi-1)
                  do ix=1,Nx
                     index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta
                     call VecSetValue(rhs, index-1, solnArray(index), INSERT_VALUES, ierr)
                  end do
               end do
            end if
         end do
      end do
      call VecRestoreArrayF90(solnLeft, solnArray, ierr)

      !       call VecView(solnLeft, PETSC_VIEWER_STDOUT_SELF, ierr)

      call KSPDestroy(KSPBoundary, ierr)
      call VecDestroy(solnLeft, ierr)
      call VecDestroy(rhsLeft, ierr)
   end if

   if (procThatHandlesRightBoundary) then
      ! This process handles the right boundary, so solve the local kinetic equation there.

      call VecAssemblyBegin(rhsRight, ierr)
      call VecAssemblyEnd(rhsRight, ierr)

      call KSPCreate(MPI_COMM_SELF, KSPBoundary, ierr)
      if (useIterativeSolver) then
         ! Use iterative solver
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4:
         call KSPSetOperators(KSPBoundary, rightMatrix, rightPreconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
         ! Syntax for PETSc version 3.5 and later
         call KSPSetOperators(KSPBoundary, rightMatrix, rightPreconditionerMatrix, ierr)
#endif
         call KSPGetPC(KSPBoundary, PCBoundary, ierr)
         call PCSetType(PCBoundary, PCLU, ierr)
         call KSPSetType(KSPBoundary, KSPBCGSL, ierr)
         call KSPSetTolerances(KSPBoundary, solverTolerance, PETSC_DEFAULT_REAL, &
              PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
         call KSPSetFromOptions(KSPBoundary, ierr)
         call KSPMonitorSet(KSPBoundary, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
      else
         ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4:
         call KSPSetOperators(KSPBoundary, rightMatrix, rightMatrix, SAME_PRECONDITIONER, ierr)
#else
         ! Syntax for PETSc version 3.5 and later
         call KSPSetOperators(KSPBoundary, rightMatrix, rightMatrix, ierr)
#endif
         call KSPGetPC(KSPBoundary, PCBoundary, ierr)
         call PCSetType(PCBoundary, PCLU, ierr)
         call KSPSetType(KSPBoundary, KSPPREONLY, ierr)
         call KSPSetFromOptions(KSPBoundary, ierr)
      end if

      call VecDuplicate(rhsRight, solnRight, ierr)
      CHKERRQ(ierr)
      select case (rightBoundaryScheme)
      case (0)
         print *,"[",myCommunicatorIndex,"] Setting f_1 at right boundary to 0."
         call VecSet(solnRight, 0, ierr)
      case (1)
         print *,"[",myCommunicatorIndex,"] Proc",myRank, &
              " is solving local kinetic equation at right boundary ..."
         if (solveSystem) then
            call KSPSolve(KSPBoundary, rhsRight, solnRight, ierr)
         end if
         CHKERRQ(ierr)

         call PetscTime(time2, ierr)
         print *,"[",myCommunicatorIndex,"] Done solving for right boundary.  Time to solve: ", time2-time1, " seconds."
         call PetscTime(time1, ierr)

         if (useIterativeSolver) then
            call KSPGetConvergedReason(KSPBoundary, reason, ierr)
            if (reason>0) then
               print *,"[",myCommunicatorIndex,"] Solution at right boundary converged!  KSPConvergedReason = ", reason
               !didItConverge = integerToRepresentTrue
            else
               print *,"[",myCommunicatorIndex,"] Solution at right boundary did not converge :(   KSPConvergedReason = ", reason
               didItConverge = integerToRepresentFalse
               if (solveSystem) then
                  stop
               end if
            end if
         else
            !didItConverge = integerToRepresentTrue
         end if
      case default
         print *,"Error! Invalid setting for rightBoundaryScheme"
         stop
      end select

      ! Where trajectories enter the domain, copy solnRight to the global rhs:
      call VecGetArrayF90(solnRight, solnArray, ierr)
      ipsi = Npsi
      do ispecies = 1,numSpecies
         do itheta=1,Ntheta
            signOfPsiDot = -IHat(ipsi)*JHat(itheta,ipsi)*dBHatdtheta(itheta,ipsi) &
                 / (psiAHat*charges(ispecies))
            if (signOfPsiDot < thresh) then
               do L=0,(Nxi-1)
                  do ix=1,Nx
                     index = (ispecies-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + L*Ntheta + itheta
                     call VecSetValue(rhs, (ipsi-1)*localMatrixSize + index-1, solnArray(index), INSERT_VALUES, ierr)
                  end do
               end do
            end if
         end do
      end do
      call VecRestoreArrayF90(solnRight, solnArray, ierr)

      !call PCDestroy(PCBoundary, ierr)
      call KSPDestroy(KSPBoundary, ierr)
      call VecDestroy(solnRight, ierr)
      call VecDestroy(rhsRight, ierr)
   end if

   ! Finally, assemble the global RHS:
   call VecAssemblyBegin(rhs, ierr)
   call VecAssemblyEnd(rhs, ierr)

   ! ***********************************************************************
   ! ***********************************************************************
   ! 
   !  Permute the rows and columns of the linear system, if desired:
   !
   ! ***********************************************************************
   ! ***********************************************************************

 !!$    if (layout /= 0) then
 !!$       call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
 !!$            1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, permutationMatrix, ierr)
 !!$       CHKERRQ(ierr)
 !!$
 !!$       if (masterProcInSubComm) then
 !!$          select case (layout)
 !!$          case (1)
 !!$             ! Main section of the global matrix:
 !!$             do ipsi=1,Npsi
 !!$                do ix=1,Nx
 !!$                   do ixi=1,Nxi
 !!$                      do itheta=1,Ntheta
 !!$                         oldIndex = (ipsi-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         newIndex = (ix-1)*Npsi*Nxi*Ntheta + (ipsi-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                      end do
 !!$                   end do
 !!$                end do
 !!$             end do
 !!$
 !!$             ! Extra rows/columns for sources and constraints:
 !!$             do i=1,2
 !!$                do ipsi=1,Npsi
 !!$                   oldIndex = Npsi*Nx*Nxi*Ntheta + (i-1)*Npsi + ipsi
 !!$                   newIndex = oldIndex
 !!$                   call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                end do
 !!$             end do
 !!$
 !!$          case (2)
 !!$             ! Main section of the global matrix:
 !!$             do ipsi=1,Npsi
 !!$                do ix=1,Nx
 !!$                   do ixi=1,Nxi
 !!$                      do itheta=1,Ntheta
 !!$                         oldIndex = (ipsi-1)*Nx*Nxi*Ntheta + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         newIndex = (ipsi-1)*(Nx*Nxi*Ntheta+2) + (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta
 !!$                         call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                      end do
 !!$                   end do
 !!$                end do
 !!$             end do
 !!$
 !!$             ! Extra rows/columns for sources and constraints:
 !!$             do i=1,2
 !!$                do ipsi=1,Npsi
 !!$                   oldIndex = Npsi*Nx*Nxi*Ntheta + (i-1)*Npsi + ipsi
 !!$                   newIndex = (ipsi-1)*(Nx*Nxi*Ntheta+2) + Nx*Nxi*Ntheta + i
 !!$                   call MatSetValueSparse(permutationMatrix, oldIndex-1, newIndex-1, one, INSERT_VALUES, ierr)
 !!$                end do
 !!$             end do
 !!$
 !!$
 !!$          case default
 !!$             print *,"Error! Invalid setting for layout"
 !!$             stop
 !!$          end select
 !!$       end if
 !!$
 !!$       call MatAssemblyBegin(permutationMatrix, MAT_FINAL_ASSEMBLY, ierr)
 !!$       call MatAssemblyEnd(permutationMatrix, MAT_FINAL_ASSEMBLY, ierr)
 !!$
 !!$       call MatPtAP(matrix, permutationMatrix, MAT_INITIAL_MATRIX, one, tempMat, ierr)
 !!$       call MatDestroy(matrix, ierr)
 !!$       matrix = tempMat
 !!$
 !!$       CHKERRQ(ierr)
 !!$       call MatPtAP(preconditionerMatrix, permutationMatrix, MAT_INITIAL_MATRIX, one, tempMat, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call MatDestroy(preconditionerMatrix, ierr)
 !!$       CHKERRQ(ierr)
 !!$       preconditionerMatrix = tempMat
 !!$       CHKERRQ(ierr)
 !!$
 !!$       ! I get an error if I don't "create" tempVec using the statement below, even though
 !!$       ! it was not necessary to initialize tempMat in the same manner just above. This seems
 !!$       ! like an inconsistency in PETSc...
 !!$       call VecDuplicate(rhs, tempVec, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call MatMultTranspose(permutationMatrix, rhs, tempVec, ierr)
 !!$       CHKERRQ(ierr)
 !!$       call VecDestroy(rhs, ierr)
 !!$       CHKERRQ(ierr)
 !!$       rhs = tempVec
 !!$       CHKERRQ(ierr)
 !!$    end if

   ! ***********************************************************************
   ! ***********************************************************************
   ! 
   !  Solve the main linear system:
   !
   ! ***********************************************************************
   ! ***********************************************************************

   call KSPCreate(MPIComm, KSPInstance, ierr)
   CHKERRQ(ierr)

   if (useIterativeSolver) then
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4:
      call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
#else
         ! Syntax for PETSc version 3.5 and later
      call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, ierr)
#endif
      CHKERRQ(ierr)
      call KSPGetPC(KSPInstance, preconditionerContext, ierr)
      CHKERRQ(ierr)
      call PCSetType(preconditionerContext, PCLU, ierr)
      CHKERRQ(ierr)
      call KSPSetType(KSPInstance, KSPBCGSL, ierr)
      CHKERRQ(ierr)
      call KSPSetTolerances(KSPInstance, solverTolerance, PETSC_DEFAULT_REAL, &
           PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      CHKERRQ(ierr)
      call KSPSetFromOptions(KSPInstance, ierr)
      CHKERRQ(ierr)
      call KSPMonitorSet(KSPInstance, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
   else
      ! Direct solver:
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
         ! Syntax for PETSc versions up through 3.4:
      call KSPSetOperators(KSPInstance, matrix, matrix, SAME_PRECONDITIONER, ierr)
#else
         ! Syntax for PETSc version 3.5 and later
      call KSPSetOperators(KSPInstance, matrix, matrix, ierr)
#endif
      CHKERRQ(ierr)
      call KSPGetPC(KSPInstance, preconditionerContext, ierr)
      CHKERRQ(ierr)
      call PCSetType(preconditionerContext, PCLU, ierr)
      CHKERRQ(ierr)
      call KSPSetType(KSPInstance, KSPPREONLY, ierr)
      CHKERRQ(ierr)
      call KSPSetFromOptions(KSPInstance, ierr)
      !        call KSPSetTolerances(KSPInstance, 1d-5, PETSC_DEFAULT_REAL, &
      !          PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      CHKERRQ(ierr)
   end if

   ! The line below is remarked out to force use of either superlu or mumps,
   ! since PETSc's built-in solver doesn't seem to work when there are zeros on the diagonal.
   ! I may eventually want to handle this problem more elegantly...
   !    if (numProcsInSubComm > 1) then
   if (.true.) then
      select case (whichParallelSolverToFactorPreconditioner)
      case (1)
         call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERMUMPS, ierr)
         if (masterProcInSubComm) then
            print *,"[",myCommunicatorIndex,"] Using mumps to factorize the preconditioner."
         end if
      case (2)
         call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
         if (masterProcInSubComm) then
            print *,"[",myCommunicatorIndex,"] Using superlu_dist to factorize the preconditioner."
         end if
      case default
         if (masterProcInSubComm) then
            print *,"Error: Invalid setting for whichParallelSolverToFactorPreconditioner"
            stop
         end if
      end select
   else
      if (masterProcInSubComm) then
         print *,"[",myCommunicatorIndex,"] Using PETSc's serial sparse direct solver to factorize the preconditioner."
      end if
   end if

   call MatGetInfo(matrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
   NNZMain = nint(myMatInfo(MAT_INFO_NZ_USED))
   if (useIterativeSolver) then
      call MatGetInfo(preconditionerMatrix, MAT_GLOBAL_SUM, myMatInfo, ierr)
      NNZPreconditioner = nint(myMatInfo(MAT_INFO_NZ_USED))
   end if
   if (masterProcInSubComm) then
      if (useIterativeSolver) then
         print *,"[",myCommunicatorIndex,"] # of nonzeros in matrix:",NNZMain,", # of nonzeros in preconditioner:", &
              NNZPreconditioner,", ratio:",((one*NNZPreconditioner)/NNZMain)
      else
         print *,"[",myCommunicatorIndex,"] # of nonzeros in matrix:",NNZMain
      end if
   end if

   call VecDuplicate(rhs, soln, ierr)
   CHKERRQ(ierr)
   if (masterProcInSubComm) then
      print *,"[",myCommunicatorIndex,"] Beginning the main solve.  This could take a while ..."
   end if

   if (solveSystem) then
      call KSPSolve(KSPInstance, rhs, soln, ierr)
   end if
   CHKERRQ(ierr)

   call PetscTime(time2, ierr)
   if (masterProcInSubComm) then
      print *,"[",myCommunicatorIndex,"] Done with the main solve.  Time to solve: ", time2-time1, " seconds."
   end if
   call PetscTime(time1, ierr)

   if (useIterativeSolver) then
      call KSPGetConvergedReason(KSPInstance, reason, ierr)
      if (reason>0) then
         if (masterProcInSubComm) then
            print *,"[",myCommunicatorIndex,"] Converged!  KSPConvergedReason = ", reason
         end if
         !didItConverge = integerToRepresentTrue
      else
         if (masterProcInSubComm) then
            print *,"[",myCommunicatorIndex,"] Did not converge :(   KSPConvergedReason = ", reason
         end if
         didItConverge = integerToRepresentFalse
      end if
   else
      !didItConverge = integerToRepresentTrue
   end if

   if (layout /= 0) then
      call MatMult(permutationMatrix, soln, tempVec, ierr)
      call VecDestroy(soln, ierr)
      soln = tempVec
   end if

   !**************************************************************************
   !**************************************************************************
   ! 
   !  Calculate moments of the solution:
   !
   !**************************************************************************
   !**************************************************************************

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

 !!$          if (psiDerivativeScheme == 0) then
 !!$             ddpsiForKTheta = ddpsiLeft
 !!$          else
 !!$             ! centered finite differences, no upwinding, 3-point stencil
 !!$             scheme = 2
 !!$             call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForKTheta, d2dpsi2)
 !!$          end if
 !!$          !       pPerpTermInKThetaWith3PointStencil = transpose(matmul(ddpsiForKTheta, &
 !!$          !            transpose(pPerpTermInKThetaBeforePsiDerivative)))
 !!$
 !!$          if (psiDerivativeScheme == 0) then
 !!$             ddpsiForKTheta = ddpsiLeft
 !!$          else
 !!$             ! centered finite differences, no upwinding, 5-point stencil
 !!$             scheme = 12
 !!$             call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, psiWeights, ddpsiForKTheta, d2dpsi2)
 !!$          end if
 !!$          !       pPerpTermInKThetaWith5PointStencil = transpose(matmul(ddpsiForKTheta, &
 !!$          !            transpose(pPerpTermInKThetaBeforePsiDerivative)))

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
  end if


  ! *********************************************************
  ! Create a PETSc viewer to record output
  ! *********************************************************

  if (saveMatlabOutput) then
     call PetscViewerASCIIOpen(MPIComm, &
          & MatlabOutputFilename,&
          & MatlabOutput, ierr)
     CHKERRQ(ierr)
     call PetscViewerSetFormat(MatlabOutput, PETSC_VIEWER_ASCII_MATLAB, ierr)
     CHKERRQ(ierr)

     call PetscObjectSetName(rhs, "rhs", ierr)
     CHKERRQ(ierr)
     call VecView(rhs, MatlabOutput, ierr)
     CHKERRQ(ierr)
     call PetscObjectSetName(soln, "soln", ierr)
     CHKERRQ(ierr)
     call VecView(soln, MatlabOutput, ierr)
     CHKERRQ(ierr)

     if (useIterativeSolver) then
        call PetscObjectSetName(preconditionerMatrix, "preconditionerMatrix", ierr)
        CHKERRQ(ierr)
        call MatView(preconditionerMatrix, MatlabOutput, ierr)
        CHKERRQ(ierr)
     end if
     call PetscObjectSetName(matrix, "matrix", ierr)
     CHKERRQ(ierr)
     call MatView(matrix, MatlabOutput, ierr)
     CHKERRQ(ierr)

     call PetscTime(time2, ierr)
     if (masterProcInSubComm) then
        print *,"[",myCommunicatorIndex,"] Time to write output: ", time2-time1, " seconds."
     end if
     call PetscTime(time1, ierr)

     call PetscViewerDestroy(MatlabOutput, ierr)
     CHKERRQ(ierr)
  end if

  !    call VecDestroy(rhs, ierr)
  call VecDestroy(soln, ierr)
  CHKERRQ(ierr)
  call MatDestroy(matrix, ierr)
  if (useIterativeSolver) then
     call MatDestroy(preconditionerMatrix, ierr)
  end if
  call KSPDestroy(KSPInstance, ierr)
  CHKERRQ(ierr)


  call PetscTime(time2, ierr)
  elapsedTime = time2 - startTime

  call printOutputs()

end subroutine solveDKE

