! For compilers that do not include the error function erf(x), the line
! below should be un-commented:
!#define USE_GSL_ERF

module profiles

  use globalVariables
  use grids
  use readHDF5Input
  use sourcesConstraints

  implicit none

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

contains

  subroutine initializeProfiles()
    ! This is the central subroutine of this module.
    ! Here we must fill the following arrays:
    ! PhiHat(Npsi)
    ! dPhiHatdpsi(Npsi)
    ! THats(Nspecies,Npsi)
    ! dTHatdpsis(Nspecies,Npsi)
    ! nHats(Nspecies,Npsi)
    ! dNHatdpsis(Nspecies,NPsi)
    ! etaHats(Nspecies,Npsi)
    ! detaHatdpsis(Nspecies,Npsi)

    implicit none

    integer :: i, ipsi, ispecies, LAPACKInfo, numIterations, nPsiFine, scheme
    PetscScalar :: temp1, temp2, FSABHat2Fine, IHatFine, analyticHeatFluxAtPsiMid
    PetscScalar :: ss, r0, rampAmplitude
    PetscScalar, dimension(:), allocatable :: nHat, THat, dnHatdpsi, dTHatdpsi, etaHat, detaHatdpsi
    PetscScalar, dimension(:), allocatable :: psiFine, PhiHatFine, rhsForT, THatFine, dTHatdpsiFine
    PetscScalar, dimension(:), allocatable :: nHatFine, dnHatdpsiFine, psiWeightsFine
    PetscScalar, dimension(:), allocatable :: etaHatFine, rhsForTPrime
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
    PetscScalar, dimension(:), allocatable :: UFine, UFactor, dPhiHatdpsiFine
    ! PetscScalar, dimension(:,:), allocatable :: ddpsi_accurate
    PetscScalar, dimension(:,:), allocatable :: ddpsiForT, ddpsiForTCopy
    PetscScalar, dimension(:,:), allocatable :: psiInterpolationMatrix
    PetscScalar, dimension(:,:), allocatable :: tempMatrix
    PetscScalar, dimension(:,:), allocatable :: matrixForT, matrixForTPrime
    PetscScalar, dimension(:), allocatable :: temp_psiWeights
    PetscScalar, dimension(:,:), allocatable :: temp_d2dpsi2
    PetscScalar, allocatable, dimension(:,:) :: TInterpolationRowVec
    PetscScalar, dimension(1) :: TAtPsiMid
    character(len=100) :: HDF5Groupname
    PetscScalar, dimension(:), allocatable :: temp_RHS
    PetscScalar, dimension(:,:), allocatable :: temp_constantSourceProfile
    
    allocate(PhiHat(Npsi))
    allocate(dPhiHatdpsi(Npsi))
    allocate(THats(Nspecies,Npsi))
    allocate(dTHatdpsis(Nspecies,Npsi))
    allocate(etaHats(Nspecies,Npsi))
    allocate(detaHatdpsis(Nspecies,Npsi))
    allocate(nHats(Nspecies,Npsi))
    allocate(dnHatdpsis(Nspecies,Npsi))
    ! multiplies global terms
    allocate(globalTermMultiplier(Npsi))
    
    select case (psiGridType)
    case(0)
       ! uniform grid
       do i=1,Npsi
          psiAHatArray(i) = psiAHat
       end do
    case(1)
       ! Create groupname to be read
       write (HDF5Groupname,"(A4,I0)") "Npsi", Npsi
       
       ! Open input file
       call openInputFile(psiAHatFilename,HDF5Groupname)
       
       call readVariable(psiAHatArray, "psiAHatArray")

       call closeInputFile() 
       
    case default
       print *,"Error! Invalid setting for psiGridType"
       stop
         
    end select

    select case (useGlobalTermMultiplier)
    case(0)
       ! no multiplier
       globalTermMultiplier = one
    case(1)
       ! Create groupname to be read
        write (HDF5Groupname,"(A4,I0)") "Npsi", Npsi
       ! Open input file
         call openInputFile(globalTermMultiplierFilename,HDF5Groupname)

         call readVariable(globalTermMultiplier, "globalTermMultiplier")

         call closeInputFile() 
       
      case default
         print *,"Error! Invalid setting for useGlobalTermMultiplier" 
         stop

      end select

      !!!!!!! TODO MOVE TO SOURCE INITIALIZATION
      allocate(sourceThetaPart(Nsources,Ntheta))
      allocate(sourceThetaPartFSA(Nsources,NEnforcedPsi))
      allocate(extraSourceThetaPart(NextraSources,Ntheta))
      allocate(extraSourceThetaPartFSA(NextraSources,NEnforcedPsi))
      allocate(extraSourceSpeciesPart(NextraSources,Nspecies))

      allocate(sourceConstraintsRHS(NextraSources,NEnforcedPsi))
      
      allocate(constantSourceProfile(NconstantSources,Nspecies,Npsi))
      
      do i = 1,NextraSources 
         select case (RHSFromFile(i))
         case(0)
            ! no charge source
            sourceConstraintsRHS(i,:) = zero
         case(1)
            ! read charge source from file
            ! Create groupname to be read
            allocate(temp_RHS(NEnforcedPsi))
            write (HDF5Groupname,"(A4,I0)") "Npsi", Npsi
            ! Open input file
            call openInputFile(sourceConstraintsFilenames(i),HDF5Groupname)
            
            call readVariable(temp_RHS, "chargeSource")
            sourceConstraintsRHS(i,:) = temp_RHS
            call closeInputFile()
            deallocate(temp_RHS)
      
         case default
            print *,"Error! Invalid setting for RHSFromFile. Suppoted values: 0,1." 
            stop
            
         end select
      end do

      do i = 1,NconstantSources 
         ! read charge source from file
         allocate(temp_constantSourceProfile(Nspecies,Npsi))
         ! Create groupname to be read
         write (HDF5Groupname,"(A4,I0)") "Npsi", Npsi
         ! Open input file
         call openInputFile(constantSourcesFilenames(i),HDF5Groupname)
         ! Read constant source
         call readVariable(temp_constantSourceProfile, "constantSource")
         constantSourceProfile(i,:,:) = temp_constantSourceProfile
         call closeInputFile()
         deallocate(temp_constantSourceProfile)
      end do
      
    if (profilesScheme .eq. 7) then
       ! Read in profiles from file

       ! Check profilesFilename is set
       if (.not. len(profilesFilename)>=0) then
          print *,"If profilesScheme==7 then profilesFilename must be set."
          stop
       end if
       
       ! Create groupname to be read
       write (HDF5Groupname,"(A4,I0)") "Npsi", Npsi

       ! Open input file
       call openInputFile(profilesFilename,HDF5Groupname)

       ! Read in profiles
       call readVariable(PhiHat, "PhiHat")
       call readVariable(dPhiHatdpsi, "dPhiHatdpsi")
       call readVariable(THats, "THats")
       call readVariable(dTHatdpsis, "dTHatdpsis")
       call readVariable(nHats, "nHats")
       call readVariable(dnHatdpsis, "dnHatdpsis")
       call readVariable(etaHats, "etaHats")
       call readVariable(detaHatdpsis, "detaHatdpsis")
       
       ! Close input file
       call closeInputFile()

    else

       allocate(THat(Npsi))
       allocate(dTHatdpsi(Npsi))
       allocate(etaHat(Npsi))
       allocate(detaHatdpsi(Npsi))
       allocate(nHat(Npsi))
       allocate(dnHatdpsi(Npsi))
       allocate(temp_psiWeights(Npsi))
       ! allocate(ddpsi_accurate(Npsi,Npsi))
       allocate(temp_d2dpsi2(Npsi,Npsi))

       select case (profilesScheme)
       case (0)
          ! Simplest radial variation
          do i=1,Npsi
             temp1 = (psi(i)-psiMid)/(0.03d+0)
#ifdef USE_GSL_ERF
             call erf(temp1, temp2)
#else
             temp2 = erf(temp1)
#endif
             nHat(i) = 1-oneHalf*temp2
             dnHatdpsi(i) = -1d0/sqrt(pi)/0.03d0*exp(-(psi(i)-psiMid)**2/9d-4)
          end do

          THat = 1 - (psi-1)
          dTHatdpsi = -1
          etaHat = 1
          detaHatdpsi = 0d0
          do i=1,Npsi
             phiHat(i) = delta/(2*omega)*THat(i) * log(etaHat(i)/nHat(i))
             dphiHatdpsi(i) = delta/2d0/omega*( dThatdpsi(i)*log(etaHat(i)/nHat(i)) &
                                                - That(i)*dnHatdpsi(i)/nHat(i)/etaHat(i) )
          end do

       case (1)
          ! More complicated radial variation
          do i=1,Npsi
             temp1 = (psi(i)-psiMid)/(0.03d+0)
#ifdef USE_GSL_ERF
             call erf(temp1, temp2)
#else
             temp2 = erf(temp1)
#endif
             nHat(i) = 1-oneHalf*temp2
             dnHatdpsi(i) = -1d0/sqrt(pi)/0.03d0*exp(-(psi(i)-psiMid)**2/9d-4)
          end do

          do i=1,Npsi
             THat(i) = exp(1-psi(i))
             dTHatdpsi(i) = -exp(1-psi(i))
             etaHat(i) = exp((psi(i)-1)*(1.7d+0))
             detaHatdpsi(i) = 1.7d0*etaHat(i)
          end do
          do i=1,Npsi
             phiHat(i) = delta/(2*omega)*THat(i) * log(etaHat(i)/nHat(i))
             dphiHatdpsi(i) = delta/2d0/omega*( dThatdpsi(i)*log(etaHat(i)/nHat(i)) &
                                                - That(i)*dnHatdpsi(i)/nHat(i)/etaHat(i) &
                                                + That(i)*detaHatdpsi(i)*nHat(i)/etaHat(i) )
          end do
       case (2)
          ! Radial profiles with a large region of constant U,
          ! for comparison with analytic theory
          delta = (0.0006d+0)*((epsil/(0.01d+0)) ** 2)/(desiredFWHMInRhoTheta/(0.24d+0))/30
          omega = delta
          psiAHat = epsil * epsil / Miller_q
          !nuPrime = nuPrimeMeaningful / Miller_q

          etaHat = 1

          ss = 150d+0
          r0 = 0.07d+0 + widthExtender
          rampAmplitude = desiredU * psiAHat / omega

          if (setTPrimeToBalanceHeatFlux) then
             numIterations = 10
             NpsiFine = 200

             allocate(psiFine(NpsiFine))
             allocate(PhiHatFine(NpsiFine))
             allocate(nHatFine(NpsiFine))
             allocate(THatFine(NpsiFine))
             allocate(dTHatdpsiFine(NpsiFine))
             allocate(dnHatdpsiFine(NpsiFine))
             allocate(etaHatFine(NpsiFine))
             allocate(ddpsiForT(NpsiFine, NpsiFine))
             allocate(tempMatrix(NpsiFine, NpsiFine))
             allocate(ddpsiForTCopy(NpsiFine, NpsiFine))
             allocate(psiWeightsFine(NpsiFine))
             allocate(rhsForT(NpsiFine))

             scheme = 12
             call uniformDiffMatrices(NpsiFine, psiMin-(1d-10), psiMax+(1d-10), scheme, psiFine, &
                  psiWeightsFine, ddpsiForT, tempMatrix)
             deallocate(tempMatrix)
             !call ChebyshevGrid(NpsiFine, psiMin, psiMax, psiFine, psiWeightsFine, ddpsiForT)

             etaHatFine = 1 + (psiFine-psiMid)*detaHatdpsiScalar
             call computePhiHat_flatURegion(NpsiFine, psiFine, psiMid, ss, r0, rampAmplitude, PhiHatFine)
             THatFine = 1
             ddpsiForT(1,:)=0
             ddpsiForT(1,1)=1

             allocate(IPIV(NpsiFine))
             do i=1,numIterations
                nHatFine = etaHatFine * exp(-2*omega/delta*PhiHatFine/THatFine)
                dTHatdpsiFine = dTHatdpsiScalar / (nHatFine ** exponent)
                rhsForT = dTHatdpsiFine
                rhsForT(1) = 0
                ddpsiForTCopy = ddpsiForT

                ! The command below overwrites both ddpsiForTCopy and rhsForT:
#if defined(PETSC_USE_REAL_SINGLE)
                call SGESV(NpsiFine, 1, ddpsiForTCopy, NpsiFine, IPIV, rhsForT, NpsiFine, LAPACKInfo)
#else
                call DGESV(NpsiFine, 1, ddpsiForTCopy, NpsiFine, IPIV, rhsForT, NpsiFine, LAPACKInfo)
#endif
                if (LAPACKInfo /= 0) then
                   print *,"LAPACK error!!  Info = ",LAPACKInfo
                   stop
                end if

                THatFine = rhsForT
                if (mod(NpsiFine,2)==1) then
                   THatFine = THatFine - THatFine((NpsiFine+1)/2)
                else
                   THatFine = THatFine - (THatFine(NpsiFine/2) + THatFine(NpsiFine/2+1))/2
                end if
                THatFine = THatFine + 1

             end do
             deallocate(IPIV)

             !call ChebyshevInterpolation(NpsiFine, Npsi, THatFine, psiMin, psiMax, psi, THat)
             !call ChebyshevInterpolation(NpsiFine, Npsi, dTHatdpsiFine, psiMin, psiMax, psi, dTHatdpsi)
             !call ChebyshevInterpolation(NpsiFine, Npsi, nHatFine, psiMin, psiMax, psi, nHat)
             !THat = THatFine
             !dTHatdpsi = dTHatdpsiFine
             !nHat = nHatFine
             allocate(psiInterpolationMatrix(Npsi,NpsiFine))
             call interpolationMatrix(NpsiFine, Npsi, psiFine, psi, psiInterpolationMatrix, -1, 0)
             THat = matmul(psiInterpolationMatrix, THatFine)
             dTHatdpsi = matmul(psiInterpolationMatrix, dTHatdpsiFine)
             nHat = matmul(psiInterpolationMatrix, nHatFine)
             ! Re-compute the derivative matrices since they have been over-written above
             call uniformDiffMatrices(NpsiFine, psiMin-(1d-10), psiMax+(1d-10), scheme, psiFine, &
                  psiWeightsFine, ddpsiForT, tempMatrix)
             dnHatdpsiFine = matmul(ddpsiForT, nHatFine)
             dnHatdpsi = matmul(psiInterpolationMatrix, dnHatdpsiFine)
             deallocate(psiInterpolationMatrix)

             deallocate(psiFine)
             deallocate(PhiHatFine)
             deallocate(nHatFine)
             deallocate(THatFine)
             deallocate(dTHatdpsiFine)
             deallocate(etaHatFine)
             deallocate(ddpsiForT)
             deallocate(ddpsiForTCopy)
             deallocate(psiWeightsFine)
             deallocate(rhsForT)

          else
             THat = 1 + dTHatdpsiScalar * (psi-psiMid)
             dTHatdpsi = dTHatdpsiScalar
          end if

          etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar
          detaHatdpsi = detaHatdpsiScalar
          call computePhiHat_flatURegion(Npsi, psi, psiMid, ss, r0, rampAmplitude, PhiHat)
          call computedPhiHatdpsi_flatURegion(Npsi, psi, psiMid, ss, r0, rampAmplitude, dPhiHatdpsi)

          if (.not. setTPrimeToBalanceHeatFlux) then
             nHat = etaHat * exp(-2*omega/delta*PhiHat/THat)
             dnHatdpsi = ( detaHatdpsi &
                           -etaHat*2d0*omega/delta*dPhiHatdpsi/THat &
                           + etaHat*2d0*omega/delta*PhiHat*dTHatdpsi/THat**2 ) * exp(-2*omega/delta*PhiHat/THat)
          end if

       case (3)
          ! Radial profiles with a Gaussian U(psi) which should have a constant
          ! radial heat flux, according to the incorrect analytic theory.

          delta = (0.0006d+0)*((epsil/(0.01d+0))**1)/(desiredFWHMInRhoTheta/(0.24d+0))/30
          omega = delta
          psiAHat = epsil * epsil / Miller_q
          !nuPrime = nuPrimeMeaningful / Miller_q

!!$       ! For profilesScheme=3, we need <B^2> computed earlier than for other profile schemes:
!!$       do i=1,Npsi
!!$          VPrimeHat(i) = dot_product(thetaWeights, 1/JHat(:,i))
!!$          FSABHat2(i) = dot_product(thetaWeights, BHat(:,i) * BHat(:,i) / JHat(:,i)) / VPrimeHat(i)
!!$       end do
          FSABHat2Fine = FSABHat2(1) ! Assumes no radial variation in <B^2>

          ss = 25
          rampAmplitude = desiredU * psiAHat / omega

          if (setTPrimeToBalanceHeatFlux) then
             numIterations = 10
             NpsiFine = 200

             allocate(psiFine(NpsiFine))
             allocate(PhiHatFine(NpsiFine))
             allocate(dPhiHatdpsiFine(NpsiFine))
             allocate(nHatFine(NpsiFine))
             allocate(dnHatdpsiFine(NpsiFine))
             allocate(THatFine(NpsiFine))
             allocate(dTHatdpsiFine(NpsiFine))
             allocate(etaHatFine(NpsiFine))
             allocate(ddpsiForT(NpsiFine, NpsiFine))
             allocate(matrixForT(NpsiFine, NpsiFine))
             allocate(matrixForTPrime(NpsiFine, NpsiFine))
             allocate(tempMatrix(NpsiFine, NpsiFine))
             allocate(psiWeightsFine(NpsiFine))
             allocate(rhsForT(NpsiFine))
             allocate(rhsForTPrime(NpsiFine))
             allocate(UFine(NpsiFine))
             allocate(UFactor(NpsiFine))

             scheme = 12
             call uniformDiffMatrices(NpsiFine, psiMin-(1d-10), psiMax+(1d-10), scheme, psiFine, &
                  psiWeightsFine, ddpsiForT, tempMatrix)
             !call ChebyshevGrid(NpsiFine, psiMin, psiMax, psiFine, psiWeightsFine, ddpsiForT)

             etaHatFine = 1 + (psiFine-psiMid)*detaHatdpsiScalar
             !call computePhiHat(NpsiFine, psiFine, psiMid, ss, r0, rampAmplitude, PhiHatFine)
             do ipsi=1,NpsiFine
                temp1 = (psiFine(ipsi)-psiMid)*ss
#ifdef USE_GSL_ERF
                call erf(temp1, temp2)
#else
                temp2 = erf(temp1)
#endif
                PhiHatFine(ipsi) = rampAmplitude * temp2 * sqrtpi/(2*ss)
                dPhiHatdpsiFine(ipsi) = rampAmplitude * exp(-((psiFine(ipsi)-psiMid)*ss) ** 2)
             end do
             THatFine = 1
             IHatFine = 1

             matrixForT = ddpsiForT
             matrixForT(1,:)=0
             matrixForT(1,1)=1

             allocate(IPIV(NpsiFine))
             do i=1,numIterations
                nHatFine = etaHatFine * exp(-2*omega/delta*PhiHatFine/THatFine)

                UFine = omega*IHatFine*dPhiHatdpsiFine/(psiAHat*sqrt(FSABHat2Fine*THatFine))
                UFactor = (one/three)*(4*(UFine**8) + 16*(UFine**6) + 24*(UFine**4) + 12*UFine**2 + 3) &
                     / (2*(UFine**4) + 2*(UFine**2) + 1) * exp(-UFine**2) * nHatFine * THatFine * sqrt(THatFine)
                rhsForTPrime = dTHatdpsiScalar
                do ipsi=1,NpsiFine
                   matrixForTPrime(ipsi,:) = -UFactor(ipsi)*delta/psiAHat*UFine(ipsi) &
                        * sqrt(THatFine(ipsi)/FSABHat2Fine)*IHatFine &
                        * ddpsiForT(ipsi,:)
                   matrixForTPrime(ipsi,ipsi) = matrixForTPrime(ipsi,ipsi) + UFactor(ipsi)
                end do

                ! The command below overwrites both matrixForTPrime and rhsForTPrime:
#if defined(PETSC_USE_REAL_SINGLE)
                call SGESV(NpsiFine, 1, matrixForTPrime, NpsiFine, IPIV, rhsForTPrime, NpsiFine, LAPACKInfo)
#else
                call DGESV(NpsiFine, 1, matrixForTPrime, NpsiFine, IPIV, rhsForTPrime, NpsiFine, LAPACKInfo)
#endif
                if (LAPACKInfo /= 0) then
                   print *,"LAPACK error 1!!  Info = ",LAPACKInfo
                   stop
                end if

                dTHatdpsiFine = rhsForTPrime

                rhsForT = dTHatdpsiFine
                tempMatrix = matrixForT
                ! The command below overwrites both matrixForTPrime and rhsForTPrime:
#if defined(PETSC_USE_REAL_SINGLE)
                call SGESV(NpsiFine, 1, matrixForT, NpsiFine, IPIV, rhsForT, NpsiFine, LAPACKInfo)
#else
                call DGESV(NpsiFine, 1, matrixForT, NpsiFine, IPIV, rhsForT, NpsiFine, LAPACKInfo)
#endif
                if (LAPACKInfo /= 0) then
                   print *,"LAPACK error 2!!  Info = ",LAPACKInfo
                   stop
                end if
                matrixForT = tempMatrix
                THatFine = rhsForT

                ! Shift THatFine so its value is 1 at psiMid:
                if (mod(NpsiFine,2)==1) then
                   THatFine = THatFine - THatFine((NpsiFine+1)/2)
                else
                   THatFine = THatFine - (THatFine(NpsiFine/2) + THatFine(NpsiFine/2+1))/2
                end if
                THatFine = THatFine + 1

             end do
             deallocate(IPIV)

             !call ChebyshevInterpolation(NpsiFine, Npsi, THatFine, psiMin, psiMax, psi, THat)
             !call ChebyshevInterpolation(NpsiFine, Npsi, dTHatdpsiFine, psiMin, psiMax, psi, dTHatdpsi)
             !call ChebyshevInterpolation(NpsiFine, Npsi, nHatFine, psiMin, psiMax, psi, nHat)
             !THat = THatFine
             !dTHatdpsi = dTHatdpsiFine
             !nHat = nHatFine
             allocate(psiInterpolationMatrix(Npsi,NpsiFine))
             call interpolationMatrix(NpsiFine, Npsi, psiFine, psi, psiInterpolationMatrix, -1, 0)
             THat = matmul(psiInterpolationMatrix, THatFine)
             dTHatdpsi = matmul(psiInterpolationMatrix, dTHatdpsiFine)
             nHat = matmul(psiInterpolationMatrix, nHatFine)
             dnHatdpsiFine = matmul(ddpsiForT, nHat)
             dnHatdpsi = matmul(psiInterpolationMatrix, dnHatdpsiFine)
             deallocate(psiInterpolationMatrix)

             deallocate(tempMatrix)
             deallocate(psiFine)
             deallocate(PhiHatFine)
             deallocate(dPhiHatdpsiFine)
             deallocate(nHatFine)
             deallocate(dnHatdpsiFine)
             deallocate(THatFine)
             deallocate(dTHatdpsiFine)
             deallocate(etaHatFine)
             deallocate(ddpsiForT)
             deallocate(matrixForT)
             deallocate(matrixForTPrime)
             deallocate(psiWeightsFine)
             deallocate(rhsForT)
             deallocate(rhsForTPrime)
             deallocate(UFine)
             deallocate(UFactor)

          else
             THat = 1 + dTHatdpsiScalar * (psi-psiMid)
             dTHatdpsi = dTHatdpsiScalar
          end if

          etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar
          detaHatdpsi = detaHatdpsiScalar
          do i=1,Npsi
             temp1 = (psi(i)-psiMid)*ss
#ifdef USE_GSL_ERF
             call erf(temp1, temp2)
#else
             temp2 = erf(temp1)
#endif
             PhiHat(i) = rampAmplitude*temp2*sqrtpi/(2*ss)
             dPhiHatdpsi(i) = rampAmplitude*exp(-((psi(i)-psiMid)*ss)**2)
          end do

          if (.not. setTPrimeToBalanceHeatFlux) then
             nHat = etaHat * exp(-2*omega/delta*PhiHat/THat)
             dnHatdpsi = ( detaHatdpsi &
                           - etaHat*2d0*omega/delta*dPhiHatdpsi/THat &
                           + etaHat*2d0*omega/delta*PhiHat*dTHatdpsi/THat**2 ) * exp(-2*omega/delta*PhiHat/THat)
          end if

       case (4,5,6)

          delta = (0.0006d+0)*((epsil/(0.01d+0))**2)/(desiredFWHMInRhoTheta/(0.24d+0))/30
          omega = delta
          psiAHat = epsil * epsil / Miller_q
          !nuPrime = nuPrimeMeaningful / Miller_q

!!$       ! For profilesScheme=3, we need <B^2> computed earlier than for other profile schemes:
!!$       do i=1,Npsi
!!$          VPrimeHat(i) = dot_product(thetaWeights, 1/JHat(:,i))
!!$          FSABHat2(i) = dot_product(thetaWeights, BHat(:,i) * BHat(:,i) / JHat(:,i)) / VPrimeHat(i)
!!$       end do
          FSABHat2Fine = FSABHat2(1) ! Assumes no radial variation in <B^2>

          if (profilesScheme==4 .or. profilesScheme==6) then
             ss = 25
          else
             ss = 150d+0
             r0 = 0.07d+0 + widthExtender
          end if
          rampAmplitude = desiredU * psiAHat / omega


          if (setTPrimeToBalanceHeatFlux) then
             numIterations = 10
             NpsiFine = 200

             allocate(psiFine(NpsiFine))
             allocate(PhiHatFine(NpsiFine))
             allocate(dPhiHatdpsiFine(NpsiFine))
             allocate(nHatFine(NpsiFine))
             allocate(dnHatdpsiFine(NpsiFine))
             allocate(THatFine(NpsiFine))
             allocate(dTHatdpsiFine(NpsiFine))
             allocate(etaHatFine(NpsiFine))
             allocate(ddpsiForT(NpsiFine, NpsiFine))
             allocate(ddpsiForTCopy(NpsiFine, NpsiFine))
             allocate(psiWeightsFine(NpsiFine))
             allocate(rhsForT(NpsiFine))
             allocate(UFine(NpsiFine))
             allocate(UFactor(NpsiFine))

             scheme = 12
             call uniformDiffMatrices(NpsiFine, psiMin-(1d-10), psiMax+(1d-10), scheme, psiFine, &
                  psiWeightsFine, ddpsiForT, ddpsiForTCopy)
             !call ChebyshevGrid(NpsiFine, psiMin, psiMax, psiFine, psiWeightsFine, ddpsiForT)

             allocate(TInterpolationRowVec(1, NpsiFine))
             call interpolationMatrix(NpsiFine, 1, psiFine, psiMid, TInterpolationRowVec, -1, 0)

             etaHatFine = 1 + (psiFine-psiMid)*detaHatdpsiScalar
             !call computePhiHat(NpsiFine, psiFine, psiMid, ss, r0, rampAmplitude, PhiHatFine)
             if (profilesScheme==4 .or. profilesScheme==6) then
                do ipsi=1,NpsiFine
                   temp1 = (psiFine(ipsi)-psiMid)*ss
#ifdef USE_GSL_ERF
                   call erf(temp1, temp2)
#else
                   temp2 = erf(temp1)
#endif
                   PhiHatFine(ipsi) = rampAmplitude * temp2 * sqrtpi/(2*ss)
                   dPhiHatdpsiFine(ipsi) = rampAmplitude * exp(-((psiFine(ipsi)-psiMid)*ss) ** 2)
                end do
             else
                call computePhiHat_flatURegion(NpsiFine, psiFine, psiMid, ss, r0, rampAmplitude, PhiHatFine)
                call computedPhiHatdpsi_flatURegion(NpsiFine, psiFine, psiMid, ss, r0, rampAmplitude, dPhiHatdpsiFine)
             end if

             THatFine = 1
             IHatFine = 1

             ddpsiForT(1,:)=0
             ddpsiForT(1,1)=1

             analyticHeatFluxAtPsiMid = (one/three)*(4*(desiredU**8) + 16*(desiredU**6) + 24*(desiredU**4) + 12*desiredU**2 + 3) &
                  / (2*(desiredU**4) + 2*(desiredU**2) + 1) * exp(-desiredU**2)

             allocate(IPIV(NpsiFine))
             do i=1,numIterations
                nHatFine = etaHatFine * exp(-2*omega/delta*PhiHatFine/THatFine)

                UFine = omega*IHatFine*dPhiHatdpsiFine/(psiAHat*sqrt(FSABHat2Fine*THatFine))
                UFactor = (one/three)*(4*(UFine**8) + 16*(UFine**6) + 24*(UFine**4) + 12*UFine**2 + 3) &
                     / (2*(UFine**4) + 2*(UFine**2) + 1) * exp(-UFine**2) * nHatFine * THatFine * sqrt(THatFine)

                if (profilesScheme==6) then
                   dTHatdpsiFine = dTHatdpsiScalar / (nHatFine ** exponent)
                else
                   dTHatdpsiFine = dTHatdpsiScalar * analyticHeatFluxAtPsiMid / UFactor
                end if
                rhsForT = dTHatdpsiFine
                rhsForT(1) = 0

                ddpsiForTCopy = ddpsiForT

                ! The command below overwrites both ddpsiForTCopy and rhsForT:
#if defined(PETSC_USE_REAL_SINGLE)
                call SGESV(NpsiFine, 1, ddpsiForTCopy, NpsiFine, IPIV, rhsForT, NpsiFine, LAPACKInfo)
#else
                call DGESV(NpsiFine, 1, ddpsiForTCopy, NpsiFine, IPIV, rhsForT, NpsiFine, LAPACKInfo)
#endif
                if (LAPACKInfo /= 0) then
                   print *,"LAPACK error!!  Info = ",LAPACKInfo
                   stop
                end if

                THatFine = rhsForT

!                if (mod(NpsiFine,2)==1) then
!                   THatFine = THatFine - THatFine((NpsiFine+1)/2)
!                else
!                   THatFine = THatFine - (THatFine(NpsiFine/2) + THatFine(NpsiFine/2+1))/2
!                end if
!                THatFine = THatFine + 1

                TAtPsiMid = matmul(TInterpolationRowVec, THatFine)
                THatFine = THatFine - TAtPsiMid(1) + 1

             end do

             deallocate(IPIV)

             !call ChebyshevInterpolation(NpsiFine, Npsi, THatFine, psiMin, psiMax, psi, THat)
             !call ChebyshevInterpolation(NpsiFine, Npsi, dTHatdpsiFine, psiMin, psiMax, psi, dTHatdpsi)
             !call ChebyshevInterpolation(NpsiFine, Npsi, nHatFine, psiMin, psiMax, psi, nHat)
             !THat = THatFine
             !dTHatdpsi = dTHatdpsiFine
             !nHat = nHatFine
             allocate(psiInterpolationMatrix(Npsi,NpsiFine))
             call interpolationMatrix(NpsiFine, Npsi, psiFine, psi, psiInterpolationMatrix, -1, 0)
             THat = matmul(psiInterpolationMatrix, THatFine)
             dTHatdpsi = matmul(psiInterpolationMatrix, dTHatdpsiFine)
             nHat = matmul(psiInterpolationMatrix, nHatFine)
             ! Re-compute the derivative matrices since they have been over-written above
             call uniformDiffMatrices(NpsiFine, psiMin-(1d-10), psiMax+(1d-10), scheme, psiFine, &
                  psiWeightsFine, ddpsiForT, ddpsiForTCopy)
             dnHatdpsiFine = matmul(ddpsiForT, nHatFine)
             dnHatdpsi = matmul(psiInterpolationMatrix, dnHatdpsiFine)
             deallocate(psiInterpolationMatrix)

             deallocate(ddpsiForTCopy)
             deallocate(psiFine)
             deallocate(PhiHatFine)
             deallocate(dPhiHatdpsiFine)
             deallocate(nHatFine)
             deallocate(dnHatdpsiFine)
             deallocate(THatFine)
             deallocate(dTHatdpsiFine)
             deallocate(etaHatFine)
             deallocate(ddpsiForT)
             deallocate(psiWeightsFine)
             deallocate(rhsForT)
             deallocate(UFine)
             deallocate(UFactor)

          else
             THat = 1 + dTHatdpsiScalar * (psi-psiMid)
             dTHatdpsi = dTHatdpsiScalar
          end if

          etaHat = 1 + (psi-psiMid)*detaHatdpsiScalar
          detaHatdpsi = detaHatdpsiScalar
          if (profilesScheme==4 .or. profilesScheme==6) then
             do i=1,Npsi
                temp1 = (psi(i)-psiMid)*ss
#ifdef USE_GSL_ERF
                call erf(temp1, temp2)
#else
                temp2 = erf(temp1)
#endif
                PhiHat(i) = rampAmplitude*temp2*sqrtpi/(2*ss)
                dPhiHatdpsi(i) = rampAmplitude*exp(-((psi(i)-psiMid)*ss)**2)
             end do
          else
             call computePhiHat_flatURegion(Npsi, psi, psiMid, ss, r0, rampAmplitude, PhiHat)
             call computedPhiHatdpsi_flatURegion(Npsi, psi, psiMid, ss, r0, rampAmplitude, dPhiHatdpsi)
          end if

          if (.not. setTPrimeToBalanceHeatFlux) then
             nHat = etaHat * exp(-2*omega/delta*PhiHat/THat)
             dnHatdpsi = ( detaHatdpsi &
                           - etaHat*2d0*omega/delta*dPhiHatdpsi/THat &
                           + etaHat*2d0*omega/delta*PhiHat*dTHatdpsi/THat**2 ) * exp(-2*omega/delta*PhiHat/THat)
          end if

       case default
          if (masterProcInSubComm) then
             print *,"Error! Invalid setting for profilesScheme"
          end if
          stop
       end select

       ! At this point, we have set BHat, dBHatdtheta, dBHatdpsi, IHat, dIHatdpsi, JHat,
       ! phiHat, dPhiHatdpsi,
       ! and, for one species, THat, dTHatdpsi, and etaHat.

       ! allocate(dnHatdpsi(Npsi))
       ! allocate(detaHatdpsi(Npsi))

       ! if (psiDerivativeScheme == 0) then
       !    !ddpsi_accurate = ddpsiLeft
       !    print *,"Error! This profileScheme is not set up yet to work with psiDerivativeScheme==0"
       !    stop
       ! else
       !    ! centered finite differences, no upwinding, 5-point stencil
       !    scheme = 12
       !    call uniformDiffMatrices(Npsi, psiMin, psiMax, scheme, psi, temp_psiWeights, ddpsi_accurate, temp_d2dpsi2)
       ! end if

       ! if (profilesScheme < 2) then
       !    dPhiHatdpsi = matmul(ddpsi_accurate, PhiHat)
       ! end if

       if (abs(omega) < 1d-13) then
          dPhiHatdpsi = 0
       end if

       ! For now, just use the same profiles for all species, scaled by the appropriate factors.
       ! Later we will want to set up experimental profiles for each species here.
       do ispecies = 1,Nspecies
          THats(ispecies,:) = THat(:)*scalarTHats(ispecies)
          dTHatdpsis(ispecies,:) = dTHatdpsi(:)*scalarTHats(ispecies)
          nHats(ispecies,:) = nHat(:)*scalarNHats(ispecies)
          etaHats(ispecies,:) = nHats(ispecies,:) * exp(charges(ispecies)*2*omega/delta*PhiHat/THats(ispecies,:))
          !dnHatdpsis(ispecies,:) = matmul(ddpsi_accurate, nHats(ispecies,:))
          !detaHatdpsis(ispecies,:) = matmul(ddpsi_accurate, etaHats(ispecies,:))
          dnHatdpsis(ispecies,:) = dnHatdpsi(:)*scalarNHats(ispecies)
          !detaHatdpsis(ispecies,:) = detaHatdpsi(:)*scalarNHats(ispecies)
          detaHatdpsis(ispecies,:) = ( dnHatdpsis(ispecies,:) &
                                       + nHats(ispecies,:)*charges(ispecies) &
                                         *2d0*omega/delta*dPhiHatdpsi/THats(ispecies,:) &
                                       - nHats(ispecies,:)*charges(ispecies) &
                                         *2d0*omega/delta*PhiHat*dTHatdpsis(ispecies,:)/THats(ispecies,:)**2 &
                                     ) * exp(charges(ispecies)*2d0*omega/delta*PhiHat/THats(ispecies,:))
       end do

    end if

  end subroutine initializeProfiles

  ! Fill some arrays that can be computed from the radial physics profiles.
  subroutine computeDerivedProfileQuantities()

    PetscScalar, dimension(:), allocatable :: rSingleSpecies
    integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
    integer :: i, iSpecies
    integer :: scheme
    integer :: LAPACKInfo

    allocate(nuPrimeProfile(Nspecies,Npsi))
    allocate(nuStarProfile(Nspecies,Npsi))
    allocate(deltaN(Nspecies,Npsi))
    allocate(deltaT(Nspecies,Npsi))
    allocate(deltaEta(Nspecies,Npsi))
    allocate(U(Nspecies,Npsi))
    allocate(r(Nspecies,Npsi))
    allocate(rSingleSpecies(Npsi))

    do ispecies = 1,Nspecies

       do i=1,Npsi
          deltaT(ispecies,i) = abs(delta*sqrt(masses(ispecies))*IHat(i)/(psiAHatArray(i)*typicalB(i) &
               *charges(ispecies)*sqrt(THats(ispecies,i)))*dTHatdpsis(ispecies,i))
          deltaN(ispecies,i) = abs(delta*sqrt(masses(ispecies)*THats(ispecies,i))*IHat(i) &
               / (psiAHatArray(i)*charges(ispecies)*typicalB(i)*nHats(ispecies,i)) * dnHatdpsis(ispecies,i))
          deltaEta(ispecies,i) = abs(delta*sqrt(masses(ispecies)*THats(ispecies,i))*IHat(i) &
               / (psiAHatArray(i)*charges(ispecies)*typicalB(i)*etaHats(ispecies,i)) * detaHatdpsis(ispecies,i))
       end do

       do i=1,Npsi
          nuPrimeProfile(ispecies,i) = nu_r * charges(ispecies)**4 &
               * Miller_q * nHats(ispecies,i) / (THats(ispecies,i)*THats(ispecies,i))
          nuStarProfile(ispecies,i) = nuPrimeProfile(ispecies,i) / (epsil*sqrt(epsil))
       end do

       U(ispecies,:) = omega*IHat*dPhiHatdpsi/psiAHatArray(:)*sqrt(masses(ispecies)/(FSABHat2*THats(ispecies,:)))

       ! Next, compute r.
       if (Npsi<5) then
         ! Cannot use 5-point stencil, and r is probably not interesting anyway
         rSingleSpecies = 0d0
       else
         ! Store dr/dpsi in the variable r, since LAPACK will over-write dr/dpsi with r in a few lines:
         rSingleSpecies = psiAHatArray(:) / delta * sqrt(FSABHat2 / THats(ispecies,:)) / IHat

         ! Re-create ddpsi_accurate, since it is over-written in the loop.
         ! centered finite differences, no upwinding, 5-point stencil
         if(leftBoundaryScheme /= 3) then
            ! non-periodic
            scheme = 12
         else
            ! periodic
            scheme = 10
         end if
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
       end if

       r(ispecies,:) = rSingleSpecies
    end do

    deallocate(rSingleSpecies)
    allocate(sqrtTHats(Nspecies,Npsi))
    sqrtTHats = sqrt(THats)

  end subroutine computeDerivedProfileQuantities

  !----------------------------------------------------------------------
  ! The subroutines below are used for simplisitic model profiles.
  !----------------------------------------------------------------------


  subroutine computePhiHat_flatURegion(Npsi, psi, psiMid, ss, r0, rampAmplitude, PhiHat)

    implicit none

    integer, intent(in) :: Npsi
    PetscScalar, intent(in) :: ss, r0, rampAmplitude, psiMid
    PetscScalar, intent(in) :: psi(Npsi)
    PetscScalar, intent(out) :: PhiHat(Npsi)

    integer :: i
    PetscScalar :: valueAt0

    valueAt0 = rampAmplitude * (-exp(r0*ss)*log(-(sqrt(exp(2*r0*ss) - 4)*exp(r0*ss) &
         - 2*exp(r0*ss) - exp(2*r0*ss))/(sqrt(exp(2*r0*ss) - 4)*exp(r0*ss) &
         + 2*exp(r0*ss) + exp(2*r0*ss)))/(sqrt(exp(2*r0*ss) - 4)*ss))

    do i=1,Npsi
       PhiHat(i) = rampAmplitude * (-exp(r0*ss)*log(-(sqrt(exp(2*r0*ss) - 4)*exp(r0*ss) &
            - 2*exp(-(psi(i)-psiMid)*ss + r0*ss) - exp(2*r0*ss))/(sqrt(exp(2*r0*ss) - 4)*exp(r0*ss) &
            + 2*exp(-(psi(i)-psiMid)*ss + r0*ss) + exp(2*r0*ss)))/(sqrt(exp(2*r0*ss) - 4)*ss)) &
            - valueAt0
    end do

  end subroutine computePhiHat_flatURegion

  subroutine computedPhiHatdpsi_flatURegion(Npsi, psi, psiMid, ss, r0, rampAmplitude, dPhiHatdpsi)

    implicit none

    integer, intent(in) :: Npsi
    PetscScalar, intent(in) :: ss, r0, rampAmplitude, psiMid
    PetscScalar, intent(in) :: psi(Npsi)
    PetscScalar, intent(out) :: dPhiHatdpsi(Npsi)

    integer :: i

    do i=1,Npsi
       dPhiHatdpsi(i) = rampAmplitude /(1 + exp(ss*(psi(i)-psiMid-r0)) + exp(-ss*(psi(i)-psiMid+r0)))
    end do

  end subroutine computedPhiHatdpsi_flatURegion

end module
