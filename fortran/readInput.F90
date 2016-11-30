module readInput

  use globalVariables

  implicit none

contains

  subroutine readNamelistInput(filename)

    implicit none

    character(len=*), intent(in) :: filename
    integer :: fileUnit, didFileAccessWork
    integer :: numMasses, numCharges, numDensities, numTemperatures, i

    namelist / flowControl / programMode, saveMatlabOutput, outputScheme, MatlabOutputFilename, &
         outputFilename, parallelizeOverScan, solveSystem

    namelist / geometryParameters / geometryToUse, geometryFilename, epsil, &
         Miller_kappa, Miller_delta, Miller_s_delta, Miller_s_kappa, Miller_dRdr, Miller_q

    namelist / speciesParameters / masses, charges, scalarNHats, scalarTHats
    
    namelist / physicsParameters / nu_r, profilesScheme, profilesFilename, &
         makeLocalApproximation, desiredU, desiredUMin, desiredUMax, desiredUNumRuns, &
         desiredFWHMInRhoTheta, dTHatdpsiScalar, &
         detaHatdpsiScalar, sourcePoloidalVariation, &
         sourcePoloidalVariationStrength, sourcePoloidalVariationPhase, &
         noChargeSource, chargeSourceFilename, noChargeSourceOption, &
         Delta, omega, psiAHat, psiMid,  &
         exponent, setTPrimeToBalanceHeatFlux, &
         includeCollisionOperator, includeddpsiTerm, &
         leftBoundaryShift, rightBoundaryShift, leftBoundaryScheme, rightBoundaryScheme

    namelist / resolutionParameters / forceOddNtheta, &
         NpsiPerDiameter, NpsiMaxFactor, NpsiMinFactor, NpsiNumRuns, &
         psiDiameter, psiDiameterMinFactor, psiDiameterMaxFactor, psiDiameterNumRuns, &
         widthExtender, widthExtenderMin, widthExtenderMax, widthExtenderNumRuns, &
         Ntheta, NthetaMaxFactor, NthetaMinFactor, NthetaNumRuns, &
         Nxi, NxiMaxFactor, NxiMinFactor, NxiNumRuns, &
         NL, NLMaxFactor, NLMinFactor, NLNumRuns, &
         Nx, NxMaxFactor, NxMinFactor, NxNumRuns, &
         xMax, xMaxMaxFactor, xMaxMinFactor, xMaxNumRuns, &
         solverTolerance, solverToleranceMinFactor, solverToleranceMaxFactor, solverToleranceNumRuns, &
         NxPotentialsPerVth, NxPotentialsPerVthMaxFactor, NxPotentialsPerVthMinFactor, NxPotentialsPerVthNumRuns, &
         xMaxForDistribution, NxUniform, NxiUniform, xUniformMax, &
         thetaGridShift

    namelist / otherNumericalParameters / thresh, xScaleFactor, &
         useIterativeSolver, useIterativeBoundarySolver, &
         psiDerivativeScheme, thetaDerivativeScheme, xDerivativeScheme, &
         psiGridType, psiAHatFilename, NpsiSourcelessRight, NpsiSourcelessLeft, &
	 useGlobalTermMultiplier, globalTermMultiplierFilename, &
         whichParallelSolverToFactorPreconditioner, PETSCPreallocationStrategy, Nxi_for_x_option

    namelist / preconditionerOptions / preconditioner_species, &
         preconditioner_x, preconditioner_psi, &
         preconditioner_theta, preconditioner_xi, preconditioner_x_min_L

    masses = speciesNotInitialized
    charges = speciesNotInitialized
    scalarNHats = speciesNotInitialized
    scalarTHats = speciesNotInitialized

    fileUnit=11
    open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
       print *,"Proc ",myRank,": Error opening ", filename
       stop
    else
       read(fileUnit, nml=flowControl, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the flowControl namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from flowControl namelist in ", filename, "."
       end if

       read(fileUnit, nml=geometryParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the geometryParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from geometryParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=speciesParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the speciesParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from speciesParameters namelist in ", filename, "."
       end if

       ! Default phase and poloidal variation
       sourcePoloidalVariationStrength = 1.0
       sourcePoloidalVariationPhase = 0.0

       ! Default option for enforcing no charge sources
       noChargeSource = 0
       noChargeSourceOption = 0

       ! include collision operator by default
       includeCollisionOperator= .true.
       
       read(fileUnit, nml=physicsParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the physicsParameters namelist in it."
          stop
       end if
       if (((leftBoundaryScheme .eq. 3) .or. (rightBoundaryScheme .eq. 3)) .and.  (leftBoundaryScheme /= rightBoundaryScheme)) then
          print *,"Error: If left or right boundary uses periodic boundary conditions , both boundaries must us periodic."
          stop
       end if
       
       if (masterProc) then
          print *,"Successfully read parameters from physicsParameters namelist in ", filename, "."
       end if

       ! default shift of theta grid
       thetaGridShift = 0.0
       
       read(fileUnit, nml=resolutionParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the resolutionParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from resolutionParameters namelist in ", filename, "."
       end if

       NpsiSourcelessLeft = 0
       NpsiSourcelessRight = 0
       useGlobalTermMultiplier = 0       

       read(fileUnit, nml=otherNumericalParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the otherNumericalParameters namelist in it."
          stop
       end if
       
       if (masterProc) then
          print *,"Successfully read parameters from otherNumericalParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=preconditionerOptions, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the preconditionerOptions namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from preconditionerOptions namelist in ", filename, "."
       end if

    end if
    close(unit = fileUnit)

    ! Validate species parameters

    numCharges = maxNspecies
    numDensities = maxNspecies
    numMasses = maxNspecies
    numTemperatures = maxNspecies

    do i=1,maxNspecies
       if (charges(i) == speciesNotInitialized) then
          numCharges = i-1
          exit
       end if
    end do

    do i=1,maxNspecies
       if (masses(i) == speciesNotInitialized) then
          numMasses = i-1
          exit
       end if
    end do

    do i=1,maxNspecies
       if (scalarNHats(i) == speciesNotInitialized) then
          numDensities = i-1
          exit
       end if
    end do

    do i=1,maxNspecies
       if (scalarTHats(i) == speciesNotInitialized) then
          numTemperatures = i-1
          exit
       end if
    end do

    if (numCharges /= numMasses) then
       print *,"Error: number of species charges differs from the number of species masses."
       stop
    end if

    if (numCharges /= numDensities) then
       print *,"Error: number of species charges differs from the number of species densities."
       stop
    end if

    if (numCharges /= numTemperatures) then
       print *,"Error: number of species charges differs from the number of species temperatures."
       stop
    end if

    Nspecies = numCharges

    if (.not. useIterativeSolver) then
      useIterativeBoundarySolver = .false.
    end if

    ! Other input validation

    if (outputScheme < 0 .or. outputScheme > 2) then
       print *,"Error: outputScheme must be 0, 1, or 2."
       stop
    end if

    if (thetaGridShift < 0 .or. thetaGridShift >= 1) then
       print *,"Error: thetaGridShift should be in [0,1),"
       print *,"where 0 corresponds to no shift (default),"
       print *,"and 1 corresponds to shifting the grid an entire gridpoint (2pi/(Ntheta+1))."
       stop
    end if
    
    if (psiGridType .eq. 1) then
       ! Check psiAHatFilename is set
       if (.not. len(psiAHatFilename)>=0) then
          print *,"If psiGridType==1 then psiAHatFilename must be set."
          stop
       end if       
    end if

    if (useGlobalTermMultiplier .eq. 1) then
       ! Check if globalTermMultiplierFilename is set
       if (.not. len(globalTermMultiplierFilename)>=0) then
          print *,"If useGlobalTermMultiplier==1 then globalTermMultiplierFilename must be set."
          stop
       end if       
    end if

    if (noChargeSource == 2 .or. noChargeSource == 3) then
       ! Check if chargeSourceFilename is set
       if (.not. len(chargeSourceFilename)>=0) then
          print *,"If noChargeSource==2 or 3 then chargeSourceFilename must be set."
          stop
       end if       
    end if

    
  end subroutine readNamelistInput

end module readInput

