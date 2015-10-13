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

    namelist / geometryParameters / geometryToUse, epsil, &
         Miller_kappa, Miller_delta, Miller_s_delta, Miller_s_kappa, Miller_dRdr, Miller_q

    namelist / speciesParameters / masses, charges, scalarNHats, scalarTHats

    namelist / physicsParameters / nu_r, profilesScheme, profilesFilename, &
         makeLocalApproximation, desiredU, desiredUMin, desiredUMax, desiredUNumRuns, &
         desiredFWHMInRhoTheta, dTHatdpsiScalar, &
         detaHatdpsiScalar, sourcePoloidalVariation, &
         Delta, omega, psiAHat, psiMid,  &
         exponent, setTPrimeToBalanceHeatFlux, includeddpsiTerm, &
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
         xMaxForDistribution, NxUniform, NxiUniform, xUniformMax

    namelist / otherNumericalParameters / thresh, xScaleFactor, &
         useIterativeSolver, &
         psiDerivativeScheme, thetaDerivativeScheme, xDerivativeScheme, &
         whichParallelSolverToFactorPreconditioner, PETSCPreallocationStrategy

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

       read(fileUnit, nml=physicsParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the physicsParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from physicsParameters namelist in ", filename, "."
       end if

       read(fileUnit, nml=resolutionParameters, iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
               " but not read data from the resolutionParameters namelist in it."
          stop
       end if
       if (masterProc) then
          print *,"Successfully read parameters from resolutionParameters namelist in ", filename, "."
       end if

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

    numCharges = maxNumSpecies
    numDensities = maxNumSpecies
    numMasses = maxNumSpecies
    numTemperatures = maxNumSpecies

    do i=1,maxNumSpecies
       if (charges(i) == speciesNotInitialized) then
          numCharges = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (masses(i) == speciesNotInitialized) then
          numMasses = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
       if (scalarNHats(i) == speciesNotInitialized) then
          numDensities = i-1
          exit
       end if
    end do

    do i=1,maxNumSpecies
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

    numSpecies = numCharges

    ! Other input validation

    if (outputScheme < 0 .or. outputScheme > 2) then
       print *,"Error: outputScheme must be 0, 1, or 2."
       stop
    end if

  end subroutine readNamelistInput

end module readInput

