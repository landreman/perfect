module readInput

  use globalVariables

  implicit none

contains

  subroutine readNamelistInput(filename)

    implicit none

    character(len=*), intent(in) :: filename
    integer :: fileUnit, didFileAccessWork
    integer :: numMasses, numCharges, numDensities, numTemperatures, i
    integer :: NsourcesTemp

    namelist / flowControl / programMode, saveMatlabOutput, outputScheme, MatlabOutputFilename, &
         outputFilename, parallelizeOverScan, solveSystem

    namelist / geometryParameters / geometryToUse, geometryFilename, epsil, &
         Miller_kappa, Miller_delta, Miller_s_delta, Miller_s_kappa, Miller_dRdr, Miller_q

    namelist / speciesParameters / masses, charges, scalarNHats, scalarTHats

    namelist / physicsParameters / nu_r, profilesScheme, profilesFilename, &
         makeLocalApproximation, desiredU, desiredUMin, desiredUMax, desiredUNumRuns, &
         desiredFWHMInRhoTheta, dTHatdpsiScalar, &
         detaHatdpsiScalar, &
         sourcePoloidalVariationStrength, sourcePoloidalVariationPhase, &
         Delta, omega, psiAHat, psiMid,  &
         exponent, setTPrimeToBalanceHeatFlux, &
         includeCollisionOperator, includeddpsiTerm, &
         leftBoundaryShift, rightBoundaryShift, leftBoundaryScheme, rightBoundaryScheme, &
         gConstraints, sourcesVStructure, sourcesThetaStructure, &
         sourceConstraints, RHSFromFile, sourceConstraintsFilenames, &
         extraSourcesVStructure, extraSourcesThetaStructure,extraSourcesSpeciesStructure, &
         miscSources, miscSourcesStrength, &
         constantSourcesVStructure, constantSourcesThetaStructure, constantSourcesFilenames

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

    namelist / otherNumericalParameters / thresh, threshholdForInclusion, xScaleFactor, &
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

    gConstraints = sourcesNotInitialized
    sourcesVStructure = sourcesNotInitialized
    sourcesThetaStructure = sourcesNotInitialized
    sourceConstraints = sourcesNotInitialized
    RHSFromFile = sourcesNotInitialized
    
    ! used to enforce constraints on particular source type
    iparticleSource = sourcesNotInitialized
    imomentumSource = sourcesNotInitialized
    iheatSource = sourcesNotInitialized
    
    extraSourcesVStructure = sourcesNotInitialized
    extraSourcesThetaStructure = sourcesNotInitialized
    extraSourcesSpeciesStructure = sourcesNotInitialized
    miscSources = sourcesNotInitialized
    miscSourcesStrength = sourcesNotInitialized
    constantSourcesVStructure = sourcesNotInitialized
    constantSourcesThetaStructure = sourcesNotInitialized

    do i=1,maxNsources
       sourceConstraintsFilenames = filenameNotInitialized
       constantSourcesFilenames = filenameNotInitialized
    end do

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

       threshholdForInclusion = 1d-12 ! old default

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

    ! Validate sources input
    Nsources = maxNsources
    do i=1,maxNsources
       if (gConstraints(i) == sourcesNotInitialized) then
          Nsources = i-1
          exit
       end if
    end do

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (sourcesVStructure(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= Nsources) then
       print *,"Error: number of sourcesVStructure differs from the number of constraints."
       stop
    end if

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (sourcesThetaStructure(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= Nsources) then
       print *,"Error: number of sourcesThetaStructure differs from the number of constraints."
       stop
    end if

    NextraSources = maxNsources
    do i=1,maxNsources
       if (sourceConstraints(i) == sourcesNotInitialized) then
          NextraSources = i-1
          exit
       end if
    end do

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (RHSFromFile(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NextraSources) then
       print *,"Error: number of RHSFromFile differs from the number of extra sources."
       stop
    end if

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (sourceConstraintsFilenames(i) == filenameNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NextraSources) then
       print *,"Error: number of sourceConstraintsFilenames differs from the number of extra sources."
       stop
    end if

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (extraSourcesVStructure(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NextraSources) then
       print *,"Error: number of extraSourcesVStructure differs from the number of extra sources."
       stop
    end if

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (extraSourcesThetaStructure(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NextraSources) then
       print *,"Error: number of extraSourcesThetaStructure differs from the number of extra sources."
       stop
    end if

    
    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (extraSourcesSpeciesStructure(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NextraSources) then
       print *,"Error: number of extraSourcesSpeciesStructure differs from the number of extra sources."
       stop
    end if


    NmiscSources = maxNsources
    do i=1,maxNsources
       if (miscSources(i) == sourcesNotInitialized) then
          NmiscSources = i-1
          exit
       end if
    end do

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (miscSourcesStrength(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NmiscSources) then
       print *,"Error: number of miscSourcesStrength differs from the number of misc sources."
       stop
    end if

    NconstantSources = maxNsources
    do i=1,maxNsources
       if (constantSourcesVStructure(i) == sourcesNotInitialized) then
          NconstantSources = i-1
          exit
       end if
    end do

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (constantSourcesThetaStructure(i) == sourcesNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NconstantSources) then
       print *,"Error: number of constantSourcesThetaStructure differs from the number of constant sources."
       stop
    end if

    NsourcesTemp = maxNsources
    do i=1,maxNsources
       if (constantSourcesFilenames(i) == filenameNotInitialized) then
          NsourcesTemp = i-1
          exit
       end if
    end do

    if (NsourcesTemp /= NconstantSources) then
       print *,"Error: number of constantSourcesFilenames differs from the number of constant sources."
       stop
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

    if (printReadInDebug == 1) then
       print *,"# sources: ",Nsources
       print *, gConstraints(1:Nsources)
       print *, sourcesVStructure(1:Nsources)
       print *, sourcesThetaStructure(1:Nsources)

       print *,"# extra sources: ",NextraSources 
       print *, sourceConstraints(1:NextraSources)
       print *, RHSFromFile(1:NextraSources)
       print *, sourceConstraintsFilenames(1:NextraSources)
       print *, extraSourcesVStructure(1:NextraSources)
       print *, extraSourcesThetaStructure(1:NextraSources)
       print *, extraSourcesSpeciesStructure(1:NextraSources)

       print *,"# misc sources: ",NmiscSources
       print *, miscSources(1:NmiscSources)
       print *, miscSourcesStrength(1:NmiscSources)

       print *,"# constant sources: ",NconstantSources
       print *, constantSourcesVStructure(1:NconstantSources)
       print *, constantSourcesFilenames(1:NconstantSources)
    end if
  end subroutine readNamelistInput

end module readInput

