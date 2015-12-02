How to use generateGeometryFromEFIT_script.m:
i) Make sure $PERFECT_DIR/fortran/generateGeometryFromEFIT.sh is in your $PATH
  somehow.
ii) Be in a directory containing:
  a) an input.namelist (that specifies a range of psi values between 0 and 1)
  b) an EFIT g-file
  c) a file EFITOptions.m containing the variables 'EFITFilename' (which should 
    be a string with the name of the g-file), 'topCropZ' and 'bottomCropZ' (which are
    numbers specifying the maximum and minimum heights (in meters) to use while 
    processing the EFIT file), 'polynomialFitDegreeForSmoothingEFITInPsi' (an
    integer; 4 may be appropriate) and
    'numFourierModesInThetaToKeepInEFITGeometry' (an integer; 5 may be appropriate).
iii) run generateGeometryFromEFIT.sh. It will output a file named according to
the 'geometryFilename' option in input.namelist
iv) Make sure geometryToUse is set to 4 in input.namelist.
v) If you change any of the options 'psiDerivativeScheme', 'psiMid',
  'psiDiameter', 'widthExtender', 'leftBoundaryShift', 'rightBoundaryShift',
  'NpsiPerDiameter', 'Ntheta', 'thetaDerivativeScheme', or 'geometryFilename' in
  input.namelist, then re-run generateGeometryFromEFIT.sh since the geometry may
  have changed.

