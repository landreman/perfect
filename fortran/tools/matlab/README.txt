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


Extrapolating geometry for outer buffer region
----------------------------------------------

PERFECT needs closed flux surfaces. When modelling an experiment the buffer region may need to
extend beyond psiN=1, so we need to extrapolate in a way that does not have an X-point.

Chosen algorithm: 
- along lines of constant poloidal angle (from the magnetic axis) extrapolate psiN
  linearly outside of some chosen value, extrapolateBeyondPsiN in the EFITOptions.m file.
- get the gradient for the linear extrapolation from some psiN range, extrapolatePsiNInterval in the
  EFITOptions.m file.

Implementation:
- create a fine r-theta grid (1024 points in each direction) centred on the magnetic axis, with a
  large enough maximum r to include the entire cartesian grid in the geqdsk file.
- interpolate psiN onto the r-theta grid using interp2
- for each theta:
  - find the largest psiN below extrapolateBeyondPsiN
  - extrapolate linearly for larger r
- interpolate back to the cartesian grid using interp2
- recalculate psi from psiN so that BR, BZ, etc. get calculated using the extrapolated field.
- set I to constant from extrapolateBeyondPsiN outwards.
- continue with the usual flux surface reconstruction, but now psiMax can be greater than 1


Non-uniform psiN grids
----------------------

It is useful, especially for global PERFECT runs, to have a non-uniform grid in psiN. This is now
implemented in the matlab geometry generation scripts.
The array of psiN values on which to output the geometry is specified in an HDF5 file, whose name is
input in the psiNFile option in EFITOptions.m. The array is called psiNArray in a group labelled by
the value of Npsi, e.g. 'Npsi100'.
