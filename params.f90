module params
  !Iterations to run
  integer :: NSTEPS=10000
  integer :: ISTEPS=0
  
  !Resolution
  integer :: NX = 64
  integer :: NY = 64
  integer :: NZ = 64
  double precision :: DSPACE = 0.2d0
  double precision :: DTSIZE = 0.01d0

  !Dump frequency: Wavefunction - Misc
  integer :: DUMPWF = 100
  integer :: DUMPUTIL = 100

  !GPE Type: 0 Natural Units - 1 Hamonic Oscillator Units - 2 Natural Units GPE with DDI (Split-step method)
  integer :: RHSType = 1
  double precision :: harm_osc_C = 300.0d0
  double precision :: harm_osc_mu = 10.136d0
  double precision :: ENERV = 0.75
  double precision :: NV = 0.75
  double precision :: EDD = -0.49d0

  !Linearly/Rotating moving frame
  double precision :: VELX = 0.0d0
  double precision :: VELY = 0.0d0
  double precision :: VELZ = 0.0d0
  double precision :: OMEGA = 0.0d0
  double precision :: DOMEGADT = 0.000d0
  
  !Dissipation
  double precision :: GAMMAC = 0.0d0

  !Boundary Conditions: 0 - reflective, 1 - periodic
  integer :: BCX = 1
  integer :: BCY = 1
  integer :: BCZ = 1

  !Potentials
  logical :: recalculatePot = .false.

  logical :: enableTrap = .false.
  double precision :: TX=0.0d0
  double precision :: TY=0.0d0
  double precision :: TZ=0.0d0
  double precision :: TXSCALE = 1.0d0
  double precision :: TYSCALE = 1.0d0
  double precision :: TZSCALE = 1.0d0

  logical :: enablePins = .false.
  integer :: numPotPins = 12000
  integer :: randSeedPins = 0
  double precision :: pinPotDepth = 100.0d0
  double precision :: pinSigma = 0.2d0

  !Initial condition: 0 - homg, 1 - TF, 2 - restart
  integer :: initialCondType = 0
  character(2048) :: ICRfilename
  integer :: INITSSTEP = 0
  integer :: RESUMESTEP = 0
  double precision :: RESUMETIME = 0

  !GLOBALS----------------------------------------------------------------------
  double precision,parameter :: PI = 4.0d0*ATAN(1.0d0)
  complex*16 :: DT,EYE = (0.0d0,1.0d0)

  contains
  subroutine init_params
    IMPLICIT NONE
    ICRfilename = repeat(" ", 2048) !Clear memory so entire string is blank
    include 'params.in'
  END subroutine
end module