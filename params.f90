module params
  !Iterations to run
  integer :: NSTEPS=10000
  integer :: ISTEPS=2000
  
  !Resolution
  integer :: NX = 64
  integer :: NY = 64
  integer :: NZ = 64
  double precision :: DSPACE = 0.2d0
  double precision :: DTSIZE = 0.01d0

  !Dump frequency: Wavefunction - Misc
  integer :: dumpwf = 100
  integer :: dumputil = 100

  !GPE Type: 0 Natural Units - 1 Hamonic Oscillator Units
  integer :: RHSType = 1
  double precision :: harm_osc_C = 300.0d0
  double precision :: harm_osc_mu = 10.136d0

  !Linearly/Rotating moving frame
  double precision :: VELX = 0.0d0
  double precision :: VELY = 0.0d0
  double precision :: VELZ = 0.0d0
  double precision :: OMEGA = 0.0d0
  
  !Dissipation
  double precision :: GAMMAC = 0.0d0

  !Boundary Conditions: 0 - reflective, 1 - periodic
  integer :: BCX = 1
  integer :: BCY = 1
  integer :: BCZ = 1

  !Potential
  logical :: enableTrap = .true.
  double precision :: TX=0.0d0
  double precision :: TY=0.0d0
  double precision :: TZ=0.0d0
  double precision :: TXSCALE = 1.0d0
  double precision :: TYSCALE = 1.0d0
  double precision :: TZSCALE = 1.0d0

  !Initial condition: 0 - homg, 1 - TF, 2 - restart
  integer :: initialCondType = 0
  character(2048) :: ICRfilename

  !GLOBALS----------------------------------------------------------------------
  double precision,parameter :: PI = 4.0d0*ATAN(1.0d0)
  complex*16 :: DT,EYE = (0.0d0,1.0d0)
  complex*16, dimension(:,:,:), ALLOCATABLE :: GRID,POT
  double precision, dimension(:), ALLOCATABLE :: GX,GY,GZ
  double precision, dimension(:,:,:), ALLOCATABLE :: MGX,MGY,MGZ
  double precision :: TIME

  !Parallel local grid sizes (plus ghost points for halo swapping)
  integer :: PSX,PEX,PSY,PEY,PSZ,PEZ

  contains
  SUBROUTINE init_params
    IMPLICIT NONE
    ICRfilename = repeat(" ", 2048) !Clear memory so entire string is blank
    include 'params.in'
  END SUBROUTINE

  subroutine init_arrays
    implicit none
    ALLOCATE(GRID(PSX:PEX,PSY:PEY,PSZ:PEZ))
    ALLOCATE(POT(PSX:PEX,PSY:PEY,PSZ:PEZ))
    ALLOCATE(GX(PSX:PEX))
    ALLOCATE(GY(PSY:PEY))
    ALLOCATE(GZ(PSZ:PEZ))
    ALLOCATE(MGX(PSX:PEX,PSY:PEY,PSZ:PEZ))
    ALLOCATE(MGY(PSX:PEX,PSY:PEY,PSZ:PEZ))
    ALLOCATE(MGZ(PSX:PEX,PSY:PEY,PSZ:PEZ))
  end subroutine
end module