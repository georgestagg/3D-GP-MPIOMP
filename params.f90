module params
  !Iterations to run
  integer :: NSTEPS=10000
  integer :: ISTEPS=0
  
  !Resolution
  integer :: NX = 64
  integer :: NY = 64
  integer :: NZ = 64
  double precision :: XSHIFT = 0.0d0
  double precision :: YSHIFT = 0.0d0
  double precision :: ZSHIFT = 0.0d0

  double precision :: DSPACE = 0.2d0
  double precision :: DTSIZE = 0.01d0

  !Dump frequency: Wavefunction - Misc
  integer :: DUMPWF = 100
  integer :: DUMPUTIL = 100

  !GPE Type: 0 Natural Units - 1 Hamonic Oscillator Units - 2 Natural Units GPE with DDI (Split-step method)
  integer :: RHSType = 1
  double precision :: harm_osc_C = 300.0d0
  double precision :: harm_osc_mu = 10.136d0
  double precision :: ENERV = 0.75d0
  double precision :: NV = 0.75d0
  double precision :: EDD = 0.0d0
  double precision :: KD = 1.0d0
  double precision :: H = 0.0d0

  !Multi-component parameters
  integer :: FLUIDS = 1
  double precision, dimension(:,:),ALLOCATABLE :: GG
  logical :: INC_MAG_FIELDS = .false.

  !Linearly/Rotating moving frame
  double precision :: VELX = 0.0d0
  double precision :: VELY = 0.0d0
  double precision :: VELZ = 0.0d0
  double precision :: OMEGA = 0.0d0
  double precision :: DOMEGADT = 0.000d0
  double precision :: DHDT = 0.000d0
  
  !Dissipation
  double precision :: GAMMAC = 0.0d0

  !Boundary Conditions: 0 - reflective, 1 - periodic
  integer :: BCX = 1
  integer :: BCY = 1
  integer :: BCZ = 1
  double precision, dimension(:,:),ALLOCATABLE :: NVORTALL
  double precision, dimension(3) :: NVORT = (/ 0.0d0, 0.0d0, 0.0d0 /)

  !Potentials
  logical :: recalculatePot = .false.

  logical :: enableTrap = .false.
  logical :: enablePot = .false.
  integer :: potType = 0
  integer :: trapType = 0

  double precision :: OBJHEIGHT = 0.0d0
  double precision :: OBJX = 0.0d0
  double precision :: OBJY = 0.0d0
  double precision :: OBJZ = 0.0d0
  double precision :: OBJXSCALE = 2.0d0
  double precision :: OBJYSCALE = 2.0d0
  double precision :: OBJZSCALE = 2.0d0

  double precision :: TRAPHEIGHT = 0.0d0
  double precision :: TRAPR = 0.0d0
  double precision :: TRAPBETA = 0.5d0
  double precision :: TX=0.0d0
  double precision :: TY=0.0d0
  double precision :: TZ=0.0d0
  double precision :: TXSCALE = 1.0d0
  double precision :: TYSCALE = 1.0d0
  double precision :: TZSCALE = 1.0d0

  integer :: initialCondType = 0


  character(2048) :: ICRfilename
  character(2048) :: SURFfilename
  integer :: INITSTEP = 0

  !GLOBALS----------------------------------------------------------------------
  double precision,parameter :: PI = 4.0d0*ATAN(1.0d0)
  complex*16 :: DT,EYE = (0.0d0,1.0d0)
  integer :: NGHOST = 1
  integer :: da_stat

  contains
  subroutine init_params
    IMPLICIT NONE
    ICRfilename = repeat(" ", 2048) !Clear memory so entire string is blank
    SURFfilename = repeat(" ", 2048)
    call set_fluid_interactions(1,reshape((/ 1.0d0 /), (/ 1, 1 /)))
    include 'params.in'
    if (FLUIDS == 1) then
      ALLOCATE(NVORTALL(3,1))
      NVORTALL(:,1) = NVORT
    end if
  END subroutine

  subroutine set_fluid_interactions(N,FG)
    IMPLICIT NONE
    integer :: N,stat
    double precision :: FG(:,:)
    DEALLOCATE(GG,STAT=stat)
    ALLOCATE(GG(N,N))
    GG = FG
    FLUIDS = N
  END subroutine

end module