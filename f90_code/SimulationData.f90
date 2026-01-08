
MODULE SimulationData
    IMPLICIT NONE

    ! BLOCK1: Main Simulation Variables
    
    REAL, allocatable :: RX(:,:,:), RY(:,:,:), RZ(:,:,:)
    INTEGER :: NMOLEC, NATOM(10), N(10)
    REAL, allocatable :: EPSI(:), SIGM(:), Q(:)
    REAL, allocatable :: RX0(:,:), RY0(:,:), RZ0(:,:), NATOMKIND(:,:)
    REAL, allocatable :: tmp_RX0(:,:), tmp_RY0(:,:), tmp_RZ0(:,:)
    REAL, allocatable :: tmp_NATOM(:,:)
    REAL, allocatable :: RX1(:), RY1(:), RZ1(:)
    REAL, allocatable :: RXC(:), RYC(:), RZC(:)
    REAL, allocatable :: QAC(:),EPSAC(:), SGC(:)
    REAL, allocatable :: UADS(:,:,:,:)
    REAL :: ACEL, ACELX, ACELY, ACELZ
    REAL, allocatable :: USS(:,:,:)
    LOGICAL :: FLAG(1000)
    REAL :: BCX, BCY, BCZ
    INTEGER :: MAT
    REAL :: ANX(5000,10), ANGY(5000,10), ANZ(5000,10)
    REAL :: EXNEW, EYNEW, EZNEW
    REAL, PARAMETER:: TOLERANCIA = 1.0E-5

    ! BLOCK2: Stores location information
    INTEGER :: LOCATE(5000,10)

    ! BLOCK3: Stores min/max values
    INTEGER :: NMAXI(1000), NMIN(1000)
    real, allocatable :: CNF(:,:,:,:,:)
    integer, allocatable :: SYMBOL(:)
END MODULE SimulationData
