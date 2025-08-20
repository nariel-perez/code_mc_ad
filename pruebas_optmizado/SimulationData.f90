!---------------------------------------------------------------------
! File: SimulationData.f90 (Refactorizado)
!
! Módulo para la gestión de datos de estado de la simulación.
! Almacena coordenadas, propiedades atómicas y arrays de control.
!---------------------------------------------------------------------
module SimulationData
    use PBC_Mod, only: rk
    
    implicit none
    private

    !======================= VARIABLES PÚBLICAS ======================
    ! Coordenadas de las partículas (allocatable para flexibilidad)
    real(rk), allocatable, public :: RX(:,:,:), RY(:,:,:), RZ(:,:,:)
    
    ! Propiedades de las partículas
    real(rk), allocatable, public :: EPSI(:), SIGM(:), Q(:)
    real(rk), allocatable, public :: RXC(:), RYC(:), RZC(:)
    real, allocatable, public :: QAC(:), EPSAC(:), SGC(:)

    ! Arrays de interacciones y energías
    real(rk), allocatable, public :: UADS(:,:,:,:)
    real(rk), allocatable, public :: USS(:,:,:)
    real, allocatable, public :: CNF(:,:,:,:,:)
    real, allocatable, public :: RX1(:), RY1(:), RZ1(:)

    ! Arrays de tamaño fijo (manteniendo la declaración original)
    ! Se mantuvieron las dimensiones fijas para compatibilidad,
    ! aunque esto puede limitar la escalabilidad.
    real, public :: ANX(5000,10), ANGY(5000,10), ANZ(5000,10)
    integer, public :: LOCATE(5000,10)
    logical, public :: FLAG(1000)
    integer, public :: NMAXI(1000), NMIN(1000)

    ! Variables para el campo de fuerza y control
    real, public :: EXNEW, EYNEW, EZNEW
    real, parameter, public :: TOLERANCIA = 1.0E-5
    integer, allocatable, public :: SYMBOL(:), N(:)
    
end module SimulationData
