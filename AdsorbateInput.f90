!─────────────────────────────────────────────────────────────────
! File: AdsorbateInput.f90
! Lee MOLEC.DAT y cada archivo de molécula; deja todo listo
! para el resto del código.
!─────────────────────────────────────────────────────────────────
module AdsorbateInput
  use PBC_Mod, only : rk
  use InputParams, only : sigmetano          ! conversión de unidades
   implicit none
   private
   !====== API pública ==========================================
   public :: read_adsorbates
   !public :: NMOLEC, maxAtoms, XT, X, NATOM, NTOTAL, NMIN, NMAXI, &
   !          RX0, RY0, RZ0, NATOMKIND, NSYM
   !====== variables globales ==================================
   integer, public           :: NMOLEC   = 0      ! nº de especies
   integer, public           :: maxAtoms = 0      ! átomos de la especie + grande
   real(rk),    public           :: XT       = 0.0_rk    ! suma de fracciones
   real(rk),    allocatable, public :: X(:)           ! fracción de cada especie
   integer, allocatable, public :: NATOM(:), NTOTAL(:), NMIN2(:), NMAXI2(:)
   real(rk),    allocatable, public :: RX0(:,:), RY0(:,:), RZ0(:,:)
   integer, allocatable, public :: NATOMKIND(:,:), NSYM(:,:)

contains
!=================================================================
subroutine read_adsorbates(molec_file, sigma, tol)
   !--------------------------------------------------------------
   ! molec_file : nombre de MOLEC.DAT
   ! tol        : tolerancia para comprobar sum(X)=1
   !--------------------------------------------------------------
   character(len=*), intent(in) :: molec_file
   real(rk),            intent(in) :: tol, sigma
   !–– auxiliares ––
   integer, parameter :: iu = 90, iu2 = 92
   integer :: ios, i, j
   ! lectura de cabeceras
   character(len=256), allocatable :: molname(:)
   integer,  allocatable :: nAtomFile(:), confMin(:), confMax(:)
   real(rk),     allocatable :: fracTmp(:)
   ! datos temporales
   integer :: atomsInFile , nconfmin, nconfmax, ikind, ns
   real(rk)    :: ax, x1, y1, z1
   !--------------------------------------------------------------

   !–– PASO 1: leer cabeceras para saber maxAtoms ––––––––––––––––
   open(unit=iu, file=trim(molec_file), status='old', action='read', iostat=ios)
   if (ios /= 0) stop 'ERROR: no se pudo abrir '//trim(molec_file)
   read(iu,*) NMOLEC

   allocate(molname(NMOLEC), nAtomFile(NMOLEC), confMin(NMOLEC), &
            confMax(NMOLEC), fracTmp(NMOLEC))

   maxAtoms = 0
   do i = 1, NMOLEC
      read(iu,'(A)') molname(i)
      open(unit=iu2, file=trim(molname(i)), status='old', action='read', iostat=ios)
      if (ios /= 0) stop 'ERROR: no se pudo abrir '//trim(molname(i))

      read(iu2,*) atomsInFile, ax, nconfmin, nconfmax
      nAtomFile(i) = atomsInFile
      fracTmp(i)   = ax
      confMin(i)   = nconfmin
      confMax(i)   = nconfmax
      if (atomsInFile > maxAtoms) maxAtoms = atomsInFile
      close(iu2)
   end do
   close(iu)

   !–– PASO 2: reservar arrays definitivos –––––––––––––––––––––––
   allocate(X(NMOLEC));               X      = fracTmp
   allocate(NATOM(NMOLEC));           NATOM  = nAtomFile
   allocate(NMIN2(NMOLEC));            NMIN2   = confMin
   allocate(NMAXI2(NMOLEC));           NMAXI2  = confMax
   allocate(NTOTAL(NMOLEC));          NTOTAL = 0

   allocate(RX0(maxAtoms,NMOLEC));    RX0 = 0.0
   allocate(RY0(maxAtoms,NMOLEC));    RY0 = 0.0
   allocate(RZ0(maxAtoms,NMOLEC));    RZ0 = 0.0
   allocate(NATOMKIND(maxAtoms,NMOLEC)); NATOMKIND = 0
   allocate(NSYM(maxAtoms,NMOLEC));      NSYM      = 0

   !–– PASO 3: cargar coordenadas y tipos atómicos –––––––––––––––
   do i = 1, NMOLEC
      open(unit=iu2, file=trim(molname(i)), status='old', action='read', iostat=ios)
      read(iu2,*) atomsInFile, ax, nconfmin, nconfmax   ! cabecera (ya conocida)

      do j = 1, NATOM(i)
         read(iu2,*) x1, y1, z1, ikind, ns
         RX0(j,i)       =  (x1 / sigmetano) * sigma
         RY0(j,i)       = (y1 / sigmetano) * sigma
         RZ0(j,i)       =  (z1 / sigmetano) * sigma
         NATOMKIND(j,i) =  ikind
         NSYM(j,i)      =  ns
      end do
      close(iu2)
   end do

   XT = sum(X)
   if (abs(XT - 1.0) > tol) then
      write(*,'(A,ES12.5)') 'WARNING: sum(X) ≠ 1  → ', XT
   end if

   write(*,*) ''
   write(*,*) '----------- ADSORBATES LOADED ---------------'
   write(*,'(A,I6)')    'NMOLEC          :', NMOLEC
   write(*,'(A,I6)')    'max atoms/specie:', maxAtoms
   write(*,'(A,ES12.5)')'Sum fractions XT:', XT
   write(*,*) '---------------------------------------------'

   ! limpia aux
   deallocate(molname, nAtomFile, confMin, confMax, fracTmp)
end subroutine read_adsorbates
!=================================================================
end module AdsorbateInput
