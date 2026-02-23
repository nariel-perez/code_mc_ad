!---------------------------------------------------------------
! File: InputParams.f90             (Fortran 90 / 2003)
! Lee input.txt y expone los parámetros de control.
!---------------------------------------------------------------
module InputParams
   implicit none
   !private                        ! Todo privado por defecto
   !------------------------------
   !  PARÁMETROS PÚBLICOS
   !------------------------------
   real, public :: P, dp, sigmetano, eps, ACEL, acelx, acely, acelz, diel, T
   integer, public ::  bcx, bcy, bcz, mat,isot, ijpasos, ikpasos, mult2, ensemble, NESTADO
   integer, public :: canonicalmolecules, NCELLMAT
   character(len= 16), public :: nam
   
contains
   subroutine read_input(fname)
      !----------------------------------------------------------
      ! Lee el archivo de entrada.  Todas las variables quedan
      ! inicializadas y “SAVE’d” dentro del módulo.
      !----------------------------------------------------------
      implicit none
      character(len=*), intent(in) :: fname
      integer :: ios

      ! Unidad de I/O local para evitar choques
      integer, parameter :: iu = 17

      open(unit=iu, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*, '(A,I0)') '*** ERROR: no pude abrir '//trim(fname)//'  iostat=', ios
         stop 1
      end if

      read(iu, *) P            ! Presión inicial   (cm Hg)
      read(iu, *) dp           ! Paso de presión
      read(iu, *) sigmetano   ! σ de metano
      read(iu, *) eps          ! ε  (unidad de energía)
      read(iu,*) ACEL, acelx, acely, acelz
      read(iu,*) diel
      read(iu,*) BCX, BCY, BCZ
      read(iu,*) T
      read(iu,*) mat
      read(iu,*) nam
      read(iu,*) isot
      read(iu,*) ijpasos
      read(iu,*) ikpasos
      read(iu,*) mult2
      read(iu,*) ensemble
      read(iu,*) NCELLMAT
      read(iu,*) NESTADO
      read(iu,*) canonicalmolecules

      close(iu)
   end subroutine read_input


subroutine print_params()
  implicit none
  
  write(*,*) '-----------------------------------------'
  write(*,*) '-----------INPUT PARAMETERS--------------'
  write(*,*) '-----------------------------------------'
  
  write(*,'(A, ES12.5)')  'Initial pressure (cm Hg):',           P
  write(*,'(A, ES12.5)')  'Final  Pressure (dp):',               dp
  write(*,'(A, F10.4)')  'Methane sigma (length reference):',   sigmetano
  write(*,'(A, F10.4)')  'Epsilon (energy reference):',         eps
  write(*,'(A, 4F10.4)') 'Cell size (Ref,x,y,z):',              ACEL, acelx, acely, acelz
  write(*,'(A, F10.4)')  'Dielectric constant:',                diel
  write(*,'(A, 3I10)') 'PBC box lengths (BCX,BCY,BCZ):',      BCX, BCY, BCZ
  write(*,'(A, F10.4)')  'Simulation temperature (T):',         T
  write(*,'(A, I10)')    'Matrix size (mat):',                  mat
  write(*,'(A, A)')      'Surface file name (nam):',            trim(nam)
  write(*,'(A, I10)')    'Points per isotherm (isot):',         isot
  write(*,'(A, I10)')    'ij-steps (averaging):',               ijpasos
  write(*,'(A, I10)')    'ik-steps (config save):',             ikpasos
  write(*,'(A, I10)')    'Multiplier (mult2):',                 mult2
  
  select case (ensemble)
   case (0)
      write(*,'(A)') 'Ensemble: CANONICAL'
   case (1)
      write(*,'(A)') 'Ensemble: GRAND CANONICAL (restart)'
   case (2)
      write(*,'(A)') 'Ensemble: GRAND CANONICAL'
   case default
      write(*,'(A)') 'Ensemble: GRAND CANONICAL (NAMD-style)'
   end select

   if (ensemble == 0) dp = 0.0        ! mismo comportamiento que tenías

   if (NCELLMAT > 50) then
      write(*,*) 'ERROR: STATICAL MATRIX SIZE IS TOO LARGE'
      stop
   end if

   write(*,'(A, I10)') 'Statical matrix size (NCELLMAT):', NCELLMAT
   write(*,'(A, I10)') 'State number (NESTADO):',          NESTADO
   write(*,'(A, I10)') 'Canonical molecules:',             canonicalmolecules

   write(*,*) '-------------------------------------------------'
   write(*,*) '------------END OF PARAMETERS READING------------'
   write(*,*) '-------------------------------------------------'
 end subroutine print_params
end module InputParams
