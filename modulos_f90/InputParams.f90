!---------------------------------------------------------------
! File: InputParams.f90             (Fortran 90 / 2003)
! Lee input.txt (formato clave-valor) y expone los parámetros.
! Construye la celda de simulación (Cell) con soporte para
! geometría triclínica.
!---------------------------------------------------------------
module InputParams
   use PBC_Mod, only : Cell_t => Cell, cell_from_lengths_angles, rk
   implicit none
   private

   !------------------------------
   !  PARÁMETROS PÚBLICOS
   !------------------------------
   ! Parámetros numéricos
   real(rk), public :: P = 0._rk, dp = 0._rk, sigmetano = 0._rk, eps = 0._rk
   real(rk), public :: diel = 1._rk, T = 298.15_rk

   ! Longitudes de celda (target para permitir pointers)
   real(rk), public, target :: acelx = 0._rk, acely = 0._rk, acelz = 0._rk
   real(rk), public :: ACEL = 0._rk

   ! Ángulos de celda (por defecto ortorrómbico: 90°)
   real(rk), public, target :: angux = 90._rk, anguy = 90._rk, anguz = 90._rk

   ! Aliases para compatibilidad con nomenclatura cristalográfica
   real(rk), public, pointer :: a_len => acelx, b_len => acely, c_len => acelz
   real(rk), public, pointer :: alpha_deg => angux, beta_deg => anguy, gamma_deg => anguz

   ! Flags PBC y dimensionalidad
   integer, public :: BCX = 1, BCY = 1, BCZ = 1   ! 1 = PBC activado
   integer, public :: dim_flag = 3                 ! 3-D por defecto

   ! Otros controles (enteros)
   integer, public :: mat = 0, isot = 0, ijpasos = 0, ikpasos = 0, mult2 = 0
   integer, public :: ensemble = 0, NESTADO = 0, canonicalmolecules = 0, NCELLMAT = 0

   ! Nombre del archivo de superficie
   character(len=64), public :: nam = ''

   ! La celda global de simulación
   type(Cell_t), public :: cell

   ! Celda en unidades reducidas (dividida por ACEL)
   ! Para uso en módulos de movimiento que trabajan en unidades reducidas
   type(Cell_t), public :: cellR

   ! Subrutinas públicas
   public :: read_input, print_params, update_cellR

contains

!====================================================================
! Lectura de input en formato clave-valor
!====================================================================
subroutine read_input(fname)
   use, intrinsic :: iso_fortran_env, only : iostat_end
   implicit none
   character(len=*), intent(in) :: fname
   integer :: iu, ios
   character(len=256) :: line, key
   real(rk) :: rval
   integer  :: ival

   call info('Leyendo '//trim(fname))

   open(newunit=iu, file=fname, status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call fatal('No se pudo abrir '//trim(fname))
   end if

   ! Leer línea por línea
   do
      read(iu, '(A)', iostat=ios) line
      if (ios == iostat_end) exit
      if (ios /= 0) call fatal('Error leyendo línea')

      ! Eliminar comentarios (todo después de !)
      if (index(line, '!') > 0) line = line(:index(line, '!') - 1)

      ! Saltar líneas vacías
      if (len_trim(line) == 0) cycle

      ! Intentar leer como entero primero
      read(line, *, iostat=ios) key, ival
      if (ios == 0) then
         call set_int(lowercase(trim(key)), ival)
      else
         ! Intentar como real
         read(line, *, iostat=ios) key, rval
         if (ios == 0) then
            call set_real(lowercase(trim(key)), rval)
         else
            ! Intentar como cadena
            call set_char(line)
         end if
      end if
   end do

   close(iu)

   ! ACEL por defecto: máximo de las longitudes
   if (ACEL == 0._rk) ACEL = max(acelx, acely, acelz)

   ! Construir la celda con longitudes y ángulos
   call cell_from_lengths_angles(cell, acelx, acely, acelz, &
                                 angux, anguy, anguz, &
                                 dim=dim_flag, centered=.true.)

   ! Aplicar flags PBC y dimensionalidad a la celda
   cell%pbc = [ BCX /= 0, BCY /= 0, BCZ /= 0 ]
   cell%dim = dim_flag

   ! Construir celda en unidades reducidas (para módulos de movimiento)
   call update_cellR()

   ! Imprimir resumen y validar
   call print_params()

end subroutine read_input

!====================================================================
! Actualizar celda en unidades reducidas
! cellR tiene longitudes = cell / ACEL, ángulos iguales
!====================================================================
subroutine update_cellR()
   implicit none

   ! Construir cellR con longitudes reducidas
   call cell_from_lengths_angles(cellR, &
                                 acelx / ACEL, acely / ACEL, acelz / ACEL, &
                                 angux, anguy, anguz, &
                                 dim=dim_flag, centered=.true.)

   ! Copiar flags PBC y dimensionalidad
   cellR%pbc = cell%pbc
   cellR%dim = cell%dim

end subroutine update_cellR

!====================================================================
! Imprimir resumen de parámetros
!====================================================================
subroutine print_params()
   use, intrinsic :: iso_fortran_env, only : output_unit
   implicit none
   integer, parameter :: iu = output_unit
   character(len=32) :: ens_txt

   select case (ensemble)
   case (0); ens_txt = 'CANONICAL'
   case (1); ens_txt = 'GC (RESTART)'
   case (2); ens_txt = 'GRAND CANONICAL'
   case default; ens_txt = 'GRAND CANONICAL (NAMD-style)'
   end select

   write(iu,*) '-----------------------------------------'
   write(iu,*) '----------- INPUT PARAMETERS ------------'
   write(iu,*) '-----------------------------------------'

   write(iu,'(A,ES12.5)') 'Initial pressure (P):              ', P
   write(iu,'(A,ES12.5)') 'Final/step pressure (dp):          ', dp
   write(iu,'(A,F10.4)')  'Methane sigma:                     ', sigmetano
   write(iu,'(A,F10.4)')  'Epsilon:                           ', eps
   write(iu,'(A,4F10.4)') 'Cell size (Ref,x,y,z):             ', ACEL, acelx, acely, acelz
   write(iu,'(A,3F10.4)') 'Angles (a,b,g) [deg]:              ', angux, anguy, anguz
   write(iu,'(A,3I5)')    'BC flags (BCx,BCy,BCz):            ', BCX, BCY, BCZ
   write(iu,'(A,I5)')     'Dimensionality:                    ', dim_flag
   write(iu,'(A,F10.4)')  'Dielectric constant:               ', diel
   write(iu,'(A,F10.4)')  'Temperature (T):                   ', T
   write(iu,'(A,I10)')    'Matrix size (mat):                 ', mat
   write(iu,'(A,A)')      'Structure file (nam):              ', trim(nam)
   write(iu,'(A,I10)')    'Isotherm points (isot):            ', isot
   write(iu,'(A,I10)')    'ij-steps (averaging):              ', ijpasos
   write(iu,'(A,I10)')    'ik-steps (config save):            ', ikpasos
   write(iu,'(A,I10)')    'Multiplier (mult2):                ', mult2
   write(iu,'(A)')        'Ensemble:                          '//trim(ens_txt)
   write(iu,'(A,I10)')    'Stat. matrix size (NCELLMAT):      ', NCELLMAT
   write(iu,'(A,I10)')    'State number (NESTADO):            ', NESTADO
   write(iu,'(A,I10)')    'Canonical molecules:               ', canonicalmolecules

   write(iu,*) '-----------------------------------------'
   write(iu,*) '--------- END OF PARAMS READING ---------'
   write(iu,*) '-----------------------------------------'

   ! Validaciones
   if (ensemble == 0) dp = 0.0_rk

   if (NCELLMAT > 50) then
      write(*,*) 'ERROR: STATICAL MATRIX SIZE IS TOO LARGE'
      stop
   end if

end subroutine print_params

!====================================================================
! Asignar valor real según clave
!====================================================================
subroutine set_real(k, v)
   character(len=*), intent(in) :: k
   real(rk), intent(in) :: v

   select case(k)
   case('p'); P = v
   case('dp'); dp = v
   case('sigmetano'); sigmetano = v
   case('eps'); eps = v
   case('diel'); diel = v
   case('t'); T = v
   case('acelx', 'a', 'ax'); acelx = v
   case('acely', 'b', 'ay'); acely = v
   case('acelz', 'c', 'az'); acelz = v
   case('angux', 'alpha', 'angx'); angux = v
   case('anguy', 'beta', 'angy'); anguy = v
   case('anguz', 'gamma', 'angz'); anguz = v
   case('acel'); ACEL = v
   case default
      call warn('Clave real desconocida: '//k)
   end select
end subroutine set_real

!====================================================================
! Asignar valor entero según clave
!====================================================================
subroutine set_int(k, v)
   character(len=*), intent(in) :: k
   integer, intent(in) :: v

   select case(k)
   case('mat'); mat = v
   case('isot'); isot = v
   case('ijpasos'); ijpasos = v
   case('ikpasos'); ikpasos = v
   case('mult2'); mult2 = v
   case('ensemble'); ensemble = v
   case('nestado'); NESTADO = v
   case('canonicalmolecules'); canonicalmolecules = v
   case('ncellmat'); NCELLMAT = v
   case('bcx'); BCX = v
   case('bcy'); BCY = v
   case('bcz'); BCZ = v
   case('dim'); dim_flag = v
   case default
      call warn('Clave entera desconocida: '//k)
   end select
end subroutine set_int

!====================================================================
! Asignar valor cadena según clave
!====================================================================
subroutine set_char(line)
   character(len=*), intent(in) :: line
   character(len=64) :: k, v
   integer :: ios

   read(line, *, iostat=ios) k, v
   if (ios /= 0) return

   select case(lowercase(trim(k)))
   case('nam', 'nombre', 'surface')
      nam = v
   case default
      call warn('Clave string desconocida: '//k)
   end select
end subroutine set_char

!====================================================================
! Convertir cadena a minúsculas
!====================================================================
pure function lowercase(str) result(out)
   character(len=*), intent(in) :: str
   character(len=len(str)) :: out
   integer :: i, ia

   do i = 1, len(str)
      ia = iachar(str(i:i))
      if (ia >= iachar('A') .and. ia <= iachar('Z')) then
         out(i:i) = achar(ia + 32)
      else
         out(i:i) = str(i:i)
      end if
   end do
end function lowercase

!====================================================================
! Mensajes de utilidad
!====================================================================
subroutine info(msg)
   character(len=*), intent(in) :: msg
   write(*,'(A)') trim(msg)
end subroutine info

subroutine warn(msg)
   character(len=*), intent(in) :: msg
   write(*,'("WARNING: ",A)') trim(msg)
end subroutine warn

subroutine fatal(msg)
   character(len=*), intent(in) :: msg
   write(*,'("FATAL: ",A)') trim(msg)
   stop
end subroutine fatal

end module InputParams
