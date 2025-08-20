!--------------------------------------------------------------------
! File: InputParams.f90 (Refactorizado)
!
! Módulo para la gestión y validación de parámetros de entrada
! para simulaciones computacionales.
!--------------------------------------------------------------------
module InputParams
    use, intrinsic :: iso_fortran_env, only : output_unit
    use PBC_Mod, only : Cell_t => Cell, cell_from_lengths_angles,rk
    implicit none
    private

    !====================== PARÁMETROS NUMÉRICOS =====================
    ! Parámetros de simulación
    real(rk), public :: P = 0.0_rk, dp = 0.0_rk, sigmetano = 0.0_rk, eps = 0.0_rk
    real(rk), public :: diel = 1.0_rk, T = 298.15_rk

    ! Dimensiones de la celda de simulación
    real(rk), public :: a_len = 0.0_rk, b_len = 0.0_rk, c_len = 0.0_rk
    real(rk), public :: alpha_deg = 90.0_rk, beta_deg = 90.0_rk, gamma_deg = 90.0_rk
    real(rk), public :: ACEL = 0.0_rk

    ! Flags de condiciones de contorno y dimensionalidad
    integer, public :: BCX = 1, BCY = 1, BCZ = 1
    integer, public :: dim_flag = 3

    ! Parámetros de control de la simulación
    integer, public :: mat = 0, isot = 0, ijpasos = 0, ikpasos = 0
    integer, public :: mult2 = 0, ensemble = 0, NESTADO = 0, canonicalmolecules = 0, NCELLMAT = 0
    character(len=16), public :: nam = ''

    ! La celda global, definida en el módulo PBC_Mod
    type(Cell_t), public :: cell

    public :: read_input

contains

!====================================================================
! Subrutinas de utilidad para manejar mensajes de la aplicación
!====================================================================
private
subroutine info(msg); character(len=*),intent(in)::msg; write(output_unit,'(A)') trim(msg); end subroutine info
subroutine warn(msg); character(len=*),intent(in)::msg; write(output_unit,'("WARNING: ",A)') trim(msg); end subroutine warn
subroutine fatal(msg); character(len=*),intent(in)::msg; write(output_unit,'("FATAL: ",A)') trim(msg); stop; end subroutine fatal

!====================================================================
! Funciones y subrutinas para el manejo de parámetros
!====================================================================
pure function lowercase(str) result(out)
    ! Convierte una cadena de caracteres a minúsculas.
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i
    do i=1,len(str)
        out(i:i) = adjustl(str(i:i))
        if (iachar(out(i:i)) >= iachar('A') .and. iachar(out(i:i)) <= iachar('Z')) then
            out(i:i) = achar(iachar(out(i:i)) + (iachar('a') - iachar('A')))
        end if
    end do
end function lowercase

private
subroutine set_param(key, value)
    ! Asigna un valor a un parámetro real de forma genérica.
    character(len=*), intent(in) :: key
    real(rk), intent(in) :: value
    select case(lowercase(trim(key)))
    case('p'); P = value
    case('dp'); dp = value
    case('sigmetano'); sigmetano = value
    case('eps'); eps = value
    case('diel'); diel = value
    case('t'); T = value
    case('acelx','a','ax'); a_len = value
    case('acely','b','ay'); b_len = value
    case('acelz','c','az'); c_len = value
    case('angux','alpha','angx','alfx'); alpha_deg = value
    case('anguy','beta','angy','alfy'); beta_deg = value
    case('anguz','gamma','angz','alfz'); gamma_deg = value
    case('acel'); ACEL = value
    case default; call warn('Clave real desconocida: '//trim(key))
    end select
end subroutine set_param

private
subroutine set_param_int(key, value)
    ! Asigna un valor a un parámetro entero de forma genérica.
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    select case(lowercase(trim(key)))
    case('mat'); mat = value
    case('isot'); isot = value
    case('ijpasos'); ijpasos = value
    case('ikpasos'); ikpasos = value
    case('mult2'); mult2 = value
    case('ensemble'); ensemble = value
    case('nestado'); NESTADO = value
    case('canonicalmolecules'); canonicalmolecules = value
    case('ncellmat'); NCELLMAT = value
    case('bcx'); BCX = value
    case('bcy'); BCY = value
    case('bcz'); BCZ = value
    case('dim'); dim_flag = value
    case default; call warn('Clave entera desconocida: '//trim(key))
    end select
end subroutine set_param_int

private
subroutine set_param_char(line)
    ! Asigna un valor a un parámetro de tipo cadena de caracteres.
    character(len=*), intent(in) :: line
    character(len=64) :: k, v
    integer :: ios
    read(line, *, iostat=ios) k, v
    if (ios /= 0) return
    select case(lowercase(trim(k)))
    case('nam', 'nombre'); nam = trim(v)
    case default; call warn('Clave string desconocida: '//trim(k))
    end select
end subroutine set_param_char

!====================================================================
! Subrutinas principales
!====================================================================
subroutine read_input(fname)
    ! Lee parámetros de entrada desde un archivo y los asigna a las variables del módulo.
    !
    ! Argumentos:
    ! fname (in): Nombre del archivo de entrada.
    use, intrinsic :: iso_fortran_env, only : iostat_end
    implicit none
    character(len=*), intent(in) :: fname
    integer :: iu, ios
    character(len=128) :: line
    character(len=64) :: key
    real(rk) :: rval
    integer :: ival
    
    call info('Leyendo '//trim(fname))
    open(newunit=iu, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        call fatal('No se pudo abrir el archivo de entrada: '//trim(fname))
    end if
    
    do
        read(iu, '(A)', iostat=ios) line
        if (ios == iostat_end) exit
        if (ios /= 0) call fatal('Error leyendo línea del archivo')
        
        line = trim(line)
        ! Limpia comentarios y líneas vacías
        if (index(line, '!') > 0) line = line(:index(line, '!') - 1)
        if (len(trim(line)) == 0) cycle
        
        ! Intenta leer primero como entero, luego como real, luego como string
        read(line, *, iostat=ios) key, ival
        if (ios == 0) then
            call set_param_int(key, ival)
        else
            read(line, *, iostat=ios) key, rval
            if (ios == 0) then
                call set_param(key, rval)
            else
                call set_param_char(line)
            end if
        end if
    end do
    close(iu)

    ! Post-procesamiento de parámetros
    if (ACEL == 0.0_rk) then
        ACEL = max(a_len, b_len, c_len)
    else
        ! Si ACEL se especifica, debe ser igual en todas las dimensiones si son 1D o 2D
        a_len = ACEL
        b_len = ACEL
        c_len = ACEL
    end if
    
    call cell_from_lengths_angles(cell, a_len, b_len, c_len, &
                                  alpha_deg, beta_deg, gamma_deg, &
                                  dim=dim_flag, centered=.true.)
    
    cell%pbc = [BCX /= 0, BCY /= 0, BCZ /= 0]
    cell%dim = dim_flag
    
    call print_summary()
end subroutine read_input

private
subroutine print_summary()
    ! Imprime un resumen de los parámetros de entrada leídos.
    implicit none
    character(len=32) :: ens_txt
    select case (ensemble)
    case (0); ens_txt = 'CANONICAL'
    case (1); ens_txt = 'GC (RESTART)'
    case (2); ens_txt = 'GRAND CANONICAL'
    case default; ens_txt = 'GRAND CANONICAL (NAMD-style)'
    end select
    
    write(output_unit, *) '-----------------------------------------'
    write(output_unit, *) '----------- INPUT PARAMETERS ------------'
    write(output_unit, *) '-----------------------------------------'
    
    write(output_unit, '(A,ES12.5)') 'Initial pressure (P): ', P
    write(output_unit, '(A,ES12.5)') 'Final/step pressure (dp): ', dp
    write(output_unit, '(A,F10.4)') 'Methane sigma: ', sigmetano
    write(output_unit, '(A,F10.4)') 'Epsilon: ', eps
    write(output_unit, '(A,F10.4)') 'Cell size (ref): ', ACEL
    write(output_unit, '(A,3F10.4)') 'Cell lengths (a,b,c): ', a_len, b_len, c_len
    write(output_unit, '(A,3F10.4)') 'Angles (alpha,beta,gamma) [deg]: ', alpha_deg, beta_deg, gamma_deg
    write(output_unit, '(A,3I10)') 'PBC (x,y,z): ', BCX, BCY, BCZ
    write(output_unit, '(A,I10)') 'Dimensionality (dim): ', dim_flag
    write(output_unit, '(A,F10.4)') 'Dielectric constant: ', diel
    write(output_unit, '(A,F10.4)') 'Temperature (T): ', T
    write(output_unit, '(A,I10)') 'Matrix size (mat): ', mat
    write(output_unit, '(A,A)') 'Structure file (nam): ', trim(nam)
    write(output_unit, '(A,I10)') 'Isotherm points (isot): ', isot
    write(output_unit, '(A,I10)') 'ij-steps (averaging): ', ijpasos
    write(output_unit, '(A,I10)') 'ik-steps (config save): ', ikpasos
    write(output_unit, '(A,I10)') 'Multiplier (mult2): ', mult2
    write(output_unit, '(A)') 'Ensemble: ' // trim(ens_txt)
    write(output_unit, '(A,I10)') 'Stat. matrix size (NCELLMAT): ', NCELLMAT
    write(output_unit, '(A,I10)') 'State number (NESTADO): ', NESTADO
    write(output_unit, '(A,I10)') 'Canonical molecules: ', canonicalmolecules
    
    write(output_unit, *) '-----------------------------------------'
    write(output_unit, *) '--------- END OF PARAMS READING ---------'
    write(output_unit, *) '-----------------------------------------'
end subroutine print_summary

end module InputParams
