!---------------------------------------------------------------------
! File: Estructura.f90 (Refactorizado)
!
! Lectura y procesamiento de las coordenadas de una superficie cristalina
! (por ejemplo, grafito) para su uso en simulaciones con condiciones
! de contorno periódicas (PBC). Se encarga de filtrar los átomos
! dentro de la celda de simulación y de convertir sus propiedades
! a unidades reducidas.
!---------------------------------------------------------------------
module EstructuraModule
    use InputParams, only: cell, sigmetano, acelx, acely, acelz
    use PhysicalConstants, only: FACTORELEC, FCLEC
    use SimulationData, only: RXC, RYC, RZC, QAC, EPSAC, SGC, SYMBOL
    use PBC_mod, only: rk
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit

    implicit none
    private
    public :: estructura

contains

!====================================================================
! Subrutinas de utilidad para un manejo de errores robusto
!====================================================================
subroutine open_file_safe(unit, file, status, action)
    ! Intenta abrir un archivo de forma segura y detiene la ejecución en caso de error.
    integer, intent(out) :: unit
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: status, action
    integer :: ios
    
    open(newunit=unit, file=trim(file), status=status, action=action, iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'FATAL ERROR: No se pudo abrir el archivo: ', trim(file)
        stop
    end if
end subroutine open_file_safe

!====================================================================
! Procedimiento principal para la lectura y procesamiento de la estructura
!====================================================================
subroutine estructura(eps, nam, sigma, sigmetano_in, NC, diel)
    ! Lee una estructura atómica, la trunca a la celda de simulación,
    ! convierte a unidades reducidas y la guarda en arrays globales.
    !
    ! Args:
    ! eps, sigma, sigmetano_in, diel (in): Constantes de la simulación.
    ! nam (in): Nombre del archivo de entrada con las coordenadas.
    ! NC (inout): Número total de átomos en el archivo (entrada)
    !             y número de átomos aceptados (salida).
    
    implicit none
    
    ! Inputs y estado de la simulación
    real(rk), intent(in) :: eps, sigma, sigmetano_in, diel
    character(len=*), intent(in) :: nam
    integer, intent(inout) :: NC
    
    ! Variables locales
    integer :: iu_in, iu_dbg, iu_out, ios
    character(len=32) :: nampro
    integer :: i, count_accepted
    
    ! Arrays temporales para la lectura del archivo
    real(rk), allocatable :: RXA(:), RYA(:), RZA(:), EPSA(:), SGCA(:), QACA_phys(:)
    integer, allocatable :: SYMA(:)
    
    ! Vectores para el filtro fraccional
    real(rk) :: r(3), s(3)
    logical :: is_inside
    
    !–– PASO 1: Apertura de archivos e inicialización de la lectura –––––––
    ! Se centraliza el manejo de archivos en la subrutina open_file_safe
    call open_file_safe(iu_dbg, 'PRUEBA.TXT', 'replace', 'write')
    call open_file_safe(iu_in, trim(nam), 'old', 'read')
    
    read(iu_in, *, iostat=ios) NC
    if (ios /= 0 .or. NC <= 0) then
        write(error_unit, '(A,A)') 'Error: NC inválido en ', trim(nam)
        stop
    end if

    !–– PASO 2: Lectura completa del archivo en arrays temporales ––––––––
    allocate(RXA(NC), RYA(NC), RZA(NC), EPSA(NC), SGCA(NC), &
             QACA_phys(NC), SYMA(NC), stat=ios)
    if (ios /= 0) then
        write(error_unit, '(A)') 'Error: no se pudo reservar memoria de lectura.'
        stop
    end if
    
    do i = 1, NC
        read(iu_in, *, iostat=ios) RXA(i), RYA(i), RZA(i), &
                                    EPSA(i), SGCA(i), QACA_phys(i), SYMA(i)
        if (ios /= 0) then
            write(error_unit, '(A,I0,A,A)') 'Error: fallo leyendo línea ', i, ' en ', trim(nam)
            stop
        end if
    end do
    close(iu_in)
    
    write(output_unit, '(A,3F8.2)') 'Dimensiones celda (l/2): ', acelx/2.0, acely/2.0, acelz/2.0

    !–– PASO 3: Filtrado de átomos dentro de la celda de simulación –––––
    count_accepted = 0
    do i = 1, NC
        r = (/ RXA(i), RYA(i), RZA(i) /) ! Coordenadas físicas (Å)
        s = matmul(cell%Ainv, r)        ! Conversión a coordenadas fraccionales
        
        where (cell%pbc) s = s - nint(s) ! Envolver solo los ejes periódicos
        
        ! Comprobar si el átomo está dentro de la celda primaria
        is_inside = (abs(s(1)) < 0.5_rk) .and. (abs(s(2)) < 0.5_rk) .and. (abs(s(3)) < 0.5_rk)
        
        if (is_inside) then
            count_accepted = count_accepted + 1
        end if
    end do
    
    !–– PASO 4: Redimensionar arrays y copiar los datos aceptados –––––––
    NC = count_accepted
    if (.not. allocated(RXC) .or. size(RXC) /= NC) then
        if (allocated(RXC)) deallocate(RXC, RYC, RZC, QAC, EPSAC, SGC, SYMBOL)
        allocate(RXC(NC), RYC(NC), RZC(NC), QAC(NC), EPSAC(NC), SGC(NC), SYMBOL(NC), stat=ios)
        if (ios /= 0) then
            write(error_unit, '(A)') 'Error: No se pudo redimensionar los arrays globales.'
            stop
        end if
    end if
    
    count_accepted = 0
    do i = 1, NC
        r = (/ RXA(i), RYA(i), RZA(i) /)
        s = matmul(cell%Ainv, r)
        where (cell%pbc) s = s - nint(s)
        is_inside = (abs(s(1)) < 0.5_rk) .and. (abs(s(2)) < 0.5_rk) .and. (abs(s(3)) < 0.5_rk)
        
        if (is_inside) then
            count_accepted = count_accepted + 1
            
            ! Conversión a unidades reducidas
            RXC(count_accepted) = RXA(i) * (sigma / sigmetano_in)
            RYC(count_accepted) = RYA(i) * (sigma / sigmetano_in)
            RZC(count_accepted) = RZA(i) * (sigma / sigmetano_in)
            
            EPSAC(count_accepted) = EPSA(i)
            SGC(count_accepted) = SGCA(i)
            QAC(count_accepted) = QACA_phys(i) * FACTORELEC * FCLEC
            SYMBOL(count_accepted) = SYMA(i)
            
            write(iu_dbg, '(3F12.6)') RXC(count_accepted), RYC(count_accepted), RZC(count_accepted)
        end if
    end do
    close(iu_dbg)
    
    write(output_unit, '(A,I0,A)') 'Estructura leída, ', NC, ' segmentos dentro de la celda.'
    
    !–– PASO 5: Escribir archivo de salida 'truncado.txt' –––––––––––––––
    ! En lugar de usar arrays temporales duplicados, se puede reusar la información
    ! de los arrays RXA, RYA, etc. y filtrar de nuevo.
    call open_file_safe(iu_out, 'truncado.txt', 'replace', 'write')
    write(iu_out, *) NC
    do i = 1, NC
        r = (/ RXA(i), RYA(i), RZA(i) /)
        s = matmul(cell%Ainv, r)
        where (cell%pbc) s = s - nint(s)
        is_inside = (abs(s(1)) < 0.5_rk) .and. (abs(s(2)) < 0.5_rk) .and. (abs(s(3)) < 0.5_rk)
        
        if (is_inside) then
            write(iu_out, '(3F12.6,2F12.6,F12.6,1X,I0)') RXA(i), RYA(i), RZA(i), &
                                                        EPSA(i), SGCA(i), QACA_phys(i), SYMA(i)
        end if
    end do
    close(iu_out)
    
    !–– PASO 6: Limpieza de la memoria ––––––––––––––––––––––––––––––––––
    if (allocated(RXA)) then
        deallocate(RXA, RYA, RZA, EPSA, SGCA, QACA_phys, SYMA)
    end if
    
end subroutine estructura

end module EstructuraModule
