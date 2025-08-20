!---------------------------------------------------------------------
! File: POTENCIAL.f90 (Refactorizado)
!
! Subrutina para el cálculo del potencial de interacción entre
! un adsorbato y una superficie. Calcula la energía de Lennard-Jones
! y el potencial electrostático Ewald para una malla 3D.
!
! El código se optimizó para claridad, eficiencia y manejo de errores,
! manteniendo los nombres de variables y la lógica original.
!---------------------------------------------------------------------
subroutine POTENCIAL(EPS, sigma, sigmetano, NC, RCUT, diel)
    use InputParams, only: mat, cell
    use PBC_Mod, only: rk, Cell_t => Cell
    use PhysicalConstants, only: FCLEC, FACTORELEC
    use SimulationData, only: UADS, RXC, RYC, RZC, EPSAC, SGC, EPSI, SIGM, Q, QAC
    use GeomUtils, only: get_cell_metrics, r2_min_image_frac, wrap_by_pbc
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

    implicit none

    !------------------------- Argumentos
    integer, intent(in) :: NC
    real(rk), intent(in) :: EPS, sigma, sigmetano, RCUT, diel

    !------------------------- Variables locales
    integer :: NKIND, INKIND, KINDI, IPOT
    integer :: i, ij, k, j, ios, m2
    real(rk) :: sI1, sI2, sI3, s1, s2, s3, r2
    real(rk) :: sr2, sr6, vij, vijr, deltv, deltw
    real(rk) :: pp
    real(rk), parameter :: RCELE = 0.5_rk
    real(rk), parameter :: BIGV = 1.0e6_rk
    real(rk), parameter :: ZERO = 0.0_rk

    ! Celda local y métricas (variables temporales para las llamadas)
    type(Cell_t) :: cellR
    real(rk), dimension(3,3) :: G, Ainv
    real(rk), dimension(3) :: s_vec_diff, r_vec
    logical, dimension(3) :: pbc_flags

    ! Sólido en fraccionales
    real(rk), allocatable :: sC1(:), sC2(:), sC3(:)

    ! Precálculos por átomo
    real(rk), allocatable :: sigma1(:), sigsq(:), rmin(:), r2min(:), r2cut(:), factor(:)

    !------------------------- Paso 1: Lectura de LJ.dat
    open(11, file='LJ.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A)') 'Error: No se pudo abrir LJ.dat'
        return
    end if

    read(11, *, iostat=ios) NKIND
    if (ios /= 0 .or. NKIND <= 0) then
        write(error_unit, '(A)') 'Error: NKIND inválido en LJ.dat'
        close(11); return
    end if

    allocate(EPSI(NKIND), SIGM(NKIND), Q(NKIND))

    do INKIND = 1, NKIND
        read(11, *, iostat=ios) KINDI, EPSI(INKIND), SIGM(INKIND), Q(INKIND)
        if (ios /= 0) then
            write(error_unit, '(A,I0,A)') 'Error leyendo línea ', INKIND, ' de LJ.dat'
            close(11); return
        end if
        Q(INKIND) = Q(INKIND) * FACTORELEC * FCLEC
    end do
    close(11)

    !------------------------- Paso 2: Configuración de la celda de simulación
    cellR = cell
    ! Usar una única operación para la conversión, más limpia y legible.
    cellR%A = cellR%A * (sigma / sigmetano)
    call cellR%update() ! Actualizar Ainv y otras propiedades internas.

    ! Obtener las métricas de la celda
    call get_cell_metrics(cellR, G, Ainv, pbc_flags)

    !------------------------- Paso 3: Conversión de coordenadas del sólido a fraccionales
    allocate(sC1(NC), sC2(NC), sC3(NC))
    do j = 1, NC
        ! Usar matmul para una conversión vectorial
        r_vec = (/ RXC(j), RYC(j), RZC(j) /)
        s_vec_diff = matmul(Ainv, r_vec)
        sC1(j) = s_vec_diff(1); sC2(j) = s_vec_diff(2); sC3(j) = s_vec_diff(3)
    end do

    !------------------------- Paso 4: Barrido de malla y cálculo de potencial
    m2 = mat / 2
    write(output_unit, '(A)') '-----------------------'

    do IPOT = 1, NKIND
        write(output_unit, '(A,I0,A,I0)') 'POTENCIAL PARA ', IPOT, ' de ', NKIND
        
        ! Pre-cálculos por tipo de átomo (dependen del tipo IPOT)
        allocate(sigma1(NC), sigsq(NC), rmin(NC), r2min(NC), r2cut(NC), factor(NC))
        do j = 1, NC
            sigma1(j) = sigma * (SGC(j) + SIGM(IPOT)) / (2.0_rk * sigmetano)
            sigsq(j) = sigma1(j) * sigma1(j)
            rmin(j) = 0.5_rk * sigma1(j)
            r2min(j) = rmin(j) * rmin(j)
            r2cut(j) = (8.0_rk * rmin(j))**2
            factor(j) = sqrt(EPSI(IPOT) * EPSAC(j)) / EPS
        end do
        
        ! Barrido de la malla 3D
        do k = 1, mat
            do ij = 1, mat
                do i = 1, mat
                    ! Conversión de índice de malla a coordenadas fraccionales
                    sI1 = (real(i, rk) - m2 - 1.0_rk) / mat
                    sI2 = (real(ij, rk) - m2 - 1.0_rk) / mat
                    sI3 = (real(k, rk) - m2 - 1.0_rk) / mat

                    deltv = ZERO
                    
                    ! Bucle sobre los átomos del sólido
                    do j = 1, NC
                        ! Δs con 'wrap' por eje
                        s1 = sI1 - sC1(j)
                        s2 = sI2 - sC2(j)
                        s3 = sI3 - sC3(j)
                        call wrap_by_pbc(s1, s2, s3, pbc_flags(1), pbc_flags(2), pbc_flags(3))
                        
                        ! Calcular r^2
                        r2 = r2_min_image_frac(G, s1, s2, s3)

                        if (r2 < r2min(j)) then
                            deltv = BIGV
                            exit
                        else
                            sr2 = sigsq(j) / r2
                            sr6 = sr2 * sr2 * sr2
                            vij = sr6 * (sr6 - 1.0_rk) * factor(j)
                            
                            if (r2 > r2cut(j)) vij = ZERO
                            
                            vijr = ZERO
                            if (r2 < RCELE**2) then
                                vijr = Q(IPOT) * QAC(j) / (diel*sqrt(r2)) ! Se incorpora diel
                            end if
                        end if
                        
                        deltv = deltv + 4.0_rk * vij + vijr
                    end do
                    
                    ! Almacenar el resultado en el array UADS
                    UADS(i, ij, k, IPOT) = deltv
                    
                end do
            end do
        end do
        
        deallocate(sigma1, sigsq, rmin, r2min, r2cut, factor)
    end do
    
    !------------------------- Paso 5: Limpieza y finalización
    deallocate(sC1, sC2, sC3)
    write(output_unit, '(A)') 'Energia Calculada'
    return
end subroutine POTENCIAL
