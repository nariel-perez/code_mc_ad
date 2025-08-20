!---------------------------------------------------------------------
! File: POTENCIALFF.f90 (Refactorizado)
!
! Subrutina para calcular el potencial de par entre dos especies
! atómicas y almacenarlo en una matriz de potenciales precalculados.
! Combina los potenciales de Lennard-Jones y electrostático (Ewald).
!---------------------------------------------------------------------
subroutine POTENCIALFF(EPS, sigma, sigmetano, NC, RCUT, diel)
    use PhysicalConstants, only: FCLEC, FACTORELEC
    use SimulationData, only: USS, FLAG, EPSI, SIGM, Q
    use PBC_Mod, only: rk
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

    implicit none
    
    !------------------------- Argumentos
    real(rk), intent(in) :: EPS, sigma, sigmetano, RCUT, diel
    integer, intent(in) :: NC ! No se usa en esta subrutina

    !------------------------- Variables locales y parámetros
    integer, parameter :: MAX_POINTS = 5000
    integer :: NKIND, INKIND, KINDI, IPOT, JPOT, i
    integer :: ios
    
    real(rk) :: REDELEC
    real(rk) :: SIGMA1, FACTOR, RCUTSQ, SIGSQ, SIGCUB
    real(rk) :: RMIN, RMINSQ, RZI, RIJSQ
    real(rk) :: SR2, SR6, VIJ, WIJ, VIJELEC
    real(rk), parameter :: RCELE = 0.35_rk 
    real(rk), parameter :: VERY_LARGE_LJ = 1.0e10_rk
    real(rk), parameter :: VERY_LARGE_ELEC = 1.0e9_rk
    real(rk), parameter :: PI = acos(-1.0_rk)

    !------------------------- Paso 1: Lectura de LJ.dat
    open(11, file='LJ.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A)') "Error: No se pudo abrir el archivo LJ.dat"
        return
    end if

    read(11, *, iostat=ios) NKIND, REDELEC
    if (ios /= 0 .or. NKIND <= 0) then
        write(error_unit, '(A)') "Error: NKIND o REDELEC inválido en LJ.dat"
        close(11)
        return
    end if

    ! Ajustar tamaño de arrays dinámicos
    ! Se asume que EPSI, SIGM y Q son allocatable en SimulationData.
    ! El código original no los aloca explícitamente, lo cual es un riesgo.
    if (.not. allocated(EPSI) .or. size(EPSI) < NKIND) then
        if (allocated(EPSI)) deallocate(EPSI, SIGM, Q)
        allocate(EPSI(NKIND), SIGM(NKIND), Q(NKIND))
    end if
    
    write(output_unit, '(A10, A10, A10, A10)') 'TYPE', 'EPSILON', 'SIGMA', 'CHARGE'
    write(output_unit, '(A10, A10, A10, A10)') '----', '-------', '-----', '------'

    do INKIND = 1, NKIND
        read(11, *, iostat=ios) KINDI, EPSI(INKIND), SIGM(INKIND), Q(INKIND)
        if (ios /= 0) then
            write(error_unit, '(A,I0,A)') "Error al leer parámetros para el tipo de molécula ", INKIND, " en LJ.dat"
            close(11)
            return
        end if
        write(output_unit, *) KINDI, EPSI(INKIND), SIGM(INKIND), Q(INKIND)

        ! Conversión de unidades de carga
        Q(INKIND) = Q(INKIND) * FACTORELEC * FCLEC
    end do
    close(11)

    !------------------------- Paso 2: Cálculo y almacenamiento de potenciales
    ! Se asume que USS es allocatable en SimulationData y se redimensiona
    ! para evitar errores de out-of-bounds. El código original no lo hace.
    if (.not. allocated(USS) .or. size(USS, dim=3) < NKIND) then
        if (allocated(USS)) deallocate(USS)
        allocate(USS(MAX_POINTS, NKIND, NKIND))
    end if
    
    ! Bucle anidado para calcular los pares de potenciales
    do IPOT = 1, NKIND
        do JPOT = IPOT, NKIND
            ! Parámetros de la interacción
            SIGMA1 = SIGMA * (SIGM(IPOT) + SIGM(JPOT)) / (2.0_rk * sigmetano)
            FACTOR = sqrt(EPSI(IPOT) * EPSI(JPOT)) / EPS
            RCUTSQ = RCUT * RCUT
            SIGSQ = SIGMA1 * SIGMA1
            RMIN = 0.5_rk * SIGMA1
            RMINSQ = RMIN * RMIN

            ! Bucle sobre los puntos de la malla
            do i = 1, MAX_POINTS
                RZI = real(i, rk) / 1000.0_rk
                RIJSQ = RZI * RZI

                ! Potencial de Lennard-Jones
                if (RIJSQ < RMINSQ) then
                    VIJ = VERY_LARGE_LJ
                else
                    SR2 = SIGSQ / RIJSQ
                    SR6 = SR2 * SR2 * SR2
                    VIJ = SR6 * (SR6 - 1.0_rk)
                    ! Corte del potencial de LJ
                    if (RIJSQ > RCUTSQ) VIJ = 0.0_rk
                end if

                ! Potencial electrostático (corrección de Ewald)
                VIJELEC = 0.0_rk
                if (RZI < RCELE) then
                    VIJELEC = Q(IPOT) * Q(JPOT) * &
                            (1.0_rk / RZI - 1.0_rk / RCELE)
                end if
                
                ! Suma total del potencial
                USS(i, IPOT, JPOT) = 4.0_rk * VIJ * FACTOR + VIJELEC
                
                ! Aplicar simetría
                if (IPOT /= JPOT) then
                    USS(i, JPOT, IPOT) = USS(i, IPOT, JPOT)
                end if
            end do
        end do
    end do

    write(output_unit, '(A)') 'Potenciales de par calculados'
    
end subroutine POTENCIALFF
