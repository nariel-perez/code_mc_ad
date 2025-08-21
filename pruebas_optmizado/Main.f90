!---------------------------------------------------------------------
! File: main.f90 (Refactorizado)
!
! Programa principal para una simulación de Monte Carlo en la
! interfaz sólido-fluido. Coordina la lectura de parámetros,
! la inicialización de estructuras de datos y el cálculo de potenciales.
!---------------------------------------------------------------------
program main
    !-----------------------------------------------------------------
    ! Módulos de Dependencias
    !-----------------------------------------------------------------
    use InputParams
    use SimulationData
    use AdsorbateInput
    use PhysicalConstants
    use EstructuraModule
    use RotationModule
    use, intrinsic :: iso_fortran_env, only: output_unit
    
    implicit none
    
    !-----------------------------------------------------------------
    ! Declaración de variables
    ! Se agrupan lógicamente y se añaden comentarios explicativos.
    !-----------------------------------------------------------------

    ! Parámetros de la simulación
    real :: p_ratio, SIGMA, TEMP, PRED, RCUT
    integer :: i, NC, auxmat
    
    ! Arrays dinámicos para la simulación
    real, allocatable :: p_vals(:), Z(:)
    
    ! Dimensiones de la celda de simulación en unidades reducidas
    real :: XMAX, YMAX, ZMAX, VOL
    
    ! Constante física
    real, parameter :: AK_input = 8.31_rk ! Constante de los gases en J/mol·K

    !-----------------------------------------------------------------
    ! Inicio del programa
    ! Se divide en secciones lógicas para mejor legibilidad.
    !-----------------------------------------------------------------
    
    !---------------------------------------
    ! 1. LECTURA DE ARCHIVOS DE ENTRADA
    !---------------------------------------
    call read_input('input.txt')
    
    !---------------------------------------
    ! 2. INICIALIZACIÓN DE ARRAYS GLOBALES
    !---------------------------------------
    ! Se centralizan las alocaciones en un solo bloque.
    auxmat = int(mat / 2)
    if (.not. allocated(UADS)) then
        allocate(UADS(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, 50))
    end if
    
    ! Arrays logarítmicos
    if (.not. allocated(p_vals)) then
        allocate(p_vals(isot + 1))
        p_ratio = (dp / p)**(1.0_rk / real(isot - 1, kind=rk))
        p_vals(1) = p
        do i = 2, isot
            p_vals(i) = p_vals(i-1) * p_ratio
        end do
    end if
    
    ! Alocación de arrays para adsorbato y superficie
    if (.not. allocated(Z)) allocate(Z(NMOLEC))
    if (.not. allocated(N)) allocate(N(NMOLEC))
    
    ! Alocación de propiedades de átomos (se asume que maxAtoms está definido en AdsorbateInput)
    if (.not. allocated(EPSI)) allocate(EPSI(maxAtoms))
    if (.not. allocated(SIGM)) allocate(SIGM(maxAtoms))
    if (.not. allocated(Q)) allocate(Q(maxAtoms))
    
    ! Alocación de arrays de coordenadas
    if (.not. allocated(RX1)) allocate(RX1(maxAtoms))
    if (.not. allocated(RY1)) allocate(RY1(maxAtoms))
    if (.not. allocated(RZ1)) allocate(RZ1(maxAtoms))

    ! Alocación de arrays principales de datos
    ! Se corrige la asignación de dimensiones, utilizando las variables del Input
    if (.not. allocated(RX)) allocate(RX(mat_sim, maxAtoms, NMOLEC))
    if (.not. allocated(RY)) allocate(RY(mat_sim, maxAtoms, NMOLEC))
    if (.not. allocated(RZ)) allocate(RZ(mat_sim, maxAtoms, NMOLEC))
    if (.not. allocated(USS)) allocate(USS(mat_sim, maxAtoms, maxAtoms))
    
    auxmat = int(NCELLMAT / 2)
    if (.not. allocated(CNF)) then
        allocate(CNF(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, NMOLEC, maxAtoms))
    end if

    !---------------------------------------
    ! 3. CONVERSIÓN A UNIDADES REDUCIDAS
    !---------------------------------------
    SIGMA = sigmetano / acel
    TEMP = T / eps
    ! P = P ! se comenta porque no se usa el valor de la variable en esta línea.
    PRED = P * SIGMA**3 / eps
    RCUT = 10.0_rk * SIGMA
    
    XMAX = acelx / acel
    YMAX = acely / acel
    ZMAX = acelz / acel
    VOL = XMAX * YMAX * ZMAX
    
    write(output_unit, '(A)') ''
    write(output_unit, '(A)') '-------------------------------------------------'
    write(output_unit, '(A)') '-------------- REDUCED UNITS --------------------'
    write(output_unit, '(A, F10.4)') 'TEMPERATURE: ', TEMP
    write(output_unit, '(A, ES12.5)') 'PRESSURE: ', PRED
    write(output_unit, '(A, F10.4)') 'SIGMA:', SIGMA
    write(output_unit, '(A, ES12.4)') 'VOLUME:', VOL

    !---------------------------------------
    ! 4. INICIALIZACIÓN DE CONSTANTES Y TABLAS
    !---------------------------------------
    call computeconstants(sigmetano, SIGMA, eps, AK_input, diel)
    call initrotationtables()

    !---------------------------------------
    ! 5. LECTURA DE ESTRUCTURAS MOLECULARES
    !---------------------------------------
    write(output_unit, '(A)') ''
    write(output_unit, '(A)') '-------------------------------------------------'
    write(output_unit, '(A)') '----------------- ADSORBATES --------------------'
    call read_adsorbates('MOLEC.DAT', sigma, 1.0e-7)

    !---------------------------------------
    ! 6. CÁLCULO DE POTENCIALES
    !---------------------------------------
    write(output_unit, '(A)') ''
    write(output_unit, '(A)') '-------------------------------------------------'
    write(output_unit, '(A)') '------------------ SURFACE ----------------------'
    open(newunit=50, file='SALIDAACTIVADO-100.TXT', status='replace')
    open(newunit=97, file='PERFILES.TXT', status='replace')
    call estructura(eps, nam, sigma, sigmetano, NC, diel)
    
    write(output_unit, '(A)') ''
    write(output_unit, '(A)') '-------------------------------------------------'
    write(output_unit, '(A)') '------------------ POTENTIALS -------------------'
    write(output_unit, '(A)') ''
    
    ! Se asume que las subrutinas POTENCIALFF y POTENCIAL están en otros módulos
    call potencialff(eps, sigma, sigmetano, NC, RCUT, diel)
    call potencial(eps, sigma, sigmetano, NC, RCUT, diel)

    !---------------------------------------
    ! 7. Limpieza y finalización (opcional)
    !---------------------------------------
    ! Se recomienda cerrar los archivos al final del programa.
    ! close(50)
    ! close(97)
    
end program main
