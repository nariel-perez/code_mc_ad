!---------------------------------------------------------------------
! File: AdsorbateInput.f90 (Refactorizado)
!
! Módulo para la lectura de estructuras de moléculas (adsorbatos)
! desde archivos de entrada. Procesa MOLEC.DAT y los archivos de
! cada especie molecular para cargar las coordenadas y propiedades.
!---------------------------------------------------------------------
module AdsorbateInput
    use PBC_Mod, only : rk
    use InputParams, only : sigmetano
    use, intrinsic :: iso_fortran_env, only : output_unit, iostat_end

    implicit none
    private

    !====== API pública ==========================================
    ! Se declara explícitamente la interfaz pública del módulo.
    public :: read_adsorbates
    public :: NMOLEC, maxAtoms, XT, X
    public :: NATOM, NTOTAL, NMIN2, NMAXI2
    public :: RX0, RY0, RZ0, NATOMKIND, NSYM

    !====== variables globales ==================================
    ! Se mantienen los nombres de variables originales
    integer, public :: NMOLEC = 0
    integer, public :: maxAtoms = 0
    real(rk), public :: XT = 0.0_rk
    real(rk), allocatable, public :: X(:)
    integer, allocatable, public :: NATOM(:), NTOTAL(:)
    integer, allocatable, public :: NMIN2(:), NMAXI2(:)
    real(rk), allocatable, public :: RX0(:,:), RY0(:,:), RZ0(:,:)
    integer, allocatable, public :: NATOMKIND(:,:)
    integer, allocatable, public :: NSYM(:,:)

contains

!=================================================================
subroutine read_adsorbates(molec_file, sigma, tol)
    ! Lee el archivo MOLEC.DAT y los archivos de cada molécula para
    ! inicializar las estructuras de datos de la simulación.
    !
    ! Args:
    ! molec_file (in): Nombre del archivo principal (ej. MOLEC.DAT).
    ! sigma (in): Parámetro sigma de referencia para la conversión
    !             de unidades.
    ! tol (in): Tolerancia para la comprobación de la suma de fracciones.
    
    character(len=*), intent(in) :: molec_file
    real(rk), intent(in) :: sigma, tol

    ! Subrutinas de utilidad para el manejo de archivos
    ! Encapsuladas para un mejor control de errores
    interface
        subroutine open_file_safe(unit, file)
            integer, intent(out) :: unit
            character(len=*), intent(in) :: file
        end subroutine
    end interface
    
    ! Variables auxiliares para la lectura
    integer :: iu_main, iu_mol, ios, i, j
    character(len=256), allocatable :: molname(:)
    integer, allocatable :: n_atoms_in_file(:), n_conf_min(:), n_conf_max(:)
    real(rk), allocatable :: fraction_tmp(:)
    integer :: atoms_read, n_conf_min_tmp, n_conf_max_tmp, ikind, ns
    real(rk) :: x1, y1, z1, ax
    
    !–– PASO 1: Leer el archivo de cabecera (MOLEC.DAT) para obtener metadatos
    call open_file_safe(iu_main, trim(molec_file))
    read(iu_main, *) NMOLEC
    
    allocate(molname(NMOLEC), n_atoms_in_file(NMOLEC), &
             n_conf_min(NMOLEC), n_conf_max(NMOLEC), fraction_tmp(NMOLEC))

    maxAtoms = 0
    do i = 1, NMOLEC
        read(iu_main, '(A)') molname(i)
        call open_file_safe(iu_mol, trim(molname(i)))
        read(iu_mol, *) atoms_read, ax, n_conf_min_tmp, n_conf_max_tmp
        
        n_atoms_in_file(i) = atoms_read
        fraction_tmp(i) = ax
        n_conf_min(i) = n_conf_min_tmp
        n_conf_max(i) = n_conf_max_tmp
        
        if (atoms_read > maxAtoms) maxAtoms = atoms_read
        close(iu_mol)
    end do
    close(iu_main)

    !–– PASO 2: Asignar y redimensionar los arrays definitivos
    ! Se realiza una sola llamada a allocate para cada grupo de arrays.
    ! La asignación se realiza después de la lectura inicial para evitar
    ! múltiples alocaciones y deslocalizaciones.
    allocate(X(NMOLEC), NATOM(NMOLEC), NMIN2(NMOLEC), NMAXI2(NMOLEC), NTOTAL(NMOLEC))
    X = fraction_tmp
    NATOM = n_atoms_in_file
    NMIN2 = n_conf_min
    NMAXI2 = n_conf_max
    NTOTAL = 0 ! Inicializar a 0
    
    allocate(RX0(maxAtoms, NMOLEC), RY0(maxAtoms, NMOLEC), RZ0(maxAtoms, NMOLEC))
    RX0 = 0.0_rk; RY0 = 0.0_rk; RZ0 = 0.0_rk
    
    allocate(NATOMKIND(maxAtoms, NMOLEC), NSYM(maxAtoms, NMOLEC))
    NATOMKIND = 0; NSYM = 0

    !–– PASO 3: Leer las coordenadas y propiedades de cada molécula
    do i = 1, NMOLEC
        call open_file_safe(iu_mol, trim(molname(i)))
        read(iu_mol, *) ! Leer la cabecera (se descarta)
        
        do j = 1, NATOM(i)
            read(iu_mol, *) x1, y1, z1, ikind, ns
            ! Conversión de unidades
            RX0(j,i) = x1 * (sigma / sigmetano)
            RY0(j,i) = y1 * (sigma / sigmetano)
            RZ0(j,i) = z1 * (sigma / sigmetano)
            NATOMKIND(j,i) = ikind
            NSYM(j,i) = ns
        end do
        close(iu_mol)
    end do

    !–– PASO 4: Realizar comprobaciones y mostrar un resumen
    XT = sum(X)
    if (abs(XT - 1.0_rk) > tol) then
        write(output_unit, '(A,ES12.5)') 'WARNING: sum(X) ≠ 1 -> ', XT
    end if

    write(output_unit, *) ''
    write(output_unit, *) '----------- ADSORBATES LOADED ---------------'
    write(output_unit, '(A,I6)') 'NMOLEC            :', NMOLEC
    write(output_unit, '(A,I6)') 'max atoms/specie  :', maxAtoms
    write(output_unit, '(A,ES12.5)') 'Sum fractions XT:', XT
    write(output_unit, *) '---------------------------------------------'

    ! Limpiar arrays temporales
    if (allocated(molname)) deallocate(molname, n_atoms_in_file, n_conf_min, n_conf_max, fraction_tmp)

end subroutine read_adsorbates

!-----------------------------------------------------------------
! Subrutinas de utilidad para un manejo de errores robusto
!-----------------------------------------------------------------
subroutine open_file_safe(unit, file)
    ! Intenta abrir un archivo de forma segura y detiene la ejecución en caso de error.
    integer, intent(out) :: unit
    character(len=*), intent(in) :: file
    integer :: ios
    
    open(newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'FATAL ERROR: No se pudo abrir el archivo: ', trim(file)
        stop
    end if
end subroutine open_file_safe

end module AdsorbateInput
