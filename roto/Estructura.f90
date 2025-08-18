!─────────────────────────────────────────────────────────────────
! File: Estructura.f90   (Fortran 90/2003)
! Lectura de la superficie (grafito) y armado de arreglos internos.
! GENERALIZADO a celda triclinic:
!   - r[Å] -> s = Ainv*r  (fraccional)
!   - wrap SOLO en ejes con PBC
!   - aceptar |s_i| < 0.5  (celda primaria, evita doble conteo)
!   - guardar en reducidas: RXC = RXA/sigmetano * sigma  (= RXA/ACEL)
!─────────────────────────────────────────────────────────────────
module EstructuraModule
  use InputParams,       only: cell, sigmetano, acelx, acely, acelz
  use PhysicalConstants, only: FACTORELEC, FCLEC
  use SimulationData,    only: RXC, RYC, RZC, QAC, EPSAC, SGC, SYMBOL
  use PBC_mod, only: rk
  implicit none
contains

  subroutine estructura(eps, nam, sigma, sigmetano_in, NC, diel)
    implicit none
    ! Inputs
    real(rk),               intent(in)    :: eps, sigma, sigmetano_in, diel
    character(len=*),   intent(in)    :: nam
    ! In/Out
    integer,            intent(inout) :: NC

    ! Locals
    integer :: ios, i, imax
    character(len=32) :: nampro

    ! Datos leídos (físicos, en Å)
    real(rk),    allocatable :: RXA(:), RYA(:), RZA(:), EPSA(:), SGCA(:), QACA_phys(:)
    integer, allocatable :: SYMA(:)

    ! Copias de los aceptados (para truncado.txt)
    real(rk),    allocatable :: RXA_in(:), RYA_in(:), RZA_in(:), EPSA_in(:), SGCA_in(:), QACA_in_phys(:)
    integer, allocatable :: SYMA_in(:)

    ! Vectores para filtro fraccional
    real(rk) :: r(3), s(3), s0(3)
    logical :: inside

    ! Archivos
    integer, parameter :: iu_in=51, iu_dbg=49, iu_out=52

    !--------------------------
    ! Abrir archivos
    !--------------------------
    open(iu_dbg, file='PRUEBA.TXT', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: no se pudo abrir PRUEBA.TXT'
      return
    end if

    open(iu_in, file=trim(nam), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: no se pudo abrir ', trim(nam)
      close(iu_dbg)
      return
    end if

    !--------------------------
    ! Leer NC y reservar
    !--------------------------
    read(iu_in, *, iostat=ios) NC
    if (ios /= 0 .or. NC <= 0) then
      print *, 'Error: NC inválido en ', trim(nam)
      close(iu_in); close(iu_dbg)
      return
    end if

    allocate(RXA(NC), RYA(NC), RZA(NC), EPSA(NC), SGCA(NC), QACA_phys(NC), SYMA(NC), stat=ios)
    if (ios /= 0) then
      print *, 'Error: no se pudo reservar memoria de lectura'
      close(iu_in); close(iu_dbg)
      return
    end if

    !--------------------------
    ! Leer líneas (Å, físicas)
    !--------------------------
    do i = 1, NC
      read(iu_in, *, iostat=ios) RXA(i), RYA(i), RZA(i), EPSA(i), SGCA(i), QACA_phys(i), SYMA(i)
      if (ios /= 0) then
        print *, 'Error: fallo leyendo línea ', i, ' en ', trim(nam)
        close(iu_in); close(iu_dbg)
        return
      end if
    end do
    close(iu_in)

    print *, 'Dimensiones celda (l/2): ', acelx/2.0, acely/2.0, acelz/2.0

    !--------------------------
    ! Contar aceptados (fraccional + wrap por PBC, general)
    !--------------------------
    imax = 0
    do i = 1, NC
      r  = (/ RXA(i), RYA(i), RZA(i) /)          ! Å
      s0 = matmul(cell%Ainv, r)                  ! Å -> frac
      s  = s0
      where (cell%pbc) s = s - nint(s)           ! envolver SOLO ejes con PBC
      inside = (abs(s(1)) < 0.5) .and. (abs(s(2)) < 0.5) .and. (abs(s(3)) < 0.5)
      if (inside) imax = imax + 1
    end do

    ! Reservar globales + buffers truncados
    if (.not. allocated(RXC))   allocate(RXC(imax))
    if (.not. allocated(RYC))   allocate(RYC(imax))
    if (.not. allocated(RZC))   allocate(RZC(imax))
    if (.not. allocated(QAC))   allocate(QAC(imax))
    if (.not. allocated(EPSAC)) allocate(EPSAC(imax))
    if (.not. allocated(SGC))   allocate(SGC(imax))
    if (.not. allocated(SYMBOL))allocate(SYMBOL(imax))

    allocate(RXA_in(imax), RYA_in(imax), RZA_in(imax), EPSA_in(imax), SGCA_in(imax), QACA_in_phys(imax), SYMA_in(imax), stat=ios)
    if (ios /= 0) then
      print *, 'Error: no se pudo reservar buffers truncados'
      close(iu_dbg)
      return
    end if

    !--------------------------
    ! Copiar/convertir aceptados
    !--------------------------
    imax = 0
    do i = 1, NC
      r  = (/ RXA(i), RYA(i), RZA(i) /)
      s0 = matmul(cell%Ainv, r)
      s  = s0
      where (cell%pbc) s = s - nint(s)
      inside = (abs(s(1)) < 0.5) .and. (abs(s(2)) < 0.5) .and. (abs(s(3)) < 0.5)
      if (inside) then
        imax = imax + 1

        ! Unidades reducidas (= RXA/ACEL)
        RXC(imax) = RXA(i) / sigmetano * sigma
        RYC(imax) = RYA(i) / sigmetano * sigma
        RZC(imax) = RZA(i) / sigmetano * sigma

        EPSAC(imax) = EPSA(i)
        SGC(imax)   = SGCA(i)
        QAC(imax)   = QACA_phys(i) * FACTORELEC * FCLEC
        SYMBOL(imax)= SYMA(i)

        write(iu_dbg, *) RXC(imax), RYC(imax), RZC(imax)

        ! Guardar originales (para truncado.txt)
        RXA_in(imax)      = RXA(i)
        RYA_in(imax)      = RYA(i)
        RZA_in(imax)      = RZA(i)
        EPSA_in(imax)     = EPSA(i)
        SGCA_in(imax)     = SGCA(i)
        QACA_in_phys(imax)= QACA_phys(i)
        SYMA_in(imax)     = SYMA(i)
      end if
    end do
    close(iu_dbg)

    NC = imax
    print *, 'Estructura leída, ', NC, ' segmentos dentro de la celda.'

    !--------------------------
    ! Escribir truncado.txt en unidades físicas
    !--------------------------
    nampro = trim(nam)//'_truncado'
    open(iu_out, file='truncado.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: no se pudo abrir ', trim(nampro)
      return
    end if
    write(iu_out, *) NC
    do i = 1, NC
      write(iu_out, '(3F12.6,2F12.6,F12.6,1X,I0)') RXA_in(i), RYA_in(i), RZA_in(i), &
                                                    EPSA_in(i), SGCA_in(i), QACA_in_phys(i), SYMA_in(i)
    end do
    close(iu_out)

    ! Limpieza
    deallocate(RXA, RYA, RZA, EPSA, SGCA, QACA_phys, SYMA)
    deallocate(RXA_in, RYA_in, RZA_in, EPSA_in, SGCA_in, QACA_in_phys, SYMA_in)

  end subroutine estructura

end module EstructuraModule

