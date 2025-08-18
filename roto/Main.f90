!─────────────────────────────────────────────────────────────────────
!  File: Main.f90   (Fortran 2003, triclinic PBC + paso MC legacy)
!─────────────────────────────────────────────────────────────────────
program main
   use PBC_Mod,           only : rk
   use InputParams,       only : read_input, P, dp, sigmetano, eps, diel, T, &
                                 ACEL, acelx, acely, acelz,                   &
                                 mat, nam, isot, ijpasos, ikpasos, mult2,     &
                                 ensemble, canonicalmolecules, NESTADO,       &
                                 ncellmat, cell
   use PhysicalConstants, only : ComputeConstants
   use RotationModule,    only : InitRotationTables
   use EstructuraModule,  only : estructura
   use AdsorbateInput     ! exporta: NMOLEC, maxAtoms, X(:), NATOM(:), ...
   use SimulationData
   implicit none

   !──────── tipos/funciones externas legacy ──────────────────────────
   real :: ranf, dummy
   external :: ranf

   !──────── variables principales ────────────────────────────────────
   real(rk)            :: p_ratio, SIGMA, TEMPSTAR, PRED, RCUT
   real(rk), allocatable :: p_vals(:), Z(:)
   real(rk)            :: VOL, XMAX, YMAX, ZMAX, detA
   integer             :: NC
   integer             :: i, j

   !──────── acumuladores/estadística ────────────────────────────────
   real(rk) :: V, VG, VA, W
   real(rk) :: Uacc, UGacc, UAacc, U2acc
   real(rk) :: UNacc, UNGacc, UNAacc, ANacc
   real(rk) :: N2acc, AN1, U1, UN1, UG1, UNA1, UA1, UNG1, AN2, ANN
   real(rk) :: CALOR, CALORG, CALORA, CALORESP1, CALORESP2, CALORESP3
   real(rk) :: auvol
   integer  :: IPASOS, JPAS, KPAS, MULT
   integer  :: IJ, ntotalB, ntotalGraf
   integer  :: natomki, icnf, jcnf, kcnf, ipull, nmolec2, natom2, molkind, ncantmol
   real(rk) :: X1, Y1, Z1
   real(rk) :: RXAI, RYAI, RZAI, EPSAI, SGCI, QACI
   integer  :: SYMBOL2
   real(rk) :: RXN, RYN, RZN
   real(rk) :: starttime, endtime, elapsedtime
   character(len=16) :: CONFIG
   character(len=16) :: CONFAT, CONFNAT
   logical           :: CREATE, GHOST
   real(rk)          :: CR, DELTVA, DELTV
   real(rk), allocatable :: ANPROM(:)

   integer :: auxmat, ensemble2

   !──────────────────────────────────────────────────────────────────
   call cpu_time(starttime)

   ! 1) Entrada (imprime resumen dentro de read_input)
   call read_input('input.txt')

   ! 2) Arrays base (UADS/CNF) con límites simétricos como el legacy
   auxmat = int(mat/2)
   if (.not.allocated(UADS)) allocate(UADS(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, 50))

   auxmat = int(NCELLMAT/2)
   if (.not.allocated(CNF)) then
      allocate(CNF(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat,  NMOLEC, maxAtoms))
   end if

   ! 3) Pasos de presión logarítmicos
   if (.not.allocated(p_vals)) allocate(p_vals(isot+1))
   p_vals = 0.0_rk
   if (isot > 1) then
      p_ratio = (dp / P)**(1.0_rk / real(isot-1, rk))
   else
      p_ratio = 1.0_rk
   end if
   p_vals(1) = P
   do i = 2, isot
      p_vals(i) = p_vals(i-1) * p_ratio
   end do
   p_vals(isot+1) = p_vals(isot) * p_ratio   ! para la línea p = p_vals(ipasos+1)

   ! 4) Unidades reducidas
   SIGMA    = sigmetano / ACEL
   TEMPSTAR = T / eps
   PRED     = P * SIGMA**3 / eps
   RCUT     = 10.0_rk * SIGMA

   ! volumen reducido (triclinic): det(A)/ACEL^3
   detA =  cell%A(1,1)*(cell%A(2,2)*cell%A(3,3) - cell%A(2,3)*cell%A(3,2)) &
         - cell%A(2,1)*(cell%A(1,2)*cell%A(3,3) - cell%A(1,3)*cell%A(3,2)) &
         + cell%A(3,1)*(cell%A(1,2)*cell%A(2,3) - cell%A(1,3)*cell%A(2,2))
   VOL  = abs(detA) / ACEL**3

   XMAX = sqrt(sum(cell%A(:,1)**2)) / ACEL
   YMAX = sqrt(sum(cell%A(:,2)**2)) / ACEL
   ZMAX = sqrt(sum(cell%A(:,3)**2)) / ACEL

   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '--------------REDUCED UNITS ---------------------'
   write(*,'(A, F10.4)') 'TEMPERATURE: ', TEMPSTAR
   write(*,'(A, ES12.5)') 'PRESSURE: ',    PRED
   write(*,'(A, F10.4)') 'SIGMA:',         SIGMA
   write(*,'(A, ES12.4)') 'VOLUME:',       VOL

   ! 5) Constantes eléctricas + tablas rotación
   call ComputeConstants(sigmetano, SIGMA, eps, 8.31_rk, diel)
   call InitRotationTables()

   ! 6) Adsorbatos y buffers (como legacy)
   call read_adsorbates('MOLEC.DAT', SIGMA, 1.0e-7_rk)
   if (.not.allocated(Z))    allocate(Z(NMOLEC)); Z = 0.0_rk
   if (.not.allocated(N))    allocate(N(NMOLEC)); N = 0
   if (.not.allocated(ANPROM)) allocate(ANPROM(NMOLEC)); ANPROM = 0.0_rk

   if (.not. allocated(EPSI)) allocate(EPSI(maxAtoms))
   if (.not. allocated(SIGM)) allocate(SIGM(maxAtoms))
   if (.not. allocated(Q))    allocate(Q(maxAtoms))
   if (.not. allocated(RX1))  allocate(RX1(maxAtoms))
   if (.not. allocated(RY1))  allocate(RY1(maxAtoms))
   if (.not. allocated(RZ1))  allocate(RZ1(maxAtoms))
   if (.not. allocated(RX))   allocate(RX(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RY))   allocate(RY(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RZ))   allocate(RZ(5000, maxAtoms, NMOLEC))
   if (.not. allocated(USS))  allocate(USS(5000, maxAtoms, maxAtoms))

   ! 7) Estructura (superficie)
   open(unit=50, file='SALIDAACTIVADO-100.TXT')
   open(unit=97, file='PERFILES.TXT')
   call estructura(eps, nam, SIGMA, sigmetano, NC, diel)

   ! 8) Potenciales
   call potencialff(eps, SIGMA, sigmetano, NC, RCUT, diel)
   call potencial  (eps, SIGMA, sigmetano, NC, RCUT, diel)

   ! 9) Reinicio si ensemble = 2
   if (ensemble == 2) then
      do i = 1, NMOLEC
         N(i)      = 0
         NTOTAL(i) = 0
      end do
      V  = 0.0_rk
      VG = 0.0_rk
      VA = 0.0_rk
   end if

   ! 10) Inicialización si ensemble < 2 (CANONICAL): leer initconf.txt + energías
   if (ensemble < 2) then
      open(unit=67, file='initconf.txt', status='old', action='read')
      read(67,*) V, VG, VA
      read(67,*) nmolec2
      if (nmolec2 /= NMOLEC) then
         write(*,*) 'ERROR: NMOLEC difiere entre initconf y entrada.'
         stop
      end if

      do i = 1, NMOLEC
         molkind = i
         read(67,*) natom2
         if (natom2 /= NATOM(i)) then
            write(*,*) 'ERROR: NATOM difiere en especie ', i
            stop
         end if
         read(67,*) ncantmol
         do j = 1, ncantmol
            N(i) = j - 1
            call lee_conf_una_molecula(j, i, NATOM(i))
            call ADD(molkind)
            N(i) = j
         end do
      end do
      close(67)

      ! energías iniciales por molécula
      do i = 1, NMOLEC
         do j = 1, N(i)
            ipull = LOCATE(j,i)
            call ADPOTOUT(ipull, i, DELTVA)
            call POTOUT  (ipull, i, DELTV)
         end do
      end do
   end if

   ! 11) Bucle principal de isoterma
   do IPASOS = 1, isot

      ! limpiar CNF
      do i = 1, NMOLEC
         do natomki = 1, NATOM(i)
            do icnf = -NCELLMAT/2, NCELLMAT/2
               do jcnf = -NCELLMAT/2, NCELLMAT/2
                  do kcnf = -NCELLMAT/2, NCELLMAT/2
                     CNF(icnf, jcnf, kcnf, i, natomki) = 0.0_rk
                  end do
               end do
            end do
         end do
      end do

      ! etiquetas de archivos como legacy
      write(CONFIG,'(I3)') IPASOS
      open(unit=40, file='CONFIG'//config//'.xyz')
      open(unit=41, file='CONFIG'//config//'.TXT')

      ! actividad Z (gas/líquido) – exactamente como legacy
      if (NESTADO == 1) then
         do i = 1, NMOLEC
            auvol   = ((sigmetano/SIGMA)**3)*VOL
            Z(i)    = X(i)*P/(8.3144_rk*T)*6.023e-7_rk*auvol
            write(*,*) Z(i), 'Z ', i
         end do
      else
         do i = 1, NMOLEC
            Z(i)    = X(i)*P*6.023e-4_rk*((ACEL)**3)*VOL
            Z(NMOLEC)  = 55.55_rk*6.023e-4_rk*((ACEL)**3)*VOL
            write(*,*) Z(i), 'Z ', i
         end do
      end if

      ! acumuladores de esta pasada
      V = V; VG = VG; VA = VA
      Uacc=0; UGacc=0; UAacc=0; U2acc=0
      UNacc=0; UNGacc=0; UNAacc=0
      ANacc=0; N2acc=0

      do JPAS = 1, ijpasos
         MULT = 1
         if (JPAS == 1) MULT = mult2

         open(unit=51, file='Ener'//config//'.TXT')
         ! reescalado a ε para los pasos internos
         V  = V  / eps
         VG = VG / eps
         VA = VA / eps
         W  = 0.0_rk
         CR = 0.0_rk

         !write(*,*) 'aqui 1'
         
         do KPAS = 1, ikpasos*MULT
            ! ========== selección legacy con goto ==========
            IJ = int(ranf(dummy)*3.0) + 1
            !write(*,*) IJ, 'aqui 2'
            if (ensemble == 0) goto 30
            goto (10,20,30,35) IJ

10          call IN(TEMPSTAR, Z, SIGMA, eps, RCUT, V, VA, VG, W, CREATE, CR, JPAS, &
                 canonicalmolecules)
           ! write(*,*) 'IN'
            goto 40

20          call OUT(TEMPSTAR, Z, SIGMA, eps, RCUT, V, VA, VG, W, GHOST, JPAS, &
                     NMIN, NMAXI, canonicalmolecules)
            !write(*,*) 'out'
            goto 40

30          call MOVE(TEMPSTAR, Z, SIGMA, eps, RCUT, V, VA, VG, W, GHOST, JPAS)
            !write(*,*) 'move'
            goto 40

35          call CHANGE(TEMPSTAR, Z, SIGMA, eps, RCUT, V, VA, VG, W, GHOST, JPAS)
            !write(*,*) 'change'
            
40          continue
         end do   ! KPAS

         ! reescalar de vuelta
         V  = V  * eps
         VG = VG * eps
         VA = VA * eps

         ! acumular conteos
         do i = 1, NMOLEC
            NTOTAL(i) = NTOTAL(i) + N(i)
         end do

         call ESTADISTICA(NCELLMAT)

         Uacc  = Uacc  + V
         UGacc = UGacc + VG
         UAacc = UAacc + VA
         ntotalB = 0
         do i = 1, NMOLEC
            ntotalB = ntotalB + N(i)
         end do
         UNacc  = UNacc  + V * ntotalB
         UNGacc = UNGacc + VG * ntotalB
         UNAacc = UNAacc + VA * ntotalB
         U2acc  = U2acc  + V*V
         ANacc  = ANacc  + ntotalB
         N2acc  = N2acc  + ntotalB*ntotalB

         close(51)
      end do   ! JPAS

      ! guardar initconf de salida (igual que legacy)
      open(unit=67, file='initconf.txt')
      write(67,*) V, VG, VA, ' energias'
      write(67,*) NMOLEC,    ' tipos de molec'
      do i = 1, NMOLEC
         molkind = i
         write(67,*) NATOM(i)
         write(67,*) N(i), ' natom ', i
         do j = 1, N(i)
            ipull = LOCATE(j,i)
           
            do natomki = 1, NATOM(i)
               X1 = RX(ipull,natomki,i)
               Y1 = RY(ipull,natomki,i)
               Z1 = RZ(ipull,natomki,i)
               write(67,*) X1, Y1, Z1, 'x y z molec'
            end do
         end do
      end do
      close(67)

      ! promedios por especie
      do i = 1, NMOLEC
         ANPROM(i) = real(NTOTAL(i), rk) / real(ijpasos, rk)
      end do
      write(*,*) ANPROM(1:NMOLEC)

      ! calores isostéricos (idéntico a legacy)
      AN1  = ANacc / real(ijpasos-1, rk)
      U1   = Uacc  / real(ijpasos-1, rk)
      UNG1 = UNGacc/ real(ijpasos-1, rk)
      UG1  = UGacc / real(ijpasos-1, rk)
      UNA1 = UNAacc/ real(ijpasos-1, rk)
      UA1  = UAacc / real(ijpasos-1, rk)
      UN1  = UNacc / real(ijpasos-1, rk)
      U2acc= U2acc / real(ijpasos-1, rk)
      AN2  = real(N2acc, rk) / real(ijpasos-1, rk)

      ANN  = AN2 - AN1**2

      CALOR  = 8.3144_rk*T - ((UN1 - U1*AN1) / max(ANN,1.0e-30_rk))*8.31_rk
      CALORG = (UNG1 - UG1*AN1) / max(ANN,1.0e-30_rk) * 8.31_rk
      CALORA = -((UNA1 - UA1*AN1)/ max(ANN,1.0e-30_rk)) * 8.31_rk

      CALORESP1 = (U2acc - U1*U1)*(8.31_rk**2) - (8.31_rk*U1*AN1)**2 / max(ANN,1.0e-30_rk)
      CALORESP2 = (UN1 - AN1*8.31_rk*T**2)
      CALORESP3 = CALORESP1 / max(CALORESP2,1.0e-30_rk)

      ! salida gráfica
      open(unit=49, file='truncado.txt', status='old', action='read')
      read(49,*) NC
      ntotalGraf = NC
      do i = 1, NMOLEC
         ntotalGraf = ntotalGraf + N(i)*NATOM(i)
      end do

      write(40,*) ntotalGraf
      write(40,*) ' '

      do i = 1, NC
         read(49,*) RXAI, RYAI, RZAI, EPSAI, SGCI, QACI, SYMBOL2
         write(40,*) SYMBOL2, RXAI, RYAI, RZAI
      end do
      close(49)

      do i = 1, NMOLEC
         do j = 1, N(i)
            ipull = LOCATE(j,i)
            do natomki = 1, NATOM(i)
               RXN = RX(ipull,natomki,i)*ACEL
               RYN = RY(ipull,natomki,i)*ACEL
               RZN = RZ(ipull,natomki,i)*ACEL
               write(40,*) NSYM(natomki,i), RXN, RYN, RZN
            end do
         end do
      end do

      do i = 1, NMOLEC
         write(*,*) i, N(i), NATOM(i), ' i natom main'
      end do

      ! volcados CNF (igual legacy)
      do i = 1, NMOLEC
         write(CONFAT,'(I0)') i
         do natomki = 1, NATOM(i)
            write(CONFNAT,'(I0)') natomki
            open(unit=101, file='CNF'//trim(CONFIG)//'-'//trim(adjustl(CONFAT))//'-'// &
                 trim(adjustl(CONFNAT))//'.TXT')
            write(101,*) NCELLMAT
            do icnf = -NCELLMAT/2, NCELLMAT/2
               do jcnf = -NCELLMAT/2, NCELLMAT/2
                  do kcnf = -NCELLMAT/2, NCELLMAT/2
                     write(101,*) icnf, jcnf, kcnf, CNF(icnf,jcnf,kcnf,i,natomki)
                  end do
               end do
            end do
            close(101)
         end do
      end do

      do i = 1, NMOLEC
         NTOTAL(i) = 0
      end do

      write(*,*) AN1, P, CALOR, ' CALOR'
      write(50,*) P, ANPROM(1:NMOLEC), CALOR, CALORA, CALORESP3
      write(*,*)  P, ANPROM(1:NMOLEC), CALOR, 8.31_rk*T - CALORG, CALORA

      close(40)
      close(41)

      ! siguiente presión
      P = p_vals(IPASOS+1)
   end do  ! IPASOS

   close(50)
   close(97)

   ! limpieza
   if (allocated(Z))      deallocate(Z)
   if (allocated(ANPROM)) deallocate(ANPROM)
   if (allocated(p_vals)) deallocate(p_vals)
   if (allocated(SYMBOL)) deallocate(SYMBOL)
   if (allocated(UADS))   deallocate(UADS)
   if (allocated(CNF))    deallocate(CNF)
   if (allocated(RX))     deallocate(RX)
   if (allocated(RY))     deallocate(RY)
   if (allocated(RZ))     deallocate(RZ)
   if (allocated(RX1))    deallocate(RX1)
   if (allocated(RY1))    deallocate(RY1)
   if (allocated(RZ1))    deallocate(RZ1)
   if (allocated(USS))    deallocate(USS)
   if (allocated(QAC))    deallocate(QAC)
   if (allocated(RXC))    deallocate(RXC)
   if (allocated(RYC))    deallocate(RYC)
   if (allocated(RZC))    deallocate(RZC)
   if (allocated(EPSAC))  deallocate(EPSAC)
   if (allocated(SGC))    deallocate(SGC)
   if (allocated(EPSI))   deallocate(EPSI)
   if (allocated(SIGM))   deallocate(SIGM)
   if (allocated(Q))      deallocate(Q)
   if (allocated(N))      deallocate(N)
   if (allocated(X))      deallocate(X)
   if (allocated(NATOM))  deallocate(NATOM)

   call cpu_time(endtime)
   elapsedtime = endtime - starttime
   write(*,*) 'Total execution time (seconds): ', elapsedtime

contains
   subroutine lee_conf_una_molecula(jmol, imol, natm)
      integer, intent(in) :: jmol, imol, natm
      integer :: k
      do k = 1, natm
         read(67,*) X1, Y1, Z1
         RX1(k) = X1
         RY1(k) = Y1
         RZ1(k) = Z1
      end do
   end subroutine lee_conf_una_molecula
end program main

