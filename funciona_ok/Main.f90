!******************************************************************************
!* FICHE F.13.  THE HEART OF A CONSTANT MU VT MONTE CARLO PROGRAM             **
!* This FORTRAN code is intended to illustrate points made in the text.       **
!* To our knowledge it works correctly.  However it is the responsibility of  **
!* the user to test it, if it is to be used in a research application.        **
!******************************************************************************

!C    *******************************************************************
!C    ** ATTEMPTED CREATIONS AND DESTRUCTIONS IN GRAND CANONICAL MC.   **
!C    **                                                               **
!C    ** THESE ROUTINES ALLOW FOR A TRIAL DESTRUCTION OR CREATION IN A **
!C    ** GRAND CANONICAL MONTE CARLO PROGRAM.                          **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                   NUMBER OF ATOMS BEFORE THE TRIAL  **
!C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING THE TRIAL  **
!C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS ALLOWED   **
!C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
!C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF ATOM     **
!C    ** REAL    RX(NMAX) ETC.       POSITIONS OF CURRENT ATOMS        **
!C    ** REAL    V                   POTENTIAL ENERGY + LRC            **
!C    ** REAL    W                   VIRIAL + LRC                      **
!C    ** REAL    DELTV               CHANGE IN ENERGY                  **
!C    ** REAL    DELTW               CHANGE IN VIRIAL                  **
!C    ** REAL    TEMP                REDUCED TEMPERATURE               **
!C    ** REAL    Z                   ABSOLUTE ACTIVITY COEFFICIENT     **
!C    ** REAL    SIGMA               LENNARD JONES DIAMETER            **
!C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
!C    ** REAL    RMIN                REDUCED MINIMUM SEPARATION        **
!C    ** LOGICAL OVRLAP              TRUE FOR SUBSTANTIAL ATOM OVERLAP **
!C    ** LOGICAL CREATE              TRUE FOR AN ACCEPTED CREATION     **
!C    ** LOGICAL GHOST               TRUE FOR AN ACCEPTED DESTRUCTION  **
!C    **                                                               **
!C    ** ROUTINES SUPPLIED:                                            **
!C    **                                                               **
!C    ** SUBROUTINE IN ( TEMP, Z, SIGMA, RCUT, N, V, W, CREATE )       **
!C    **    PERFORMS A TRIAL CREATION                                  **
!C    ** SUBROUTINE OUT ( TEMP, Z, SIGMA, RCUT, N, V, W, GHOST )       **
!C    **    PERFORMS A TRIAL DESTRUCTION                               **
!C    ** SUBROUTINE POTIN ( RXNEW, RYNEW, RZNEW, N, SIGMA, RCUT, RMIN, **
!C    ** :                  DELTV, DELTW, OVRLAP )                     **
!C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON CREATION         **
!C    ** SUBROUTINE POTOUT ( IPULL, N, SIGMA, RCUT, DELTV, DELTW )     **
!C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON DESTRUCTION      **
!C    **                                                               **
!C    ** ROUTINES REFERENCED:                                          **
!C    **                                                               **
!C    **    REAL FUNCTION RANF ( DUMMY )  (GIVEN IN F.11)              **
!C    **    RETURNS A UNIFORM RANDOM VARIATE ON ZERO TO ONE            **
!C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
!C    **    UPDATES LOCATE AFTER ADDITION (GIVEN IN F.14)              **
!C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
!C    **    UPDATES LOCATE AFTER REMOVAL (GIVEN IN F.14)               **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** ROUTINES IN AND OUT SHOULD BE CALLED WITH EQUAL PROBABILITY   **
!C    ** IN A GRAND CANONICAL MONTE CARLO SIMULATION. IF A TRIAL       **
!C    ** CREATION IS ACCEPTED THEN CREATE IS SET TO TRUE. IF A TRIAL   **
!C    ** DESTRUCTION IS ACCEPTED THEN GHOST IS SET TO TRUE. THE        **
!C    ** ROUTINES ARE WRITTEN FOR LENNARD-JONES ATOMS. THE BOX IS OF   **
!C    ** UNIT LENGTH, ALL DISTANCES ARE SCALED TO THE BOX LENGTH.      **
!C    ** TRIAL INPUTS WHICH RESULT IN A SEPARATION OF LESS THAN        **
!C    ** 0.5*SIGMA ARE REJECTED. THE LONG-RANGE CORRECTIONS ARE        **
!C    ** INCLUDED IN V AND W. ALL ACCUMULATORS ARE UPDATED IN THE MAIN **
!C    ** PART OF THE PROGRAM WHICH IS NOT GIVEN HERE.                  **
!C    *******************************************************************

program main
  use InputParams
  use simulationdata
  use AdsorbateInput
  use physicalconstants
  use estructuramodule
  use rotationmodule
  
   implicit none
   
   ! ================================================================
   ! ============== Declaración de variables ========================
   ! ================================================================

   
   real               p_ratio
   real, allocatable:: p_vals(:), Z(:)
   integer           NTOTALP, NTOTALB
   real              V, VA, VG, u2
   real              W
   real              DELTV
   real              TEMP
   real              SIGMA
   real              RCUT
   real              ANPROM(10)
   logical           CREATE
   logical           GHOST
   character(len=3)  CONFIG
   character(len=3)  CONFAT
   character(len=3)  CONFNAT
   character(len=32) archivo_truncado
   real              auvol
   integer           IPASOS, JPASOS, KPASOS
   integer           MULT
   real              PRED
   real              VOL, XMAX, YMAX, ZMAX
   integer           NC
   integer           I, J
   integer           NKIND, KATOM
   real              X1, Y1, Z1
   integer           INMOLEC
   real              CR
   integer           IJ
   real              DUMMY
   real              RANF
   real              U, UG, UA, UN, UNG, UNA, AN, ANN
   real              N2
   real              AN1, U1, UNG1, UG1, UNA1, UA1, UN1, AN2
   real              CALOR, CALORG, CALORA
   real              RXAI, RYAI, RZAI, EPSAI, SGCI, QACI
   integer           SYMBOL2
   integer           K, jin
   real              RXN, RYN, RZN
   real              starttime, endtime, elapsedtime
   integer           nmolec2
   integer           natom2
   integer           molkind
   integer           ncantmol
   real              deltva
   integer           ipull
   integer           ICNF, JCNF, KCNF
   integer           ipos_dot
   integer           NATOMKINDI
   real              CALORESP1, CALORESP2, CALORESP3
   integer           ntotalGRAF
   real, parameter   :: AK_input = 8.31   ! constante de los gases en j/mol·K
   integer           auxmat

   
   
   ! ================================================================
   ! =================== Inicio del programa =========================
   ! ================================================================

   call cpu_time(starttime)  ! tiempo inicial

   ! ----------------------------------------------------------------
   ! LECTURA DE ENTRADA
   ! ----------------------------------------------------------------

   call read_input('input.txt')  

   call print_params()

   auxmat = int(mat/2)
   allocate(uads(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, 50))
   
   !---------------------------------------
   !  pasos logaritmicos
   !--------------------------------------
   
   if (.not. allocated(p_vals)) allocate(p_vals(isot +1 ))
   p_vals = 0 
   
   p_ratio = real((dp/p)**(1.0d0/ real(isot -1, kind = 8)), kind=4)

   p_vals(1) = p
   do i = 2, isot
      p_vals(i) = p_vals(i-1)*p_ratio
   end do
   
   
   ! ----------------------------------------------------------------
   ! TRANSFORMACIÓN A UNIDADES REDUCIDAS
   ! ----------------------------------------------------------------
   SIGMA = sigmetano / acel
   TEMP  = T / eps
   PRED  = P * SIGMA**3 / eps
   RCUT  = 10*SIGMA
   
   XMAX = acelx / acel
   YMAX = acely / acel
   ZMAX = acelz / acel
   VOL  = XMAX * YMAX * ZMAX
   
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '--------------REDUCED UNITS ---------------------'
   write(*,'(A, F10.4)') 'TEMPERATURE: ', TEMP
   write(*,'(A, ES12.5)') 'PRESSURE: ', PRED
   
   write(*,'(A, F10.4)') 'SIGMA:', SIGMA
   write(*,'(A, ES12.4)') 'VOLUME:', VOL

   
   ! ----------------------------------------------------------------
   ! CÁLCULO DE CONSTANTES ELÉCTRICAS Y TABLAS DE ROTACIÓN
   ! ----------------------------------------------------------------
   call computeconstants(sigmetano, SIGMA, eps, AK_input, diel)
   call initrotationtables()

   ! ----------------------------------------------------------------
   ! LECTURA DE ESTRUCTURAS MOLECULARES
   ! ----------------------------------------------------------------
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '----------------ADSORBATES-----------------------'

   call read_adsorbates('MOLEC.DAT',sigma,  1.0e-7)
   
   ! ----------------------------------------------------------------
   ! LECTURA DE NKIND DESDE LJ.dat (para dimensionar tablas por tipo)
   ! ----------------------------------------------------------------
   block
      integer :: ios, iu
      real :: redelec
      character(len=256) :: linebuf
      iu = 12
      open(unit=iu, file='LJ.dat', status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'ERROR: no se pudo abrir LJ.dat'
         stop
      end if
      read(iu, '(A)', iostat=ios) linebuf
      if (ios /= 0) then
         write(*,*) 'ERROR: no se pudo leer la primera linea de LJ.dat'
         close(iu)
         stop
      end if
      redelec = 0.0
      read(linebuf, *, iostat=ios) NKIND, redelec
      if (ios /= 0) then
         read(linebuf, *, iostat=ios) NKIND
      end if
      close(iu)
      if (ios /= 0 .or. NKIND <= 0) then
         write(*,*) 'ERROR: NKIND invalido en LJ.dat'
         stop
      end if
   end block
   
   ! Validar que los tipos atomicos de MOLEC.DAT existan en LJ.dat
   do I = 1, NMOLEC
      do KATOM = 1, NATOM(I)
         if (NATOMKIND(KATOM, I) < 1 .or. NATOMKIND(KATOM, I) > NKIND) then
            write(*,*) 'ERROR: tipo atomico fuera de rango en MOLEC.DAT'
            write(*,*) 'Molecula:', I, 'Atomo:', KATOM, 'Tipo:', NATOMKIND(KATOM, I), 'NKIND:', NKIND
            stop
         end if
      end do
   end do
   
   if (.not. allocated(Z))    allocate(Z(NMOLEC))
   if (.not. allocated(N))    allocate(N(NMOLEC))
   
   if (.not. allocated(EPSI)) allocate(EPSI(NKIND))
   if (.not. allocated(SIGM)) allocate(SIGM(NKIND))
   if (.not. allocated(Q))    allocate(Q(NKIND))
   
   if (.not. allocated(RX1))  allocate(RX1(maxAtoms))
   if (.not. allocated(RY1))  allocate(RY1(maxAtoms))
   if (.not. allocated(RZ1))  allocate(RZ1(maxAtoms))
   if (.not. allocated(RX))   allocate(RX(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RY))   allocate(RY(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RZ))   allocate(RZ(5000, maxAtoms, NMOLEC))
   if (.not. allocated(USS))  allocate(USS(5000, NKIND, NKIND))
   
   auxmat = int(NCELLMAT/2)
   
   if (.not. allocated(CNF)) then
      allocate(CNF(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, NMOLEC, maxAtoms))
   end if
   
   
   
   ! ----------------------------------------------------------------
   ! LLAMADA A SUBRUTINAS DE POTENCIALES
   ! ----------------------------------------------------------------
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '------------------SURFACE------------------------'
   
   open(unit=50, file='SALIDAACTIVADO-100.TXT')
   open(unit=97, file='PERFILES.TXT')
   call estructura(eps, nam, sigma, sigmetano, NC, diel)
   ipos_dot = index(trim(nam), '.', back=.true.)
   if (ipos_dot > 0) then
      archivo_truncado = trim(nam(1:ipos_dot-1)) // '_truncado.txt'
   else
      archivo_truncado = trim(nam) // '_truncado.txt'
   end if
   
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '------------POTENTIALS--------------------'
   write(*,*)
   
   call potencialff(eps, sigma, sigmetano, NC, RCUT, diel)
   call potencial(eps, sigma, sigmetano, NC, RCUT, diel)
   
   ! -------------------------------------------------------------
   ! REINICIO DE VALORES SI ensemble = 2
   ! -------------------------------------------------------------
   if (ensemble.eq.2) then
      do I = 1, NMOLEC
         N(I)      = 0
         NTOTAL(I) = 0
         NTOTALP   = 0
         do J = 1, 5000
            LOCATE(J,I) = 0
         end do
      end do
      V  = 0.
      VG = 0.
      VA = 0.
   end if
   
   ! -------------------------------------------------------------
   ! INICIALIZACIÓN SI ensemble < 2 (CANONICAL)
   ! -------------------------------------------------------------
   if (ensemble .lt. 2) then
      open(unit=67, file='initconf.txt')
      write(*,*) 'archivo abierto'
      read(67,*) V, VG, VA
      write(*,*) V, VG, VA
      read(67,*) nmolec2
      write(*,*) nmolec2
      
      if (nmolec2 .ne. NMOLEC) then
         write(*,*) 'El número de moleculas no es el mismo'
         write(*,*) nmolec2, NMOLEC
         stop
      end if
      
      do i = 1, nmolec2
         molkind = i
         read(67,*) natom2
         write(*,*) natom2

         if (natom2 .ne. NATOM(i)) then
            write(*,*) 'El número de átomos no es el mismo'
            write(*,*) natom2, NATOM(i), i
            stop
         end if
         
         read(67,*) ncantmol
         write(*,*) ncantmol, i
         do j = 1, ncantmol
            n(i) = j - 1
            do k = 1, natom2
               read(67,*) X1, Y1, Z1
               write(*,*) X1, Y1, Z1, j, k
               RX1(k) = X1
               RY1(k) = Y1
               RZ1(k) = Z1
            end do
            call add(molkind)
            n(i) = j
         end do
      end do
      close(67)
      
      ! ----------------------------------------------------------
      ! CÁLCULO DE LA ENERGÍA INICIAL
      do I = 1, NMOLEC
         do J = 1, N(I)
            ipull = locate(J,I)
            call adpotout(ipull, I, deltva)
            call potout(ipull, I, deltv)
         end do
      end do
      
      write(*,*) VA, VG, V, ' Energias Iniciales'
      do I = 1, NMOLEC
         anprom(I) = 0
         ntotal(I) = 0
      end do
   end if

   ! -------------------------------------------------------------
   ! BUCLE PRINCIPAL DE LA SIMULACIÓN (ISOT)
   ! -------------------------------------------------------------
   do IPASOS = 1, isot
      
      ! Limpieza de estadísticas espaciales en CNF
      do I = 1, NMOLEC
         do NATOMKINDI = 1, NATOM(I)
            do ICNF = -NCELLMAT/2, NCELLMAT/2
               do JCNF = -NCELLMAT/2, NCELLMAT/2
                  do KCNF = -NCELLMAT/2, NCELLMAT/2
                     CNF(ICNF, JCNF, KCNF, I, NATOMKINDI) = 0
                  end do
               end do
            end do
         end do
      end do
      
      write(CONFIG, '(I3)') IPASOS
      open(unit=40, file='CONFIG'//CONFIG//'.xyz')
      open(unit=41, file='CONFIG'//CONFIG//'.TXT')
      
      ! Conversión de presión en actividad
      if (NESTADO.eq.1) then
         do INMOLEC = 1, NMOLEC
            auvol = ((sigmetano/SIGMA)**3)*VOL
            Z(INMOLEC) = X(INMOLEC)*P/(8.3144*T)*6.023E-7*auvol
            write(*,*) Z(INMOLEC), 'Z ', INMOLEC
         end do
      else
         do INMOLEC = 1, NMOLEC
            Z(INMOLEC) = X(INMOLEC)*P*6.023E-4*((ACEL)**3)*VOL
            write(*,*) Z(INMOLEC), 'Z ', INMOLEC
         end do
         ! Forzar actividad de la última especie (agua) si aplica
         Z(NMOLEC) = 55.55*6.023E-4*((ACEL)**3)*VOL
      end if
      
      ! Calores isostéricos (acumuladores en cero)
      U   = 0
      UG  = 0
      UA  = 0
      UN  = 0
      UNG = 0
      UNA = 0
      AN  = 0
      N2  = 0
      u2  = 0
      
      ! Sub-bucle JPASOS
      do JPASOS = 1, ijpasos
         if (mod(JPASOS, 500) == 0) then
            write(*,*) JPASOS, N(1:NMOLEC)
         end if
         
         V  = V / eps
         VG = VG / eps
         VA = VA / eps
         W  = 0.
         CR = 0.
         MULT = 1
         
         if (JPASOS.eq.1) then
            MULT = mult2
         end if
         
         do KPASOS = 1, ikpasos*MULT
            ! Elección del paso
            if (ensemble == 0) then
               call move(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos)
            else
               IJ = int(ranf(DUMMY)*4) + 1
               select case (IJ)
               case (1)
                  call in(temp, z, sigma, eps, rcut, v, va, vg, w, create, cr, jpasos, &
                       canonicalmolecules)
               case (2)
                  call out(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos, &
                       canonicalmolecules)
               case (3)
                  call move(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos)
               case (4)
                  call change(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos)
               end select
            end if
         end do

         ! Reescalado de energías
         V  = V * eps
         VG = VG * eps
         VA = VA * eps
         
         ! Acumulación de números de partículas
         do I = 1, NMOLEC
            NTOTAL(I) = NTOTAL(I) + N(I)
            NTOTALP   = NTOTALP + N(I)
         end do
         
         call estadistica(NCELLMAT)
         
         U  = U  + V
         UG = UG + VG
         UA = UA + VA
         
         ! estadística para un solo tipo o total
         NTOTALB = 0
         do I = 1, NMOLEC
            NTOTALB = NTOTALB + N(I)
         end do
         
         UN   = UN   + V*NTOTALB
         UNG  = UNG  + VG*NTOTALB
         UNA  = UNA  + VA*NTOTALB
         u2   = u2   + V**2
         AN   = AN   + NTOTALB
         N2   = N2   + NTOTALB**2

      end do  ! fin do JPASOS
      
      write(*,*) '-----'
      
      ! Guardar configuración final en initconf.txt
      open(unit=67, file='initconf.txt')
      write(67,*) V, VG, VA, ' energias'
      write(67,*) NMOLEC, ' tipos de molec'
      do i = 1, NMOLEC
         molkind = i
         write(67,*) NATOM(i)
         write(67,*) N(i), ' natom ', i
         do j = 1, N(i)
            jin = locate(j,i)
            do k = 1, NATOM(i)
               X1 = RX(jin,k,i)
               Y1 = RY(jin,k,i)
               Z1 = RZ(jin,k,i)
               write(67,*) X1, Y1, Z1, 'x y z molec'
            end do
         end do
      end do
      close(67)
      
      ! Promedios
      do I = 1, NMOLEC
         ANPROM(I) = real(NTOTAL(I)) / real(ijpasos)
      end do
      
      write(*,*) ANPROM(1:NMOLEC)
      
      ! Cálculo de calores isostéricos
      AN1 = AN / real(ijpasos - 1)
      U1  = U  / real(ijpasos - 1)
      write(*,*) U1, ' U1'
      
      UNG1 = UNG / real(ijpasos - 1)
      write(*,*) UNG1,' UNG1'
      
      UG1 = UG / real(ijpasos - 1)
      write(*,*) UG1, ' UG1'
      
      UNA1 = UNA / real(ijpasos - 1)
      write(*,*) UNA1, ' UNA1'
      
      UA1 = UA / real(ijpasos - 1)
      write(*,*) UA1,'UA1'
      
      UN1 = UN / real(ijpasos - 1)
      write(*,*) UN1,' UN1'
      
      u2 = u2 / real(ijpasos - 1)

      AN2 = real(N2) / real(ijpasos - 1)
      write(*,*) AN2, ' AN2'
      
      ANN = AN2 - AN1**2
      write(*,*) ANN, ' ANN'

      CALOR  = 8.3144*T - ((UN1 - U1*AN1)/ (ANN))*8.31
      CALORG = (UNG1 - UG1*AN1)/ (ANN)*8.31
      CALORA = -((UNA1 - UA1*AN1)/ (ANN))*8.31

      CALORESP1 = (u2 - U1**2)*(8.31**2) - (8.31*U1*AN1)**2 / ANN
      CALORESP2 = (UN1 - AN1*8.31*T**2)
      CALORESP3 = CALORESP1 / CALORESP2
      
      open(unit=49, file=archivo_truncado)
      read(49,*) NC
      ntotalGRAF = NC
      do I = 1, NMOLEC
         ntotalGRAF = ntotalGRAF + N(I)*NATOM(I)
      end do

      write(40,*) ntotalGRAF
      write(40,*) ' '
      write(*,*) ntotalGRAF
      write(*,*) ' '
      
      write(*,*) acelx, acely, acelz
      
      do I = 1, NC
         read(49,*) RXAI, RYAI, RZAI, EPSAI, SGCI, QACI, SYMBOL2
         write(40,*) SYMBOL2, RXAI, RYAI, RZAI
      end do
      close(49)

      do I = 1, NMOLEC
         do j = 1, N(I)
            jin = locate(j,I)
            do k = 1, NATOM(I)
               RXN = RX(jin,k,I)*acel
               RYN = RY(jin,k,I)*acel
               RZN = RZ(jin,k,I)*acel
               write(40,*) NSYM(k,I), RXN, RYN, RZN
            end do
         end do
      end do
      
      do i = 1, NMOLEC
         write(*,*) i, N(i), NATOM(i), ' i natom main'
      end do
      
      ! Estadística espacial CNF
      do I = 1, NMOLEC
         write(CONFAT, '(I0)') I
         CONFAT = adjustl(CONFAT)
         
         do NATOMKINDI = 1, NATOM(I)
            write(CONFNAT, '(I0)') NATOMKINDI
            CONFNAT = adjustl(CONFNAT)
            
            open(unit=101, file='CNF'//trim(CONFIG)//'-'//trim(CONFAT)//'-'// &
                 trim(CONFNAT)//'.TXT')
            write(101,*) NCELLMAT
            do ICNF = -NCELLMAT/2, NCELLMAT/2
               do JCNF = -NCELLMAT/2, NCELLMAT/2
                  do KCNF = -NCELLMAT/2, NCELLMAT/2
                     write(101,*) ICNF, JCNF, KCNF, &
                          CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
                  end do
               end do
            end do
            close(101)
         end do
      end do
      
      do I = 1, NMOLEC
         NTOTAL(I) = 0
      end do

      write(*,*) AN1, P, CALOR, ' CALOR'
      write(50,*) P, ANPROM(1:NMOLEC), CALOR, CALORA, CALORESP3
      
      write(*,*) P, ANPROM(1:NMOLEC), CALOR, 8.31*T - CALORG, CALORA
      close(40)
      close(41)
      
      AN  = 0
      AN1 = 0
      !P   = P + dp*1333.22
      p = p_vals(ipasos +1 )
   end do  ! fin IPASOS
   
   close(50)
   close(97)

   deallocate(qac)
   
   if (allocated(SYMBOL)) deallocate(SYMBOL)
   if (allocated(UADS))   deallocate(UADS)
   if (allocated(CNF))   deallocate(CNF)

   deallocate(Z)
   deallocate(X)
   deallocate(N)
   deallocate(NATOM)
   deallocate(rxc)
   deallocate(ryc)
   deallocate(rzc)
   deallocate(epsac)
   deallocate(sgc)
   deallocate(epsi, sigm, q)
   deallocate(rx1, ry1, rz1)
   deallocate(rx, ry, rz)
   deallocate(uss)
   deallocate(p_vals)
   call cpu_time(endtime)
   elapsedtime = endtime - starttime
   write(*,*) 'Total execution time (seconds): ', elapsedtime
   
end program main
