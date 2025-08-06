!─────────────────────────────────────────────────────────────────────
!  File: Main.f90                (Fortran 2003)
!  Compilación típica:
!     gfortran -c PBC_Mod.f90
!     gfortran -c InputParams.f90
!     gfortran -c Main.f90
!     gfortran -o test_pbc Main.o InputParams.o PBC_Mod.o
!─────────────────────────────────────────────────────────────────────
program test_pbc
   use InputParams, only : read_input,            &
                           a_len, b_len, c_len,   &
                           alpha_deg, beta_deg, gamma_deg, cell

   use InputParams, only : P, dp, sigmetano, eps, ACEL, acelx, acely, acelz, &
                        T, isot, cell,diel,nam,ncellmat,mat

   use PhysicalConstants, only : ComputeConstants
   use RotationModule,   only : InitRotationTables
   use EstructuraModule
   use SimulationData
   use AdsorbateInput
   
   use PBC_Mod,   only : rk, cart_to_frac, frac_to_cart, min_image
   implicit none

   !--- Variables para un test rápido --------------------------------
   real(rk) :: r1(3), r2(3), dr(3), s(3)
   !-----------------------------------
   real(rk) :: p_ratio, SIGMA, TEMP, PRED, RCUT
   real(rk), allocatable :: p_vals(:), Z(:)
   real(rk) :: detA, XMAX, YMAX, ZMAX, VOL
   integer  :: i, NC,auxmat
   ! --- constantes eléctricas -------------------------------------
   real(rk), parameter :: AK_input = 8.31_rk   ! Avogadro (valor original)

   
   !------------------------------------------------------------------
   call read_input('input.txt')

   print *, '================  CELDA LEÍDA  ================'
   write(*,'(A,3F12.6)') 'Longitudes (a,b,c)  : ', a_len, b_len, c_len
   write(*,'(A,3F12.6)') 'Ángulos    (α,β,γ)  : ', alpha_deg, beta_deg, gamma_deg
   print *, 'Matriz A (vectores columna):'
   write(*,'(3F14.6)') cell%A(1,:)
   write(*,'(3F14.6)') cell%A(2,:)
   write(*,'(3F14.6)') cell%A(3,:)
   print *, 'Matriz Ainv:'
   write(*,'(3F14.6)') cell%Ainv(1,:)
   write(*,'(3F14.6)') cell%Ainv(2,:)
   write(*,'(3F14.6)') cell%Ainv(3,:)

   !------------------------------------------------------------------
   ! Test de conversión cart ↔ frac y mínima-imagen
   !------------------------------------------------------------------
   r1 = [ 0.0_rk, 0.0_rk, 0.0_rk ]
   r2 = [ 0.6_rk*a_len, 0.6_rk*b_len, 0.6_rk*c_len ]   ! fuera de la celda

   s   = cart_to_frac(cell, r2)
   print *, 'Frac (sin wrap) de r2: ', s
   print *, 'Frac (wrap)           : ', s - nint(s)

   dr  = min_image(cell, r1, r2)
   write(*,'(A,3F12.6)') 'Vector mínima-imagen : ', dr
   write(*,'(A,F12.6)')  'Distancia            : ', sqrt(sum(dr*dr))

   !======================================================================
   ! 1. PASOS LOGARÍTMICOS DE PRESIÓN (array p_vals)
   !======================================================================
   if (.not. allocated(p_vals)) allocate(p_vals(isot))
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
   
   !======================================================================
   ! 2. TRANSFORMACIÓN A UNIDADES REDUCIDAS
   !======================================================================
   SIGMA = sigmetano / ACEL        ! longitud de referencia
   TEMP  = T         / eps         ! T* = kT/ε
   PRED  = P * SIGMA**3 / eps      ! P* = P σ³ / ε
   RCUT  = 10.0_rk * SIGMA         ! radio de corte en unidades absolutas
   
   !======================================================================
   ! 3. DIMENSIONES REDUCIDAS – versión general (funciona para orto y triclinic)
   !======================================================================
  
   
   ! --- volumen reducido (V* = V/ACEL³) ---
   detA =  cell%A(1,1)*(cell%A(2,2)*cell%A(3,3) - cell%A(2,3)*cell%A(3,2)) &
        - cell%A(2,1)*(cell%A(1,2)*cell%A(3,3) - cell%A(1,3)*cell%A(3,2)) &
        + cell%A(3,1)*(cell%A(1,2)*cell%A(2,3) - cell%A(1,3)*cell%A(2,2))
   
   VOL  = abs(detA) / ACEL**3        ! V*  (adimensional)
   
   ! --- longitudes reducidas de cada vector de red (opcional) ----------
   XMAX = sqrt(sum(cell%A(:,1)**2)) / ACEL
   YMAX = sqrt(sum(cell%A(:,2)**2)) / ACEL
   ZMAX = sqrt(sum(cell%A(:,3)**2)) / ACEL
   !======================================================================

   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '--------------REDUCED UNITS ---------------------'
   write(*,'(A, F10.4)') 'TEMPERATURE: ', TEMP
   write(*,'(A, ES12.5)') 'PRESSURE: ', PRED
   
   write(*,'(A, F10.4)') 'SIGMA:', SIGMA
   write(*,'(A, ES12.4)') 'VOLUME:', VOL
   
   
   !---------------------------------------------------------------
   !  CONSTANTES ELÉCTRICAS  +  TABLAS DE ROTACIÓN
   !---------------------------------------------------------------
   call ComputeConstants(sigmetano, SIGMA, eps, AK_input, diel)
   call InitRotationTables()

   !------------------------------------------------------------------
   ! LECTURA DE ESTRUCTURAS MOLECULARES (adsorbatos)
   !------------------------------------------------------------------
   write(*,*) ''
   write(*,*) '-------------------------------------------------'
   write(*,*) '----------------  ADSORBATES  -------------------'

   call read_adsorbates('MOLEC.DAT', sigma, 1.0e-7_rk)
   ! 1D arrays por especie
   if (.not. allocated(Z))    allocate(Z(NMOLEC))
   
   ! 1D arrays por átomo
   if (.not. allocated(EPSI)) allocate(EPSI(maxAtoms))
   if (.not. allocated(SIGM)) allocate(SIGM(maxAtoms))
   if (.not. allocated(Q))    allocate(Q(maxAtoms))
   if (.not. allocated(RX1))  allocate(RX1(maxAtoms))
   if (.not. allocated(RY1))  allocate(RY1(maxAtoms))
   if (.not. allocated(RZ1))  allocate(RZ1(maxAtoms))
   
   ! 3D por configuración / átomo / especie
   if (.not. allocated(RX))   allocate(RX(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RY))   allocate(RY(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RZ))   allocate(RZ(5000, maxAtoms, NMOLEC))
   if (.not. allocated(USS))  allocate(USS(5000, maxAtoms, maxAtoms))
   
   auxmat = int(NCELLMAT/2)
   if (.not. allocated(CNF)) then
      allocate(CNF(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, NMOLEC, maxAtoms))
   end if

   auxmat = int(mat/2)
   allocate(uads(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, 50))
   
   
  
   
   
   
   call estructura(eps, nam, sigma, sigmetano, NC, diel)

   call potencialff(eps, sigma, sigmetano, NC, RCUT, diel)
   call potencial(eps, sigma, sigmetano, NC, RCUT, diel)
   
end program test_pbc

