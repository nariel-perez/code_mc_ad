!--------------------------------------------------------------------
! test_pbc.f90 - Verificación de PBC_Mod y GeomUtils
!
! Comprueba que:
! 1. Celda ortorrómbica se construye correctamente (A diagonal)
! 2. min_image coincide con fórmula ortorrómbica dr - L*nint(dr/L)
! 3. wrap_by_pbc funciona correctamente
! 4. Celda triclínica simple se construye sin errores
! 5. cellR (celda reducida) se construye correctamente
!--------------------------------------------------------------------
program test_pbc
   use PBC_Mod, only: rk, Cell_t => Cell, cell_from_lengths_angles, min_image
   use GeomUtils
   use InputParams, only: ip_cell => cell, ip_cellR => cellR, &
                          acelx, acely, acelz, ACEL, &
                          angux, anguy, anguz, BCX, BCY, BCZ, dim_flag, &
                          update_cellR
   implicit none

   type(Cell_t) :: c
   real(rk) :: r1(3), r2(3), dr_min(3), dr_orto(3)
   real(rk) :: Lx, Ly, Lz
   real(rk) :: s1, s2, s3
   real(rk) :: G(3,3), Ainv(3,3)
   logical :: pbcx, pbcy, pbcz
   real(rk) :: tol, diff
   integer :: i, ntest, npass

   tol = 1.0e-10_rk
   ntest = 0
   npass = 0

   write(*,*) "=================================================="
   write(*,*) "    TEST DE PBC_Mod y GeomUtils"
   write(*,*) "=================================================="
   write(*,*)

   !-----------------------------------------------------------------
   ! TEST 1: Celda ortorrómbica (90, 90, 90)
   !-----------------------------------------------------------------
   write(*,*) "TEST 1: Construcción de celda ortorrómbica"
   write(*,*) "------------------------------------------"

   Lx = 40.0_rk
   Ly = 40.0_rk
   Lz = 20.0_rk

   call cell_from_lengths_angles(c, Lx, Ly, Lz, 90._rk, 90._rk, 90._rk, &
                                 dim=3, centered=.true.)
   c%pbc = [.true., .true., .false.]  ! PBC en x,y pero no en z (slab)

   write(*,'(A,3F10.4)') "  Longitudes (a,b,c):  ", c%a_len, c%b_len, c%c_len
   write(*,'(A,3F10.4)') "  Angulos (a,b,g):     ", c%alpha_deg, c%beta_deg, c%gamma_deg
   write(*,'(A,3L3)')    "  PBC (x,y,z):         ", c%pbc
   write(*,*)
   write(*,*) "  Matriz A (vectores de red en columnas):"
   do i = 1, 3
      write(*,'(A,3F12.6)') "    ", c%A(i,:)
   end do

   ! Verificar que A es diagonal para ortorrómbico
   ntest = ntest + 1
   if (abs(c%A(1,1) - Lx) < tol .and. abs(c%A(2,2) - Ly) < tol .and. &
       abs(c%A(3,3) - Lz) < tol .and. &
       abs(c%A(1,2)) < tol .and. abs(c%A(1,3)) < tol .and. &
       abs(c%A(2,1)) < tol .and. abs(c%A(2,3)) < tol .and. &
       abs(c%A(3,1)) < tol .and. abs(c%A(3,2)) < tol) then
      write(*,*) "  [PASS] Matriz A es diagonal correctamente"
      npass = npass + 1
   else
      write(*,*) "  [FAIL] Matriz A NO es diagonal como se esperaba"
   end if
   write(*,*)

   !-----------------------------------------------------------------
   ! TEST 2: min_image vs fórmula ortorrómbica (PBC en x,y,z)
   !-----------------------------------------------------------------
   write(*,*) "TEST 2: min_image vs formula ortorombica (PBC x,y,z)"
   write(*,*) "----------------------------------------------------"

   c%pbc = [.true., .true., .true.]  ! PBC en todos los ejes

   r1 = [5.0_rk, 5.0_rk, 5.0_rk]
   r2 = [35.0_rk, -15.0_rk, 18.0_rk]  ! Diferencia que cruza bordes

   ! min_image del módulo
   dr_min = min_image(c, r1, r2)

   ! Fórmula ortorrómbica directa: dr - L*nint(dr/L)
   dr_orto(1) = (r2(1) - r1(1)) - Lx * nint((r2(1) - r1(1)) / Lx)
   dr_orto(2) = (r2(2) - r1(2)) - Ly * nint((r2(2) - r1(2)) / Ly)
   dr_orto(3) = (r2(3) - r1(3)) - Lz * nint((r2(3) - r1(3)) / Lz)

   write(*,'(A,3F12.6)') "  r1:           ", r1
   write(*,'(A,3F12.6)') "  r2:           ", r2
   write(*,'(A,3F12.6)') "  dr = r2 - r1: ", r2 - r1
   write(*,'(A,3F12.6)') "  min_image:    ", dr_min
   write(*,'(A,3F12.6)') "  dr_orto:      ", dr_orto

   ntest = ntest + 1
   diff = sqrt(sum((dr_min - dr_orto)**2))
   if (diff < tol) then
      write(*,'(A,ES12.4)') "  [PASS] Diferencia: ", diff
      npass = npass + 1
   else
      write(*,'(A,ES12.4)') "  [FAIL] Diferencia: ", diff
   end if
   write(*,*)

   !-----------------------------------------------------------------
   ! TEST 3: min_image con slab (PBC solo en x,y)
   !-----------------------------------------------------------------
   write(*,*) "TEST 3: min_image con slab (PBC solo x,y)"
   write(*,*) "-----------------------------------------"

   c%pbc = [.true., .true., .false.]  ! Slab: no PBC en z

   r1 = [5.0_rk, 5.0_rk, 2.0_rk]
   r2 = [35.0_rk, -15.0_rk, 18.0_rk]

   dr_min = min_image(c, r1, r2)

   ! Para slab: wrap en x,y pero NO en z
   dr_orto(1) = (r2(1) - r1(1)) - Lx * nint((r2(1) - r1(1)) / Lx)
   dr_orto(2) = (r2(2) - r1(2)) - Ly * nint((r2(2) - r1(2)) / Ly)
   dr_orto(3) = r2(3) - r1(3)  ! Sin wrap en z

   write(*,'(A,3F12.6)') "  r1:           ", r1
   write(*,'(A,3F12.6)') "  r2:           ", r2
   write(*,'(A,3F12.6)') "  min_image:    ", dr_min
   write(*,'(A,3F12.6)') "  dr_orto:      ", dr_orto

   ntest = ntest + 1
   diff = sqrt(sum((dr_min - dr_orto)**2))
   if (diff < tol) then
      write(*,'(A,ES12.4)') "  [PASS] Diferencia: ", diff
      npass = npass + 1
   else
      write(*,'(A,ES12.4)') "  [FAIL] Diferencia: ", diff
   end if
   write(*,*)

   !-----------------------------------------------------------------
   ! TEST 4: wrap_by_pbc
   !-----------------------------------------------------------------
   write(*,*) "TEST 4: wrap_by_pbc"
   write(*,*) "-------------------"

   call cell_to_metric(c, G, Ainv, pbcx, pbcy, pbcz)

   ! Coordenadas fraccionales fuera de rango
   s1 = 1.3_rk   ! Debería ir a 0.3 - 1 = -0.7 -> pero nint(1.3)=1, s1-1=0.3
   s2 = -0.8_rk  ! Debería ir a -0.8 + 1 = 0.2 -> nint(-0.8)=-1, s2-(-1)=0.2
   s3 = 2.5_rk   ! NO debería cambiar (pbcz = false)

   write(*,'(A,3F12.6)') "  Antes wrap:  s1,s2,s3 = ", s1, s2, s3
   call wrap_by_pbc(s1, s2, s3, pbcx, pbcy, pbcz)
   write(*,'(A,3F12.6)') "  Despues wrap:s1,s2,s3 = ", s1, s2, s3

   ntest = ntest + 1
   ! s1 debe ser 0.3, s2 debe ser 0.2, s3 debe seguir siendo 2.5
   if (abs(s1 - 0.3_rk) < tol .and. abs(s2 - 0.2_rk) < tol .and. &
       abs(s3 - 2.5_rk) < tol) then
      write(*,*) "  [PASS] wrap_by_pbc funciona correctamente"
      npass = npass + 1
   else
      write(*,*) "  [FAIL] wrap_by_pbc no dio los valores esperados"
      write(*,'(A,3F12.6)') "  Esperado: ", 0.3_rk, 0.2_rk, 2.5_rk
   end if
   write(*,*)

   !-----------------------------------------------------------------
   ! TEST 5: Celda triclínica simple (verificar que no falla)
   !-----------------------------------------------------------------
   write(*,*) "TEST 5: Construccion de celda triclinica"
   write(*,*) "-----------------------------------------"

   call cell_from_lengths_angles(c, 10._rk, 10._rk, 10._rk, &
                                 80._rk, 85._rk, 75._rk, &
                                 dim=3, centered=.true.)
   c%pbc = [.true., .true., .true.]

   write(*,'(A,3F10.4)') "  Longitudes (a,b,c):  ", c%a_len, c%b_len, c%c_len
   write(*,'(A,3F10.4)') "  Angulos (a,b,g):     ", c%alpha_deg, c%beta_deg, c%gamma_deg
   write(*,*)
   write(*,*) "  Matriz A:"
   do i = 1, 3
      write(*,'(A,3F12.6)') "    ", c%A(i,:)
   end do
   write(*,*)
   write(*,*) "  Matriz Ainv:"
   do i = 1, 3
      write(*,'(A,3F12.6)') "    ", c%Ainv(i,:)
   end do

   ! Verificar que A * Ainv = I
   ntest = ntest + 1
   diff = 0._rk
   do i = 1, 3
      diff = diff + abs(sum(c%A(i,:) * c%Ainv(:,i)) - 1._rk)
   end do
   if (diff < tol) then
      write(*,*) "  [PASS] A * Ainv = I (diagonal = 1)"
      npass = npass + 1
   else
      write(*,'(A,ES12.4)') "  [FAIL] A * Ainv != I, error: ", diff
   end if
   write(*,*)

   !-----------------------------------------------------------------
   ! TEST 6: cellR (celda reducida) desde InputParams
   !-----------------------------------------------------------------
   write(*,*) "TEST 6: cellR (celda reducida) desde InputParams"
   write(*,*) "-------------------------------------------------"

   ! Configurar parámetros como si se hubiera leído un input
   acelx = 40.0_rk
   acely = 40.0_rk
   acelz = 20.0_rk
   ACEL = 40.0_rk
   angux = 90.0_rk
   anguy = 90.0_rk
   anguz = 90.0_rk
   BCX = 1
   BCY = 1
   BCZ = 0
   dim_flag = 3

   ! Construir ip_cell manualmente (como lo hace read_input)
   call cell_from_lengths_angles(ip_cell, acelx, acely, acelz, &
                                 angux, anguy, anguz, &
                                 dim=dim_flag, centered=.true.)
   ip_cell%pbc = [ BCX /= 0, BCY /= 0, BCZ /= 0 ]
   ip_cell%dim = dim_flag

   ! Construir ip_cellR
   call update_cellR()

   write(*,'(A,3F10.4)') "  cell longitudes:  ", ip_cell%a_len, ip_cell%b_len, ip_cell%c_len
   write(*,'(A,3F10.4)') "  cellR longitudes: ", ip_cellR%a_len, ip_cellR%b_len, ip_cellR%c_len
   write(*,'(A,F10.4)')  "  ACEL:             ", ACEL
   write(*,'(A,3F10.4)') "  Esperado cellR:   ", acelx/ACEL, acely/ACEL, acelz/ACEL

   ntest = ntest + 1
   ! ip_cellR debe tener longitudes = ip_cell / ACEL
   if (abs(ip_cellR%a_len - acelx/ACEL) < tol .and. &
       abs(ip_cellR%b_len - acely/ACEL) < tol .and. &
       abs(ip_cellR%c_len - acelz/ACEL) < tol) then
      write(*,*) "  [PASS] cellR tiene longitudes correctas"
      npass = npass + 1
   else
      write(*,*) "  [FAIL] cellR NO tiene longitudes correctas"
   end if

   ! Verificar que ip_cellR tiene los mismos PBC que ip_cell
   ntest = ntest + 1
   if (all(ip_cellR%pbc .eqv. ip_cell%pbc)) then
      write(*,*) "  [PASS] cellR tiene mismos PBC que cell"
      npass = npass + 1
   else
      write(*,*) "  [FAIL] cellR NO tiene mismos PBC que cell"
   end if
   write(*,*)

   !-----------------------------------------------------------------
   ! RESUMEN
   !-----------------------------------------------------------------
   write(*,*) "=================================================="
   write(*,'(A,I2,A,I2,A)') "    RESULTADO: ", npass, " / ", ntest, " tests pasaron"
   write(*,*) "=================================================="

   if (npass == ntest) then
      write(*,*) "    TODOS LOS TESTS PASARON"
   else
      write(*,*) "    ALGUNOS TESTS FALLARON"
   end if

end program test_pbc
