!─────────────────────────────────────────────────────────────────
module EstructuraModule
   use InputParams,      only : cell             ! celda activa
   use PBC_Mod,          only : rk, cart_to_frac, wrap_frac_centered, frac_to_cart
   use PhysicalConstants 
   use SimulationData    
   implicit none
   private
   public :: estructura
contains
!=================================================================
subroutine estructura(eps, nam, sigma, sigmetano, NC, diel)
   implicit none
   !------------------ argumentos --------------------------------
   character(len=*), intent(in) :: nam
   real(rk),        intent(in)  :: eps, sigma, sigmetano, diel
   integer,         intent(inout) :: NC     ! se actualiza al número final

   !------------------ variables internas ------------------------
   integer :: i, ios, imax
   character(len=32) :: nampro
   real(rk) :: r_cart(3), s(3)
   real(rk), allocatable :: RXA(:), RYA(:), RZA(:)

   !------------------ abrir archivo de estructura ---------------
   open(unit=15, file=trim(nam), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: no se pudo abrir ', trim(nam)
      return
   end if

   read(15,*,iostat=ios) NC
   if (ios /= 0 .or. NC<=0) then
      write(*,*) 'ERROR al leer número de átomos en ', trim(nam)
      close(15); return
   end if

   !------------------ alocar arrays temporales ------------------
   allocate(RXA(NC), RYA(NC), RZA(NC))
   if (.not. allocated(RXC))    allocate(RXC(NC))
   if (.not. allocated(RYC))    allocate(RYC(NC))
   if (.not. allocated(RZC))    allocate(RZC(NC))
   if (.not. allocated(QAC))    allocate(QAC(NC))
   if (.not. allocated(EPSAC))  allocate(EPSAC(NC))
   if (.not. allocated(SGC))    allocate(SGC(NC))
   if (.not. allocated(SYMBOL)) allocate(SYMBOL(NC))

   !------------------ lectura y envoltura -----------------------
   imax = 0
   do i = 1, NC
      read(15,*,iostat=ios) r_cart(1), r_cart(2), r_cart(3), &
                            EPSAC(i), SGC(i), QAC(i), SYMBOL(i)
      if (ios /= 0) then
         write(*,*) 'ERROR leyendo línea ', i, ' en ', trim(nam)
         close(15); return
      end if

      ! → fraccional + wrap a [-0.5,0.5)
      s  = cart_to_frac(cell, r_cart)
      s  = wrap_frac_centered(cell, s)

      ! → de vuelta a cartesiano (ya dentro de la celda)
      r_cart = frac_to_cart(cell, s)

      imax = imax + 1
      RXC(imax) = r_cart(1) / sigmetano * sigma
      RYC(imax) = r_cart(2) / sigmetano * sigma
      RZC(imax) = r_cart(3) / sigmetano * sigma
      QAC(imax) = QAC(i) * FACTORELEC * FCLEC
   end do
   close(15)

   NC = imax               ! número real de átomos almacenados
   write(*,'(A,I0)') 'Estructura leída, ', NC, ' segmentos dentro de la celda.'

   !------------------ escribir archivo truncado -----------------
   nampro = 'truncado.txt'
   open(unit=16, file=nampro, status='replace', iostat=ios)
   if (ios == 0) then
      write(16,*) NC
      do i = 1, NC
         write(16,'(3F12.6,3F10.4,I4)') &
              RXC(i)*sigmetano/sigma, RYC(i)*sigmetano/sigma, RZC(i)*sigmetano/sigma, &
              EPSAC(i), SGC(i), QAC(i)/(FCLEC*FACTORELEC), SYMBOL(i)
      end do
      close(16)
   else
      write(*,*) 'AVISO: no se pudo escribir ', nampro
   end if

   !------------------ liberar temporales ------------------------
   deallocate(RXA, RYA, RZA)
end subroutine estructura
!=================================================================
end module EstructuraModule

