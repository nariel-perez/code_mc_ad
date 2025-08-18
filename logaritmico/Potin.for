C------------------------------------------------------------------------------


      SUBROUTINE POTIN ( MOLKIND, SIGMA,RCUT, RMIN, DELTV,
     :            OVRLAP )
      implicit none
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS, ACEL, acelx,acely,acelz,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE ON ADDING AN ATOM.        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE ADDITION **
C    ** REAL    RXI,RYI,RZI       THE COORDINATES OF THE ADDED ATOM   **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    ** REAL    RMIN              MINIMUM ALLOWED APPROACH OF ATOMS   **
C    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL ADDITION OF AN ATOM TO THE FLUID. THE LONG     **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************
      REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW

      REAL RXC(9000),RYC(9000),RZC(9000)
      REAL EPSAC(9000),SGC(9000),QAC(9000)
        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, RMIN, SIGMA, RXI, RYI, RZI, DELTV, DELTW
        INTEGER     NTRIAL
        LOGICAL     OVRLAP

        REAL        RCUTSQ, RMINSQ, SIGSQ, SR2, SR6, SR3, SR9
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )
        INTEGER I1
        INTEGER NATOM(10)
        REAL RX1(50),RY1(50),RZ1(50)
        INTEGER IPOT,JPOT
        INTEGER NATOMKIND(50,10)
        INTEGER I
        INTEGER NMOLEC
        INTEGER N(10)
        INTEGER LOCATE(5000,10)
        INTEGER MOLKIND
        INTEGER K
        REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)
        REAL BCX,BCY,BCZ ,ACEL
 
	real acelx,acely,acelz

        REAL RIJ
        INTEGER IDIST 
        REAL UADS(-100:100,-100:100,-100:100,50)
        REAL USS(5000,50,50)
        REAL EPSI(50),SIGM(50),Q(50)
        REAL RX0(50,10),RY0(50,10),RZ0(50,10)
        INTEGER MAT
        LOGICAL FLAG(1000)
        real xmax,ymax,zmax
!        WRITE(*,*) 'MAT', MAT
        
C     ******************************************************************
!	write(*,*)ACEL, acelx,acely,acelz,' aceles******5'

!	pause	
       xmax=acelx/acel
       ymax=acely/acel
       zmax=acelz/acel

        OVRLAP = .FALSE.
        RCUTSQ = RCUT * RCUT
        RMINSQ = RMIN * RMIN
        SIGSQ  = SIGMA * SIGMA


C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0

C    ** LOOP OVER ALL ATOMS  **
!      write(*,*)'------------------------'
!      write(*,*) 'potin'
!      write(*,*) molkind,' molkind'
!      write(*,*)natom(molkind),' natom'
      
      DO I1=1,NATOM(MOLKIND)
          
          RXI=RX1(I1)
          RYI=RY1(I1)
          RZI=RZ1(I1)
          IPOT=NATOMKIND(I1,MOLKIND)
!          write(*,*) rxi,ryi,rzi, 'RX TARGET'
!          write(*,*) ipot,' ipot', i1,' atomo de la molecula'
          DO  I = 1,NMOLEC
!              write(*,*) nmolec,' nmolec'
!              write(*,*)i,n(i),' n(i)'
              DO J=1,N(I)
              
                  !     ** PICK ACTIVE ATOMS FROM THE ARRAY LOCATE **

                  JIN   = LOCATE(J,I)
!                  write(*,*)jin,' jin'
!                  write(*,*) i,natom(i),' natom i'
                  DO K=1,NATOM(I)
!                      write(*,*)'atom k  ',k
                      JPOT=NATOMKIND(K,I)
!                      write(*,*) jpot, ' JPOT'
                      RXIJ  = RXI - RX(JIN,K,I)
                      RYIJ  = RYI - RY(JIN,K,I)
                      RZIJ  = RZI - RZ(JIN,K,I)
 !                     WRITE(*,*)JIN,K,I,' JIN K I'
  !            write(*,*)RX(JIN,K,I),RY(JIN,K,I),RZ(JIN,K,I),'RX STAND'
                      RXIJ  = RXIJ - BCX*xmax*ANINT ( RXIJ/xmax )
                      RYIJ  = RYIJ - BCY*ymax*ANINT ( RYIJ/ymax )
                      RZIJ  = RZIJ - BCZ*zmax*ANINT ( RZIJ/zmax )
   !                   write(*,*) '1***********'
		      !write(*,*) xmax,ymax,zmax,' x y z'

                      RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
                      
                      RIJ=SQRT(RIJSQ)
                      
                      IDIST=RIJ*1000+1
                      
                      IF ( RIJ .LT. RMIN) THEN
                          OVRLAP = .TRUE.
                          DELTV=1E10
    !                      WRITE(*,*) IDIST,'IDIST'
                          !PAUSE
                          RETURN
                      ENDIF
			!write(*,*) '2***********'
 
			!write(*,*) IDIST,JPOT,IPOT

                      VIJ   = USS(IDIST,JPOT,IPOT)
			 !write(*,*) vij, ' vij*****'

                      !IF (VIJ.GT.0) then
                      !write(*,*) vij,idist,jpot,ipot, 'vij parcial'
                       !   end if
                      !write(*,*) deltv,' deltv antes'
                      DELTV = DELTV + VIJ
                      !write(*,*) deltv, ' deltv despues'
                      
                      DELTW = DELTW + WIJ
                  ENDDO
                  ENDDO
                  
              ENDDO
      ENDDO
      
        DELTV = DELTV
!        write(*,*) deltv,' deltv final'
!        DELTW = 48.0 * DELTW / 3.0

C    ** ADD CHANGE IN LONG RANGE CORRECTION **

!        DELTV = DELTV + ( 2.0 * REAL ( N ) + 1.0 ) * VLRC0
!        DELTW = DELTW + ( 2.0 * REAL ( N ) + 1.0 ) * WLRC0

        RETURN
        END
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
