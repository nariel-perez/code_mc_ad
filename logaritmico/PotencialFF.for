C------------------------------------------------------------------------
C------------------------------------------------------------------------
C    ** SUBRUTINA POTENCIAL2
C    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZ
C    ** USS(1000),FLAG(1000)
C------------------------------------------------------------------------
	SUBROUTINE POTENCIALFF(EPS,sigma,sigmetano,NC,RCUT,diel)
      implicit none
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS, acel,acelx,acely,acelz,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE
!----------------------------------------------------
	real tamacel
            REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW

        INTEGER N(10)
       
        INTEGER MOLKIND
        REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)
        
        INTEGER NATOMKIND(50,10)
        
      REAL UADS(-100:100,-100:100,-100:100,50)
      REAL USS(5000,50,50)
      LOGICAL FLAG(1000)
      INTEGER LOCATE(5000,10)
      REAL RXC(9000),RYC(9000),RZC(9000)
      REAL EPSAC(9000),SGC(9000),QAC(9000)
      REAL RX1(50),RY1(50),RZ1(50)
      REAL EPSI(50),SIGM(50),Q(50)
      REAL RX0(50,10),RY0(50,10),RZ0(50,10)
      REAL BCX,BCY,BCZ   ,acel,acelx,acely,acelz
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
!-----------------------------------------------------------------------
      REAL PI
      	PARAMETER ( PI = 3.14159265 )
      REAL RCELE
      REAL FACTORELEC
      REAL AK
      INTEGER NKIND, INKIND
      INTEGER KINDI
      REAL AE0
      REAL FCLEC1AUX,FCLEC2AUX,FCLEC3AUX,FCLECAUX,FCLEC
      REAL SIGMETANO,SIGMA,EPS
      INTEGER IPOT,JPOT
      REAL SIGMA1,FACTOR
      REAL RCUTSQ,RCUT
      REAL SIGSQ,SIGCUB
      REAL RMIN,RMINSQ
      REAL SR3,SR9,SR2,SR6
      REAL VLRC0,WLRC0
      REAL DELTV,DELTW
      INTEGER I
      REAL RZI,RIJSQ
      REAL VIJ,WIJ
      REAL VIJE1
      INTEGER NC
      real diel
      real REDELEC
      
      
      
      
      tamacel=7/sigmetano/sigma
	
	write(*,*) tamacel, ' TAMACEL'
	

      RCELE=0.35 !RADIO DE CORTE POTENCIAL ELECTROSTÁTICO
C-----------------------------------------------------------------------
C	UNIDADES ELECTROSTÁTICAS EN TÉRMINOS DE EPS1 Y SIGMA
C-----------------------------------------------------------------------
	FACTORELEC=96500/6.023E23
      AK=8.31/6.023E23
      write(*,*)'POTENCIAL FLUIDO-FLUIDO'
      WRITE(*,*)'-----------------------'
      WRITE(*,*)MAT,BCX,BCY,BCZ
      
      open(87,file='testpot.txt')
C------------------------------------------------------------------------
      open(11,file='LJ.dat')
      READ(11,*) NKIND, REDELEC
	DO INKIND=1,NKIND
          READ(11,*)KINDI,EPSI(INKIND),SIGM(INKIND),Q(INKIND)
          write(*,*)KINDI,EPSI(INKIND),SIGM(INKIND),Q(INKIND)
          Q(INKIND)=Q(INKIND)*FACTORELEC
          
 	AE0=8.85E-12 !PERMITIVITY IN FREE SPACE



      FCLEC1AUX=SQRT(4*3.14*8.85E-12)					 !PERIMITIVIDAD POR 4 PI
      FCLEC2AUX=SQRT(SIGMETANO/SIGMA*1E-10)			 !UNIDADES DE LA CAJA
      FCLEC3AUX=SQRT((EPS*AK)*diel)						 !EPS EN UNIDADES DE PASADAS A J
      FCLECAUX=(FCLEC1AUX)*(FCLEC2AUX)*(FCLEC3AUX)
      WRITE(*,*)FCLECAUX,'FCLECAUX'
     	FCLEC=1/FCLECAUX								 !FACTOR DE REDUCCION

	WRITE(*,*)FCLEC, SIGMETANO,SIGMA,EPS,' FCLEC'
	WRITE(*,*)SIGMETANO/SIGMA*1E-10,(EPS*8.31/6.023E23)

	Q(INKIND)=REAL(Q(INKIND))*FCLEC
      WRITE(*,*)Q(INKIND), INKIND
      !PAUSE
	ENDDO
      
      CLOSE(11)
      DO IPOT=1,NKIND
          DO JPOT=1,NKIND      !ES MEDIO CHAPUCERO, YA QUE CALCULO DOS VECES LO MISMO, PERO BUEH, ES LO QUE HAY...
              SIGMA1=SIGMA*(SIGM(IPOT)+SIGM(JPOT))/(2*SIGMETANO)
              FACTOR=SQRT(EPSI(IPOT)*EPSI(JPOT))/EPS
              RCUTSQ = RCUT * RCUT
              SIGSQ  = SIGMA1 * SIGMA1
              SIGCUB = SIGSQ * SIGMA1
              RMIN=0.5*SIGMA1
              RMINSQ=RMIN*RMIN    
              !C    ** CALCULATE LONG RANGE CORRECTIONS **
              !C    ** NOTE: SPECIFIC TO LENNARD-JONES  **
              SR3    = ( SIGMA1 / RCUT ) ** 3
              SR9    = SR3 ** 3
          VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
          WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )
              !C    ** ZERO ACCUMULATORS **
              DELTV  = 0.0
              DELTW  = 0.0
              DO I=1,5000
                  RZI    = REAL(I)/1000.
                  RIJSQ = RZI*RZI
                  IF ( RIJSQ .LT. RMINSQ) THEN
                      !FLAG(I) = .TRUE.
                      VIJ=1E10*SIGM(IPOT)*SIGM(JPOT)
                      ELSE
                      SR2   = SIGSQ / RIJSQ
                      SR6   = SR2 * SR2 * SR2
                      VIJ   = SR6 * ( SR6 - 1.0 )
		      IF(RZI.GT.8*RMIN) VIJ=0
                      WIJ   = SR6 * ( SR6 - 0.5 )
                  ENDIF

                      DELTV =  4.0 * VIJ *FACTOR
          VIJE1=Q(IPOT)*Q(JPOT)*(1/RZI-1/RCELE+(1/RCELE**2)*(RZI-RCELE)) #vije1 es una variable no vij*e1
	
                      IF (RZI.GT.RCELE) VIJE1=0
		       !IF (RZI.LT.tamacel) VIJE1=VIJE1
		       IF (RZI.LT.rmin) VIJE1=1e9*abs(Q(IPOT)*Q(JPOT))
                      !WRITE(*,*)VIJE1,DELTV,I,RZI
                      
                      USS(I,IPOT,JPOT)=DELTV+VIJE1
              !WRITE(87,*)deltv,vije1,I,IPOT,JPOT
              !WRITE(*,*),deltv,vije1,I,IPOT,JPOT
                      !PAUSE
                      
		   ENDDO
		ENDDO
	     ENDDO
          WRITE(*,*)'-----------------------'
              close(87)    
	RETURN
	END


