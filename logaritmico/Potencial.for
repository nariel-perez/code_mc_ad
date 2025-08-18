C------------------------------------------------------------------------
C------------------------------------------------------------------------
C    ** SUBRUTINA POTENCIAL
C    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZ
C    ** UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(1000)
C------------------------------------------------------------------------
	SUBROUTINE POTENCIAL(EPS,sigma,sigmetano,NC,RCUT,diel)
      implicit none
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS,acel,acelx,acely,acelz,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE
!----------------------------------------------------
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
      REAL BCX,BCY,BCZ,acel,acelx,acely,acelz
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
        

C    *******************************************************************

      INTEGER NC
      REAL PI
	PARAMETER ( PI = 3.14159265 )
      REAL RCELE
      REAL FACTORELEC,AK
      INTEGER NKIND, INKIND
      INTEGER KINDI
      REAL EPSIINKIND,SIGMINKIND,QINKIND
      REAL AE0

       REAL FCLEC1AUX, FCLEC2AUX, FCLEC3AUX,FCLECAUX,FCLEC
       REAL SIGMETANO
       REAL SIGMA,EPS
       INTEGER IPOT
       INTEGER I, IJ,K
       REAL RXI,RYI,RZI
       REAL DELTV
       INTEGER J
       REAL SIGMA1, FACTOR
       REAL RCUTSQ
       REAL RCUT
       REAL SIGSQ, SIGCUB
       REAL RMIN, RMINSQ
       REAL SR3,SR9,SR2,SR6
       REAL VLRC0, WLRC0,DELTW
       REAL RXIJ,RYIJ,RZIJ
       REAL RIJSQ, RIJ
       REAL VIJ,WIJ,VIJR
       REAL ZESACT,ZI
       real diel
       real xmax,ymax,zmax
       REAL REDELEC
       
       xmax=acelx/acel
       ymax=acely/acel
       zmax=acelz/acel
       
!------------------------------------------------------------------------      
      RCELE=0.5 !RADIO DE CORTE POTENCIAL ELECTROSTÁTICO
C-----------------------------------------------------------------------
C	UNIDADES ELECTROSTÁTICAS EN TÉRMINOS DE EPS1 Y SIGMA
C-----------------------------------------------------------------------
	FACTORELEC=96500/6.023E23
      AK=8.31/6.023E23
      WRITE(*,*)'---------------------------'
      write(*,*)'POTENCIAL SOLIDO-FLUIDO'
      
C------------------------------------------------------------------------
!	OPEN(47,FILE='ENER.DAT')
      open(11,file='LJ.dat')
      READ(11,*) NKIND
      DO INKIND=1,NKIND
          READ(11,*)KINDI,EPSIINKIND,SIGMINKIND,QINKIND
          write(*,*)KINDI,EPSIINKIND,SIGMINKIND,QINKIND
          EPSI(INKIND)=EPSIINKIND
          SIGM(INKIND)=SIGMINKIND
          Q(INKIND)=QINKIND
          Q(INKIND)=Q(INKIND)*FACTORELEC
          write(*,*)q(inkind),' q ',inkind
	!pause
          AE0=8.85E-12 !PERMITIVITY IN FREE SPACE
          FCLEC1AUX=SQRT(4*3.14*8.85E-12)					 !PERIMITIVIDAD POR 4 PI
          FCLEC2AUX=SQRT(SIGMETANO/SIGMA*1E-10)			 !UNIDADES DE LA CAJA
          FCLEC3AUX=SQRT((EPS*AK)*diel)						 !EPS EN UNIDADES DE PASADAS A J
          FCLECAUX=(FCLEC1AUX)*(FCLEC2AUX)*(FCLEC3AUX)
          !WRITE(*,*)FCLECAUX,'FCLECAUX'
          FCLEC=1/FCLECAUX								 !FACTOR DE REDUCCION
          !WRITE(*,*)FCLEC, SIGMETANO,SIGMA,EPS,' FCLEC'
          !WRITE(*,*)SIGMETANO/SIGMA*1E-10,(EPS*8.31/6.023E23)
          Q(INKIND)=(Q(INKIND))*FCLEC
          write(*,*)q(inkind),' q ',inkind
      ENDDO
      write(*,*)NC, ' NUMERO DE ATOMOS EN EL ADSORBENTE'
      WRITE(*,*)MAT, 'TAMAÑO DE LA MATRIZ'

      DO IPOT=1,NKIND
          WRITE(*,*) 'POTENCIAL PARA ', IPOT,' DE ', NKIND
          DO I=-mat/2,mat/2
              RZI    = REAL(I)/real(mat)
              DO IJ=-mat/2,mat/2
                  RYI    = REAL(IJ)/real(mat)
                  DO K=-mat/2,mat/2
                      RXI    = REAL(K)/real(mat)
                      !    ** LOOP OVER ALL ATOMS  EXCEPT IPULL **
                      DELTV=0.    
                      vij=0
                      vijr=0
                      
                      
                      DO J = 1, NC
                          !*********************************************************************************************************
                          SIGMA1=SIGMA*(SGC(J)+SIGM(IPOT))/(2*SIGMETANO)
                          
                          FACTOR=SQRT(EPSI(IPOT)*EPSAC(J))/EPS
 !                     WRITE(*,*)sigma1, FACTOR, ' NFACTOR',EPSAC(J), J
                          !-----------------------------------------------------------------------------------------------
                          RCUTSQ = RCUT * RCUT
                          SIGSQ  = SIGMA1 * SIGMA1
                          SIGCUB = SIGSQ * SIGMA1
                          RMIN=0.5*SIGMA1
                          RMINSQ=RMIN*RMIN
 !                         write(*,*)rcut, 'rcut'
 !                         pause
                          !C    ** CALCULATE LONG RANGE CORRECTIONS **
                          !C    ** NOTE: SPECIFIC TO LENNARD-JONES  **
                          SR3    = ( SIGMA1 / RCUT ) ** 3
                          SR9    = SR3 ** 3
C-----------------------------------------------------------------------------------------------
C-----------------------------------------------------------------------------------------------


              RXIJ  = RXI - RXC(J)
              RYIJ  = RYI - RYC(J)
              RZIJ  = RZI - RZC(J)
              
!              write(*,*)rxij,ryij,rzij, ' RXXXXXX'
!              pause

              RXIJ  = RXIJ -BCX* xmax*ANINT ( RXIJ/xmax )
              RYIJ  = RYIJ -BCY* ymax*ANINT ( RYIJ/ymax )
	      RZIJ  = RZIJ -BCZ* zmax*ANINT ( RZIJ/zmax )
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
!              write(*,*) bcx,bcy,bcz, ' pbc'
!              write(*,*)rxij,ryij,rzij, ' RXXXXXX2'
!              pause

C	WRITE(*,*) RIJSQ

	IF ( RIJSQ .LT. RMINSQ) THEN
		VIJ=1E6
          DELTV=1e6
          exit

	ELSE
	
                 SR2   = SIGSQ / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ   = SR6 * ( SR6 - 1.0 )*FACTOR
		IF (RIJSQ.GT.(8*RMIN)**2) VIJ=0
!                 write(*,*) rxc(j),ryc(j),rzc(j),j, 'coord sust'
!                 write(*,*) rxi,ryi,rzi, vij
!                 pause
                 !write(*,*)vij, 'vij',i,ij,k
                 !pause
                 WIJ   = SR6 * ( SR6 - 0.5 )*FACTOR
                     IF(RIJSQ.LT.RCELE**2) THEN
                         RIJ=SQRT(RIJSQ)
                         ZESACT=QAC(J)
                         ZI=Q(IPOT)
              VIJR=ZI *ZESACT*(1/RIJ-1/RCELE+(1/RCELE**2)*(RIJ-RCELE))
                     else
                         vijr=0.
                 ENDIF
                 
                 
                 
	ENDIF
                 DELTV = DELTV +4* VIJ+VIJR
                ! write(*,*) deltv,vij,j,rijsq,vijr,4*vij
                 !pause
                 ! pause
                 DELTW = DELTW + WIJ


100     ENDDO
	UADS(K,IJ,I,IPOT)=DELTV
	IF (DELTV.GT.0) DELTV=0
	deltv=deltv/eps*8.31
!	WRITE(47,*)K,IJ,I,DELTV 
!      WRITE(*,*) DELTV,K,IJ,I,IPOT 
!      pause
	ENDDO
	ENDDO
	WRITE(*,*) DELTV,K,IJ,I,IPOT 
      ENDDO
      ENDDO
         
8000	write(*,*)'Energia Calculada'
!	CLOSE(47)
	RETURN
	END
C------------------------------------------------------------------------
