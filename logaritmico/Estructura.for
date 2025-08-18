C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    ** SUBRUTINA ESTRUCTURA
C	** LEE LAS COORDENADAS DE LOS ATOMOS DE CARBONO
C	** Y REALIZA UNA TRANSFORMACION PARA NORMALIZARLAS
C-----------------------------------------------------------------------
	SUBROUTINE estructura(eps,nam,sigma,sigmetano,NC,diel)
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
C    *******************************************************************
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
      REAL BCX,BCY,BCZ
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
        

C    *******************************************************************
      
      CHARACTER NAM*16
	REAL RXA(9000),RYA(9000),RZA(9000)
      REAL FCLEC1AUX, FCLEC2AUX, FCLEC3AUX,FCLECAUX,FCLEC,FACTORELEC
      REAL SIGMETANO,EPS
      REAL AK
      REAL SIGMA
      REAL RXAMAX,RYAMAX,RZYAMAX
      INTEGER NC,IMAX
      INTEGER I
      INTEGER SYMBOL(9000)
      REAL AE0
      real qnuevo
      real acel,acelx,acely,acelz
      real diel
	OPEN (49,FILE='PRUEBA.TXT')
!---------------------------------------------------------------------------
      
      AK=8.31/6.023E23  

      FCLEC1AUX=SQRT(4*3.14*8.85E-12)					 !PERIMITIVIDAD POR 4 PI
      FCLEC2AUX=SQRT(SIGMETANO/SIGMA*1E-10)			 !UNIDADES DE LA CAJA
      FCLEC3AUX=SQRT((EPS*AK))						 !EPS EN UNIDADES DE PASADAS A J
      FCLECAUX=(FCLEC1AUX)*(FCLEC2AUX)*(FCLEC3AUX)
      !WRITE(*,*)FCLECAUX,'FCLECAUX'
     	FCLEC=1/FCLECAUX								 !FACTOR DE REDUCCION
	write(*,*) nam
	!pause
!---------------------------------------------------------------------------

	OPEN (51,FILE=NAM)
	rxamax=acelx/2
	write(*,*)rxamax,' rxmax'
	ryamax=acely/2
	write(*,*)ryamax,' rymax'
      rzyamax=acelz/2
      write(*,*)rzyamax,' rzmax'
      
	READ (51,*) NC
	
	imax=0
	DO I=1,NC
      READ (51,*)RXA(I),RYA(I),RZA(I),EPSAC(I),SGC(I),QAC(I),SYMBOL(I)
      !write(*,*)RXA(I),RYA(I),RZA(I)
            
!---------------------------------------------------------------
C-----------------------------------------------------------------------
C	UNIDADES ELECTROSTÁTICAS EN TÉRMINOS DE EPS1 Y SIGMA
C-----------------------------------------------------------------------
	FACTORELEC=96500/6.023E23

C------------------------------------------------------------------------
C	CARGAS EN COULOMBS
C------------------------------------------------------------------------

		QAC(I)=QAC(I)*FACTORELEC
C----------------------------------------------------------------------------------

	AE0=8.85E-12 !PERMITIVITY IN FREE SPACE
      FCLEC1AUX=SQRT(4*3.14*8.85E-12)					 !PERIMITIVIDAD POR 4 PI
      FCLEC2AUX=SQRT(SIGMETANO/SIGMA*1E-10)			 !UNIDADES DE LA CAJA
      FCLEC3AUX=SQRT((EPS*AK)*diel)						 !EPS EN UNIDADES DE PASADAS A J
      FCLECAUX=(FCLEC1AUX)*(FCLEC2AUX)*(FCLEC3AUX)
     	FCLEC=1/FCLECAUX								 !FACTOR DE REDUCCION
	QAC(I)=REAL(QAC(I))*FCLEC

        if(rxa(i).lt.rxamax) then
	if(rxa(i).gt.-rxamax) then
        if(rya(i).lt.ryamax) then
	if(rya(i).gt.-ryamax) then
        if(rza(i).lt.rzyamax) then
	if(rza(i).gt.-rzyamax) then
                              
				        imax=imax+1
				        RXC(imax)=RXA(I)/SIGMETANO*SIGMA 
				        RYC(imax)=RYA(I)/SIGMETANO*SIGMA 
				        RZC(imax)=RZA(I)/SIGMETANO*SIGMA
                          QAC(IMAX)=QAC(I)
				        WRITE(49,*) RXC(Imax),RYC(Imax),RZC(Imax)
      !                    WRITE(*,*) rxa(i),Rya(I),Rza(I), imax,i
	endif
	endif
	endif
	endif
	endif
	endif
	

	ENDDO
	NC=imax
	write(*,*)'Estructura leida',NC
	CLOSE(51)
	OPEN (51,FILE='TRUNCADO.TXT')
	write (51,*) NC
	DO I=1,NC
	qnuevo=QAC(I)/FCLEC/FACTORELEC
      write (51,*)RXA(I),RYA(I),RZA(I),EPSAC(I),SGC(I),qnuevo,SYMBOL(I)
	ENDDO
	close(51)
	CLOSE(49)
	RETURN
	END
C------------------------------------------------------------------------
