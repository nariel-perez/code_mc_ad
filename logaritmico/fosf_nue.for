********************************************************************************
** FICHE F.13.  THE HEART OF A CONSTANT MU VT MONTE CARLO PROGRAM             **
** THIS FORTRAN CODE IS INTENDED TO ILLUSTRATE POINTS MADE IN THE TEXT.       **
** TO OUR KNOWLEDGE IT WORKS CORRECTALY.  HOWEVER IT IS THE RESPONSIBILITY OF  **
** THE USER TO TEST IT, IF IT IS TO BE USED IN A RESEARCH APPLICATION.        **
********************************************************************************

C    *******************************************************************
C    ** ATTEMPTED CREATIONS AND DESTRUCTIONS IN GRAND CANONICAL MC.   **
C    **                                                               **
C    ** THESE ROUTINES ALLOW FOR A TRIAL DESTRUCTION OR CREATION IN A **
C    ** GRAND CANONICAL MONTE CARLO PROGRAM.                          **
C    **                                                               **							    
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF ATOMS BEFORE THE TRIAL  **
C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING THE TRIAL  **
C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS ALLOWED   **									 
C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF ATOM     **
C    ** REAL    RX(NMAX) ETC.       POSITIONS OF CURRENT ATOMS        **
C    ** REAL    V                   POTENTIAL ENERGY + LRC            **
C    ** REAL    W                   VIRIAL + LRC                      **
C    ** REAL    DELTV               CHANGE IN ENERGY                  **
C    ** REAL    DELTW               CHANGE IN VIRIAL                  **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    Z                   ABSOLUTE ACTIVITY COEFFICIENT     **
C    ** REAL    SIGMA               LENNARD JONES DIAMETER            **
C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
C    ** REAL    RMIN                REDUCED MINIMUM SEPARATION        **
C    ** LOGICAL OVRLAP              TRUE FOR SUBSTANTIAL ATOM OVERLAP **
C    ** LOGICAL CREATE              TRUE FOR AN ACCEPTED CREATION     **
C    ** LOGICAL GHOST               TRUE FOR AN ACCEPTED DESTRUCTION  **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE IN ( TEMP, Z, SIGMA, RCUT, N, V, W, CREATE )       **
C    **    PERFORMS A TRIAL CREATION                                  **
C    ** SUBROUTINE OUT ( TEMP, Z, SIGMA, RCUT, N, V, W, GHOST )       **
C    **    PERFORMS A TRIAL DESTRUCTION                               **
C    ** SUBROUTINE POTIN ( RXNEW, RYNEW, RZNEW, N, SIGMA, RCUT, RMIN, **
C    ** :                  DELTV, DELTW, OVRLAP )                     **
C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON CREATION         **
C    ** SUBROUTINE POTOUT ( IPULL, N, SIGMA, RCUT, DELTV, DELTW )     **
C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON DESTRUCTION      **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** REAL FUNCTION RANF ( DUMMY )  (GIVEN IN F.11)                 **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON ZERO TO ONE            **
C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
C    **    UPDATES LOCATE AFTER ADDITION (GIVEN IN F.14)              **
C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
C    **    UPDATES LOCATE AFTER REMOVAL (GIVEN IN F.14)               **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** ROUTINES IN AND OUT SHOULD BE CALLED WITH EQUAL PROBABILITY   **
C    ** IN A GRAND CANONICAL MONTE CARLO SIMULATION. IF A TRIAL       **
C    ** CREATION IS ACCEPTED THEN CREATE IS SET TO TRUE. IF A TRIAL   **
C    ** DESTRUCTION IS ACCEPTED THEN GHOST IS SET TO TRUE. THE        **
C    ** ROUTINES ARE WRITTEN FOR LENNARD-JONES ATOMS. THE BOX IS OF   **
C    ** UNIT LENGTH, ALL DISTANCES ARE SCALED TO THE BOX LENGTH.      **
C    ** TRIAL INPUTS WHICH RESULT IN A SEPARATION OF LESS THAN        **
C    ** 0.5*SIGMA ARE REJECTED. THE LONG-RANGE CORRECTIONS ARE        **
C    ** INCLUDED IN V AND W. ALL ACCUMULATORS ARE UPDATED IN THE MAIN **
C    ** PART OF THE PROGRAM WHICH IS NOT GIVEN HERE.                  **
C    *******************************************************************
      PROGRAM MAIN
    ! ! implicit none
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

      REAL CNF1(-225:225,-225:225),CNF2(-225:225,-225:225)
	REAL CNF3(-225:225,-225:225),CNF4(-225:225,-225:225)
	REAL CNFQ(-225:225,-225:225),CNFQCI(-225:225,-225:225)
      REAL CNFQOH(-225:225,-225:225)
      REAL DCI1XY(-225:225,-225:225),DCI2XY(-225:225,-225:225)
	REAL NCICON01,NCICON02,KAPPA
	REAL P,DP,KA,KA2,DATOSQR(-9000:9000)
      DOUBLE PRECISION FCLECAUX
	
      INTEGER N,NCHAIN,JPASOS    
      INTEGER A(6)
   !   INTEGER NTRIAL              
   !   INTEGER NMAX                
   !   REAL    RXNEW,RYNEW,RZNEW   
      REAL    V                   
      REAL    W                   
      !REAL    DELTV,
	REAL DATOS(-9000:9000),DATOSCI(-9000:9000),DATOSQ(-9000:9000)
	INTEGER DATOSCI1(-9000:9000),DATOSCI2(-9000:9000)
 !     REAL    DELTW               
      REAL    TEMP                
      REAL    SIGMA               
      REAL    RCUT                
!


!     REAL    RMIN                
!      LOGICAL OVRLAP              
      LOGICAL CREATE              
      LOGICAL GHOST
	CHARACTER NAM*16
	CHARACTER NAMPOALY*16
!	CHARACTER NAMRELL*16
!	CHARACTER NMBR*2
	CHARACTER IADS*3
	CHARACTER CONFIG*3
	REAL CICON(10)
!	REAL CENTRO_MAS(100)
      CALL CPU_TIME(AINICIO)
	OPEN(10,FILE='INPUT.TXT')
	READ(10,*)P
	READ(10,*)DP
	READ(10,*)SIGMETANO
	READ(10,*)SIGMA
      READ(10,*)RCELE
      RCELE=RCELE/SIGMETANO
	WRITE(*,*) SIGMA,' TAMAÑO DE LA CELDA'
      
	AK=8.31/6.023E23
	READ(10,*)EPS
	READ(10,*)T
	TEMP=T/EPS
	WRITE(*,*) TEMP,'TEMPERATURE'
	PRED=P*SIGMA**3/EPS
	!WRITE(*,*) PRED,SIGMA,P,' '
	RCUT=10*SIGMA
	XMAX=1.
	YMAX=1.
	ZMAX=1.
	READ(10,*)MAT
	WRITE(*,*) MAT, ' MATRIX SIZE'
	VOL=XMAX*YMAX*ZMAX
	READ(10,*)NAM
	WRITE(*,*)NAM,' STRUCTURE TO BE READ'
	READ(10,*)NAMPOALY	
	WRITE(*,*)NAMPOALY, ' AUXILIAR STRUCTURE TO BE READ'
	READ(10,*)ALREF
	WRITE(*,*)ALREF, ' REFERENCE SIZE BETWEEN THE PARTICLES IN  CHAIN'
	READ(10,*)ITETH
	WRITE(*,*)ITETH, '  THETERED?'
	READ(10,*)DAT
	WRITE(*,*)DAT, ' ENERGY SAVED'
	READ(10,*)ISOT
	WRITE(*,*)ISOT,' ISOTHERMS POINT'
	READ(10,*)IJPASOS
	WRITE(*,*)IJPASOS, 'BLOCK AVERAGING POINTS'
	READ(10,*)IKPASOS
	WRITE(*,*)IKPASOS, ' INNER MONTECARLO STEPS'
	READ(10,*)MULT2
	WRITE(*,*)MULT2,' MULT2 FOR FIRST POINT'
	READ(10,*) PBX
	WRITE(*,*)PBX,'PBX'
	READ(10,*) PBY
	WRITE(*,*)PBY,'PBY'
	READ(10,*) PBZ
	WRITE(*,*)PBZ,'PBZ'
	READ(10,*) APOT
	WRITE(*,*) APOT,' APOT'
	READ(10,*) APOT2
	WRITE(*,*) APOT2, ' APOT2'
	READ(10,*) KA
	WRITE(*,*) KA,' KA'

	READ(10,*) KA2
	WRITE(*,*) KA2,' KA2'

	READ(10,*) QLEC
	WRITE(*,*)QLEC, ' QLEC'
	READ(10,*) DIELEC
	WRITE(*,*) DIELEC,' DIELEC'
	READ(10,*)ALX
	READ(10,*)ALY
	READ(10,*)ALZ
	WRITE(*,*) ALX,ALY,ALZ, ' ALX,ALY,ALZ'
	READ(10,*)IKAPPA
	WRITE(*,*) IKAPPA, ' IKAPPA'
      
	READ(10,*)CICON(1)
	READ(10,*)CICON(2)
	READ(10,*)NZI(1)
	READ(10,*)NZI(2)
	READ(10,*)NCHARGEMONOMER1
	WRITE(*,*)'NCHARGEMONOMER1',NCHARGEMONOMER1
	READ(10,*)NCHARGEMONOMER2
	WRITE(*,*) CICON(1),CICON(2),NZI(1),NZI(2),NH, ' COUNTERIONS'
	READ(10,*)IEXPLICIT
	WRITE(*,*)IEXPLICIT,' EXPLICIT_COUNTERION'
	READ(10,*)POLCIL
	READ(10,*)ZCIL
	READ(10,*)RCIL
	WRITE(*,*) POLCIL,ZCIL,RCIL, ' CILINDRO'
	READ(10,*) TAMCELDISC
      READ(10,*) TAMDEF
      NCELLMAT=INT(SIGMA/TAMDEF)
      IF(POLCIL.GT.0) NRCELLMAT=2*INT(RCIL/TAMDEF)
      IF(POLCIL.GT.0) NCELLMAT=2*INT(RCIL/TAMDEF)
      NCELLMAT=100
      NRCELLMAT=100
      WRITE(*,*)NCELLMAT,' NCELLMAT'
      
      CLOSE(10)
      
	CONOH=P
	CONCIA=P
!-----------------------------------------------------------------------------------------------------
!      
!      CONCENTRACIONES DE LOS IONES
!----------------------------------------------------------------------------------------------------
      OPEN(10,FILE='IONES.TXT')
      READ (10,*) IONES
      DO I=1,IONES
      READ (10,*)NQCI(I)
      READ (10,*)CICON(I)
      
      ENDDO
      ciconaux3=cicon(3)
      ciconaux4=cicon(4)
      
      CLOSE(10)
!-----------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
	
	NDISC=SIGMA/TAMCELDISC
      
	PHBULK=-LOG10(P)
      PHBULK2=-LOG10(CICON(1))
      
	WRITE(*,*) PHBULK,PHBULK2
	!PAUSE
      DO J=1,4
	DO I=1,9000
	LOCATE(I,J)=0
      ENDDO
      ENDDO
	VOLPH=1/TAMCELDISC**3/6.023E23*1E24*1000
	VOLNORMAL=1/6.023E23/(SIGMA**3*ALX*ALY*ALZ)*1E27
	ZCIL=ZCIL/SIGMA
	RCIL=RCIL/SIGMA
	SIGMA=SIGMETANO/SIGMA
      WRITE(*,*)VOLPH,VOLNORMAL,' VOL'
	WRITE(*,*)SIGMA
      FCLEC1AUX=SQRT(8.85E-12*12.56637061)
      FCLEC2AUX=SQRT(SIGMETANO/SIGMA*1E-10)
      FCLEC3AUX=SQRT((EPS*8.31/6.023E23)*DIELEC)
      FCLECAUX=(FCLEC1AUX)*(FCLEC2AUX)*(FCLEC3AUX)
      !WRITE(*,*)FCLEC1AUX, FCLEC2AUX,      FCLEC3AUX
      !!PAUSE
      !FCLECAUX=SQRT(FCLEC1AUX*FCLEC2AUX*FCLEC3AUX)
      
      WRITE(*,*)FCLECAUX,'FCLECAUX'
      !!PAUSE
      	FCLEC=1.6E-19/FCLECAUX

	WRITE(*,*)FCLEC, DIELEC,SIGMETANO,SIGMA,EPS,' FCLEC'
	!!PAUSE
	WRITE(*,*)SIGMETANO/SIGMA*1E-10,(EPS*8.31/6.023E23),DIELEC
     !!!!!!!PAUSE

	!FCLEC=FCLEC*SIGMA/SIGMETANO
	AKAPPA=0
	DO I=1,IONES
	
	AKAPPA=AKAPPA+REAL(NQCI(I)**2)*CICON(I)/(8.85E-12
     +*DIELEC*(1.3807E-23)*T)
	WRITE(*,*)AKAPPA,AKAPPA
      ENDDO
      
	IF(IEXPLICIT.EQ.1) IKAPPA=0

	AKAPPA=REAL(IKAPPA)*AKAPPA*1000*(1.6E-38)*6.023E23 *
     +(1E-10)**2
	AKAPPA12=SQRT(AKAPPA) *SIGMETANO/SIGMA
	!WRITE(*,*) FCLEC,AKAPPA,AKAPPA12,1/AKAPPA12,EXP(-AKAPPA12)
		WRITE(*,*) MAT, ' MAT ENTRADA'

	!!!!!!!PAUSE

      

	CALL STRUCPOALY(NAMPOALY,SIGMA,SIGMETANO,N,NCHAIN,EPS)
	WRITE(*,*) 'INICIANDO1'
		WRITE(*,*) MAT, ' MAT ENTRADA'

		OPEN(50,FILE='SALIDAACTIVADO-100.TXT')
	OPEN(97,FILE='PERFILES.TXT')
	CALL ESTRUCTURA(NAM,SIGMA,SIGMETANO)
		WRITE(*,*) MAT, ' MAT ENTRADA'

	WRITE(*,*) 'INICIANDO2'
	CALL POTENCIAL2(SIGMA,SIGMETANO,RCUT,APOT2,ALREF)
		WRITE(*,*) 'INICIANDO3'
	WRITE(*,*) MAT, ' MAT ENTRADA'
	CALL POTENCIAL(EPS,SIGMA,SIGMETANO,RCUT,DAT,APOT)
	DO J=1,4
      DO I=1,9000
	LOCATE(I,J)=0
      ENDDO
      ENDDO
	WRITE(*,*) 'INICIANDO'
	!!!!!!!!PAUSE
C-----------------------------------------------------------------
	IPH=0
	IPOH=0
	DO I=-NDISC/2,NDISC/2
	DO J=-NDISC/2,NDISC/2
	DO K=-NDISC/2,NDISC/2
	PHPARCIAL(I,J,K)=0
	ENDDO
	ENDDO
	ENDDO

	DO I=-225,225
	DO J=-225,225
	CNF1(I,J)=0
	CNF2(I,J)=0
	CNF3(I,J)=0
      DCI1XY(I,J)=0
      DCI2XY(I,J)=0
	CNF4(I,J)=0
	CNFQ(I,J)=0
	CNFQCI(I,J)=0
	DATOSQ(I)=0
	DATOS(I)=0
      DO ION=1,4
	DCI(I,ION)=0
      ENDDO

	DATOSCI1(I)=0
	DATOSQR(I)=0

	ENDDO
      ENDDO
      DO I=1,9000
     	DATOSCI2(I)=0
	DATOSCIR1(I)=0
	DATOSCIR2(I)=0
	DATOSR(I)=0

      ENDDO

	WRITE(*,*) 'INICIANDO5'

	DO I=-NDISC/2,NDISC/2
	DO J=-NDISC/2,NDISC/2
	DO K=-NDISC/2,NDISC/2
	PHPARCIAL(I,J,K)=0
	PH(I,J,K)=PHBULK
	ENDDO
	ENDDO
	ENDDO

	DO IPASOS=1,ISOT
	PHBULK=-LOG10(P)
	WRITE(*,*) PHBULK,' PHBULK'
!	!!!PAUSE
	IF(QLEC.GT.0) THEN
	DO ION=1,4
          DO I=1,9000
          
	LOCATE(I,ION)=0
      ENDDO
	ENDDO

	NK(1)=0
	NK(2)=0
      DO ION=1,4
	NCI(ION)=0
	ENDDO
	
	CALL STRUCPOALY(NAMPOALY,SIGMA,SIGMETANO,N,NCHAIN,EPS)
	WRITE(*,*) 'INICIANDO6'

	IF(ABS(DP).GT.0) THEN



	DO I=-225,225
	DO J=-225,225
	CNF1(I,J)=0
	CNF2(I,J)=0
	CNF3(I,J)=0
      DCI1XY(I,J)=0
      DCI2XY(I,J)=0
	CNF4(I,J)=0
	CNFQ(I,J)=0
	CNFQCI(I,J)=0
	DATOSQ(I)=0
	DATOS(I)=0
      DO ION=1,4
	DCI(I,ION)=0
      ENDDO
	DATOSQR(I)=0
	ENDDO
      ENDDO
      DO I=1,9000
      DO ION=1,4
          	DATOSN(I,ION)=0
      ENDDO
      
	DATOSCIR1(I)=0
	DATOSCIR2(I)=0
	DATOSR(I)=0

      ENDDO
	ENDIF

	ELSE

	DO I=-225,225
	DO J=-225,225
	CNF1(I,J)=0
	CNF2(I,J)=0
	CNF3(I,J)=0
      DCI1XY(I,J)=0
      DCI2XY(I,J)=0
	CNF4(I,J)=0
	CNFQ(I,J)=0
	CNFQCI(I,J)=0
	DATOSQ(I)=0
	DATOS(I)=0
      DO ION=1,4
	DCI(I,ION)=0
      ENDDO
      
	ENDDO
	ENDDO

	ENDIF
	WRITE(*,*) 'INICIANDO6'


	CONFIG='CONFIG'
	WRITE(*,*) CONFIG,IPASOS
	WRITE(CONFIG,'(I3)') IPASOS
	WRITE(*,*) CONFIG,IPASOS
!	!!!!!!PAUSE
	IADS='IADS'
	IPASOS1=1
	WRITE(IADS,'(I3)')IPASOS
	WRITE(*,*)IADS
	OPEN(39,FILE='CONFIGQ'//IADS//'.TXT')
	OPEN(40,FILE='CONFIG'//IADS//'.XYZ')
	OPEN(41,FILE='CONFIGZ' //IADS//'.TXT')
	OPEN(42,FILE='CONFIGCIZ' //IADS//'.TXT')
	OPEN(43,FILE='CONFIGXY' //IADS//'.TXT')
	OPEN(44,FILE='CONFIGXZ' //IADS//'.TXT')
      
      OPEN(67,FILE='CONFIGXYCI1'//IADS//'.TXT')
      OPEN(68,FILE='CONFIGXYCI2'//IADS//'.TXT')
      OPEN(69,FILE='CONFIGXYCI3'//IADS//'.TXT')
      OPEN(70,FILE='CONFIGXYCI4'//IADS//'.TXT')

      OPEN(167,FILE='CONFIGRCI1'//IADS//'.TXT')
      OPEN(168,FILE='CONFIGRCI2'//IADS//'.TXT')
      OPEN(169,FILE='CONFIGRCI3'//IADS//'.TXT')
      OPEN(170,FILE='CONFIGRCI4'//IADS//'.TXT')
      
      
	OPEN(45,FILE='CONFIGNXY' //IADS//'.TXT')
	OPEN(46,FILE='CONFIGNXZ' //IADS//'.TXT')
	OPEN(47,FILE='CONFIGQXY' //IADS//'.TXT')
	OPEN(48,FILE='CONFIGQCIXY' //IADS//'.TXT')

	OPEN(72,FILE='CONFIGCIZ1' //IADS//'.TXT')
	OPEN(73,FILE='CONFIGCIZ2' //IADS//'.TXT')

	OPEN (49,FILE='PRUEBA.TXT')


* CONVERSION OF PRESSURE INTO ACTIVITY ZP: 1 CM HG=1333.22 PA; NA=6.0220E+23;
* R=8.3144 J/(MOL*K); Z IN ATOM/A**3
C	WRITE(*,*) AP ,DP
!	WRITE(*,*)P,KA
!	!!!!!!!PAUSE
       VOLREAL=ALX*ALY*ALZ
      IF(POLCIL.NE.0) VOLREAL=3.14*ZCIL*RCIL**2
      AL1=ALX
      AL2=ALY
      AL3=ALZ
      
C--------------------------------------------------------------------
C	IF IEXPLICIT=1 -->ACTIVITIES FOR IONS
      
      CICON(3)=P+1E-14/P+ciconaux3
      CICON(4)=1E-14/P+P+ciconaux4
      WRITE(*,*) CICON(3),CICON(4)
      !PAUSE
	IF(IEXPLICIT.EQ.1) THEN
      DO I=1,IONES      
      Z(I)=(CICON(I))*(SIGMETANO/SIGMA)**3*VOLREAL*6.023E23*1E-27
      ENDDO


      WRITE(*,*)Z(3),Z(4),Z(1),Z(2),' Z3, Z4, Z1, Z2'
      !PAUSE
	ENDIF
C--------------------------------------------------------------------
      P1=(KA/P)
	
!-----------------------------------------------------------------------      
	IF(KA.NE.0) THEN
	WRITE(*,*) KA
	WRITE(*,*) P
	FR=1/(1+P/KA)
	
	PRED=P*SIGMA**3/EPS
	ESC=SIGMETANO/SIGMA
	
	Z(5)=-LOG10(P)+LOG10(KA)
	WRITE(*,*)Z(1),' Z1'
		
      ENDIF
      IF(KA2.NE.0) THEN
 	WRITE(*,*) KA2
	WRITE(*,*) P
	ESC=SIGMETANO/SIGMA
	Z(6)=-LOG10(P)+LOG10(KA2)
	WRITE(*,*)Z(6),' Z6'
	
      ENDIF
!------------------------------------------------------------------------

		ESC=SIGMETANO/SIGMA
		WRITE(*,*) 'INICIANDO7'
      WRITE(*,*)IPH,IPOH, 'iph, ipoh'
C---------------------------------
	QT=0.
	QT1=0.
	QT2=0.
	QTOH=0.
	NCICON01=0.
	NCICON02=0.

	U=0
	UN=0
	UG=0
 	UNG=0
	UA=0
 	UNA=0
	AN=0
	N2=0
	AN1=0
	NTOTAL=0
 	U=0
	UG=0		   
	UA=0
	UN=0
	UNG=0
	UNA=0
	AN=0
	N2=0

C	V=0
C	N=0
* SECOND VIRIAL CORRECTION SHOULD BE ADDED
!	WRITE(*,*) NCI
!	!!!!!!!PAUSE
	JNKP=0
	DO JPASOS=1,IJPASOS
!	JNKP=JNKP+1
	!WRITE(*,*)JPASOS,JNKP
!	OPEN(51,FILE='ENER'//CONFIG//'.TXT')
	MULT=1
	IF (JPASOS.EQ.1) THEN
	MULT=MULT2
	ENDIF
	IF (IPASOS.EQ.1.AND.JPASOS.EQ.1) THEN
	 MULT=MULT2
	ENDIF
!		WRITE(*,*) NCI
		!!!!!!!PAUSE		  

!	WRITE(*,*)NCICON(1)
!	WRITE(*,*)NCICON(2)
!	WRITE(*,*)'----------------------ARRANCA-----------------'
!	!!!PAUSE
!	OPEN(93,FILE='LOG.TXT')
!	!!!!PAUSE
	!WRITE(*,*)'START ',JPASOS
	DO KPASOS=1,IKPASOS*MULT    !*N*NCHAIN

C------------------------------------------------------------------
C	ELECCION DEL PASO
C------------------------------------------------------------------
1887	IJ=INT(RANF(DUMMY)*((IEXPLICIT*4)))+1
      
	!IF(KA.EQ.0) THEN
	!IJ=INT(RANF(DUMMY)*(1+(IEXPLICIT*)))+1
	!ENDIF
	GOTO (10,20,36,37) IJ
	!IF(KA.EQ.0)GOTO (30,35,34,36,37) IJ

  10	 CALL IN(TEMP,SIGMA,EPS,RCUT,N,NCHAIN,
     +CREATE,NCHARGEMONOMER1,NCHARGEMONOMER2,PHBULK,PHBULK2)
!	IJ=INT(RANF(DUMMY)*4)+1
!	GOTO (34,36,46,39) IJ
	GOTO 30

  20	CALL OUT(TEMP,  SIGMA,EPS, RCUT, N,NCHAIN,
     +NCHARGEMONOMER1,NCHARGEMONOMER2,phbulk,PHBULK2)
!	IJ=INT(RANF(DUMMY)*6)+1
!	GOTO (35,37,47,38,10,40) IJ


	GOTO 30
  30  cALL MOVE(SIGMETANO,TEMP,SIGMA,EPS,RCUT,N,NCHAIN,V,VA,VG,
     +W,GHOST,ITETH,ALREF)

	GOTO 40

  36	CALL INCION(SIGMETANO,TEMP,SIGMA,EPS,RCUT,N,NCHAIN,V,VA,VG,W,
     +CREATE,CR)
      GOTO 30


  37	CALL OUTION(SIGMETANO,TEMP,  SIGMA,EPS, 
     +RCUT, N,NCHAIN, 
     +V,VA,VG, W, GHOST )
!	IJ=INT(RANF(DUMMY)*6)+1
!	GOTO (34,36,46,39,20,40) IJ
      goto 30
      
  40	ENDDO !WRITE(*,*)NCI,NCIOH,NCICON(1),NCICON(2),NK(1),NK(2)
  	


C---------------------------------------------------------------------------
C----------------------------------------------------------------------
	QT=QT+REAL(NCI(1))
	QTOH=QTOH+REAL(NCI(2))

	QT1=QT1+REAL(NK(1))
	QT2=QT2+REAL(NK(2))
	NCICON01=NCICON01+REAL(NCI(3))
	NCICON02=NCICON02+REAL(NCI(4))
	NTOTAL=1

	V=1

	VG=1
!	WRITE(*,*)1

	VA=1
	!WRITE(*,*)QTPARCIAL, QTOHPARCIAL,VOLPH,IPH,IPOH,VOLNORMAL 
      GOTO 678
C----------------------------------------------------------------------
678	RXESC=RX(N,NCHAIN)*SIGMETANO/SIGMA
	RYESC=RY(N,NCHAIN)*SIGMETANO/SIGMA
	RZESC=RZ(N,NCHAIN)*SIGMETANO/SIGMA

!	WRITE(*,*)JPASOS,RZESC,NCI,NCICON(1),NCICON(2),PH(1,1,1)
C----------------------------------------------------------------------
CÁLCULO DE LA DISTANCIA EXTREMO A EXTREMO
C----------------------------------------------------------------------
	ETER=0.
	DO I=1,NCHAIN
	ETEX=(RX(1,I)-RX(N,I))**2
	ETEY=(RY(1,I)-RY(N,I))**2
	ETEZ=(RZ(1,I)-RZ(N,I))**2
	ETEN=SQRT(ETEX+ETEY+ETEZ)
      ETER=ETER+ETEN
	ENDDO
	ETE=ETE+ETER/REAL(NCHAIN)
	!WRITE(*,*)'P1'
	!!!PAUSE
C-------------------------------------------------------------------------------
C PERFILES DE DENSIDAD EN FUNCIÓN DE Z
C----------------------------------------------------------------------

	DO ICHAIN=1,NCHAIN
	DO I=2,N
	DO IK=0,200
C	DX=SQRT(RX(I)**2+RY(I)**2)
	DX=RZ(I,NCHAIN)+0.5
	AIK=(REAL(IK)*0.005)
	IF(DX.LE.AIK) THEN
	!DATOS(IK)=DATOS(IK)+1 

	ENDIF
	ENDDO
	ENDDO
	ENDDO

	!WRITE(*,*)'P2'
	!!!PAUSE

	DO ICHAIN=1,NCHAIN
 	DO I=2,N


	IJCN=INT(RZ(I,ICHAIN)*NCELLMAT)
	DATOS(IJCN)=DATOS(IJCN)+1
	ENDDO
	ENDDO

	!WRITE(*,*)'P3'
	!!!PAUSE
      IF(POLCIL.NE.1) GOTO 365
C----------------------------------------------------------------------
C	PERFILES DE DENSIDAD RADIALES
C----------------------------------------------------------------------
	DO ICHAIN=1,NCHAIN
 	DO I=2,N


	IJCN=INT(SQRT((RX(I,ICHAIN))**2+(RY(I,ICHAIN))**2)*100)
      !WRITE(*,*)IJCN
	DATOSR(IJCN)=DATOSR(IJCN)+1
	ENDDO
	ENDDO

	!WRITE(*,*)'P4'
	!!!PAUSE

	DO ICHAIN=1,NCHAIN
 	DO I=2,N


	IJCN=INT(SQRT((RX(I,ICHAIN))**2+(RY(I,ICHAIN))**2)*100)
      

	IF(Q(I,ICHAIN).NE.0) THEN
	DATOSQR(IJCN)=DATOSQR(IJCN)+1
	ENDIF

	ENDDO
	ENDDO

	!WRITE(*,*)'P5'
	!!!PAUSE
C----------------------------------------------------------------------
C	COUNTERIONS
C----------------------------------------------------------------------

!	DO I=1,NCI
!
!	IPULL3=LOCATE(I)
!	IJCN=INT(RZCI(IPULL3)*100)
!	DATOSCI(IJCN)=DATOSCI(IJCN)+1
!	ENDDO

	DO I=1,NCI(3)

	IPULL3=LOCATE(I,3)
	IJCN=INT(RZCI(IPULL3,3)*NCELLMAT)
	DATOSCI1(IJCN)=DATOSCI1(IJCN)+1


	ENDDO

	!WRITE(*,*)'P8'
	!!!PAUSE

	DO I=1,NCI(4)

	IPULL3=LOCATE(I,4)
	IJCN=INT(RZCI(IPULL3,4)*NCELLMAT)
	DATOSCI2(IJCN)=DATOSCI2(IJCN)+1
	ENDDO

	!WRITE(*,*)'P9'
	!!!PAUSE

	DO ICHAIN=1,NCHAIN
 	DO I=2,N


	IJCN=INT(RZ(I,ICHAIN)*NCELLMAT)
	IF(Q(I,ICHAIN).NE.0) THEN
	DATOSQ(IJCN)=DATOSQ(IJCN)+1
	ENDIF

	ENDDO
	ENDDO

	!WRITE(*,*)'P10'
	!!!PAUSE

C----------------------------------------------------------------------
  !-----------------------------------------------------------
  !PERFIL DE DENSIDADES EN EL PLANO XY
  !-----------------------------------------------------------
365	DO ICHAIN=1,NCHAIN
	DO I=1,N
!	WRITE(*,*) I,ICHAIN
	ICNF=INT(RX(I,ICHAIN)/ALX*NCELLMAT/2)
      JCNF=INT(RY(I,ICHAIN)/ALY*NCELLMAT/2)
!	WRITE(*,*)ICNF,JCNF
	CNF1(ICNF,JCNF)=CNF1(ICNF,JCNF)+1
	ENDDO
	ENDDO

	!WRITE(*,*)'P11'
	!!!PAUSE
	DO ICHAIN=1,NCHAIN
	DO I=1,N
	ICNF=INT(RX(I,ICHAIN)/ALX*NCELLMAT/2)
      JCNF=INT(RY(I,ICHAIN)/ALY*NCELLMAT/2)
	CNFQ(ICNF,JCNF)=CNFQ(ICNF,JCNF)+ABS(Q(I,ICHAIN)/
     +(Q(I,ICHAIN)+1E-23))

	ENDDO
	ENDDO
	!WRITE(*,*)'P12'
	!!!PAUSE

	DO ICHAIN=1,NCHAIN
	ICNF=INT(RX(N,ICHAIN)/ALX*NCELLMAT/2)
      JCNF=INT(RY(N,INT(ICHAIN/ALY))*NCELLMAT/2)
	CNF3(ICNF,JCNF)=CNF3(ICNF,JCNF)+1
	ENDDO

	!WRITE(*,*)'P13'
	!!!PAUSE

	DO I=1,NCI(1)
	IPULL3=LOCATE(I,1)
	ICNF=INT(RXCI(IPULL3,1)/ALX*NCELLMAT/2)
      JCNF=INT(RYCI(IPULL3,1)/ALY*NCELLMAT/2)
	CNFQCI(ICNF,JCNF)=CNFQCI(ICNF,JCNF)+1
      ENDDO

      DO I=1,NCI(2)
	IPULL3=LOCATE(I,2)
	ICNF=INT(RXCI(IPULL3,2)/ALX*NCELLMAT/2)
      JCNF=INT(RYCI(IPULL3,2)/ALY*NCELLMAT/2)
	CNFQOH(ICNF,JCNF)=CNFQOH(ICNF,JCNF)+1
	ENDDO

      DO I=1,NCI(3)
	IPULL3=LOCATE(I,3)
	ICNF=INT(RXCI(IPULL3,3)/ALX*NCELLMAT/2)
      JCNF=INT(RYCI(IPULL3,3)/ALY*NCELLMAT/2)
	DCI1XY(ICNF,JCNF)=DCI1XY(ICNF,JCNF)+1
	ENDDO
	DO I=1,NCI(4)
	IPULL3=LOCATE(I,4)
	ICNF=INT(RXCI(IPULL3,4)/ALX*NCELLMAT/2)
      JCNF=INT(RYCI(IPULL3,4)/ALY*NCELLMAT/2)
	DCI2XY(ICNF,JCNF)=DCI2XY(ICNF,JCNF)+1
	ENDDO

  !-----------------------------------------------------------
  !PERFIL DE DENSIDADES EN EL PLANO XZ
  !-----------------------------------------------------------
	DO ICHAIN=1,NCHAIN
	DO I=1,N
	ICNF=INT(RX(I,ICHAIN)/ALX*NCELLMAT)
      JCNF=INT(RZ(I,ICHAIN)/ALZ*NCELLMAT)
      


	CNF2(ICNF,JCNF)=CNF2(ICNF,JCNF)+1

	ENDDO
	ENDDO
	!WRITE(*,*)'P15'
	!!!PAUSE

	DO ICHAIN=1,NCHAIN
	ICNF=INT(RX(N,ICHAIN)/ALX*NCELLMAT)
      JCNF=INT(RZ(N,ICHAIN)/ALZ*NCELLMAT)
	CNF4(ICNF,JCNF)=CNF4(ICNF,JCNF)+1
	ENDDO
C----------------------------------------------------------------
	!WRITE(*,*)'P16'
	!!!PAUSE

	
	ENDDO
C-----------------------------------------------------------------------
C	ESCRIBE ARCHIVO DE SALIDA PARA VISUALIZAR
C----------------------------------------------------------------------
	!WRITE(40,*)NCI
	DO I=1,NC
          write(*,*) nci
	READ(49,*) RXC(I),RYC(I),RZC(I)
      !write(*,*) RXC(I),RYC(I),RZC(I)

	WRITE(40,*) 16,RXC(I)*ESC,RYC(I)*ESC,RZC(I)*ESC

      ENDDO
      !pause
	A(1)=1
      A(2)=8
      A(3)=3
      A(4)=7
      DO ION=1,4
	DO I=1,NCI(ION)
	J=LOCATE(I,ION)
	WRITE(40,*)A(ION) ,RXCI(J,ION)*ESC,RYCI(J,ION)*ESC,RZCI(J,ION)*ESC
      ENDDO
      ENDDO

	


	DO JCHAIN=1,NCHAIN
	DO J=1,N
	J2=J+NC
	PERF=(DATOS(J)-DATOS(J-1))/REAL(NTOTAL)
	IF(Q(J,JCHAIN).LT.0) IAT=37
	IF(Q(J,JCHAIN).GT.0) IAT=55
	IF(Q(J,JCHAIN).EQ.0) IAT=56

	!IAT=6+1*(Q(J,JCHAIN))/(ABS(Q(J,JCHAIN))+0.0000001)

	WRITE(40,*) IAT,RX(J,JCHAIN)*ESC,RY(J,JCHAIN)*ESC,
     +(RZ(J,JCHAIN)*ESC)
	ENDDO
	ENDDO

!	WRITE(*,*)NCI,NCICON(1),NCICON(2)
!	!!!!!!!PAUSE
!	CLOSE(40)
C------------------------------------------------------------------------
C	ESCRIBE EL PERFIL EN Z
C------------------------------------------------------------------------
	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		PERF=DATOS(JCON2)
		IF(POLCIL.GT.0) THEN
              if(jcon2.gt.0) then
              PERF=DATOSR(JCON2)
              WRITE(41,*) JCON2,PERF
              ENDIF
              endif
	IF(POLCIL.LT.1) THEN
          PERF=DATOS(JCON2)
          WRITE(41,*) JCON2+NCELLMAT/2,PERF
              ENDIF
		

	ENDDO

	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		
	IF(POLCIL.GT.0) THEN
          IF(JCON2.GT.0) THEN
          PERF=DATOSCIR1(JCON2)
          WRITE(72,*) JCON2,PERF
          ENDIF
          ENDIF
	IF(POLCIL.LT.1) THEN
          PERF=DCI(JCON2+NCELLMAT/2,1)
          WRITE(72,*) JCON2+NCELLMAT/2,PERF
          ENDIF
		
	ENDDO

	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		
	IF(POLCIL.GT.0) THEN
          IF(JCON2.GT.0) THEN
          PERF=DATOSCIR2(JCON2)
          WRITE(73,*) JCON2,PERF
          ENDIF
          ENDIF
	IF(POLCIL.LT.1) THEN
          PERF=DCI2(JCON2+NCELLMAT/2)
		WRITE(73,*) JCON2+NCELLMAT/2,PERF
          ENDIF

	ENDDO


	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		PERF=DATOSQ(JCON2)
	IF(POLCIL.GT.0) THEN
          PERF=DATOSQR(JCON2)
          WRITE(39,*) JCON2,PERF
          ENDIF
 	IF(POLCIL.LT.1) THEN
      PERF=DATOSQ(JCON2+NCELLMAT/2)

		WRITE(39,*) JCON2+NCELLMAT/2,PERF
          ENDIF

	ENDDO

	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		
		IF(POLCIL.GT.0) THEN
              IF(JCON2.GT.0) THEN 
              PERF=DATOSN(JCON2,1)
              WRITE(42,*) JCON2,PERF
              ENDIF
          ENDIF
          
		IF(POLCIL.LT.1) THEN
              PERF=DCI(JCON2+NCELLMAT/2,1)
		WRITE(42,*) JCON2+NCELLMAT/2,PERF
          ENDIF

	ENDDO


C-----------------------------------------------
!--------------------------------------------------------------------------------------
!ESCRIBE EL PERFIL DE DENSIDADES
!--------------------------------------------------------------------------------------
	DO JCON1=-NCELLMAT/4,NCELLMAT/4
	DO JCON2=-NCELLMAT/4,NCELLMAT/4

		PERF=CNF1(JCON1,JCON2)
		WRITE(43,*) JCON1,JCON2,PERF

	ENDDO

	ENDDO

	DO JCON1=-NCELLMAT/4,NCELLMAT/4
	DO JCON2=-NCELLMAT/4,NCELLMAT/4

		PERF=CNFQ(JCON1,JCON2)
		WRITE(47,*) JCON1,JCON2,PERF

	ENDDO
	ENDDO

	DO JCON1=-NCELLMAT/4,NCELLMAT/4
	DO JCON2=-NCELLMAT/4,NCELLMAT/4

		PERF=CNFQCI(JCON1,JCON2)
		WRITE(48,*) JCON1,JCON2,PERF

	ENDDO
	ENDDO

	DO JCON1=-NCELLMAT/4,NCELLMAT/4
	DO JCON2=-NCELLMAT/4,NCELLMAT/4

		PERF=CNF2(JCON1,JCON2)

	WRITE(44,*) JCON1,JCON2,PERF
	ENDDO

	ENDDO

	DO JCON1=-NCELLMAT/4,NCELLMAT/4
	DO JCON2=-NCELLMAT/4,NCELLMAT/4

		PERF=CNF3(JCON1,JCON2)
		WRITE(45,*) JCON1,JCON2,PERF

	ENDDO

      ENDDO
!-----------------------------------------------------------------------------------
!     PERFILES X-Y DE LOS IONES
!-----------------------------------------------------------------------------------
      
     	DO JCON1=-NCELLMAT/2,NCELLMAT/2
	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		PERF=DXY(JCON1,JCON2,1)
		WRITE(67,*) JCON1,JCON2,PERF

	ENDDO
      ENDDO

      DO JCON1=-NCELLMAT/2,NCELLMAT/2
	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		PERF=DXY(JCON1,JCON2,2)
		WRITE(68,*) JCON1,JCON2,PERF

	ENDDO
      ENDDO
     	DO JCON1=-NCELLMAT/2,NCELLMAT/2
	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		PERF=DXY(JCON1,JCON2,3)
		WRITE(69,*) JCON1,JCON2,PERF

	ENDDO
      ENDDO

      
     	DO JCON1=-NCELLMAT/2,NCELLMAT/2
	DO JCON2=-NCELLMAT/2,NCELLMAT/2

		PERF=DXY(JCON1,JCON2,4)
		WRITE(70,*) JCON1,JCON2,PERF

	ENDDO
      ENDDO
!-----------------------------------------------------------------------------------
!     PERFILES RADIALES DE LOS IONES
!-----------------------------------------------------------------------------------
      
	DO JCON2=1,101

		PERF=DATOSN(JCON2,1)
		WRITE(167,*) JCON2,PERF

	ENDDO

	DO JCON2=1,101

		PERF=DATOSN(JCON2,2)
		WRITE(168,*) JCON2,PERF

      ENDDO

	DO JCON2=1,101

		PERF=DATOSN(JCON2,3)
		WRITE(169,*) JCON2,PERF

	ENDDO

      
	DO JCON2=1,101

		PERF=DATOSN(JCON2,4)
		WRITE(170,*)JCON2,PERF

	ENDDO


!-------------------------------------------------------------------------------------

	DO JCON1=-NCELLMAT/4,NCELLMAT/4
	DO JCON2=-NCELLMAT/4,NCELLMAT/4

		PERF=CNF4(JCON1,JCON2)
		WRITE(46,*) JCON1,JCON2,PERF

	ENDDO

	ENDDO

C---------------------------------------------------------------------------------------------
	CLOSE(39)
      CLOSE(40)
	CLOSE(41)
	CLOSE(42)
	CLOSE(43)
	CLOSE(44)
	CLOSE(45)
	CLOSE(46)
	CLOSE(47)
	CLOSE(48)
	CLOSE(49)
      CLOSE(67)
      CLOSE(68)
      CLOSE(69)
      CLOSE(70)
      CLOSE(167)
      CLOSE(168)
      CLOSE(169)
      CLOSE(170)
      CLOSE(72)
	CLOSE(73)

!	!!!!!!!PAUSE

C-----
!	CLOSE(51)
	WRITE(*,*)'-----'
 27	FORMAT('CONECT',1X,I4,1X,I4)

 25	FORMAT('HETATM',1X,I4,2X,'C',7X,I4,6X,F7.3,1X,F7.3,1X,F7.3)
 26	FORMAT('HETATM',1X,I5,2X,'S',5X,I5,6X,F7.3,1X,F7.3,1X,F7.3)
	ETE=ETE/REAL(JPASOS)
	QT=QT/REAL(JPASOS)
	QTOH=QTOH/REAL(JPASOS)
	QT1=QT1/REAL(JPASOS)
	QT2=QT2/REAL(JPASOS)
	NCICON01=NCICON01/REAL(JPASOS)
	NCICON02=NCICON02/REAL(JPASOS)

88	FORMAT(E9.4,1X,E10.4,1X,E10.4,1X,E10.4,1X,E10.4,1X,E10.4,1X,
     +E10.4,1X,E10.4,1X,E10.4)
	WRITE(50,88)P,QT,QTOH,QT1,QT2,NCICON01,NCICON02,PH(0,0,-1),
     +PH(0,0,1)
	ETE=0.
	WRITE(*,*) P,NCI,NCIOH,QT
	CLOSE(40)
	CLOSE(21)
	CLOSE(22)
	FP=10.**DP
	P=P*FP
      ENDDO
      CALL CPU_TIME(AFINAL)
      WRITE(50,*)AFINAL-AINICIO
	CLOSE(50)
	CLOSE(97)
	END
C--------------------------------------------------------------------
!-----------------------------
      SUBROUTINE INCION(SIGMETANO,TEMP,SIGMA,EPS, RCUT, N,NCHAIN,
     +V,VA,VG,W,
     :CREATE,CR)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)


C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL CREATION                           **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP            TEMPERATURE                           **
C    ** REAL    Z               ABSOLUTE ACTIVITY                     **
C    ** REAL    SIGMA           LENNARD-JONES DIAMETER                **
C    ** REAL    RCUT            CUT-OFF DISTANCE                      **
C    ** REAL    V               POTENTIAL ENERGY                      **
C    ** REAL    W               VIRIAL                                **
C    ** INTEGER N               NUMBER OF ATOMS BEFORE TRIAL CREATION **
C    ** LOGICAL CREATE          TRUE FOR A SUCCESSFUL CREATION        **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 9000 )


        REAL        TEMP, Z, SIGMA, RCUT, V, W

        INTEGER     N
        LOGICAL     CREATE

        REAL        BETA, RXNEW, RYNEW, RZNEW, DELTV, DELTW, DELTCB
        REAL        RANF, DUMMY, RMIN
        INTEGER     NTRIAL
        LOGICAL     OVRLAP
		XMAX=1.
		YMAX=1.
		ZMAX=1.


C    *******************************************************************

        CREATE = .FALSE.
        BETA   = 1.0 / (TEMP)
        RMIN   = 0.5 * SIGMA
!---------------------------------------------------
        !SELECCIONA TIPO DE ION
        !
        NION=INT(RANF(DUMMY)*REAL(IONES))+1
        
        NCITRIAL = NCI(NION)+1
      
        IF ( NCITRIAL .GE. NMAX ) STOP 'MAXIMUM NUMBER OF ATOMS IN BOX'

C    ** GENERATE THE POSITION OF THE TRIAL ATOM **


C----------------------------------------------------------------------
C	CREA LAS COORDENADAS DEL CONTRAION
C-----------------------------------------------------------------------

C    ** GENERATE THE POSITION OF THE TRIAL ATOM **

        RXNEW  = (RANF(DUMMY)-0.5)*XMAX*ALX
        RYNEW  = (RANF(DUMMY)-0.5)*YMAX*ALY
        RZNEW  = (RANF(DUMMY)-0.5)*ZMAX*ALZ

        IF(POLCIL.EQ.0) THEN       
        RXNEWIM=RXNEW
        RYNEWIM=RYNEW
        RZNEWIM=-RZNEW
      ENDIF
      

	IF(POLCIL.EQ.1) THEN
	  RADIO=(RANF(DUMMY))*RCIL
        DELTARADIO=RCIL-RADIO
        RADIOEXTRA=RCIL+DELTARADIO
	TITA2=(RANF(DUMMY))*6.28
	  RXNEW  = RADIO*COS(TITA2)
        RYNEW  = RADIO*SIN(TITA2)
!----------------------------------------------------------------------
        RXNEWIM  = RADIOEXTRA*COS(TITA2)
        RYNEWIM  = RADIOEXTRA*SIN(TITA2)
        RZNEWIM  = RZNEW
        
	ENDIF
	IPULL=0
	IPULLCHAIN=0
C------------------------------------------------------------------


C------------------------------------------------------------------
	RXI=RXNEW
	RYI=RYNEW
	RZI=RZNEW

C--------------------------------------------------------------------------
C	ENERGÍA DE LOS CONTRAIONES
C-----------------------------------------------------------------------
	QCICON=REAL(NQCI(NION)*FCLEC)
      CALL ENERGYCHARGECION(RXI, RYI, RZI,QCICON,N,IPULL,IPULLCHAIN,
     +NCHAIN,SIGMA,EPS, RCUT,NION, 
     +DELTV3,
     + DELTW )

C--------------------------------------------------------------------
	DELTATOTAL2=DELTV3
C------------------------------------------------------------------------
C	CAMBIO EN EL LA INTERACCION SOLIDO-FLUIDO

	OVRLAP=.FALSE. 
      
        IF ( .NOT. OVRLAP ) THEN
      RANDOM=RANF ( DUMMY )
           DELTCB = BETA * DELTATOTAL2 - LOG(Z(NION)/REAL(NCITRIAL))
      !WRITE(*,*)'----------------------'
      !WRITE(*,*)' INCI '
      !WRITE(*,*)DELTCB,' DELTCB'
      !WRITE(*,*)DELTV3,' DELTV3'
      !WRITE(*,*)Z(NION),NION,' Z'
      !	WRITE(*,*)NCI(1),' NC1'
      !WRITE(*,*)NCI(2),' NC2'
      !WRITE(*,*)NCI(3),' NC3'
      !WRITE(*,*)NCI(4),' NC4'
      !WRITE(*,*)BETA, 'BETA '
      !WRITE(*,*)LOG(Z(NION)/REAL(NCITRIAL)),' LOG'
      !

           IF ( DELTCB .LT. 75.0 ) THEN

				IF ( DELTCB .LE. 0.0 ) THEN

      CALL ADDION ( RXNEW, RYNEW, RZNEW,RXNEWIM,RYNEWIM,RZNEWIM,NION)

                 NCI(NION)    = NCITRIAL
                 NCIM(NION)=NCITRIAL
              !WRITE(*,*)' SUCCES 1'
              !PAUSE
				ELSE IF ( EXP ( - DELTCB ) .GT. RANDOM ) THEN
      CALL ADDION ( RXNEW, RYNEW, RZNEW,RXNEWIM,RYNEWIM,RZNEWIM,NION)
                !  PAUSE
                 NCI(NION)    = NCITRIAL
                 NCIM(NION)=NCITRIAL
                  ELSE
                 !     WRITE(*,*)' NADA 1'
                 !     PAUSE
				RETURN
				ENDIF
           ELSE
                  !    WRITE(*,*)' NADA 2'
                   !   PAUSE

               RETURN
           ENDIF

        ENDIF
	

        RETURN
        END
                
      SUBROUTINE ADDION (RXNEW,RYNEW,RZNEW,RXNEWIM,RYNEWIM,RZNEWIM,NION)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
C    *******************************************************************
C    ** SUBROUTINE TO ADD AN ATOM TO THE ARRAY LOCATE.                **
C    **                                                               **
C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE NEW ADDITION   **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )




        REAL        RXNEW, RYNEW, RZNEW
        INTEGER     N, IPULL

        INTEGER     INEW, NTRIAL

C    *******************************************************************

        NTRIAL = NCI(NION) + 1

        INEW = LOCATE(NTRIAL,NION)
        NDISC=100
	I=NINT(RXNEW*REAL(NDISC))
	J=NINT(RYNEW*REAL(NDISC))
	K=NINT(RZNEW*REAL(NDISC))
	IJK=NINT(RZNEW*100)+50
      
      DXY(I,J,NION)=DXY(I,J,NION)+1
      
	IJKR=INT(SQRT(RXNEW**2+RYNEW**2)*100)+1
	DATOSN(IJKR,NION)=DATOSN(IJKR,NION)+1
	DCI(IJK,NION)=DCI(IJK,NION)+1.
        IF ( INEW .EQ. 0 ) THEN

C       ** ATOM REQUIRES A NEW NUMBER **

           LOCATE(NTRIAL,NION) = NTRIAL
           LOCATEIM(NTRIAL,NION) = NTRIAL
           INEW           = NTRIAL
        ENDIF

C    ** FIT NEW ATOM INTO THE ARRAY **

        RXCI(INEW,NION) = RXNEW
        RYCI(INEW,NION) = RYNEW
        RZCI(INEW,NION) = RZNEW

        RXCIM(INEW,NION) = RXNEWIM
        RYCIM(INEW,NION) = RYNEWIM
        RZCIM(INEW,NION) = RZNEWIM

!	!!!WRITE(*,*) INEW,RXCI(INEW),RYCI(INEW),RZCI(INEW),NTRIAL
!	 !!!!!!!!PAUSE
       RETURN
      END
      

!-----------------------------------------------------------------------
      SUBROUTINE OUTION (SIGMETANO, TEMP, SIGMA,EPS, RCUT,N,NCHAIN,
     +V,
     :VA,VG,W,GHOST )
C--------

C---------
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL DESTRUCTION                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP         TEMPERATURE                              **
C    ** REAL    Z            ABSOLUTE ACTIVITY                        **
C    ** REAL    SIGMA        LENNARD-JONES DIAMETER                   **
C    ** REAL    RCUT         CUT-OFF DISTANCE                         **
C    ** REAL    V            POTENTIAL ENERGY                         **
C    ** REAL    W            VIRIAL                                   **
C    ** INTEGER N            NUMBER OF ATOMS BEFORE TRIAL DESTRUCTION **
C    ** LOGICAL GHOST        TRUE FOR A SUCCESSFUL DESTRUCTION        **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )


        REAL        TEMP, Z, SIGMA, RCUT, V, W

        INTEGER     N
        LOGICAL     GHOST

        REAL        BETA, DELTV, DELTW, DELTDB, RANF, DUMMY
        INTEGER     NTRIAL, NLOC, IPULL
        LOGICAL     OVRLAP,CREATE
C    *******************************************************************
	BETA   = 1.0 / (TEMP)
              NION  = INT ( REAL ( IONES ) * RANF ( DUMMY ) ) + 1
        NCITRIAL = NCI(NION)-1
	
      IF(NCITRIAL.LT.0) THEN
	RETURN
	ENDIF
	
C------------------------------------------------------------------------
C    ** COUNTERION
      NLOC  = INT ( REAL ( NCITRIAL ) * RANF ( DUMMY ) ) + 1

	
        IPULL2 = LOCATE(NLOC,NION)
	

C------------------------------------------------------------------------


C	ENERGY CHANGE ON COUNTERION
	CALL POTOUTION (IPULL2,NION, IPULL, IPULLCHAIN,
     +  N,NCHAIN, SIGMA,EPS, RCUT,DELTV2, 
     +  DELTW)

C--------------------------------------------------------------------------
	DELTATOTAL2=(DELTV2)
      
	RANDOM=RANF ( DUMMY )
C    ** CHECK FOR ACCEPTANCE **
        DELTDB = BETA * (DELTATOTAL2)  -LOG(REAL(NCITRIAL+1)/Z(NION))
      !WRITE(*,*)'----------------------'
      !WRITE(*,*)' OUT '
      !WRITE(*,*)DELTDB,' DELTDB'
      !WRITE(*,*)DELTV2,' DELTV3'
      !WRITE(*,*)Z(NION),NION,' Z'
      !WRITE(*,*)NCI(1),' NC1'
      !WRITE(*,*)NCI(2),' NC2'
      !WRITE(*,*)NCI(3),' NC3'
      !WRITE(*,*)NCI(4),' NC4'
      !WRITE(*,*)BETA, 'BETA '
      !WRITE(*,*)LOG(REAL(NCITRIAL+1)/Z(NION)),' LOG'	
      
        IF ( DELTDB .LT. 75.0 ) THEN

           IF ( DELTDB .LT. 0.0 ) THEN
	                    
	
			CALL REMOVEION ( NLOC, IPULL2,NION )
	NCI(NION) = NCITRIAL
      NCIM(NION)=NCITRIAL
           ELSE IF ( EXP( -DELTDB ) .GT.RANF ( DUMMY )  ) THEN
			
			CALL REMOVEION ( NLOC, IPULL2,NION ) 
	NCI(NION) = NCITRIAL
      NCIM(NION)=NCITRIAL
      

		ELSE
		RETURN


           ENDIF

        ENDIF
!	!!!WRITE(*,*)DELTATOTAL2,DELTDB, Z(3), NCICON(1),' OUT '
 	!!!!!PAUSE


        RETURN
        END
!----------------------------------------------------------------------------------------------
      SUBROUTINE POTOUTION (IPULL2,NION, IPULL, IPULLCHAIN,
     +  N,NCHAIN, SIGMA,EPS, RCUT,DELTV, 
     +  DELTW)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA, DELTV, DELTW,KAPPA
        INTEGER     N, IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************


C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0

	  QNAT	=REAL(NQCI(NION))*FCLEC
        
        
C-------------------------------------------------------------
C----------------------------------------------------------------------------------------------------
      NES=1

        RXES(NES)    = RXCI(IPULL2,NION)
        RYES(NES)    = RYCI(IPULL2,NION)
        RZES(NES)    = RZCI(IPULL2,NION)

        ZES(1)=QNAT
          DO I=1,NC
            IF(QAC(I).NE.0) THEN
              NES=NES+1
              RXES(NES)  = RXC(I)
              RYES(NES)  = RYC(I)
              RZES(NES)  = RZC(I)

              ZES(NES)=QAC(I)
            ENDIF
        END DO
C-------------------------------------------------------------
C    ** LOOP OVER ALL ATOMS  IN ALL THE  CHAINS  **
	  DO JCHAIN=1,NCHAIN
        DO J = 1, N
	!
C       ** PICK ACTIVE ATOMS FROM LOCATE **

	!IF(JCHAIN.NE.IPULLCHAIN. AND .J.NE.IPULL) THEN
      IF(Q(J,JCHAIN).NE.0) THEN
		NES=NES+1
              RXES(NES)  = RX(J,JCHAIN)
              RYES(NES)  = RY(J,JCHAIN)
              RZES(NES)  = RZ(J,JCHAIN)

              ZES(NES)=Q(J,JCHAIN)
		NES=NES+1
              RXES(NES)  = XIM(J,JCHAIN)
              RYES(NES)  = YIM(J,JCHAIN)
              RZES(NES)  = ZIM(J,JCHAIN)

              ZES(NES)=QIM(J,JCHAIN)
          
              
              ENDIF


           

100     ENDDO
      ENDDO

C    ** LOOP OVER ALL COUNTERIONS  EXCEPT IPULL **
        DO J = 1, 4
          IF (J.NE.NION ) THEN
          DO JION=1,NCI(J)
C       ** PICK ACTIVE ATOMS FROM LOCATE **

           
              JIN = LOCATE(JION,j)
           !
               NES=NES+1

              RXES(NES)  = RXCI(JIN,J)
              RYES(NES)  = RYCI(JIN,J)
              RZES(NES)  = RZCI(JIN,J)

              ZES(NES)=REAL(NQCI(J))*FCLEC
              NES=NES+1
              RXES(NES)  = RXCIM(JIN,J)
              RYES(NES)  = RYCIM(JIN,J)
              RZES(NES)  = RZCIM(JIN,J)

              ZES(NES)=REAL(NQCI(J))*FCLEC
              
              ENDDO

          ENDIF
          ENDDO
          DO JION=1,NCI(NION)
              JIN = LOCATE(JION,NION)
              IF ( JIN .NE. IPULL2) THEN
               NES=NES+1

              RXES(NES)  = RXCI(JIN,NION)
              RYES(NES)  = RYCI(JIN,NION)
              RZES(NES)  = RZCI(JIN,NION)

              ZES(NES)=REAL(NQCI(NION))*FCLEC

              NES=NES+1

              RXES(NES)  = RXCIM(JIN,NION)
              RYES(NES)  = RYCIM(JIN,NION)
              RZES(NES)  = RZCIM(JIN,NION)

              ZES(NES)=REAL(NQCI(NION))*FCLEC

              ENDIF
          ENDDO
C-------------------------------------------------------------
      LTYPE=1
	CALL RWALD (VR,NES,LTYPE)

      
      DELTV=VR


        DELTV =   DELTV


C    ** CHANGE SIGN OF DELTV AND DELTW FOR A REMOVAL **

        DELTV = - DELTV
        DELTW = - DELTW
        RETURN
        END
        
      SUBROUTINE REMOVEION ( NLOC, IPULL2,NION )

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
C    *******************************************************************
C    ** SUBROUTINE TO REMOVE AN ATOM FROM THE ARRAY LOCATE.           **
C    **                                                               **
C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE REMOVAL.       **
C    ** ELEMENT IPULL OF LOCATE IS TO BE DESTROYED.                   **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )



        INTEGER     NCI, IPULL2, NLOC

        INTEGER     K
	
	!!!!WRITE(*,*) NLOC,NCICON(1),IPULL2, ' NLOC,NCI, IPULL EN REMOVE1'
!	!!!WRITE(*,*) ' BAZINGA!!'
C    *******************************************************************

        IF ( NLOC .LT. NCI(NION) ) THEN

C       ** CLOSE UP THE ARRAY LOCATE AFTER THE REMOVAL **

           DO K = NLOC + 1, NCI(NION)

              LOCATE(K - 1,NION) = LOCATE(K,NION)
              LOCATEIM(K - 1,NION) = LOCATEIM(K,NION)

         ENDDO

C       ** PLACE THE GHOST ATOM IPULL JUST OUTSIDE THE ACTIVE **
C       ** RANGE OF THE ARRAY LOCATE FOR FUTURE USE           **

           LOCATE(NCI(NION),NION) = IPULL2
           LOCATEIM(NCI(NION),NION) = IPULL2
!	!!!WRITE(*,*) LOCATE(NCI),NCI, ' LOCATE EN REMOVE!!!!!!!!'
!	!!!!!!!!PAUSE

        ENDIF
!	!!!WRITE(*,*) LOCATE(NCI),NCI, ' LOCATE EN REMOVE'
!	!!!!!!!!PAUSE

        RETURN
        END
      
C--------------------------------------------------------------------
        SUBROUTINE MOVE(SIGMETANO,TEMP,SIGMA,EPS,RCUT,N,NCHAIN,V,VA,VG,
     +W,GHOST,ITETH,ALREF)
C--------------------------------------------------------------------
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)


        INTEGER     NMAX
        PARAMETER ( NMAX = 1500 )


        REAL        TEMP, Z, SIGMA, RCUT, V, W

        INTEGER     N
        LOGICAL     GHOST

        REAL        BETA, DELTV, DELTW,  RANF, DUMMY
        INTEGER     NTRIAL, NLOC, IPULL
        LOGICAL     OVRLAP,CREATE
!	DOUBLE PRECISION COSTITA,SENTITA,RXNEW,RYNEW,RZNEW

C    *******************************************************************
	CREATE=.FALSE.
        GHOST  = .FALSE.
        BETA   = 1.0 / (TEMP)
        NTRIAL = N 
        IF ( NTRIAL .EQ. 0 ) STOP 'ONaly ONE ATOM REMAINS'

C    ** PICK A RANDOM ELEMENT FROM THE ACTIVE PART OF LOCATE **

        NLOC  = INT ( REAL ( NTRIAL-iteth ) * RANF ( DUMMY ) ) + iteth+1
        IPULL = NLOC
	  IPULLCHAIN=INT ( REAL ( NCHAIN ) * RANF ( DUMMY ) ) + 1
	
      !ipullchain=nchain
      
C----------------------------------------------------------------
C	CONSERVO LAS POSICIONES INICIALES
C----------------------------------------------------------------
	RXI=RX(IPULL,IPULLCHAIN)
	RYI=RY(IPULL,IPULLCHAIN)
	RZI=RZ(IPULL,IPULLCHAIN)
      
      RXIIM=XIM(IPULL,IPULLCHAIN)
      RYIIM=YIM(IPULL,IPULLCHAIN)
      RZIIM=ZIM(IPULL,IPULLCHAIN)

	DELTA=(alref/10)/(sigmetano/sigma)

	VN1=V
	VG1=VG
	VA1=VA
      !!write(*,*)vn1
      !!!pause
C    ** CALCULATE ENERGY CHANGE ON REMOVAL OF ATOM IPULL **

	CALL ENERGY ( IPULL, IPULLCHAIN,N,NCHAIN, SIGMA,EPS, RCUT, DELTV,
     + DELTW )
      !!write(*,*) deltv
      !!!pause
	CALL energycharge(IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS, RCUT, 
     +DELTV1,
     + DELTW )

    	CALL ADPOTIN(SIGMETANO,IPULL, IPULLCHAIN, SIGMA,EPS, RCUT, 
     +DELTVA,DELTW)

      !write(*,*) deltv,deltv1,deltva, 'move1'
	DELTATOTAL1=(DELTV+DELTVA+deltv1)

	VN=V+DELTATOTAL1
	VGN=VG+DELTV
	VAN=VA+DELTVA

	NTIPO=INT(RANF(DUMMY)*3)+1
!	!write(*,*)ntipo,ipull
	if (ntipo.ge.2) ntipo=2
              if(iteth.ne.1) then
        if(ipull.eq.int(n/2)) return
        endif

C----------------------------------------------------

C    ** GENERATE THE POSITION OF THE TRIAL ATOM **
C-----------------------------------------------------------------
c	Transformacion de rotacion
!	ntipo=2
C-----------------------------------------------------------------
	GOTO(60,59) NTIPO

59	IF (IPULL.LT.N.and.ipull.gt.1) THEN
C------------------------------------------------------------------
C	BUSCO EL VECTOR UNITARIO (A,B,C)
C------------------------------------------------------------------
	distplane1=sqrt((RX(IPULL,IPULLCHAIN)-RX(IPULL+1,IPULLCHAIN))**2+
     +(RY(IPULL,IPULLCHAIN)-RY(IPULL+1,IPULLCHAIN))**2+
     +(Rz(IPULL,IPULLCHAIN)-Rz(IPULL+1,IPULLCHAIN))**2)


	distplane2=sqrt((RX(IPULL,IPULLCHAIN)-RX(IPULL-1,IPULLCHAIN))**2+
     +(RY(IPULL,IPULLCHAIN)-RY(IPULL-1,IPULLCHAIN))**2+
     +(Rz(IPULL,IPULLCHAIN)-Rz(IPULL-1,IPULLCHAIN))**2)
	
      !!write(*,*) distplane1,distplane2, ' dist1 dist2'
	RXM1=RX(IPULL-1,IPULLCHAIN)
	RYM1=RY(IPULL-1,IPULLCHAIN)
	RZM1=RZ(IPULL-1,IPULLCHAIN)

	RXP1=RX(IPULL+1,IPULLCHAIN)
 	RYP1=RY(IPULL+1,IPULLCHAIN)
	RZP1=RZ(IPULL+1,IPULLCHAIN)

	RXCHAIN1=RX(IPULL,IPULLCHAIN)
	RYCHAIN1=RY(IPULL,IPULLCHAIN)
	RZCHAIN1=RZ(IPULL,IPULLCHAIN)


      IF(DISTPLANE1.GT.0.1.and.distplane2.lt.0.1) THEN

	RTESTX=RXCHAIN1-RXM1
	RTESTY=RYCHAIN1-RYM1
	RTESTZ=RYCHAIN1-RZM1


		IF(ABS(RTESTX).GT.0.1)  RXM1=RXM1+SIGN(ALX,RTESTX)
		IF(ABS(RTESTY).GT.0.1)  RYM1=RYM1+SIGN(ALY,RTESTY)
		IF(ABS(RTESTZ).GT.0.1)  RZM1=RZM1+SIGN(ALZ,RTESTZ)
		
      elseIF(DISTPLANE2.GT.0.1.and.disTplane1.lt.0.1) THEN

	RTESTX=RXCHAIN1-RXP1
	RTESTY=RYCHAIN1-RYP1
	RTESTZ=RYCHAIN1-RZP1

		IF(ABS(RTESTX).GT.0.1)  RXP1=RXP1+SIGN(ALX,RTESTX)
		IF(ABS(RTESTY).GT.0.1)  RYP1=RYP1+SIGN(ALY,RTESTY)
		IF(ABS(RTESTZ).GT.0.1)  RZP1=RZP1+SIGN(ALZ,RTESTZ)
      elseIF(DISTPLANE2.GT.0.1.and.disTplane1.Gt.0.1) THEN
	

	!GOTO 60
	RTESTX=RXCHAIN1-RXP1
	RTESTY=RYCHAIN1-RYP1
	RTESTZ=RYCHAIN1-RZP1

		IF(ABS(RTESTX).GT.0.1)  RXCHAIN1=RXCHAIN1-SIGN(ALX,RTESTX)
		IF(ABS(RTESTY).GT.0.1)  RYCHAIN1=RYCHAIN1-SIGN(ALY,RTESTY)
		IF(ABS(RTESTZ).GT.0.1)  RZCHAIN1=RZCHAIN1-SIGN(ALZ,RTESTZ)


	ENDIF



	A=RXM1-RXP1
	B=RYM1-RYP1
	C=RZM1-RZP1
!	!!!write(*,*) A,B,C,' A, B, C'
!	A=RX(IPULL-1,IPULLCHAIN)-RX(IPULL+1,IPULLCHAIN)
!	B=RY(IPULL-1,IPULLCHAIN)-RY(IPULL+1,IPULLCHAIN)
!	C=RZ(IPULL-1,IPULLCHAIN)-RZ(IPULL+1,IPULLCHAIN)

!	if(distplane1.gt.0.1. AND .distplane2.gt.0.1) goto 60




!	A=A-alx*anint(A/alx)*pbx
!	B=B-aly*anint(B/aly)*pby

!      A= A - ANINT ( A )
!      B= B- ANINT ( B )



	AMODULO=SQRT(A**2+B**2+C**2)

	A=A/AMODULO

	B=B/AMODULO

	C=C/AMODULO

	D=SQRT(B**2+C**2)

	if(d.eq.0) goto 60

	X1=RXCHAIN1 !rX(IPULL-1,IPULLCHAIN)
	Y1=RYCHAIN1 !RY(IPULL-1,IPULLCHAIN)
	Z1=RZCHAIN1 !RZ(IPULL-1,IPULLCHAIN)

	X2=RX(IPULL,IPULLCHAIN)
 	Y2=RY(IPULL,IPULLCHAIN)
	Z2=RZ(IPULL,IPULLCHAIN)

	TITA=(real((RANF(DUMMY))))*6.28
	z12=z1-z2
	z12=z12-anint(z12)
	z21=-z12
	y12=y1-y2
	y12=y12-anint(y12)
	y21=-y12
	x12=x1-x2
	x12=x12-anint(x12)
	x21=-x12

	COSTITA=COS(TITA)
	SENTITA=SIN(TITA)


      RXNEW=(a*(costita-1)*(b*(y12)+c*(z12))+
     +b*sentita*(z21)+c*sentita*(y12)+costita*d**2*(x21)+d**2*
     +(x12))/(a**2+d**2)+x2
      
	RYNEW=b*(a*(b*sentita*(z12)+c*sentita*(y21)+
     +d**2*(costita-1)*(x12))+d**2*(costita-1)*(b*(y12)+
     +c*(z12)))/((a**2+d**2)*(b**2+c**2))+(a*c*sentita
     $*(b*(y12)+c*(z12))+b**2*(y1-costita*(y12))-c*(c*(costita*(y12)
     $-y1)+d**2*sentita*(x12)))/(b**2+c**2)

      
	RZNEW=c*(a*(b*sentita*(z12)+c*sentita*(y21)+
     +d**2*(costita-1)*(x12))+d**2*(costita-1)*(b*(y12)+
     +c*(z12)))/((a**2+d**2)*(b**2+c**2))-(a*b*sentita*(b*(y12)+
     +c*(z12))+b**2*(costita*(z12)-z1)+b*d**2*sentita*(x21)+
     +c**2*(costita*(z12)-z1))/(b**2+c**2)
   
C--------------------------------------------------------------------
      ELSE

	
          RXNEW  =RXI+ (RANF(DUMMY)-0.5)*DELTA
		RYNEW  =RYI+(RANF(DUMMY)-0.5)*DELTA
		RZNEW  =RZI+ (ranf(dummy)-0.5)*DELTA
 
	ENDIF

	rxnew=rxnew-alx*anint(rxnew/alx)*pbx
	rynew=rynew-aly*anint(rynew/aly)*pby
	rznew=rznew-alz*anint(rznew/alz)*pbz
    
	if(abs(rznew).gt.alz/2)return
	if(abs(rynew).gt.aly/2)return
	if(abs(rxnew).gt.alx/2)return
      
      rpol2=rxnew**2+rynew**2
      rcil2=rcil**2
      if(rpol2.gt.rcil2.and.polcil.eq.1)return		

	RX(IPULL,IPULLCHAIN)=rxnew
	RY(IPULL,IPULLCHAIN)=RYnew
	RZ(IPULL,IPULLCHAIN)=RZnew
!	!!!write(*,*)'move 2'	
      
	RX(IPULL,IPULLCHAIN)=RXNEW
	RY(IPULL,IPULLCHAIN)=RYNEW
	RZ(IPULL,IPULLCHAIN)=RZNEW

      if(polcil.eq.0) then
      !Posiciones en el plano
      
      XIM(IPULL,IPULLCHAIN)=RX(IPULL,IPULLCHAIN)
      YIM(IPULL,IPULLCHAIN)=RY(IPULL,IPULLCHAIN)
      ZIM(IPULL,IPULLCHAIN)=-RZ(IPULL,IPULLCHAIN)
      endif
!------------------------------------------------------------------------------------
     	IF(POLCIL.eq.1) THEN
      !Posiciones en el cilindro
          
	  RADIO=SQRT(RXNEW**2+RYNEW**2)
        DELTARADIO=RCIL-RADIO
        RADIOEXTRA=RCIL+DELTARADIO
!----------------------------------------------------------------------
        XIM(IPULL,IPULLCHAIN)  = RADIOEXTRA*RXNEW/RADIO
        YIM(IPULL,IPULLCHAIN)  = RADIOEXTRA*RYNEW/RADIO
        ZIM(IPULL,IPULLCHAIN)  = RZNEW
        
	ENDIF

 
	GOTO 50
!	!!!WRITE(*,*) RXI,RYI,RZI
      
!**********************************************************************************
60	  RXNEW  =RXI+ (RANF(DUMMY)-0.5)*DELTA
		RYNEW  =RYI+(RANF(DUMMY)-0.5)*DELTA
		RZNEW  =RZI+ (ranf(dummy)-0.5)*DELTA

	rxnew=rxnew-alx*anint(rxnew/alx)*pbx
	rynew=rynew-aly*anint(rynew/aly)*pby
	rznew=rznew-alz*anint(rznew/alz)*pbz
!	!!!write(*,*)rxnew,rynew,rznew, ipull,ipullchain
            rpol2=rxnew**2+rynew**2
      rcil2=rcil**2
!--------------------------------------------------------
!     CONTROLES
      if(rpol2.gt.rcil2.and.polcil.ne.0)return
	if(abs(rznew).gt.alz/2)return
	if(abs(rynew).gt.aly/2)return
	if(abs(rxnew).gt.alx/2)return

	
 
	RX(IPULL,IPULLCHAIN)=rxnew
	RY(IPULL,IPULLCHAIN)=RYnew
	RZ(IPULL,IPULLCHAIN)=RZnew
	
!-----------------------------------------------------------------

      

      IF(POLCIL.EQ.0) THEN	
      XIM(IPULL,IPULLCHAIN)=RX(IPULL,IPULLCHAIN)
      YIM(IPULL,IPULLCHAIN)=RY(IPULL,IPULLCHAIN)
      ZIM(IPULL,IPULLCHAIN)=-RZ(IPULL,IPULLCHAIN)
      ENDIF
       
!------------------------------------------------------------------------------------
     	IF(POLCIL.eq.1) THEN
      !Posiciones en el cilindro
          
	  RADIO=SQRT(RXNEW**2+RYNEW**2)
        DELTARADIO=RCIL-RADIO
        RADIOEXTRA=RCIL+DELTARADIO
!----------------------------------------------------------------------
        XIM(IPULL,IPULLCHAIN)  = RADIOEXTRA*RXNEW/RADIO
        YIM(IPULL,IPULLCHAIN)  = RADIOEXTRA*RYNEW/RADIO
        ZIM(IPULL,IPULLCHAIN)  = RZNEW
        
	ENDIF

C--------------------------------------------------------
C---------------------------------------------------------
c--------------------------------------------------------
C---------------------------------------------------------
50	CALL ENERGY( IPULL, IPULLCHAIN,N,NCHAIN, SIGMA,EPS, RCUT, DELTV,
     + DELTW )
	CALL energycharge(IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS, RCUT, 
     +DELTV1,
     + DELTW )
    	CALL ADPOTIN(SIGMETANO,IPULL, IPULLCHAIN, SIGMA,EPS, RCUT, 
     +DELTVA,DELTW)
      !write(*,*) deltv,deltv1,deltva, 'move2'
      !pause
      
	if(deltv.gt.25) then
	RX(IPULL,IPULLCHAIN)=RXI
	RY(IPULL,IPULLCHAIN)=RYI
	RZ(IPULL,IPULLCHAIN)=RZI
      XIM(IPULL,IPULLCHAIN)=RXIIM
      YIM(IPULL,IPULLCHAIN)=RYIIM
      ZIM(IPULL,IPULLCHAIN)=RZIIM
	return
	endif
C------------------------------------------------------------------------
C	CAMBIO EN EL LA INTERACCION SOLIDO-FLUIDO

C    ** CHECK FOR ACCEPTANCE **
!	!!!write(*,*) deltva,' 2'
!	!!!!!!!!pause

	DELTATOTAL=(DELTV+DELTVA+deltv1)
	DELTATOTAL2=deltatotal-deltatotal1
      !write(*,*) deltatotal2, 'deltatotal'
	!pause
	VN2=VN+DELTATOTAL
	VGN2=VGN+DELTV
 	VAN2=VAN+DELTVA
	OVRLAP=.FALSE.
      IF ( .NOT. OVRLAP ) THEN
		DELTCB = (deltatotal2)*BETA

          IF ( DELTCB .LT. 75.0 ) THEN
               IF ( DELTCB .LE. 0.0 ) THEN
                 CREATE = .TRUE.
                 V    = VN2 
	           VG	=VGN2
			   VA	=VAN2
                 W    = W + DELTW
                 N    = NTRIAL
!	!!!write(*,*)'move 7'
				ELSE IF ( EXP ( - DELTCB ) .GT. RANF ( DUMMY ) ) THEN

                 CREATE = .TRUE.

                 V    = VN2
	           VG	=VGN2
 			   VA	=VAN2

                 W    = W + DELTW
                 N    = NTRIAL
!	!!!write(*,*)'move 8'
	ELSE

	RX(IPULL,IPULLCHAIN)=RXI
	RY(IPULL,IPULLCHAIN)=RYI
	RZ(IPULL,IPULLCHAIN)=RZI
	XIM(IPULL,IPULLCHAIN)=RXIIM
      YIM(IPULL,IPULLCHAIN)=RYIIM
      ZIM(IPULL,IPULLCHAIN)=RZIIM
              ENDIF
	ELSE
	RX(IPULL,IPULLCHAIN)=RXI
	RY(IPULL,IPULLCHAIN)=RYI
	RZ(IPULL,IPULLCHAIN)=RZI
	XIM(IPULL,IPULLCHAIN)=RXIIM
      YIM(IPULL,IPULLCHAIN)=RYIIM
      ZIM(IPULL,IPULLCHAIN)=RZIIM
           ENDIF
	ELSE
	RX(IPULL,IPULLCHAIN)=RXI
	RY(IPULL,IPULLCHAIN)=RYI
	RZ(IPULL,IPULLCHAIN)=RZI
	XIM(IPULL,IPULLCHAIN)=RXIIM
      YIM(IPULL,IPULLCHAIN)=RYIIM
      ZIM(IPULL,IPULLCHAIN)=RZIIM

        ENDIF

        RETURN
        END


c-------------------------------------------------------------------------------------------------
C-------------------------------------------------------------------------------------------------
      SUBROUTINE IN (TEMP,SIGMA,EPS, RCUT, N,NCHAIN,
     :CREATE,NCHARGEMONOMER1,NCHARGEMONOMER2,PHBULK,PHBULK2)
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)



	

C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL CREATION                           **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP            TEMPERATURE                           **
C    ** REAL    Z               ABSOLUTE ACTIVITY                     **
C    ** REAL    SIGMA           LENNARD-JONES DIAMETER                **
C    ** REAL    RCUT            CUT-OFF DISTANCE                      **
C    ** REAL    V               POTENTIAL ENERGY                      **
C    ** REAL    W               VIRIAL                                **
C    ** INTEGER N               NUMBER OF ATOMS BEFORE TRIAL CREATION **
C    ** LOGICAL CREATE          TRUE FOR A SUCCESSFUL CREATION        **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )


        REAL        TEMP, Z, SIGMA, RCUT !, V, W

        INTEGER     N !,NC
        LOGICAL     CREATE

        REAL        BETA, RXNEW, RYNEW, RZNEW, DELTV1, DELTV3,
     +DELTW, DELTCB,KA,KA2
        REAL        RANF, DUMMY, RMIN
        INTEGER     NCITRIAL
        LOGICAL     OVRLAP
		XMAX=1.
		YMAX=1.
		ZMAX=1.


C    *******************************************************************

        CREATE = .FALSE.
        BETA   = 1.0 / (TEMP)
        RMIN   = 0.5 * SIGMA
  !      NCITRIAL = NCI+1

        !  WRITE(*,*) 'IN ETHER'
  !      IF ( NCITRIAL .GE. NMAX ) STOP 'MAXIMUM NUMBER OF ATOMS IN BOX'

C    ** GENERATE THE POSITION OF THE TRIAL ATOM **

        IPULL = (INT ( REAL ( N/1) * RANF ( DUMMY ) ) + 1)*1
	  IPULLCHAIN=INT ( REAL ( NCHAIN ) * RANF ( DUMMY ) ) + 1
	
       QANT=Q(IPULL,IPULLCHAIN)
	KTYPE=KAM(IPULL,IPULLCHAIN)
!       write (*,*)KTYPE, ' KTYPE ANTES'
      IF(KTYPE.NE.0) RETURN
      
      KTYPE=INT(RANF(DUMMY)*2)+1
      
	I=INT(RX(IPULL,IPULLCHAIN)*NDISC)
	J=INT(RY(IPULL,IPULLCHAIN)*NDISC)
	K=INT(RZ(IPULL,IPULLCHAIN)*NDISC)

	Z(5)=(PHBULK+LOG10(KA)) *(REAL(-NCHARGEMONOMER1))

      IF(KA.EQ.0.AND.KTYPE.EQ.1) RETURN
	IF(KA2.EQ.0.AND.KTYPE.EQ.2) THEN
	RETURN
      ELSE
      Z(6)=(PHBULK2+LOG10(KA2)) *(REAL(-NCHARGEMONOMER1))
	ENDIF
      KAM(IPULL,IPULLCHAIN)=KTYPE
      KTYPE2=KTYPE+4
 
C    ** CALCULATE ENERGY CHANGE ON ADDITION **
C--------------------------------------------------------------
C	ENERGÍA DE CREACION DE LA CARGA DEL SEGMENTO CARGADO
C---------------------------------------------------------------
	IF(KTYPE.EQ.1) THEN
          Q(IPULL,IPULLCHAIN)=REAL(NCHARGEMONOMER1)*FCLEC
          
      ENDIF
      
	
      IF(KTYPE.EQ.2) THEN
          Q(IPULL,IPULLCHAIN)=REAL(NCHARGEMONOMER2)*FCLEC
          
      ENDIF
      
    
      QIM(IPULL,IPULLCHAIN)=Q(IPULL,IPULLCHAIN)
      
      CALL ENERGYCHARGE(IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS, RCUT, 
     +DELTV1,
     + DELTW )
   
C-----------------------------------------------------------------------

	DELTATOTAL2=DELTV1 
C------------------------------------------------------------------------
C	CAMBIO EN EL LA INTERACCION SOLIDO-FLUIDO
	
C    ** CHECK FOR ACCEPTANCE **
C	!!!WRITE(*,*) OVRLAP
     
	OVRLAP=.FALSE.

        IF ( .NOT. OVRLAP ) THEN

           DELTCB = BETA * DELTATOTAL2 - LOG(10.)*Z(KTYPE2) 
           
           IF ( DELTCB .LT. 75.0 ) THEN

				IF ( DELTCB .LE. 0.0 ) THEN
	

	IF(KTYPE.EQ.1) THEN
          NK(KTYPE)=NK(KTYPE)+1
      ELSE
          NK(KTYPE)=NK(KTYPE)+1
      ENDIF
      

				ELSE IF ( EXP ( - DELTCB ) .GT. RANF ( DUMMY ) ) THEN

	IF(KTYPE.EQ.1) THEN
          NK(KTYPE)=NK(KTYPE)+1
      ELSE
          NK(KTYPE)=NK(KTYPE)+1
      ENDIF
				ELSE
      Q(IPULL,IPULLCHAIN)=0
      QIM(IPULL,IPULLCHAIN)=Q(IPULL,IPULLCHAIN)
      KAM(IPULL,IPULLCHAIN)=0
      !KAM(IPULL,IPULLCHAIN)=0
                 
				ENDIF
		ELSE
      Q(IPULL,IPULLCHAIN)=0
      QIM(IPULL,IPULLCHAIN)=Q(IPULL,IPULLCHAIN)
      KAM(IPULL,IPULLCHAIN)=0
           ENDIF

        ENDIF


        RETURN
        END

C------------------------------------------------------------------------------------
C------------------------------------------------------------------------------------

        SUBROUTINE OUT ( TEMP, SIGMA,EPS, RCUT,N,NCHAIN,
     +NCHARGEMONOMER1,NCHARGEMONOMER2,PHBULK,PHBULK2)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL DESTRUCTION                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP         TEMPERATURE                              **
C    ** REAL    Z            ABSOLUTE ACTIVITY                        **
C    ** REAL    SIGMA        LENNARD-JONES DIAMETER                   **
C    ** REAL    RCUT         CUT-OFF DISTANCE                         **
C    ** REAL    V            POTENTIAL ENERGY                         **
C    ** REAL    W            VIRIAL                                   **
C    ** INTEGER N            NUMBER OF ATOMS BEFORE TRIAL DESTRUCTION **
C    ** LOGICAL GHOST        TRUE FOR A SUCCESSFUL DESTRUCTION        **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )
		REAL KA,KA2

        REAL        TEMP, Z, SIGMA, RCUT !, V, W

        INTEGER     N
        !LOGICAL     GHOST

        REAL        BETA, DELTW, DELTDB, RANF, DUMMY
        INTEGER     NCITRIAL, NLOC, IPULL
        !LOGICAL     CREATE
C    *******************************************************************
	BETA   = 1.0 / (TEMP)	
 
	
C	** SEGMENTO
	IPULL = (INT ( REAL ( N/1) * RANF ( DUMMY ) ) + 1)*1
	IPULLCHAIN=INT ( REAL ( NCHAIN ) * RANF ( DUMMY ) ) + 1
	
C------------------------------------------------------------------------
	QTEMP=0	
!      KTYPE1=ABS(Q(IPULL,IPULLCHAIN))/FCLEC
	IF(Q(IPULL,IPULLCHAIN).EQ.QTEMP) THEN 
         
	RETURN
      ENDIF
	
      KTYPE=KAM(IPULL,IPULLCHAIN)

      I=INT(RX(IPULL,IPULLCHAIN)*NDISC)
	J=INT(RY(IPULL,IPULLCHAIN)*NDISC)
	K=INT(RZ(IPULL,IPULLCHAIN)*NDISC)


	Z(5)=(PHBULK+LOG10(KA)) *(-REAL(NCHARGEMONOMER1))
	IF(KA2.GT.0) THEN
	Z(6)=(PHBULK2+LOG10(KA2)) *(-REAL(NCHARGEMONOMER1))
	ELSEIF(KA2.EQ.0.AND.KTYPE.EQ.2) THEN
	RETURN
	ENDIF
      KTYPE2=KTYPE+4
      IF(KTYPE.EQ.2) QNEW=0
      IF(KTYPE.EQ.3) QNEW=0
      
C------------------------------------------------------------------------

C    ** CALCULATE ENERGY CHANGE ON ADDITION CHARGE ON POALYMER SEGMENT

      CALL ENERGYOUT(IPULL2,IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS, 
     +RCUT, 
     +DELTV1,
     + DELTW)



C--------------------------------------------------------------------------
	DELTATOTAL2=(DELTV1)
        DELTDB = BETA * (DELTATOTAL2)  + LOG(10.)*Z(KTYPE2)
        IF ( DELTDB .LT. 75.0 ) THEN

           IF ( DELTDB .LT. 0.0 ) THEN
	      Q(IPULL,IPULLCHAIN)=QNEW 
            QIM(IPULL,IPULLCHAIN)=Q(IPULL,IPULLCHAIN)
            KAM(IPULL,IPULLCHAIN)=0
!          KAM(IPULL,IPULLCHAIN)=0
              
	IF(KTYPE.EQ.1)THEN
          NK(KTYPE)=NK(KTYPE)-1
      ELSE
          NK(KTYPE)=NK(KTYPE)-1
      ENDIF
      
           ELSE IF ( EXP( -DELTDB ) .GT. RANF ( DUMMY ) ) THEN
	      Q(IPULL,IPULLCHAIN)=0
            QIM(IPULL,IPULLCHAIN)=0
            KAM(IPULL,IPULLCHAIN)=0
	IF(KTYPE.EQ.1)THEN
          NK(KTYPE)=NK(KTYPE)-1
      ELSE
          NK(KTYPE)=NK(KTYPE)-1
          
      ENDIF

		ELSE
		RETURN


           ENDIF
	

        ENDIF
        RETURN
        END
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------

C-----------------------------------------------------------------------------


      SUBROUTINE ENERGY(IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS, RCUT, 
     +DELTV,
     + DELTW )


	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA, DELTV, DELTW
        INTEGER     N, IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA
        SIGCUB = SIGSQ * SIGMA

C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
        WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0

        RXI    = RX(IPULL,IPULLCHAIN)
        RYI    = RY(IPULL,IPULLCHAIN)
        RZI    = RZ(IPULL,IPULLCHAIN)
	  QNAT	=Q(IPULL,IPULLCHAIN)
     
C    ** LOOP OVER ALL ATOMS  IN ALL THE  CHAINS EXCEPT IPULL **
	  DO JCHAIN=1,NCHAIN
        DO J = 1, N

C       ** PICK ACTIVE ATOMS FROM LOCATE **


		

           IF  (JCHAIN.NE.IPULLCHAIN) THEN

              RXIJ  = RXI - RX(J,JCHAIN)
              RYIJ  = RYI - RY(J,JCHAIN)
              RZIJ  = RZI - RZ(J,JCHAIN)

              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
			RIJ=SQRT(RIJSQ)
			I=RIJ*1000+1
          IF(RIJ.GT.1) then
                  
          VIJ=0
          else
               VIJ   = USS(I)
          endif
          
	    	
      	IF(RIJ.LT.0.5*SIGMA) THEN
      	VIJ=100000000
        	ENDIF

              DELTV = DELTV + VIJ !+VIJE
              


           ENDIF
	 

      ENDDO
        ENDDO
        
  	IF (IPULL.GT.3) THEN
        DO J = 1, IPULL-2

C       ** PICK ACTIVE ATOMS FROM LOCATE **

              RXIJ  = RXI - RX(J,IPULLCHAIN)
              RYIJ  = RYI - RY(J,IPULLCHAIN)
              RZIJ  = RZI - RZ(J,IPULLCHAIN)

              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

			RIJ=SQRT(RIJSQ)

			I=RIJ*1000+1
              VIJ   = USS(I)

	    	IF(RIJ.GT.1) VIJ=0
              DELTV = DELTV + VIJ !+VIJE
        
	     ENDDO
      ENDIF
  	IF(IPULL.LT.N-2) THEN
        DO J = IPULL+2, N

C       ** PICK ACTIVE ATOMS FROM LOCATE **


		



              RXIJ  = RXI - RX(J,IPULLCHAIN)
              RYIJ  = RYI - RY(J,IPULLCHAIN)
              RZIJ  = RZI - RZ(J,IPULLCHAIN)

              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

			RIJ=SQRT(RIJSQ)

			I=RIJ*1000+1
              VIJ   = USS(I)
	    	IF(RIJ.GT.1) VIJ=0
              DELTV = DELTV + VIJ !+VIJE


	     ENDDO
      ENDIF
  	IPULLF=IPULL+1
	IPULLI=IPULL-1
      
	IF(IPULL.EQ.N) THEN
  
              RXIJ  = RXI - RX(IPULL-1,IPULLCHAIN)
              RYIJ  = RYI - RY(IPULL-1,IPULLCHAIN)
             RZIJ  = RZI - RZ(IPULL-1,IPULLCHAIN)
!
              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

			RIJ=SQRT(RIJSQ)
!	!!!WRITE(*,*)RIJ

			I=RIJ*1000+1
	        VIJRES=USRE(I)
!              write(*,*)vijres
!              pause
              
	    	IF(RIJ.GT.1) THEN
			!VIJE=0
			VIJRES=1E6
              ENDIF
              DELTV = DELTV + VIJRES
              ELSEIF (IPULL.EQ.1) THEN
	
                  
	        RXIJ  = RXI - RX(IPULL+1,IPULLCHAIN)
              RYIJ  = RYI - RY(IPULL+1,IPULLCHAIN)
             RZIJ  = RZI - RZ(IPULL+1,IPULLCHAIN)
!
              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

			RIJ=SQRT(RIJSQ)

			I=RIJ*1000+1
              VIJRES=USRE(I)
	    	IF(RIJ.GT.1) THEN
			!VIJE=0
			VIJRES=1E6
              ENDIF
              DELTV = DELTV +VIJRES
  	ELSE

	DO J = IPULLI, IPULLF

	IF(J.NE.IPULL) THEN

              RXIJ  = RXI - RX(J,IPULLCHAIN)
              RYIJ  = RYI - RY(J,IPULLCHAIN)
             RZIJ  = RZI - RZ(J,IPULLCHAIN)
!
              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

			RIJ=SQRT(RIJSQ)

			I=RIJ*1000+1
              !VIJE=0
          	VIJRES=USRE(I)
	    	IF(RIJ.GT.1) THEN
			!VIJE=0
			VIJRES=1E6
              ENDIF
              if(n.eq.1)VIJRES=0
              DELTV = DELTV + VIJRES
      ENDIF        
!
!
	ENDDO
	ENDIF



        DELTV =   DELTV



        RETURN
        END
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
	REAL FUNCTION RANF ( DUMMY )

C    *******************************************************************
C    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
C    **                                                               **
C    **                 ***************                               **
C    **                 **  WARNING  **                               **
C    **                 ***************                               **
C    **                                                               **
C    ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
C    ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
C    *******************************************************************

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
        REAL        DUMMY
        SAVE        SEED
        DATA        SEED / 0 /

C    *******************************************************************

        SEED = MOD ( SEED * L + C, M )
        RANF = REAL ( SEED ) / M

        RETURN
        END
C-------------------------------------------------------------------------------
C------------------------------------------------------------------------
C   SUBRUTINA ADPOTIN
C   CALCULA LA ENERGIA DE INTERACCION ENTRE LA SUPERFICIE Y EL FLUIDO
C------------------------------------------------------------------------
      SUBROUTINE ADPOTIN(SIGMETANO,IPULL, IPULLCHAIN, SIGMA,EPS, 
     +RCUT, 
     +DELTV,DELTW)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WITH SURFACE ON ADDING AN ATOM.        **
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

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, RMIN, SIGMA, RXI, RYI, RZI, DELTV, DELTW
        INTEGER     N, NTRIAL
        LOGICAL     OVRLAP

        REAL        RCUTSQ, RMINSQ, SIGSQ, SR2, SR6, SR3, SR9
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0
!	!!!WRITE(*,*) IPULL,IPULLCHAIN
	RXI=RX(IPULL,IPULLCHAIN)
	RYI=RY(IPULL,IPULLCHAIN)
	RZI=RZ(IPULL,IPULLCHAIN)
!	!!!WRITE(*,*)RXI,RYI,RZI,MAT
	I=ANINT(RXI*MAT)
	J=ANINT(RYI*MAT)
	K=ANINT(RZI*MAT)
	DELTV1=UADS(I,J,K)
!	!!!WRITE(*,*)I,J,K, DELTV1

	DELTV=DELTV1

	!DELTV=0
        RETURN
        END
C------------------------------------------------------------------------

C------------------------------------------------------------------------
C   SUBRUTINA ADPOTIN
C   CALCULA LA ENERGIA DE INTERACCION ENTRE LA SUPERFICIE Y EL FLUIDO
C------------------------------------------------------------------------
      SUBROUTINE ADPOTIN2(SIGMETANO,RXI,RYI,RZI, SIGMA,EPS, 
     +RCUT, 
     +DELTV,DELTW)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WITH SURFACE ON ADDING AN ATOM.        **
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

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, RMIN, SIGMA, RXI, RYI, RZI, DELTV, DELTW
        INTEGER     N, NTRIAL
        LOGICAL     OVRLAP

        REAL        RCUTSQ, RMINSQ, SIGSQ, SR2, SR6, SR3, SR9
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************
        DELTV  = 0.0
        DELTW  = 0.0
!	!!!WRITE(*,*) IPULL,IPULLCHAIN

	I=ANINT(RXI*MAT)
	J=ANINT(RYI*MAT)
	K=ANINT(RZI*MAT)
	DELTV1=UADS2(I,J,K)
!	!!!WRITE(*,*)I,J,K, DELTV1


	DELTV=DELTV1

	!DELTV=0
        RETURN
        END
C------------------------------------------------------------------------
C-----------------------------------------------------------------------

C    ** SUBRUTINA ESTRUCTURA
C	** LEE LAS COORDENADAS DE LOS ATOMOS DE CARBONO
C	** Y REALIZA UNA TRANSFORMACION PARA NORMALIZARLAS
C-----------------------------------------------------------------------
	SUBROUTINE ESTRUCTURA(NAM,SIGMA,SIGMETANO)
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
	CHARACTER NAM*16
	!REAL RXA(9000),RYA(9000),RZA(9000)
	OPEN (49,FILE='PRUEBA.TXT')
	OPEN (51,FILE=NAM)
	!WRITE(*,*)NAM
	WRITE(*,*) MAT, ' MAT ENTRADA ES1'


	RXAMAX=SIGMETANO/(2*SIGMA)*ALX
	WRITE(*,*)RXAMAX,' RXMAX'
	RYAMAX=SIGMETANO/(2*SIGMA)*ALY
	RZAMAX=SIGMETANO/(2*SIGMA)*ALZ
      
	WRITE(*,*)RYAMAX,' RXMAX'
	READ (51,*) NC
	!WRITE(*,*) NC
      PAUSE
	IMAX=0
	DO I=1,NC
	READ (51,*)RXA,RYA,RZA,QACX
	IF(RXA.LT.RXAMAX) THEN
		IF(RXA.GT.-RXAMAX) THEN
			IF(RYA.LT.RYAMAX) THEN
				IF(RYA.GT.-RYAMAX) THEN
				IF(RZA.LT.RZAMAX) THEN
				IF(RZA.GT.-RZAMAX) THEN
				IMAX=IMAX+1
				RXC(IMAX)=RXA/SIGMETANO*SIGMA 
				RYC(IMAX)=RYA/SIGMETANO*SIGMA 
				RZC(IMAX)=RZA/SIGMETANO*SIGMA
                  QAC(IMAX)=REAL(QACX)*FCLEC*2
				WRITE(49,*) RXC(IMAX),RYC(IMAX),RZC(IMAX)
                  !WRITE(49,*) RXC(IMAX),RYC(IMAX),RZC(IMAX)
				ENDIF
				ENDIF
				ENDIF
			ENDIF
		ENDIF
	ENDIF
	

	ENDDO
	NC=IMAX
	WRITE(*,*)'ESTRUCTURA LEIDA ',NC
	CLOSE(51)
	CLOSE(49)
	!PAUSE
	RETURN
	END
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C    ** SUBRUTINA POTENCIAL
C    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZMA
C    ** UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000)
C------------------------------------------------------------------------
	SUBROUTINE POTENCIAL(EPS,SIGMA,SIGMETANO,RCUT,DAT,APOT)

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)





	!INTEGER NC
	LOGICAL FLAG
	PARAMETER ( PI = 3.14159265 )
	WRITE(*,*)NC
	SIGMA1=SIGMA*(3.4+SIGMETANO)/(2*SIGMETANO)
	FACTOR=SQRT(EPS*28)/EPS
      RCUTSQ = RCUT * RCUT
      SIGSQ  = SIGMA1 * SIGMA1
      SIGCUB = SIGSQ * SIGMA1
	RMIN=0.5*SIGMA1
	IF(APOT.EQ.0) RMIN=SIGMA1
	RMINSQ=RMIN*RMIN
C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA1 / RCUT ) ** 3
        SR9    = SR3 ** 3
        VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
        WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0
	IF (DAT.EQ.1) THEN
	GOTO 7000
	ENDIF
	OPEN(11,FILE='ENERGIA.DAT')
	MATX=MAT/2*ALX
	MATY=MAT/2*ALY
	MATZ=MAT/2*ALZ
	WRITE(*,*) MAT,ALX,ALY,ALZ,MATX,MATY,MATZ
	!!!!!!!!PAUSE

	DO I=-MATZ,MATZ
        RZI    = REAL(I)/REAL(MAT)

	DO IJ=-MATY,MATY
        RYI    = REAL(IJ)/REAL(MAT)

	DO K=-MATX,MATX
        RXI    = REAL(K)/REAL(MAT)

C    ** LOOP OVER ALL ATOMS  EXCEPT IPULL **
       DELTV=0.  
       WRITE(*,*)I,IJ,K
			DO J = 1, NC


C       ** PICK ACTIVE ATOMS FROM LOCATE **


              RXIJ  = RXI - RXC(J)
              RYIJ  = RYI - RYC(J)
              RZIJ  = RZI - RZC(J)

              RXIJ  = RXIJ - PBX*ALX*ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
			RZIJ  = RZIJ - PBZ*ALZ*ANINT ( RZIJ/ALZ )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

	IF ( RIJSQ .LT. RMINSQ) THEN

              FLAG(J) = .TRUE.

		VIJ=1E9
	

	ELSE
	
                 SR2   = SIGSQ / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ   = SR6 * ( SR6 - 1.0 )
	
	IF (APOT.EQ.0.) THEN
	VIJ=0 !SR2
	ENDIF
              
	VIJ2=0.
	ENDIF
                 DELTV = DELTV + VIJ
			   

                 


100     ENDDO
	DELTV =  4.0 * DELTV * FACTOR

	UADS(K,IJ,I)=DELTV
	UADS2(K,IJ,I)=DELTV
	WRITE(11,*) K,IJ,I,DELTV 
	ENDDO
	ENDDO
	WRITE(*,*) DELTV,K,IJ,I 
	ENDDO
	CLOSE(11)
	GOTO 8000
7000	OPEN(10,FILE='ENERGIA.DAT')
	DO I=-MAT/2,MAT/2
	DO IJ=-MAT/2,MAT/2
	DO K=-MAT/2,MAT/2
	READ(10,*)KA,IJA,IA,DELTV
      UADS(K,IJ,I)=DELTV
	ENDDO
	ENDDO
	ENDDO
	CLOSE(10)
8000	WRITE(*,*)'ENERGIA CALCULADA'
	RETURN
	END
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C    ** SUBRUTINA POTENCIAL2
C    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZ
C    ** USS(1000),FLAG(9000)
C------------------------------------------------------------------------
	SUBROUTINE POTENCIAL2(SIGMA,SIGMETANO,RCUT,APOT2,ALREF)
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

	!INTEGER NC
	LOGICAL FLAG
	SIGMA1=SIGMA
	  FACTOR=0.434958836
       RCUTSQ = RCUT * RCUT
      SIGSQ  = SIGMA1 * SIGMA1
      SIGCUB = SIGSQ * SIGMA1
	RMIN=0.5*SIGMA1
	RMINSQ=RMIN*RMIN
	RMAXRES=1.25*(ALREF/SIGMETANO*SIGMA)
	RMINRES=0.75*(ALREF/SIGMETANO*SIGMA)

        DELTV  = 0.0
        DELTW  = 0.0
	DO I=1,1000
        RZI    = REAL(I)/1000.

       RIJSQ = RZI*RZI
	IF(RZI.LT.RMINRES) THEN
	USRE(I)=1E6
	ELSEIF(RZI.GT.RMAXRES) THEN
 	USRE(I)=1E6
	ELSE
	USRE(I)=0
	ENDIF

	
      USSCE(I)=1/RZI

	IF ( RIJSQ .LT. RMINSQ) THEN

              FLAG(I) = .TRUE.
	VIJ=1E6
              
	
	ELSE
	IF(APOT2.EQ.1.) THEN
	
	           SR2   = SIGSQ / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ   = SR6 * ( SR6 - 1.0 )
                 WIJ   = SR6 * ( SR6 - 0.5 )
	ELSEIF (APOT2.EQ.0) THEN
	IF (RIJSQ.LT.SIGMA1**2) THEN
	VIJ=1E6
	ELSE
	VIJ=0
	ENDIF
	ELSEIF (APOT2.EQ.-1) THEN
	VIJ=SR6
	ENDIF


               
	ENDIF
	IF (RIJSQ.LT.SIGMA1**2) THEN
	VIJCI=1E6
	ELSE
	VIJCI=0
	ENDIF
	

	DELTV =  4.0 * VIJ
	USS(I)=DELTV
	USSCI(I)=VIJCI
	ENDDO
	RETURN
	END

C-----------------------------------------------------------------------
C----------------------------------------------------
      SUBROUTINE STRUCPOALY(NAMPOALY,SIGMA,SIGMETANO,N,NCHAIN,EPS)
	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax     

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

	CHARACTER NAM*16
	CHARACTER NAMPOALY*16

	OPEN (51,FILE=NAMPOALY)
	WRITE(*,*)NAMPOALY
	RXAMAX=SIGMETANO/(2*SIGMA)
	!WRITE(*,*)RXAMAX,' RXMAX'
	RYAMAX=SIGMETANO/(2*SIGMA)
	!WRITE(*,*)RYAMAX,' RXMAX'
      WRITE(*,*)MAT
      !PAUSE
      
	READ (51,*) N
	READ (51,*)NCHAIN
	WRITE(*,*) N,NCHAIN
	READ(51,*)NP2
C	!!!!!!!!PAUSE
	IMAX=0
	RXAMAX=SIGMETANO/(2*SIGMA)*ALX
	RYAMAX=SIGMETANO/(2*SIGMA)*ALY
	RZAMAX=SIGMETANO/(2*SIGMA)*ALZ

	DO ICHAIN=1,NCHAIN
	DO I=1,N
	READ (51,*)QNAT,RXCHAIN,RYACHAIN,RZACHAIN ,KAM(I,ICHAIN)
      KAM(I,ICHAIN)=0

	Q(I,ICHAIN)=QNAT*FCLEC
      
!	QT=QT+QNAT
	RX(I,ICHAIN)=RXCHAIN/SIGMETANO*SIGMA 
	RY(I,ICHAIN)=RYACHAIN/SIGMETANO*SIGMA 
	RZ(I,ICHAIN)=RZACHAIN/SIGMETANO*SIGMA
      
      RXNEW=RX(I,ICHAIN)
      RYNEW=RY(I,ICHAIN)
      RZNEW=RZ(I,ICHAIN)
      
     	IF(POLCIL.eq.1) THEN
      !Posiciones en el cilindro de las cargas imagenes
          
	  RADIO=SQRT(RXNEW**2+RYNEW**2)
        DELTARADIO=RCIL-RADIO
        RADIOEXTRA=RCIL+DELTARADIO
!----------------------------------------------------------------------
        XIM(I,ICHAIN)  = RADIOEXTRA*RXNEW/RADIO
        YIM(I,ICHAIN)  = RADIOEXTRA*RYNEW/RADIO
        ZIM(I,ICHAIN)  = RZNEW
        QIM(I,ICHAIN)  = Q(I,ICHAIN)
      ELSE
        XIM(I,ICHAIN)  = RXNEW
        YIM(I,ICHAIN)  = RYNEW
        ZIM(I,ICHAIN)  = -RZNEW
        QIM(I,ICHAIN)  = Q(I,ICHAIN)
       
        
      ENDIF


	ENDDO
	ENDDO
	!!!WRITE(*,*)'ESTRUCTURA LEIDA '
	CLOSE(51)
      WRITE(*,*)MAT
            IPULL=N
      IPULLCHAIN=1
              RXIJ  = RX(IPULL,IPULLCHAIN) - RX(IPULL-1,IPULLCHAIN)
              RYIJ  = RY(IPULL,IPULLCHAIN) - RY(IPULL-1,IPULLCHAIN)
              RZIJ  = RZ(IPULL,IPULLCHAIN) - RZ(IPULL-1,IPULLCHAIN)
!
              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

			RIJ=SQRT(RIJSQ)

	RESMAX=RIJ
      WRITE(*,*) RESMAX
      
      

      !PAUSE
	N=NP2
	RETURN
	END

C-------------------------------------------------------------------------
C-------------------------------------------------------------------------

C------------------------------------------------------------------------------

C-------------------------------------------------------------------------
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------


      SUBROUTINE ENERGYCHARGE (IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS,
     +RCUT, 
     +DELTV,
     + DELTW )

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA, DELTV, DELTW,KAPPA
        INTEGER     N, IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA
        SIGCUB = SIGSQ * SIGMA

C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
        WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0
        IF(Q(IPULL,IPULLCHAIN).EQ.0) RETURN
     
      NES=1
        RXES(NES)    = RX(IPULL,IPULLCHAIN)
        RYES(NES)    = RY(IPULL,IPULLCHAIN)
        RZES(NES)    = RZ(IPULL,IPULLCHAIN)
	  ZES(NES)	=Q(IPULL,IPULLCHAIN)
       Q(IPULL,IPULLCHAIN)=0
       
C    ** LOOP OVER ALL ATOMS  IN ALL THE  CHAINS EXCEPT IPULL **
        DO I=1,NC
               
            IF(QAC(I).NE.0) THEN
              NES=NES+1
              RXES(NES)  = RXC(I)
              RYES(NES)  = RYC(I)
              RZES(NES)  = RZC(I)

              ZES(NES)=QAC(I)
            ENDIF
        END DO
            
                
 	  DO JCHAIN=1,NCHAIN
        DO J = 1, N
             
C       ** PICK ACTIVE ATOMS FROM LOCATE **
          IF(Q(J,JCHAIN).NE.0) THEN

               NES=NES+1

              RXES(NES)  = RX(J,JCHAIN)
              RYES(NES)  = RY(J,JCHAIN)
              RZES(NES)  = RZ(J,JCHAIN)

              ZES(NES)=Q(J,JCHAIN)
              
              !CARGAS IMAGENES
              
              NES=NES+1
              RXES(NES)  = XIM(J,JCHAIN)
              RYES(NES)  = YIM(J,JCHAIN)
              RZES(NES)  = ZIM(J,JCHAIN)

              ZES(NES)=QIM(J,JCHAIN)

          ENDIF
      ENDDO
	ENDDO

      Q(IPULL,IPULLCHAIN)=ZES(1)
C-----------------------------------------------------------------------
C	ENERGÍA DE CREACIÓN CON LOS CONTRAIONES
C------------------------------------------------------------------------
        DO ION=1,4
         DO J = 1, NCI(ION)
        
            NES=NES+1
C       ** PICK ACTIVE ATOMS FROM LOCATE **
              JIN = LOCATE(J,ION)


              RXES(NES)  = RXCI(JIN,ION)
              RYES(NES)  = RYCI(JIN,ION)
              RZES(NES)  = RZCI(JIN,ION)

              ZES(NES)=REAL(NQCI(ION))*FCLEC
              NES=NES+1

              RXES(NES)  = RXCIM(JIN,ION)
              RYES(NES)  = RYCIM(JIN,ION)
              RZES(NES)  = RZCIM(JIN,ION)

              ZES(NES)=REAL(NQCI(ION))*FCLEC

         ENDDO
      ENDDO

      LTYPE=0
	CALL RWALD (VR,NES,LTYPE)

      
      DELTV=VR
      


C----------------------------------------------------------------------
        DELTV =   DELTV
!	!!!WRITE(*,*) DELTV,' 6'
!	!!!!!!!!PAUSE
        RETURN
        END
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------


      SUBROUTINE ENERGYCHARGECION(RXI,RYI,RZI,QCICON,N,IPULL,IPULLCHAIN,
     +NCHAIN,SIGMA,EPS, RCUT,NION, 
     +DELTV,
     + DELTW )


	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA, DELTV, DELTW,KAPPA
        INTEGER     N, IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

        
C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA
        SIGCUB = SIGSQ * SIGMA

C    ** CALCULATE LONG RANGE CORRECTIONS **

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0


	  QNAT	=REAL(NQCI(NION))*FCLEC
        NES=1
        RXES(NES)    = RXI
        RYES(NES)    = RYI
        RZES(NES)    = RZI
	  ZES(NES)	=QNAT
        
!------------------------------------------------------------------------------
!             Estructura SOLIDO
!------------------------------------------------------------------------------
        DO I=1,NC
            
           
            IF(QAC(I).NE.0) THEN
              NES=NES+1
              RXES(NES)  = RXC(I)
              RYES(NES)  = RYC(I)
              RZES(NES)  = RZC(I)

              ZES(NES)=QAC(I)
              !write(*,*)RXES(NES),RYES(NES),RZES(NES)
              !write(*,*)nes,zes(nes),'solido'

              
            ENDIF
        END DO
!---------------------------------------------------------------------------
C    ** LOOP OVER ALL ATOMS  IN ALL THE  CHAINS EXCEPT IPULL **
!---------------------------------------------------------------------------
	  DO JCHAIN=1,NCHAIN
        DO J = 1, N
              RXIJ  = RXI - RX(J,JCHAIN)
              RYIJ  = RYI - RY(J,JCHAIN)
              RZIJ  = RZI - RZ(J,JCHAIN)

              RXIJ  = RXIJ -PBX*ALX* ANINT ( RXIJ/ALZ )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ* ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
			RIJ=SQRT(RIJSQ)
			I=RIJ*1000+1
              IF(I.LE.1000) VIJ=USSCI(I)
	        IF(I.GT.1000) VIJ=0
	        
	


              DELTV = DELTV +VIJ

C       ** PICK ACTIVE ATOMS FROM LOCATE **
          IF(Q(J,JCHAIN).NE.0) THEN

           
               NES=NES+1

              RXES(NES)  = RX(J,JCHAIN)
              RYES(NES)  = RY(J,JCHAIN)
              RZES(NES)  = RZ(J,JCHAIN)

              ZES(NES)=Q(J,JCHAIN)
              !write(*,*)RXES(NES),RYES(NES),RZES(NES)
          !write(*,*)nes,zes(nes),'polimero'
           ENDIF
          
          IF(QIM(J,JCHAIN).NE.0) THEN

           
               NES=NES+1

              RXES(NES)  = XIM(J,JCHAIN)
              RYES(NES)  = YIM(J,JCHAIN)
              RZES(NES)  = ZIM(J,JCHAIN)

              ZES(NES)=QIM(J,JCHAIN)
              !write(*,*)RXES(NES),RYES(NES),RZES(NES)
          !write(*,*)nes,zes(nes),'pol im'
           ENDIF
          
	 

100     ENDDO
      ENDDO
C-----------------------------------------------------------------------
C	ENERGÍA DE CREACIÓN CON LOS CONTRAIONES
C------------------------------------------------------------------------
C-----------------------------------------------------------------------
C	ENERGÍA DE CREACIÓN CON LOS CONTRAIONES
C------------------------------------------------------------------------
      Do JION=1,4
      DO J = 1, NCI(JION)
            NES=NES+1
C       ** PICK ACTIVE ATOMS FROM LOCATE **
              JIN = LOCATE(J,JION)
              RXIJ  = RXI - RXCI(JIN,JION)
              RYIJ  = RYI - RYCI(JIN,JION)
              RZIJ  = RZI - RZCI(JIN,JION)

              RXIJ  = RXIJ - PBX*ALX*ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ - PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ - PBZ*ALZ*ANINT ( RZIJ/ALZ )
	
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
	
			RIJ=SQRT(RIJSQ)
			I=RIJ*1000+1
			VIJ=USSCI(I)
          	IF(RIJ.GT.1) THEN
	        VIJ=0
	        ENDIF


              DELTV = DELTV + VIJ
              


              RXES(NES)  = RXCI(JIN,JION)
              RYES(NES)  = RYCI(JIN,JION)
              RZES(NES)  = RZCI(JIN,JION)

              ZES(NES)=REAL(NQCI(JION))*FCLEC
              !write(*,*)RXES(NES),RYES(NES),RZES(NES)
           !write(*,*)nes,zes(nes),'iones'
              NES=NES+1
C       ** PICK ACTIVE ATOMS FROM LOCATE **

              RXES(NES)  = RXCIM(JIN,JION)
              RYES(NES)  = RYCIM(JIN,JION)
              RZES(NES)  = RZCIM(JIN,JION)

              ZES(NES)=REAL(NQCI(JION))*FCLEC
             ! write(*,*)RXES(NES),RYES(NES),RZES(NES)
              !write(*,*)nes,zes(nes),'iones im'
      ENDDO
      ENDDO
!----------------------------------------------------------------------------------------------------------------
      LTYPE=2
      !pause
      CALL RWALD (VR,NES,LTYPE)
    
      
      DELTV=VR+DELTV

C----------------------------------------------------------------------
        DELTV =   DELTV
C    ** ADD CHANGE IN LONG RANGE CORRECTION **
        DELTV =  DELTV
        RETURN
        END
C------------------------------------------------------------------------------





      SUBROUTINE ENERGYOUT(IPULL2,IPULL,IPULLCHAIN,N,NCHAIN,SIGMA,EPS, 
     +RCUT, 
     +DELTV,
     + DELTW )

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA, DELTV, DELTW
        INTEGER     N, IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0,KAPPA
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA
        SIGCUB = SIGSQ * SIGMA

C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
        WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0


	  QNAT	=Q(IPULL,IPULLCHAIN)
!	!WRITE(*,*) QNAT, ' QNAT', IPULL, IPULLCHAIN
!	!!!PAUSE
C----------------------------------------------------------------------------------------------------
      NES=1

        RXES(1)    = RX(IPULL,IPULLCHAIN)
        RYES(1)    = RY(IPULL,IPULLCHAIN)
        RZES(1)    = RZ(IPULL,IPULLCHAIN)

        ZES(1)=QNAT
      Q(IPULL,IPULLCHAIN)=0       
        DO I=1,NC
            IF(QAC(I).NE.0) THEN
              NES=NES+1
              RXES(NES)  = RXC(I)
              RYES(NES)  = RYC(I)
              RZES(NES)  = RZC(I)

              ZES(NES)=QAC(I)
            ENDIF
        END DO
C-------------------------------------------------------------
C    ** LOOP OVER ALL ATOMS  IN ALL THE  CHAINS  **
	  DO JCHAIN=1,NCHAIN
        DO J = 1, N
	!
C       ** PICK ACTIVE ATOMS FROM LOCATE **
      IF(Q(J,JCHAIN).NE.0) THEN
		NES=NES+1
              RXES(NES)  = RX(J,JCHAIN)
              RYES(NES)  = RY(J,JCHAIN)
              RZES(NES)  = RZ(J,JCHAIN)

              ZES(NES)=Q(J,JCHAIN)

              NES=NES+1
              RXES(NES)  = XIM(J,JCHAIN)
              RYES(NES)  = YIM(J,JCHAIN)
              RZES(NES)  = ZIM(J,JCHAIN)

              ZES(NES)=QIM(J,JCHAIN)

      ENDIF
      ENDDO
      ENDDO
      Q(IPULL,IPULLCHAIN)=QNAT
C    ** LOOP OVER ALL COUNTERIONS  EXCEPT IPULL **
      DO ION=1,4
        DO J = 1, NCI(ION)

C       ** PICK ACTIVE ATOMS FROM LOCATE **

           JIN = LOCATE(J,ION)

         !  IF ( JIN .NE. IPULL2 ) THEN
               NES=NES+1

              RXES(NES)  = RXCI(JIN,ION)
              RYES(NES)  = RYCI(JIN,ION)
              RZES(NES)  = RZCI(JIN,ION)

              ZES(NES)=FCLEC*REAL(NQCI(ION))

              NES=NES+1
C       ** PICK ACTIVE ATOMS FROM LOCATE **
              JIN = LOCATEIM(J,ION)
              


              RXES(NES)  = RXCIM(JIN,ION)
              RYES(NES)  = RYCIM(JIN,ION)
              RZES(NES)  = RZCIM(JIN,ION)

              ZES(NES)=REAL(NQCI(ION))*FCLEC
          ! ENDIF

        ENDDO
      ENDDO

 
      LTYPE=0
        CALL RWALD(VR,NES,LTYPE)

C------------------------------------------------------------------------

C--------------------------------------------------------------------------
        
      DELTV =   DELTV
      DELTV=VR+DELTV
C    ** ADD CHANGE IN LONG RANGE CORRECTION **

C    ** CHANGE SIGN OF DELTV AND DELTW FOR A REMOVAL **

        DELTV =-VR

        RETURN
        END
C------------------------------------------------------------------------------

C------------------------------------------------------------------------------
C    *******************************************************************
C    ** THIS FORTRAN CODE IS INTENDED TO ILLUSTRATE POINTS MADE IN    **
C    ** THE TEXT. TO OUR KNOWLEDGE IT WORKS CORRECTLY. HOWEVER IT IS  **
C    ** THE RESPONSIBILITY OF THE USER TO TEST IT, IF IT IS USED IN A **
C    ** RESEARCH APPLICATION.                                         **
C    *******************************************************************

C    *******************************************************************
C    ** FICHE F.22                                                    **
C    ** REAL-SPACE AND RECIPROCAL-SPACE PARTS OF EWALD SUM FOR IONS.  **
C    *******************************************************************

C    *******************************************************************
C    ** REAL-SPACE AND RECIPROCAL-SPACE PARTS OF EWALD SUM FOR IONS.  **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** WOODCOCK AND SINGER, TRANS. FARADAY SOC. 67, 12, 1971.        **
C    ** DE LEEUW ET AL., PROC. ROY. SOC. A 373, 27, 1980.             **
C    ** HEYES, J. CHEM. PHYS. 74, 1924, 1981.                         **
C    ** SEE ALSO FINCHAM, MDIONS, CCP5 PROGRAM LIBRARY.               **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE SETUP ( KAPPA )                                    **
C    **    SETS UP THE WAVEVECTORS FOR USE IN THE EWALD SUM           **
C    ** SUBROUTINE RWALD ( KAPPA, VR )                                **
C    **    CALCULATES THE R-SPACE PART OF THE SUM                     **
C    ** SUBROUTINE KWALD ( KAPPA, VK )                                **
C    **    CALCULATES THE K-SPACE PART OF THE SUM                     **
C    ** REAL FUNCTION ERFC ( X )                                      **
C    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER  TOTK         THE TOTAL NUMBER OF K-VECTORS STORED    **
C    ** INTEGER  MAXK         MAXIMUM POSSIBLE NUMBER OF K-VECTORS    **
C    ** INTEGER  KMAX         MAX INTEGER COMPONENT OF THE K-VECTOR   **
C    ** INTEGER  KSQMAX       MAX SQUARE MOD OF THE K-VECTOR REQUIRED **
C    ** REAL     VR           ENERGY FROM R-SPACE SUM                 **
C    ** REAL     VK           ENERGY FROM K-SPACE SUM                 **
C    ** REAL     KVEC(MAXK)   ARRAY USED TO STORE K-VECTORS           **
C    ** REAL     KAPPA        WIDTH OF CANCELLING DISTRIBUTION        **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SETUP IS CALLED ONCE AT THE BEGINNING OF THE SIMULATION       **
C    ** TO CALCULATE ALL THE K-VECTORS REQUIRED IN THE EWALD SUM.     **
C    ** THESE VECTORS ARE USED THROUGHOUT THE SIMULATION IN THE       **
C    ** SUBROUTINE KWALD TO CALCULATE THE K-SPACE CONTRIBUTION TO THE **
C    ** POTENTIAL ENERGY AT EACH CONFIGURATION. THE SELF TERM IS      **
C    ** SUBTRACTED FROM THE K-SPACE CONTRIBUTION IN KWALD.            **
C    ** THE SURFACE TERM FOR SIMULATIONS IN VACUUM IS NOT INCLUDED.   **
C    ** ROUTINE RWALD RETURNS THE R-SPACE CONTRIBUTION TO THE EWALD   **
C    ** SUM AND IS CALLED FOR EACH CONFIGURATION IN THE SIMULATION.   **
C    ** A CUBIC BOX AND UNIT BOX LENGTH ARE ASSUMED THROUGHOUT.       **
C    *******************************************************************






        SUBROUTINE RWALD ( VR,NES,LTYPE )

	COMMON/BLOCK1/RX(9000,900),RY(9000,900),RZ(9000,900),RXC(99000),
     +RYC(99000),MAT,Q(9000,900),QT,RCELE, IONES,XIM(9000,900),
     +RXES(9000),RYES(9000),RZES(9000),ZES(9000),KVEC,YIM(9000,900),
     +UADS2(-225:225,-225:225,-225:225),USSCI(1000),USSCE(1000),FCLEC,
     +RZC(99000),UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(9000),
     +RXCA(9000),RYCA(9000),RZCA(9000),PBX,PBY,PBZ,DIELEC,ALX,ALY,ALZ,
     +AKAPPA12,POLCIL,ZCIL,RCIL,USRE(1000),Z(6),KAM(9000,900),NK(2),
     +NQCI(4),NZI(2),NH,KA,KA2,DXY(-225:225,-225:225,10),QIM(9000,900),
     +PH(-40:40,-40:40,-40:40),QAC(99000),NC,DCI(-9000:9000,4),
     +PHPARCIAL(-40:40,-40:40,-40:40),NDISC,IPH,NCELLMAT,
     +NCI(4),RXCI(9000,4),RYCI(9000,4),RZCI(9000,4),ZIM(9000,900),
     +NCIM(4),RXCIM(9000,4),RYCIM(9000,4),RZCIM(9000,4),resmax      

 

     +/BLOCK2/LOCATE(9000,4),LOCATEIM(9000,4),
     +POHPARCIAL(-40:40,-40:40,-40:40),IPOH
     +/BLOCK3/DCI2(-9000:9000),
     +DATOSR(9000),DATOSN(9000,4),DATOSCIR1(9000),DATOSCIR2(9000)
        
!        COMMON / BLOCK1 / RXES, RYES, RZES, ZES

C    *******************************************************************
C    ** CALCULATES R-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                     NUMBER OF IONS                  **
C    ** REAL    RX(N),RY(N),RZ(N)     POSITIONS OF IONS               **
C    ** REAL    Z(N)                  IONIC CHARGES                   **
C    ** REAL    VR                    R-SPACE POTENTIAL ENERGY        **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION ERFC ( X )                                      **
C    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
C    *******************************************************************

        INTEGER     NES
!        REAL        RXES(5000), RYES(5000), RZES(5000), ZES(5000)
        REAL        KAPPA, VR

        REAL        RXI, RYI, RZI, ZI, RXIJ, RYIJ, RZIJ
        REAL        RIJSQ, RIJ, KRIJ, ERFC, VIJ

        INTEGER     I, J

C    *******************************************************************
      
        VR = 0.0
!        WRITE(*,*)'---------------'
!       WRITE(*,*) 'EWALD ',LTYPE
     

           RXI = RXES(1)
           RYI = RYES(1)
           RZI = RZES(1)
           ZI  = ZES(1)
 !         WRITE(*,*)RXI,RYI,RZI,ZI,NES
          
          !WRITE(*,*)RIJ,1,RXES(1),RYES(1),RZES(1),ZES(1)
           DO 99 J = 2, NES
              !!WRITE(*,*) ZES(J),NES
              
              RXIJ  = RXI - RXES(J)
              RYIJ  = RYI - RYES(J)
              RZIJ  = RZI - RZES(J)

              RXIJ  = RXIJ -PBX*ALX*ANINT ( RXIJ/ALX )
              RYIJ  = RYIJ -PBY*ALY*ANINT ( RYIJ/ALY )
              RZIJ  = RZIJ -PBZ*ALZ*ANINT ( RZIJ/ALZ )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
              RIJ   = SQRT ( RIJSQ )
              
            IF (RIJ.GT.RCELE) THEN
              VIJ=0
            ELSE
                       
            VIJ=ZI*ZES(J)*(1/RIJ-1/RCELE+(1/RCELE**2)*(RIJ-RCELE))
            ENDIF
 !             WRITE(*,*)RIJ,RXES(J),RYES(J),RZES(J),ZES(J),J
 !             WRITE(*,*)VIJ,' VIJ EWALD'
 !             !PAUSE
              VR    = VR + VIJ

99         CONTINUE

 !         WRITE(*,*) VR,' VR'
 !         WRITE(*,*)'-----------------------'
 !         PAUSE

        RETURN
        END






        REAL FUNCTION ERFC ( X )

C    *******************************************************************
C    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
C    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
C    *******************************************************************

        REAL        A1, A2, A3, A4, A5, P

        PARAMETER ( A1 = 0.254829592, A2 = -0.284496736 )
        PARAMETER ( A3 = 1.421413741, A4 = -1.453152027 )
        PARAMETER ( A5 = 1.061405429, P  =  0.3275911   )

        REAL        T, X, XSQ, TP

C    *******************************************************************

        T  = 1.0 / ( 1.0 + P * X )
        XSQ = X * X

        TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

        ERFC = TP * EXP ( -XSQ )

        RETURN
        END
