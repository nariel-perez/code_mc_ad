********************************************************************************
** FICHE F.13.  THE HEART OF A CONSTANT MU VT MONTE CARLO PROGRAM             **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
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
      Program main
      implicit none
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS,acel,acelx,acely,acelz ,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE
     +	  /BLOCK3/nmaxi,nmin

      REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW
      REAL UADS(-100:100,-100:100,-100:100,50)
      REAL USS(5000,50,50)
      LOGICAL FLAG(1000)
      INTEGER LOCATE(5000,10)
      REAL RXC(9000),RYC(9000),RZC(9000)
      REAL EPSAC(9000),SGC(9000),QAC(9000)
      REAL RX1(50),RY1(50),RZ1(50)
      REAL EPSI(50),SIGM(50),Q(50)
	real P,dp, p_ratio
      REAL X(10),XT,AX,p_vals(31)
      INTEGER NSYM(50,10)
      INTEGER N(10)              
      INTEGER MAT
      INTEGER NTRIAL              
      INTEGER NMAX     
      INTEGER NTOTAL(10),NTOTALP,NTOTALB
      REAL    RXNEW,RYNEW,RZNEW   
      REAL    V,VA,VG,u2                   
      REAL    W                   
      REAL    DELTV,DATOS(5000)               
      REAL    DELTW               
      REAL    TEMP                
      REAL    Z(10)                   
      REAL    SIGMA               
      REAL    RCUT                
      REAL    RMIN
      REAL ANPROM(10)
      LOGICAL OVRLAP              
      LOGICAL CREATE              
      LOGICAL GHOST
	CHARACTER NAM*16
      CHARACTER MOLEC1*16
	CHARACTER nmbr*2
	CHARACTER CONFIG*3
      CHARACTER CONFAT*3
      CHARACTER CONFNAT*3
      real sigmetano
      real eps
	real auvol
      real acel, acelx,acely,acelz
      real bcx,bcy,bcz
      real T
      INTEGER ISOT, IJPASOS, IKPASOS, IPASOS,JPASOS,KPASOS
      INTEGER MULT2, MULT
      REAL AK, PRED
      REAL VOL, XMAX,YMAX, ZMAX
      INTEGER NC, NMOLEC
      INTEGER I,J
      INTEGER NMATOM
      REAL X1,Y1,Z1,RX0(50,10),RY0(50,10),RZ0(50,10)
      INTEGER IKIND,NS
      INTEGER NATOM(10),NMIN(1000),NMAXI(1000)
      INTEGER NATOMKIND (50,10)
      INTEGER INMOLEC
      REAL CR
      INTEGER IJ
      REAL DUMMY
      REAL RANF
      REAL U,UG,UA,UN,UNG,UNA,AN, ANN
      REAL N2
      REAL AN1,U1,UNG1,UG1,UNA1,UA1,UN1,AN2
      REAL CALOR,CALORG,CALORA
      REAL ESCALA
      REAL RXAI,RYAI,RZAI, EPSAI,SGCI,QACI
      INTEGER SYMBOL
      INTEGER K,jin
      REAL RXN,RYN,RZN
      REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)
      REAL  starttime, endtime, elapsedtime

      real aitest76,aitest77
      real diel
      integer ensemble
      integer nmolec2
      integer natom2
      integer molkind
      integer ncantmol

	!real deltv
	real deltva
	integer ipull
      
      REAL CNF(-25:25,-25:25,-25:25,50,10)
      INTEGER ICNF,JCNF,KCNF,NCELLMAT
      INTEGER NATOMKINDI
      REAL CALORESP1,CALORESP2,CALORESP3
      INTEGER NCONFMIN,NCONFMAX
      INTEGER NESTADO
	INTEGER ntotalGRAF
	real escalax,escalay,escalaz
	integer canonicalmolecules
	integer ensemble2

      
 
  ! Record the start time
        call CPU_TIME(starttime)     
      
!-------------------------------------------------------------------------------------
! Lectura de archivos de entrada
!--------------------------------------------------------------------------------------
	open(10,file='input.txt')
	read(10,*)P !PRESI? INICIAL DE LA SIMULACI? en cm de Hg
      write(*,*) P
	read(10,*)dp !CAMBIO EN LA PRESI? PARA LA ISOTERMA
	write(*,*) dp
      read(10,*)sigmetano ! ESTO SIRVE COMO REFERENCIA PARA LAS DEM?S UNIDADES DE LONGITUD
	write(*,*) sigmetano
      read(10,*)eps !Referencia para la energ?
      write(*,*)eps
     	read(10,*)ACEL, acelx,acely,acelz !TAMA? DE LA CELDA
      write(*,*)ACEL, acelx,acely,acelz
	!pause
      read (10,*) diel
      write(*,*) diel
      READ(10,*)BCX,BCY,BCZ !CONDICIONES DE CONTORNO
      write(*,*)BCX,BCY,BCZ
	read(10,*)T !Temperatura de la simulaci?
      write(*,*)T
	read(10,*)mat
      write(*,*)mat
	read(10,*)nam
      write(*,*) nam
	read(10,*)isot
      write(*,*)isot
	read(10,*)ijpasos
      write(*,*)ijpasos
	read(10,*)ikpasos
      write(*,*)ikpasos
	read(10,*)mult2
      write(*,*) mult2
      read(10,*) ensemble
      if(ensemble.eq.0) then 
          write(*,*) 'CANONICAL'
      ELSE IF (ENSEMBLE.EQ.1) THEN 
          WRITE(*,*) 'GRAND CANONICAL WITH INITIAL CONFIGURATION'
      ELSE   IF (ENSEMBLE.EQ.2)  then
          WRITE(*,*) 'GRAND CANONICAL'
          ELSE
          write(*,*) ' GRAND CANONICAL NAMD CONFIGURATION'
      ENDIF
      
      
      
      if(ensemble.eq.0) then 
          dp=0
      endif
      READ(10,*) NCELLMAT
      IF(NCELLMAT.GT.50) THEN
          WRITE(*,*) 'EL TAMA? DE LA MATRIZ ESTADISTICA ES MUY GRANDE'
          STOP
      ENDIF
      READ (10,*) NESTADO
	read (10,*)canonicalmolecules
	close(10)
      write(*,*)mult2
      
      
      
      if(ensemble.eq.3) then
      call NAMD1
      ensemble2=ensemble
      ensemble=1
      endif


!-----------------------------
      p_ratio = (dp/ p)**(1.0d0 / real(isot - 1, kind=8))
      
      p_vals = 0
      p_vals(1) = p
      do i = 2, isot
         p_vals(i) = p_vals(i-1) *p_ratio
      end do
      



       
!----------------------------------------------------------
!     Transformaci? a unidades reducidas
!----------------------------------------------------------
      SIGMA=SIGMETANO/ACEL !TAMA? DE CELDA REDUCIDO
	AK=8.31/6.023E23
	TEMP=T/EPS
	!P=P*1333.22
	PRED=P*SIGMA**3/EPS
	RCUT=10*SIGMA
      
!-----------------------------------------------------------
!     Volumen
!-----------------------------------------------------------
	

      XMAX=ACELX/ACEL
	YMAX=ACELY/ACEL
	ZMAX=ACELZ/ACEL
      
      VOL=XMAX*YMAX*ZMAX
      
      WRITE(*,*) TEMP,'TEMP'
	WRITE(*,*) PRED,SIGMA,P,' PRES'
	write(*,*) VOL, ' VOL'
	!pause
!-------------------------------------------------------------------
	OPEN(50,FILE='SALIDAACTIVADO-100.TXT')
	OPEN(97,FILE='PERFILES.TXT')
!-----------------------------------------------------------------      
!     LLAMADA A SUBRUTINAS DE POTENCIALES
!-----------------------------------------------------------------
      call estructura(eps,nam,sigma,sigmetano,NC,diel)
	call POTENCIALFF(EPS,sigma,sigmetano,NC,RCUT,diel)
	call POTENCIAL(EPS,sigma,sigmetano,NC,RCUT,diel)
!-----------------------------------------------------------------
!     LLAMADA A SUBRUTINAS DE LECTURAS DE ESTRUCTURAS MOLECULARES
C-----------------------------------------------------------------
!
      OPEN(90,FILE='MOLEC.DAT')
      READ(90,*) NMOLEC
      XT=0.
      DO I=1, NMOLEC
	NTOTAL(I)=0
          READ(90,*)MOLEC1
          OPEN(92,FILE=MOLEC1)
          READ(92,*)NMATOM,AX,NCONFMIN,NCONFMAX
		NMIN(I)=NCONFMIN
		NMAXI(I)=NCONFMAX
		WRITE(*,*) I,NMIN(I),NMAXI(I)
!		PAUSE
          NATOM(I)=NMATOM
          X(I)=AX
          XT=XT+X(I)
          WRITE(*,*) XT, ' XT'
          
          IF(NMATOM.GT.50) STOP
          DO J=1,NMATOM
              READ(92,*)X1,Y1,Z1,IKIND,NS
              RX0(J,I)=X1/SIGMETANO*SIGMA 
              RY0(J,I)=Y1/SIGMETANO*SIGMA 
              RZ0(J,I)=Z1/SIGMETANO*SIGMA 
              NATOMKIND(J,I)=IKIND
              NSYM(J,I)=NS
              write(*,*) j,i
              
          ENDDO
          CLOSE(92)
      ENDDO
      !PAUSE
      IF (XT.NE.1) THEN
          WRITE(*,*) 'CUIDADO!!! LA SUMA DE LAS FRACCIONES NO ES 1'
          !PAUSE
      ENDIF

      !STOP
!-------------------------------------------------------------
!     RESET DE LOS VALORES
      
      if(ensemble.eq.2) then
      DO I=1,NMOLEC
          N(I)=0
          ntotal(i)=0
          NTOTALP=0
          DO J=1,5000
              LOCATE(J,I)=0
          ENDDO
          
      ENDDO
      
      V=0.
      VG=0.
      VA=0.
      endif
      
!-------------------------------------------------------------
!     COMIENZA LA SIMULACI?
!-------------------------------------------------------------
!*************************************************************
!     Canonical ensemble
!
      if (ensemble.lt.2) then
          open(unit=67,file='initconf.txt')
		write(*,*) 'archivo abierto'
          read(67,*)V,VG,VA
          write(*,*) V,VG,VA
          read(67,*)nmolec2
          write(*,*) nmolec2
!	pause
          if(nmolec2.ne.nmolec) then
              write(*,*)'El n?ero de moleculas no es el mismo'
              write(*,*) nmolec2,nmolec
      !        pause
              stop
          endif
          do i=1,nmolec2
              molkind=i
              read(67,*)natom2
              write(*,*)natom2
		!pause
              if(natom2.ne.natom(i)) then
              write(*,*)'El n?ero de ?omos no es el mismo'
              write(*,*) natom2,natom(i),i
              !pause
              stop
          endif
          read(67,*)ncantmol
          write(*,*) ncantmol,i
          do j=1,ncantmol
              n(i)=j-1
              DO k=1,natom2
              READ(67,*)X1,Y1,Z1
              write(*,*)x1,y1,z1,j,k
              RX1(k)=X1 
              RY1(k)=Y1 
              RZ1(k)=Z1 
              enddo
              call add(molkind)
              n(i)=j
          enddo
          
      ENDDO
      close(67)
!---------------------------------------------------------------------------------------------------
!	CALCULO ENERGIA
      !V=0.
      !VG=0.
      !VA=0.

	DO I=1,NMOLEC
	DO J=1,N(I)
	IPULL = LOCATE(J,I)
	CALL ADPOTOUT(IPULL, I,DELTVA )
	CALL POTOUT ( IPULL, I,  DELTV)
	!V=V+DELTVA+DELTV
	!VG=VG+DELTV
	!VA=VA+DELTVA
	ENDDO
	ENDDO
	!V=-V/2
	!VG=-VG/2
	!VA=-VA/2
	write(*,*) VA,VG,V,' Energias Iniciales'
	ANPROM(I)=0
	ntotal(i)=0
      endif
 !--------------------------------------------------------------------------------------------------     
              
              
      
	DO IPASOS=1,isot
      !-------------------------------------------
      !RESET VALORES ESTADISTICAS ESPACIALES
          Do I=1,nmolec
              DO NATOMKINDI=1,NATOM(I)
                  DO ICNF=-NCELLMAT/2,NCELLMAT/2
                      DO JCNF=-NCELLMAT/2,NCELLMAT/2
                          DO KCNF=-NCELLMAT/2,NCELLMAT/2
                              CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)=0
                          ENDDO
                      ENDDO
                  enddo
              ENDDO
          ENDDO
          
                  
               !---------------------------------------
          
	CONFIG='CONFIG'
	write(CONFIG,'(i3)') ipasos


	open(40,file='CONFIG'//CONFIG//'.xyz')
	open(41,file='CONFIG'//CONFIG//'.TXT')

* conversion of pressure into activity zp: 1 cm Hg=1333.22 Pa; NA=6.0220E+23;
* R=8.3144 J/(mol*K); z in atom/A**3
!	write(*,*) aP ,dp

      IF (NESTADO.EQ.1) THEN
      DO INMOLEC=1,NMOLEC
	auvol=((sigmetano/SIGMA)**3)*VOL
      Z(INMOLEC)=X(INMOLEC)*P/(8.3144*T)*6.023E-7*auvol
          WRITE(*,*)Z(INMOLEC), 'Z ',INMOLEC

      ENDDO
      ELSE
      DO INMOLEC=1,NMOLEC
      Z(INMOLEC)=X(INMOLEC)*P*6.023E-4*((ACEL)**3)*VOL
	Z(nmolec)=55.55*6.023E-4*((ACEL)**3)*vol
          WRITE(*,*)Z(INMOLEC), 'Z ',INMOLEC


      ENDDO
      ENDIF
	!write(*,*) 'estoy aca'
!---------------------------------------------------------------------------------
!     VALORES EN CERO PARA CALORES ISOST?ICOS
!---------------------------------------------------------------------------------
	U=0
	UG=0
	UA=0
      UN=0
	UNG=0
	UNA=0
	AN=0
	N2=0     
      
!---------------------------------------------------------------------------------
	DO JPASOS=1,ijpasos
          	AITEST76=REAL(JPASOS)/500
	WRITE(58,*) AITEST76
	AITEST77=AITEST76-INT(AITEST76)
		MULT=1
	IF(AITEST77.EQ.0) THEN
	WRITE(*,*)JPASOS,N(1:nmolec)
 
      
	!N=0.
 	!NGAS2=0
	!V=0.
	!MULT=MULT2
      ENDIF				  

	open(51,file='Ener'//CONFIG//'.TXT')
	V=V/EPS
	VG=VG/EPS
	VA=Va/Eps
	W=0.
	cr=0.
	mult=1
	if (jpasos.eq.1) then
	mult=mult2 
	endif
  	DO KPASOS=1,ikpasos*mult
C------------------------------------------------------------------
C	ELECCION DEL PASO
C------------------------------------------------------------------
	ij=ranf(dummy)*3+1
	!write(*,*) 'estoy aca',ij
      if(ensemble.eq.0) goto 30
      !if(jpasos.ge.111)  then 
      !      write(*,*) ij, ' N ',N(1), '--', kpasos
      !      if(kpasos.ge.1781) then
      !pause      
      !endif
      !endif
      !write(*,*) V, ' V ***************************'
      !write(*,*) VG, ' VG ***************************'
      !write(*,*) VA, ' VA ***************************'
      !write(*,*)VG+VA,' suma'
      !write(*,*)'*******************'
      !pause
!      write(*,*) nmolec, 'nmolec main'
!      write(*,*) ij,jpasos,n(1:3), ' MAIN'
!	write(*,*) ntotal(1:nmolec)
	
!	write(*,*) ij,jpasos,n(1), 'elec paso'
	goto (10,20,30,35) ij
      
  10  CALL IN(TEMP,Z,SIGMA,EPS,RCUT,V,va,vg,W,CREATE,cr,jpasos,
     +canonicalmolecules)
      !write(*,*) V, ' V ***************************'
      !write(*,*) VG, ' VG ***************************'
      !write(*,*) VA, ' VA ***************************'
      !write(*,*)VG+VA,' suma'
      !pause
	GOTO 40
  20  CALL OUT( TEMP, Z,SIGMA,EPS, RCUT,V,va,vg,W,GHOST,jpasos,
     +NMIN,NMAXI,canonicalmolecules )
      !write(*,*) V, ' V ***************************'
      !write(*,*) VG, ' VG ***************************'
      !write(*,*) VA, ' VA ***************************'
      !write(*,*)VG+VA,' suma'
      !pause
      GOTO 40
  30  CALL MOVE( TEMP, Z,SIGMA,EPS, RCUT,V,va,vg,W,GHOST,jpasos )
      !write(*,*) V, ' V ***************************'
      !write(*,*) VG, ' VG ***************************'
      !write(*,*) VA, ' VA ***************************'
      !write(*,*)VG+VA,' suma'
      !pause
      goto 40
  35  CALL change( TEMP, Z,SIGMA,EPS, RCUT,V,va,vg,W,GHOST,jpasos )

  40  	ENDDO
       !     write(*,*) ij, ' N ',N(1), 'eleccion del paso', Jpasos
      !pause
      
    !------------------------------------------------------------------
!     REESCALA ENERG?S
!------------------------------------------------------------------
	V=V*EPS
	VG=VG*eps
	va=va*eps
!-------------------------------------------------------------------
!     ESTADISTICAS ESPACIALES
!-------------------------------------------------------------------
!      write(*,*) ntotal(1:nmolec), 'NTOTALLLLLLL'
!      pause
      DO I=1,NMOLEC
          NTOTAL(I)=NTOTAL(I)+N(I)
          NTOTALP=NTOTALP+N(I)
      ENDDO
      !WRITE(*,*) 'LLAMO A ESTADISTICA'
      !PAUSE
      CALL ESTADISTICA(CNF,NCELLMAT)
      !WRITE(*,*) ' VUELVO DE ESTADISTICA'
      !PAUSE
      !write(*,*) ntotal(1), 'B'
	U=U+V
	UG=UG+VG
	UA=UA+VA
!---------------------------------------------------------------------
!     ESTADISTICAS VALIDAS PARA UN SOLO TIPO DE MOLECULAS
!---------------------------------------------------------------------
      NTOTALB=0
      DO I=1,NMOLEC
          NTOTALB=NTOTALB+N(I)
      ENDDO
      
      UN=UN+V*NTOTALB
	UNG=UNG+VG*NTOTALB
	UNA=UNA+VA*NTOTALB
      U2=U2+V**2
	AN=AN+NTOTALB
	N2=N2+NTOTALB**2
!---------------------------------------------------------------------      
      ENDDO
            
      
      !write(*,*) ntotal(1), 'C'
      !pause
	close(51)

      write(*,*)'-----'
      !*************************************************************
!     Canonical ensemble
!
      
          open(unit=67,file='initconf.txt')
          write(67,*)V,VG,VA,' energias'
          write(67,*)nmolec,' tipos de molec'
          do i=1,nmolec
              molkind=i
              write(67,*)natom(i)
          write(67,*)n(i), ' natom ', i
          do j=1,n(i)
              jin=locate(j,i)
              DO k=1,NATOM(i)
                  x1=rx(jin,k,i)
                  y1=ry(jin,k,i)
                  z1=rz(jin,k,i)
              write(67,*)X1,Y1,Z1, 'x y z molec'
              enddo
          enddo
          enddo
      
       ! write(*,*) ntotal(1), 'D'
      !pause
        
          
      
      close(67)
      
 !--------------------------------------------------------------------------------------------------     


!----------------------------------------------------------------------
!     PROMEDIOS PARA MOLECULAS
!----------------------------------------------------------------------
 !              write(*,*) ntotal(1:nmolec), 'E'
 !     pause
        
      DO I=1,NMOLEC
        
          ANPROM(I)=REAL(NTOTAL(I))/REAL(JPASOS)

      ENDDO
      
      write(*,*)ANPROM(1:nmolec)
	!pause
!----------------------------------------------------------------------
!      PROMEDIOS PARA CALCULAR CALORES ISOSTERICOS
!-----------------------------------------------------------------------
	AN1=AN/REAL(JPASOS-1)
 	U1=U/REAL(JPASOS-1)
      WRITE(*,*) U1, ' U1'
	UNG1=UNG/REAL(JPASOS-1)
      WRITE(*,*) UNG1,' UNG1'
	UG1=UG/REAL(JPASOS-1)
      WRITE(*,*)UG1, ' UG1'
 	UNA1=UNA/REAL(JPASOS-1)
      WRITE(*,*)UNA1, ' UNA1'
	UA1=UA/REAL(JPASOS-1)
      WRITE(*,*) UA1,'UA1'
	UN1=UN/REAL(JPASOS-1)
      WRITE(*,*)UN1,' UN1'
      U2=U2/REAL(JPASOS-1)
      

	AN2=real(N2)/REAL(JPASOS-1)
      WRITE(*,*)AN2, ' AN2'
	ANN=an2-an1**2
      WRITE(*,*) ANN, ' ANN'
	escala=ACEL
	
	!if(an2.eq.0) an2=1e20
	CALOR=8.3144*T-((UN1-U1*AN1)/(ANn))*8.31
	CALORG=(UNG1-UG1*AN1)/(ann)*8.31
	CALORA=-((UNA1-UA1*AN1)/(ann))*8.31
      
      Caloresp1=(U2-U1**2)*8.31**2-(8.31*U1*AN1)**2/(ANN)
      CALORESP2=(UN1-AN1*8.31*T**2)
      
      CALORESP3=CALORESP1/CALORESP2
      
!---------------------------------------------------------------------
	OPEN (49,FILE='TRUNCADO.TXT')
	read(49,*) NC
	NTOTALGRAF=NC
	DO I=1,NMOLEC
	ntotalGRAF=NTOTALGRAF+n(i)*natom(i)
	ENDDO
	write(40,*) NTOTALGRAF
	WRITE(40,*) ' ' 
	write(*,*) NTOTALGRAF
	WRITE(*,*) ' ' 
	ESCALAx=acelx
	escalay=acely
	escalaz=acelz

	write(*,*)escalax,escalay,escalaz
	!pause
	DO I=1,NC
	READ (49,*)RXAI,RYAI,RZAI,EPSAI,SGCI,QACI,SYMBOL
	write(40,*) SYMBOL,RXAI,RYAI,RZAI
      NTOTALGRAF=NC

      ENDDO
      close(49)

      DO I=1,NMOLEC
          do j=1,N(I)
              jin=locate(j,i)
              DO K=1,NATOM(I)
                  RXN=RX(Jin,K,I)*ESCALA
                  RYN=RY(Jin,K,I)*ESCALA
                  RZN=RZ(Jin,K,I)*ESCALA
                  write(40,*)NSYM(K,I),RXN,RYN,RZN
              enddo
          ENDDO
      ENDDO

	!write(*,*) escala,NC,N,NMOLEC,RX,RY,RZ
	do i=1,nmolec
	write(*,*)i,n(i),natom(i), ' i natom main'
	enddo
	if(ensemble2.eq.3) then
      !call namd2(escala,NC,N,NMOLEC,RX,RY,RZ,natom,locate)
	endif
!-------------------------------------------------------------------------------------










!--------------------------------------------------------------------------------------



      
      
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i2)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
	write(CONFNAT,'(i2)') NATOMKINDI
      open(101,file='CNF'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
      WRITE(101,*)NCELLMAT
                  DO ICNF=-NCELLMAT/2,NCELLMAT/2
                      DO JCNF=-NCELLMAT/2,NCELLMAT/2
                          DO KCNF=-NCELLMAT/2,NCELLMAT/2
          write(101,*)icnf,jcnf,kcnf,CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
                          ENDDO
                      ENDDO
                      
                  enddo
                  CLOSE(101)
              ENDDO
          ENDDO
          
          
          
          
          
      
      
       
      DO I=1,NMOLEC
          ntotal(I)=0
      ENDDO
      

	WRITE(*,*) AN1,P,CALOR,' CALOR'
!88	format(E10,1x,F10.3,1x,f10.3,1x,f10.3,1x,f10.3,F10.3,1X,F10.3)
	WRITE(50,*)P,ANPROM(1:NMOLEC) ,CALOR,calora,CALORESP3
      !WRITE(50,*) P,ANPROM(1),CALOR,8.31*T-calorg,calora
      
	!WRITE(*,*) P,ANPROM(1),ANPROM(2),ANPROM(3),CALOR,calora      
	WRITE(*,*) P,ANPROM(1:nmolec),CALOR,8.31*T-calorg,calora
	CLOSE(40)
		CLOSE(41)
	CLOSE(21)
	CLOSE(22)
	AN=0
	an1=0
	!p=P+DP                  !*1333.22
        p = p_vals(ipasos +1)
	ENDDO
	CLOSE(50)
	CLOSE(97)



  ! Record the end time
        call CPU_TIME(endtime)
  
  ! Calculate the elapsed time
      
      elapsedtime= endtime -starttime
  
  ! Output the execution time
      write(*,*) 'Total execution time (seconds): ', elapsedtime



      END
C--------------------------------------------------------------------
