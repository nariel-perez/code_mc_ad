      SUBROUTINE ESTADISTICA(CNF,NCELLMAT)
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
	real P,dp
      REAL X(10),XT,AX
      INTEGER NSYM(50,10)
      INTEGER N(10)              
      INTEGER MAT
      INTEGER NTRIAL              
      INTEGER NMAX     
      INTEGER NTOTAL(10)
      REAL    RXNEW,RYNEW,RZNEW   
      REAL    V,VA,VG                   
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
      real sigmetano
      real eps
      real bcx,bcy,bcz,acel,acelx,acely,acelz
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
      INTEGER NATOM(10)
      INTEGER NATOMKIND (50,10)
      INTEGER INMOLEC
      REAL CR
      INTEGER IJ
      REAL DUMMY
      REAL RANF
      REAL U,UG,UA,UN,UNG,UNA,AN, ANN
      INTEGER N2
      REAL AN1,U1,UNG1,UG1,UNA1,UA1,UN1,AN2
      REAL CALOR,CALORG,CALORA
      REAL ESCALA
      REAL RXAI,RYAI,RZAI, EPSAI,SGCI,QACI
      INTEGER SYMBOL
      INTEGER K,jin
      REAL RXN,RYN,RZN
      REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)
      INTEGER ICNF,JCNF,KCNF,NCELLMAT
      real aitest76,aitest77
      real diel
      integer ensemble
      integer nmolec2
      integer natom2
      integer molkind
      integer ncantmol
            REAL CNF(-25:25,-25:25,-25:25,50,10)

      INTEGER NATOMKINDI, IMOL
      
      

      !CNF(i,j,k,ni)=0
      Do I=1,nmolec
         ! WRITE(*,*) I,NMOLEC, 'a'
          !pause
	DO imol=1,N(I)
          !WRITE(*,*) IMOL,N(I),I,' b'
          !pause
          JIN   = LOCATE(imol,I)
          !WRITE(*,*) jin,IMOL,N(I),I,' bprima'
          DO NATOMKINDI=1,NATOM(I)
            !  WRITE(*,*) NATOMKINDi,NATOM(I),I, ' c'
              !pause
          
	ICNF=INT(RX(jin,NATOMKINDI,I)*NCELLMAT)
      JCNF=INT(RY(jin,NATOMKINDI,I)*NCELLMAT)
      KCNF=INT(RZ(jin,NATOMKINDI,I)*NCELLMAT)
      
      !write(*,*) icnf,jcnf,kcnf, 'i,j,k'
      !pause

	CNF(icnf,jcnf,KCNF,I,NATOMKINDI)=
     +CNF(icnf,jcnf,KCNF,I,NATOMKINDI)+1
      ENDDO
	ENDDO
      enddo
      RETURN
      END
      
