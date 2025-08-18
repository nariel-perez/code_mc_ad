C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
********************************************************************************
** FICHE F.14.  ALGORITHM TO HANDLE INDICES IN CONSTANT MU VT MONTE CARLO     **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** INDEX-HANDLING IN GRAND CANONICAL MONTE CARLO SIMULATION.     **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
C    **    ADDS AN ATOM TO THE ARRAY LOCATE.                          **
C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
C    **    REMOVES AN ATOM FROM THE ARRAY LOCATE.                     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF ATOMS BEFORE TRIAL      **
C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS           **
C    ** INTEGER IPULL               INDEX OF ATOM FOR REMOVAL         **
C    ** INTEGER NLOC                POSITION OF N IN LOCATE           **
C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING TRIAL      **
C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF AN ATOM  **
C    ** REAL    RX(NMAX), ETC.      POSITIONS OF ATOMS                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** ROUTINE ADD IS CALLED AFTER A SUCCESSFUL TRIAL ADDITION.      **
C    ** ROUTINE REMOVE IS CALLED AFTER A SUCCESSFUL TRIAL REMOVAL.    **
C    ** THE ARRAY LOCATE IS UPDATED IN EACH CASE.                     **
C    *******************************************************************


C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
        SUBROUTINE ADD ( MOLKIND)
        IMPLICIT NONE
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
C    *******************************************************************
C    ** SUBROUTINE TO ADD AN ATOM TO THE ARRAY LOCATE.                **
C    **                                                               **
C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE NEW ADDITION   **
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
      real acel, acelx,acely,acelz
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
        

C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )




        REAL        RXNEW, RYNEW, RZNEW
        INTEGER     IPULL,I

        INTEGER     INEW, NTRIAL
      INTEGER POSXIN,POSZIN,POSYIN, j, jin
      real rxi,ryi,rzi
C    *******************************************************************
      !write(*,*) '----'
      !write(*,*) 'ADD'
      !write(*,*)molkind, 'molkind'
      !write(*,*)N(MOLKIND),' nmol'
      posxin=int(rx1(1)*10000)
      posyin=int(ry1(1)*1000)
      poszin=int(rZ1(1)*100)
      

        NTRIAL = N(MOLKIND) + 1
       ! write(*,*) ntrial, n(molkind), ' NTRIAL'
        INEW = LOCATE(NTRIAL,MOLKIND)
      !write(*,*) ntrial, n(molkind)
      !WRITE(*,*) INEW, ' INEW'
        IF ( INEW .EQ. 0 ) THEN

C       ** ATOM REQUIRES A NEW NUMBER **

           LOCATE(NTRIAL,MOLKIND) = NTRIAL
           INEW           = NTRIAL
       !    write(*,*) ntrial,  inew,n(molkind), 'ntrial, inew, n'
        ENDIF
!      WRITE(*,*) 'ADD'
C    ** FIT NEW ATOM INTO THE ARRAY **
        DO I=1,NATOM(MOLKIND)  
        !    write(*,*) ntrial, natom(molkind),n(molkind), 'natom'
        RX(INEW,I,MOLKIND) = RX1(I)
        RY(INEW,I,MOLKIND) = RY1(I)
        RZ(INEW,I,MOLKIND) = RZ1(I)
!        WRITE(*,*)RX1(I),RY1(I),RZ1(I),I,INEW
        ENDDO
        IF(NATOM(MOLKIND).GT.1) THEN
            ANX(INEW,MOLKIND)=EXNEW
            ANGY(INEW,MOLKIND)=EYNEW
            ANZ(INEW,MOLKIND)=EZNEW
        ENDIF
        
!            write(*,*)'---------------------'
        RETURN
        END
C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------
