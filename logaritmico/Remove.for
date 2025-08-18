C------------------------------------------------------------------------------
C-----------------------------------------------------------------------------

        SUBROUTINE REMOVE (NLOC,IPULL,MOLKIND)
      implicit none
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS,    acel,acelx,acely,acelz,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE
C    *******************************************************************
C    ** SUBROUTINE TO REMOVE AN ATOM FROM THE ARRAY LOCATE.           **
C    **                                                               **
C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE REMOVAL.       **
C    ** ELEMENT IPULL OF LOCATE IS TO BE DESTROYED.                   **
C    *******************************************************************
      REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        


        INTEGER     IPULL, NLOC

        INTEGER     K
!----------------------------------------------------
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
      REAL BCX,BCY,BCZ  ,acel,acelx,acely,acelz
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
      
      !WRITE(*,*) ' REMOVE'
        
      !write(*,*)NLOC, IPULL,MOLKIND,'NLOC, IPULL,MOLKIND'
      !return
C    *******************************************************************

        IF ( NLOC .LT. N(MOLKIND) ) THEN

C       ** CLOSE UP THE ARRAY LOCATE AFTER THE REMOVAL **

           DO  K = NLOC + 1, N(MOLKIND)

              LOCATE(K - 1,MOLKIND) = LOCATE(K,MOLKIND)
       !       write(*,*)LOCATE(K - 1,MOLKIND),LOCATE(K,MOLKIND)

         enddo

C       ** PLACE THE GHOST ATOM IPULL JUST OUTSIDE THE ACTIVE **
C       ** RANGE OF THE ARRAY LOCATE FOR FUTURE USE           **

           LOCATE(N(MOLKIND),MOLKIND) = IPULL
           
        !   WRITE(*,*) N(MOLKIND),LOCATE(N(MOLKIND),MOLKIND),IPULL

        ENDIF
C	WRITE(21,*) 'IPULL',IPULL
C	PAUSE
        RETURN
        END
C-------------------------------------------------------------------------
