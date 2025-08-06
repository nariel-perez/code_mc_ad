      subroutine namd2(escala,NC,N,NMOLEC,RX,RY,RZ,natom,locate)

!!! compilation command:  gfortran -O -ffree-form Gas_Generator-CW.f -o Gas_Generator-CW.exe
!!! running the code:  ./Gas_Generator-CW.exe
!!!
!!! generates initial configuration for natural gas mixture in the system
!!! 
!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! basic parameters

       implicit none

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      real acel
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
      INTEGER N2
      REAL AN1,U1,UNG1,UG1,UNA1,UA1,UN1,AN2
      REAL CALOR,CALORG,CALORA
      REAL ESCALA
      REAL RXAI,RYAI,RZAI, EPSAI,SGCI,QACI
      INTEGER SYMBOL
      INTEGER K,jin
      REAL RXN,RYN,RZN
      REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)

      real aitest76,aitest77
      real diel
      integer ensemble
      integer nmolec2
      integer natom2
      integer molkind
      integer ncantmol

      REAL CNF(-25:25,-25:25,-25:25,10,50)
      !REAL, allocatable :: CNF(:,:,:,:,:)
      INTEGER ICNF,JCNF,KCNF,NCELLMAT
      INTEGER NATOMKINDI
      REAL CALORESP1,CALORESP2,CALORESP3
      INTEGER NCONFMIN,NCONFMAX
      INTEGER NESTADO
!--------------------------------------------------------------------------------------------------
     
     
     
     
     
     
     
     
     

      !real, parameter :: RMIN=2.00        ! minimum possible distance between atoms
                                         ! (eventually could be substituted by Rmin vdW)
      real xij,yij,zij,rij2               ! to calculate distance between atoms

      integer, parameter :: ngasmax=10000    ! maximum number of GAS atoms to be generated
      integer, parameter :: namax=1000        ! maximum number of atoms to be read from GAS
      integer, parameter :: nstrucmax = 5000 ! maximum number of atoms to be read from STRUCTURE

1000  FORMAT (A4,2X,I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,F8.3,F8.3,F8.3,
     +F6.2,F6.2,6X,A4,A2)
1001  format (A6,   I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,F8.3,F8.3,F8.3,
     +F6.2,F6.2,6X,A4,A2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DECLARATIONS RELATED TO MISC STUFF
!!!
      integer iclock,iseed                            ! seeds for rand()
      
      integer iatom,jprev,inew,icollision
      real pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DECLARATIONS RELATED TO READING THE GAS & STRUCTURE IN
!!!
      integer imolinfile,iresmolec,iatomsmolec,imolec  ! for gas molecule in
      real xmolec(50,namax),ymolec(50,namax),zmolec(50,namax)          ! for gas molecule in
      character*6 ATOMmolec(50,namax)                            ! for gas molecule in
      character*4 ATOMNAMEmolec(50,namax)                        ! for gas molecule in
      character*3 RESNmolec(50,namax)                            ! for gas molecule in
      character*4 SEGMmolec(50,namax)                            ! for gas molecule in
      character*2 ELEMmolec(50,namax)                            ! for gas molecule in
      character*20 molinfile                                  ! for gas molecule in


      integer istructinfile,iresstruct(nstrucmax),iatomsstruct,istruct  ! for structure in
      real xstruct(nstrucmax),ystruct(nstrucmax),zstruct(nstrucmax)     ! for structure in
      character*6 ATOMstruct(nstrucmax)                                 ! for structure in
      character*4 ATOMNAMEstruct(nstrucmax)                             ! for structure in
      character*3 RESNstruct(nstrucmax)                                 ! for structure in
      character*4 SEGMstruct(nstrucmax)                                ! for structure in
      character*2 ELEMstruct(nstrucmax)                                 ! for structure in
      character*20 structureinfile                                      ! for structure in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DECLARATIONS RELATED TO CREATING THE GAS (OUT)
!!!
      integer igasfile,iresgas(ngasmax),iatomsgas,igas        ! for gas OUT
      real xgas(ngasmax),ygas(ngasmax),zgas(ngasmax)          ! for gas OUT
      character*6 ATOMgas(ngasmax)                            ! for gas OUT
      character*4 ATOMNAMEgas(ngasmax)                        ! for gas OUT
      character*3 RESNgas(ngasmax)                            ! for gas OUT
      character*4 SEGMgas(ngasmax)                            ! for gas OUT
      character*2 ELEMgas(ngasmax)                            ! for gas OUT
      character*20 gasOUTfile                                 ! for gas OUT
      integer Ngas                                            ! for gas OUT
      real LzMIN,LzMAX,Ly,Lx,deltax,deltay,deltaz             ! for gas OUT
      real eulerphi,eulerpsi,eulertheta                       ! for gas OUT
      real xtemp,ytemp,ztemp                                  ! for gas OUT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! MISC DECLARATIONS
      character*80 readline
      real occ,beta
      pi=4.*atan(1.)

	!WRITE(*,*)'-------------------------------'
	!write(*,*) escala,NC,N,NMOLEC,RX,RY,RZ
	!WRITE(*,*)'-------------------------------'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Input GAS molecule file 
!!! 
!!! reads file molinfile (GAS_MOLECULE_IN.pdb)
!!! stores:
!!!   ATOMNAMEmolec(i) - atom names (C1,H1X,H1Y,H1Z,H1W)
!!!   RESNmolec(i)     - RESIDUE (NC1)
!!!   iresmolec(i)        - molecule number (should be just 1)
!!!   xmolec(i), ymolec(i), zmolec(i) - coordinates
!!!   occ              - occupation, should be 1.00
!!!   beta             - beta, should be 0.00
!!!   SEGMmolec(i)     - segment (MET)
!!!   ELEMmolec(i)     - atom type (C, H)
!!! at the end of the read, should also provide 
!!!   iatomsmolec = number of atoms
!!!   imolec = number of molecules
!!!
!!!
      molinfile = 'GAS_MOLECULE_IN.pdb'  ! file where molecule is read from
      imolinfile = 11                    ! unit where molecule is read from 
      open(unit=imolinfile,file=molinfile,status='OLD',
     +form='formatted',access='sequential',action='read')
!!! now read input molecule 
      i = 1
      do j=1,namax
        read (imolinfile,'(A)',end=10) readline
        if ((readline(1:4) .eq. 'ATOM') .or.
     +(readline(1:6) .eq. 'HETATM')) then
            read (readline(7:11),*)  iatom
            read (readline(23:26),*) iresmolec
            read (readline(13:16),'(A4)') ATOMNAMEmolec(iresmolec,iatom)
            read (readline(18:20),'(A3)') RESNmolec(iresmolec,iatom)
            read (readline(31:38),*) xmolec(iresmolec,iatom)
            read (readline(39:46),*) ymolec(iresmolec,iatom)
            read (readline(47:54),*) zmolec(iresmolec,iatom)
            read (readline(55:60),*) occ
            read (readline(61:66),*) beta
            read (readline(73:76),'(A4)') SEGMmolec(iresmolec,iatom)
            read (readline(77:78),'(A2)') ELEMmolec(iresmolec,iatom)
            i = i + 1
        ENDiF
      enddo
 10   close(imolinfile)
      iatomsmolec = i - 1               ! contains number of atoms in molecule read
      imolec = iresmolec   ! number of molecules read
      print *, '------------------------------------'
      print *, 'finished reading file ', molinfile
      print *, 'read molecule(s) = ', imolec
      print *, 'read     atom(s) = ', iatomsmolec
      print *, 'atoms/molecule   = ',iatomsmolec/imolec
      print *, '------------------------------------'
!
!
!!!!!! JUST PRINTING THE READ GAS MOLECULE TO CHECK
!1011 format (A6,   I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,F8.3,F8.3,F8.3,F6.2,F6.2,6X,A4,A2)
!     do i = 1, iatomsmolec 
!         write (*,1011) 'ATOM  ',i,ATOMNAMEmolec(i),RESNmolec(i),iresmolec(i),xmolec(i) &
!          ,ymolec(i),zmolec(i),occ,beta,SEGMmolec(i),ELEMmolec(i)    
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!
!
!!!!!! JUST PRINTING THE READ STRUCTURE TO CHECK
!1021 format (A6,   I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,F8.3,F8.3,F8.3,F6.2,F6.2,6X,A4,A2)
!     do i = 1, iatomsstruct 
!         write (*,1021) 'ATOM  ',i,ATOMNAMEstruct(i),RESNstruct(i),iresstruct(i),xstruct(i) &
!           ,ystruct(i),zstruct(i),occ,beta,SEGMstruct(i),ELEMstruct(i)    
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Now generate Nmolec gas molecules at random locations (and orientations)
!!! while avoiding the structure and previous atoms (rij < RMIN ==> reject)
!!!
!!! Ngas = number of molecules to create

	write(*,*)NMOLEC
	!pause

            DO I=1,NMOLEC
      print *, 'GENERATING GAS...'
      print *, 'number of molecules = ', n(i)
      print *, 'number of atoms = ', NATOM(I),n(i)
      print *, '******************************************'

          do j=1,N(I)
              jin=locate(j,i)
              DO K=1,NATOM(I)
                  RXN=RX(Jin,K,I)*ESCALA
                  RYN=RY(Jin,K,I)*ESCALA
                  RZN=RZ(Jin,K,I)*ESCALA
                  write(*,*)NSYM(K,I),RXN,RYN,RZN
              enddo
          ENDDO
      ENDDO


	iatom = 0.

      do k=1,nmolec
      
      do j = 1, N(k)          ! number of molecules to generate
        jin=locate(j,k)
!666     print *, 'I am on molecule ', j 
      do i = 1, natom(k)      ! number of atoms/molecule
          iatom = iatom + 1       ! index, atom number
          ATOMNAMEgas(iatom) = ATOMNAMEmolec(k,i)
          RESNgas(iatom) = RESNmolec(k,i)
          SEGMgas(iatom) = SEGMmolec(k,i)
          ELEMgas(iatom) = ELEMmolec(k,i)
          iresstruct(iatom) = j
          !! rotating the molecule (all atoms around the origin)
          xgas(iatom) = RX(Jin,i,k)*ESCALA
          ygas(iatom) = RY(Jin,i,k)*ESCALA
          zgas(iatom) = RZ(Jin,i,k)*ESCALA
        enddo
!! COLLISION DETECTION... 
        !
!! COLLISION DETECTION... WITH EARLIER MOLECULES/ATOMS
      enddo
      enddo
!
!
      print *, 'DONE...'
      print *, 'created gas atoms = ', iatom
      print *, '******************************************'
!
!
      gasOUTfile = 'GAS.pdb'         ! file where GAS goes to
      igasfile = 31               ! unit  
      open(unit=igasfile,file=gasOUTfile,status='unknown',
     +form='formatted', access='sequential',action='write')
!!!!!! WRITING THE GAS MOLECULES TO gasOUTfile
1031  format (A6,   I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,F8.3,F8.3,
     +F8.3,F6.2,F6.2,6X,A4,A2)
      do i = 1, iatom
         write (igasfile,1031) 'ATOM  ',i,ATOMNAMEgas(i),
     +RESNgas(i),iresstruct(i),xgas(i)
     +,ygas(i),zgas(i),occ,beta,SEGMgas(i),ELEMgas(i)
      enddo
      write (igasfile,'(A3)') 'END'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





       return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


