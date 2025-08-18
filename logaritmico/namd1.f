        SUBROUTINE NAMD1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DECLARATIONS RELATED TO READING THE GAS & STRUCTURE IN
!!!
      integer imolinfile2,iresmolec2,iatomsmolec2,imolec2  ! for gas molecule in
      real xmolec2(10000,100),ymolec2(10000,100),zmolec2(10000,100)          ! for gas molecule in
      character*6 ATOMmolec2(10000)                            ! for gas molecule in
      character*4 ATOMNAMEmolec2(10000,100)                        ! for gas molecule in
      character*3 RESNmolec2(10000)                            ! for gas molecule in
      character*4 SEGMmolec2(10000,100)                            ! for gas molecule in
      character*2 ELEMmolec2(10000,100)                            ! for gas molecule in
	real V,VG,VA
      character*30 filenamegas(50)
      integer numeatomgasl(50)
!--------------------------------------------------------------------------------------
      real beta
      CHARACTER NAM*16
      integer iatom
      integer iatomsmolec
      integer ii
      integer ij
      integer ijgas
      integer ijkl
      real xmolec,ymolec,zmolec
      integer i
      integer imolinfile
      integer imolinfile3
      integer imolinfile4
      integer inj
      integer iresmolec
      integer istrucinfile
      integer j
      integer ji
      integer k
      integer natomkind
      integer ngas
      integer nj
      integer nkindresi
      integer nn
      integer ntipo
      integer numeatomico
      integer numeatomic
      integer nstrucmax
      integer istructinfile
      real occ
      real x1,y1,z1
      character *50 namdfile1
	character *50 namdfile2

	character *50 namdfile3
      character *50 namdfile4
	
	character *50 namdfile5
	character *50 namdfile6
	character *50 namdfile7
      
      
      








		integer imolinfilefix
       integer totatomfix
       character (80) :: readline
       character (30) :: DUMMYATOM1
       character (30) :: DUMMYATOM2
       character (30) :: ATOMTYPE(10000,100)
       real charge (10000,10000)
       character (30) :: TEXTAUX1
       character (30) :: TEXTAUX2
       character (30) :: TEXTAUX3
       character (4):: atommolec(10000,100)
       character (30) :: COMMAND
       character (30) :: DELIM
       character (30) :: NTYPESEGM(50)
       logical salida
       INTEGER INDEX
       integer nresi
       integer natom(1000)
       integer auxnatom
       character *50 molinfile
       character *50 molinfile2
       character *50 molinfile3

       character *50 molinfile4
!       character *50 molinfile2
        character *4 nkindseg (20)
        character *20 currentseg
        character *4 ATOMNAMEmolec
        character *4 RESNmolec
        character *4 SEGMmolec
        character *4 currentsegmolec
        character *2 ELEMmolec
        character *2 ELEMmolecfix
     	real ACEL, acelx,acely,acelz !TAMAÑO DE LA CELDA

        integer atominmolec
        integer imolec
        REAL X(5000,50,10),Y(5000,50,10),Z(5000,50,10)
        character *30 param
        integer paramfile
        INTEGER NKINDATOM
        REAL EPS(50),SIGMA(50),DUMMYNATOM
        CHARACTER *6 ATOMKIND(50)
        character *4 atomfix (50)
        character *3 atomGAS (50)
        integer natomfix
        logical latomfix(50)
          logical lauxatomfix
          integer numefix(50),numegas(50)
          logical numeatomosfix
          integer natomgasmole(50)
          character *2 numeatomicread(120)
          character *2 numeatomicreadaux
        logical lnatompot
	open(20,file='NAMD_MC.inp')
	read(20,*)namdfile1
	read(20,*)namdfile2
	read(20,*)namdfile3

	read(20,*)namdfile4
	read(20,*)namdfile5
	read(20,*)namdfile6
	read(20,*)namdfile7


	close(20)

	open(10,file='input.txt')
	read(10,*)readline
	read(10,*)readline
      read(10,*)readline
      read(10,*)readline
     	read(10,*)ACEL, acelx,acely,acelz !TAMAÑO DE LA CELDA

      read (10,*) readline      
	READ(10,*)readline
	read(10,*)readline
	read(10,*)readline
	read(10,*)nam
	close(10)










       open (file=namdfile1, unit=85)
       do i=1,118
       read (85,*,end=21) numeatomicreadaux,numeatomic
	write(*,*) numeatomicreadaux,numeatomic
       numeatomicread(numeatomic)=adjustl(trim( numeatomicreadaux))
       enddo
  21  close(85)
      write(*,*) 'Simbolos leidos'
!      pause



        salida=.false.
          istructinfile=10
          nresi=0

         open(file=namdfile2,unit=istructinfile)
      nkindresi=0
      nstrucmax=9000
      do j=1,nstrucmax
        read (istructinfile,'(A)',end=20) readline
 !                  write(*,*) readline
             if ((readline(1:4) .eq. 'RESI')) then
                nresi=nresi+1
                                     natom(nresi)=0
  !                        pause
                DELIM=' '
                INDEX= SCAN(READLINE,DELIM)
!                WRITE(*,*) INDEX
                read (readline(1:INDEX),*)  COMMAND
!                WRITE(*,*) command, ' command'

        read (readline(INDEX+1:LEN(READLINE)),*) NTYPESEGM(nresi)
                NTYPESEGM(nresi)=adjustl(trim(NTYPESEGM(nresi)))
                !write(*,*) NTYPESEGM(NRESI) , ' ntypesegm'
                nkindresi=nkindresi+1
                !PAUSE
                
 !               write(*,*)nkindresi ,' nkindresi'
                nkindseg(nkindresi)=ntypesegm(NKINDRESI)
 !                write(*,*)nkindseg(nkindresi), 'nkindseg'
                 do nn=1,nstrucmax
                read (istructinfile,'(A)',end=20) readline
!                                     write(*,*)readline
                             !        pause


                     if ((readline(1:5) .eq. 'GROUP')) then
                     salida=.false.
                        do i=1,100000
                           read (istructinfile,'(A)',end=20) readline
 !                           write(*,*)readline
                            !pause
                           if ((readline(1:4) .eq. 'ATOM')) then
                           natom(nresi)=natom(nresi)+1
                           !write(*,*)nresi,natom(nresi)
                              INDEX= SCAN(READLINE,DELIM)
                              DUMMYATOM1=READLINE(1:index)
                           !   write(*,*)DUMMYATOM1
                              textaux1=READLINE(index+1:len(readline))
                           !   write(*,*)textaux1
                              INDEX= SCAN(TEXTAUX1,DELIM)
                              DUMMYATOM2=TEXTAUX1(1:index)
                            !  write(*,*)DUMMYATOM2
                              atommolec(natom(nresi),nresi)=dummyatom2
                              textaux2=(TEXTAUX1(index+1:len(TEXTAUX1)))
                            !  write(*,*) textaux2
                              textaux2=adjustl(textaux2)
                            !  write(*,*) textaux2
                              INDEX= SCAN(TEXTAUX2,DELIM)

                              auxnatom= natom(nresi)
                              ATOMTYPE(auxnatom,nresi)=TEXTAUX2(1:index)
                              
                              textaux3=(TEXTAUX2(index+1:len(TEXTAUX2)))
                             ! INDEX= SCAN(TEXTAUX3,DELIM)
                              textaux3=adjustl(textaux3)
                              
                              

                              
                              !write(*,*)ATOMTYPE(natom(nresi),nresi),j,' atomtype'
                              read(TEXTAUX3,*)  CHARGE(auxnatom,nresi)
                              if(charge(natom(nresi),nresi).ne.0) then
                              !write(*,*) charge (j)
                             ! pause
                              endif
                              !pause
                              elseif ((readline(1:5) .eq. 'GROUP')) then
                              salida=.true.
                               exit
                           endif
                        enddo
                endif
                
                if(salida) then
                salida=.false.
                exit
                endif

                enddo
                
            endif
            
            enddo
            

  20  write(*,*) 'Termina lectura PSF'
      write(*,*) 'nresi=',nresi
      do i=1,nresi

 !     write(*,*)nkindseg(i)
      do j=1,natom(i)
!      write(*,*)atommolec(j,i),ATOMTYPE(j,i),charge(j,i)
      enddo
      !pause
      enddo

!------------------------------------------------------------------------------------------
!      Lectura del archivo de par metros
!------------------------------------------------------------------------------------------
       param= namdfile3
       paramfile=94
      open (file=param,unit=paramfile)
       open (file='LJ.dat',unit=25)
       natomkind=1
      do j=1,100000
        read (PARAMFILE,'(A)',end=1) readline
        if ((readline(1:9) .eq. 'NONBONDED')) then
           DO I=1,50
           latomfix(i)=.false.
           lauxatomfix=.false.
           read (PARAMFILE,'(A)',end=1) readline
           !WRITE(*,*)READLINE
                DO II=1,NRESI
                  DO JI=1,NATOM(II)
           !WRITE(*,*)ATOMTYPE(JI,II)
                          IF (READLINE(1:6).EQ.ATOMTYPE(JI,II)) THEN
      READ(READLINE,*)ATOMKIND(natomkind),DUMMYNATOM,EPS(Natomkind),
     +SIGMA(Natomkind)
                      do k=1,natomkind-1
                      if(ATOMKIND(natomkind).eq.atomkind(k)) then
                lauxatomfix=.true.
                latomfix(i)=latomfix(i).or.lauxatomfix
                                                      endif
                      enddo
                      

                if(.not.latomfix(i)) then
                      
          EPS(Natomkind) =EPS(Natomkind)*(-1000)/1.98
          SIGMA(Natomkind)=SIGMA(Natomkind)*2
      natomkind=natomkind+1
      
      exit
      endif
      
                         ENDIF
                  ENDDO
             ENDDO
      ENDDO

      ENDIF
      ENDDO
      

 1    WRITE(*,*)'TERMINA LECTURA PARAMETROS'
      write(25,*) natomkind-1, ' atomos distintos'

	
      natomkind=natomkind-1
      write(*,*) natomkind, ' atomos distintos'

 !     PAUSE
       DO K=1, NATOMKIND
             write(*,*) k,eps(k),sigma(k),charge(j,i),atomkind(k)
             ENDDO
!             PAUSE
             

!---------------------------------------------------------------------------------------------------
!      Separaci¢n entre estructura y mol‚culas del gas
!---------------------------------------------------------------------------------------------------
      molinfile3 = namdfile4  ! file where molecule is read from
      imolinfile3 = 13
        open(unit=imolinfile3,file=molinfile3)
        do i=1,50000
        natomfix=i
        read (imolinfile3,'(A)',end=5) atomfix(i)
        write(*,*)natomfix,atomfix(i),i
        enddo

 5    WRITE(*,*)'TERMINA LECTURA ATOMOS FIJOS'
       NATOMFIX=NATOMFIX-1
      WRITE(*,*) NATOMFIX, ' ATOMOS FIJOS'
	do i=1,natomfix
	imolinfilefix=i+151
	write(*,*)imolinfilefix,atomfix(i)
	open(unit=imolinfilefix,file=atomfix(i))
	enddo

        ijgas=0
        do i=1,nresi
        latomfix(i)=.false.
        do j=1,natomfix

        if (atomfix(j).EQ.NTYPESEGM(i)) then
        lauxatomfix=.true.
        latomfix(i)=latomfix(i).or.lauxatomfix
        lauxatomfix=.false.
        ENDIF
!        exit
        ENDDO
        ENDDO
        
        do i=1,nresi
        if(.not.latomfix(i)) then
        ijgas=ijgas+1
        atomgas(ijgas)=NTYPESEGM(i)
        natomgasmole(ijgas)=natom(i)
        endif
        enddo
           ! pause
        DO I=1,IJGAS
        WRITE(*,*)IJGAS,' MOLECULAS GASEOSAS'
        WRITE(*,*)ATOMGAS(I)
        ENDDO

       ngas=ijgas

      do i=1,natomfix
      numefix(i)=0
      enddo
      do i=1,ngas
      numegas(i)=0
      enddo

!-------------------------------------------------------------------------------------------


              !PAUSE

      molinfile = namdfile5  ! file where molecule is read from
      imolinfile = 11                    ! unit where molecule is read from
      molinfile2 = nam
      imolinfile2 = 14
      molinfile4 = namdfile6
      imolinfile4 = 15                    ! unit where molecule is read from

        open(unit=imolinfile,file=molinfile)
        open(unit=imolinfile2,file=molinfile2)
        open(unit=imolinfile4,file=molinfile4)

!-------------------------------------------------------------------------------------------------
!      COMIENZA LECTURA
      i = 1
      numefix=0
      do j=1,80000
      numeatomosfix=.false.
        read (imolinfile,'(A)',end=10) readline
        !write(*,*)readline
        if (readline(1:4).eq.'ATOM'.or.readline(1:6).eq.'HETATM') then
            read (readline(7:11),*)  iatom
            read (readline(13:16),'(A4)') ATOMNAMEmolec
            read (readline(18:22),'(A4)') RESNmolec
            read (readline(23:26),*) iresmolec
            read (readline(31:38),*) xmolec
            read (readline(39:46),*) ymolec
            read (readline(47:54),*) zmolec
			if(abs(xmolec).gt.acelx/2) then
			write(*,*) 'OUT OF LIMITS'
			write(*,*) readline,' x ', acelx/2
			 stop
			end if


			if(abs(ymolec).gt.acely/2) then 

			write(*,*) 'OUT OF LIMITS'

			write(*,*) readline,' y ', acely/2
			 stop
			end if

			if(abs(zmolec).gt.acelz/2) then 

			write(*,*) 'OUT OF LIMITS'

			write(*,*) readline,' z ', acelz/2
			 stop
			end if

            read (readline(55:60),*) occ
            read (readline(61:66),*) beta
            read (readline(73:76),'(A4)') SEGMmolec
            read (readline(77:78),'(A2)') ELEMmolec
                 do nj=1,natomfix
                 currentseg=adjustl(trim(atomfix(nj)))
                 currentsegmolec=adjustl(trim(resNmolec))
		!write(*,*)currentsegmolec,' ', currentseg
		!pause
                  if (currentsegmolec.eq. currentseg) then
                  numefix(nj)=numefix(nj)+1
                ENDIF
                enddo
                  do nj=1,ngas
                 currentseg=adjustl(trim(atomgas(nj)))
                 currentsegmolec=adjustl(trim(resNmolec))
                  if (currentsegmolec.eq. currentseg) then


                  numegas(nj)=numegas(nj)+1
                ENDIF
                enddo

        ENDiF
      enddo
 10   close(imolinfile)
!------------------------------------------------------------------------------------------------
      do k=1,natomkind
       lnatompot=.false.
      do i=1,nresi
      write(*,*) i,natom(i)
!      pause
      do j=1,natom(i)

!      write(*,*)k
!      write(*,*)i,j

      atomkind(k)=adjustl(trim(atomkind(k)))
      ATOMTYPE(j,i)=adjustl(trim(ATOMTYPE(j,i)))
 !     write(*,*)'1',atomkind(k),' 2 ',ATOMTYPE(j,i)
!      write(*,*) k,eps(k),sigma(k),charge(j,i)
 !     pause
      if(atomkind(k).eq.ATOMTYPE(j,i)) then


      write(25,*)k,eps(k),sigma(k),charge(j,i),atomkind(k),ntypesegm(i),
     + j,i
      write(*,*) k,eps(k),sigma(k),charge(j,i),atomkind(k),ntypesegm(i),
     +j,i
      lnatompot=.true.
!      pause
      exit
      endif
      enddo
      if(lnatompot)exit
      enddo
      enddo

 	close(25)
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
      molinfile = namdfile5  ! file where molecule is read from
      imolinfile = 11                    ! unit where molecule is read from
      molinfile2 = nam
      imolinfile2 = 14
      molinfile4 = namdfile6
      imolinfile4 = 15
      totatomfix=0.
      do i=1,natomfix
      write(*,*) i, numefix(i)  ,ATOMFIX(i)
      enddo
      do i=1,ngas
      numegas(i)=numegas(i)/natomgasmole(i)
      write(*,*) i,numegas(i) ,' molec. gas',natomgasmole(i)
      enddo

	!pause

      do i=1,natomfix
      totatomfix=totatomfix+numefix(i)
      enddo

                          ! unit where molecule is read from
       write( imolinfile2,*) totatomfix
       write(*,*)     totatomfix
       !pause

        open(unit=imolinfile,file=molinfile)
      do j=1,80000
	!if(iatom.gt.140) pause
        read (imolinfile,'(A)',end=80) readline
	write(*,*)'--------------------------------'
        write(*,*)readline
	!pause
        if (readline(1:4).eq.'ATOM'.or.readline(1:6).eq.'HETATM') then
            read (readline(7:12),*)  iatom
		write(*,*) iatom
            read (readline(13:17),'(A4)') ATOMNAMEmolec
		write(*,*)ATOMNAMEmolec
            read (readline(18:22),'(A4)') RESNmolec
		write(*,*) RESNmolec
            read (readline(23:26),*) iresmolec
		write(*,*)iresmolec
            read (readline(31:38),*) xmolec
		write(*,*) xmolec
            read (readline(39:46),*) ymolec
		write(*,*) ymolec
            read (readline(47:54),*) zmolec
		write(*,*) zmolec
            read (readline(55:60),*) occ
		write(*,*) occ
            read (readline(61:66),*) beta
		write(*,*)beta
            read (readline(73:76),'(A4)') SEGMmolec
		write(*,*)SEGMmolec
            read (readline(77:78),'(A2)') ELEMmolec
		write(*,*)ELEMmolec
		write(*,*)'--------------------------------'!---------------------------
            do numeatomic=1,118
            ELEMmolecfix=adjustl(trim(ELEMmolec))
            
            if  (ELEMmolecfix.eq.numeatomicread(numeatomic)) then
            numeatomico=numeatomic
            exit
            endif
            enddo
!-----------------------------------------------------------------------------

            do nj=1,nresi
            currentseg=adjustl(trim(nkindseg(nj)))
            currentsegmolec=adjustl(trim(resNmolec))
		write(*,*)'********'
		write(*,*)currentseg,len(currentseg)
		write(*,*)currentsegmolec,len(currentsegmolec),' **'
		write(*,*)'********'

           ! pause
            if (currentsegmolec.eq. currentseg) then
		write(*,*) nkindseg(nj),'//////////////////son iguales'


            IMOLEC=nJ
                write(*,*) natom(imolec), 'natom'
        !pause
               DO IJ=1,NATOM(IMOLEC)
!      write(*,*) ATOMNAMEMOLEC,currentseg,nkindseg(nj),nj ,
!     +atommolec(ij,nj)
            !pause
               atommolec(ij,nj)=adjustl(trim(atommolec(ij,nj)))
               ATOMNAMEMOLEC=adjustl(trim(ATOMNAMEMOLEC))
               IF(ATOMNAMEMOLEC.EQ.atommolec(ij,nj)) THEN
               ATOMINMOLEC=IJ
               x(iresmolec,ATOMINMOLEC,iMOLEC)=xmolec
               Y(iresmolec,ATOMINMOLEC,iMOLEC)=Ymolec
               Z(iresmolec,ATOMINMOLEC,iMOLEC)=Zmolec
               do inj=1,natomkind
               atomtype(ij,nj)=adjustl(trim(atomtype(ij,nj)))
               atomkind(inj)=adjustl(trim(atomkind(inj)))
               !write(*,*) atomkind(inj),atomtype(ij,nj)
              ! pause
               if (ATOMKIND(inj).eq.atomtype(ij,nj)) then
               ntipo=inj
               !write(*,*)natomfix
!------------------------------------------------------------------------------------               !pause
                DO IJKL=1,NATOMFIX
                resnmolec=adjustl(trim(resnmolec))
                atomfix(ijkl)=adjustl(trim(atomfix(ijkl)))
                write(*,*)RESNmolec,' ',ATOMFIX(IJKL),ijkl
                !pause
                IF (RESNmolec.EQ.ATOMFIX(IJKL)) THEN
                write(*,*)RESNmolec,' ',ATOMFIX(IJKL),ijkl
		write(ijkl+151,'(A)')readline
!                write(*,*) 'atomo fijo',RESNmolec,atomfix(ijkl)
		write(*,*) 'EQ ',RESNmolec,' ',ATOMFIX(IJKL),ijkl,iatom
               write(imolinfile2,*)x(iresmolec,ATOMINMOLEC,iMOLEC),
     +  y(iresmolec,ATOMINMOLEC,iMOLEC),z(iresmolec,ATOMINMOLEC,iMOLEC),
     +EPS(inj),
     +SIGMA(inj),CHARGE(ij,nj),  numeatomic

                  !pause
                exit
		else 
		if(iatom.ge.150) then
		write(*,*) 'NE ',RESNmolec,' ',ATOMFIX(IJKL),ijkl,iatom
		!pause
		endif

                ENDIF
                ENDDO
!-------------------------------------------------------------------------------------
                DO IJKL=1,NGAS
                resnmolec=adjustl(trim(resnmolec))
                atomgas(ijkl)=adjustl(trim(atomGAS(ijkl)))
!                write(*,*)RESNmolec,' ',ATOMFIX(IJKL),ijkl
!                pause
                IF (RESNmolec.EQ.ATOMGAS(IJKL)) THEN
        x(iresmolec,ATOMINMOLEC,ijkl)=x(iresmolec,ATOMINMOLEC,iMOLEC)
        y(iresmolec,ATOMINMOLEC,ijkl)=y(iresmolec,ATOMINMOLEC,iMOLEC)
        z(iresmolec,ATOMINMOLEC,ijkl)=z(iresmolec,ATOMINMOLEC,iMOLEC)
               write(imolinfile4,*)x(iresmolec,ATOMINMOLEC,ijkl),
     +  y(iresmolec,ATOMINMOLEC,ijkl),z(iresmolec,ATOMINMOLEC,ijkl),
     +iresmolec,atominmolec,ijkl,ntipo,CHARGE(ij,nj),  ATOMKIND(inj)
	write(*,*) readline
               write(*,*)x(iresmolec,ATOMINMOLEC,ijkl),
     +  y(iresmolec,ATOMINMOLEC,ijkl),z(iresmolec,ATOMINMOLEC,ijkl),
     +iresmolec,atominmolec,ijkl,ntipo,CHARGE(ij,nj),  ATOMKIND(inj)
	!pause

                !exit
                ENDIF
                ENDDO
!------------------------------------------------------------------------------------
      exit
                     endif
               enddo



               ENDIF
               ENDDO
            ENDIF
            ENDDO


        ENDiF
      enddo



      
      
      
 80   close(imolinfile)
           open(unit=67,file='initconf.txt')
          write(67,*)V,VG,VA,' energias'
          write(67,*)ngas,' tipos de molec'
          do i=1,ngas
              !molkind=i
              write(67,*)natomgasmole(i)
          write(67,*)numegas(i), ' natom ', i
	write(*,*)numegas(i), ' natom ', i
	!pause

          do j=1,numegas(i)
		write(*,*)j,numegas(i),natomgasmole(i)

              DO k=1,natomgasmole(i)
                  x1=x(j,k,i)/acel
                  y1=y(j,k,i)/acel
                  z1=z(j,k,i)/acel
              write(67,*)X1,Y1,Z1, 'x y z molec'
		write(*,*)X1,Y1,Z1, j,k,i
	!pause
              enddo
	!pause
          enddo
          enddo
	close(67)


!------------------------------------------------------------------------------------------
!     Lectura de la mol‚cula gaseosa
!-----------------------------------------------------------------------------------------

                           
      molinfile =namdfile7	! file where molecule is read from
      imolinfile = 11                    ! unit where molecule is read from
      open(unit=imolinfile,file=molinfile,status='OLD',
     +form='formatted',access='sequential',action='read')
!!! now read input molecule

      !do i=1,9000
      numeatomgasl=0
      !enddo
             i = 1

      do j=1,100000
        read (imolinfile,'(A)',end=110) readline
        write(*,*)readline
        if ((readline(1:4) .eq. 'ATOM') .or.
     +(readline(1:6) .eq. 'HETATM')) then
            read (readline(7:11),*)  iatom

            read (readline(23:26),*) iresmolec2
            read (readline(18:20),'(A3)') RESNmolec2(iresmolec2)
            i=numeatomgasl(iresmolec2)  +1
            read (readline(13:16),'(A4)') ATOMNAMEmolec2(i,iresmolec2)
            read (readline(31:38),*) xmolec2(i,iresmolec2)
            read (readline(39:46),*) ymolec2(i,iresmolec2)
            read (readline(47:54),*) zmolec2(i,iresmolec2)
            read (readline(55:60),*) occ
            read (readline(61:66),*) beta
            read (readline(73:76),'(A4)') SEGMmolec2(i,iresmolec2)
            read (readline(77:78),'(A2)') ELEMmolec2(i,iresmolec2)
            numeatomgasl(iresmolec2)=numeatomgasl(iresmolec2)+1

        ENDiF
      enddo
 110   close(imolinfile)
      iatomsmolec = i  -1
      write(*,*)iresmolec2,' ********'
      do i=1,iresmolec2               ! contains number of atoms in molecule read
      imolec = iresmolec2   ! number of molecules read
      print *, '------------------------------------'
      print *, 'finished reading file ', molinfile
      print *, 'read molecule(s) = ', imolec
      print *, 'read     atom(s) = ', iatomsmolec
      print *, 'atoms/molecule   = ',iatomsmolec/imolec
      print *, '------------------------------------'
      enddo
      open(file='MOLEC.DAT',unit=45)
      write(45,*) imolec
      do i=1,imolec
       write(45,*)RESNmolec2(i)
       enddo
       close(45)
       
       
       
       
      do i=1,imolec
      
      write(*,*) 'arranco molecula', i
      open(file=RESNmolec2(i),unit=45)
      write(45,*) numeatomgasl(i),1,i,i
       do j=1,numeatomgasl(i)
          write(*,*) ' atomo ',j, 'de', i
        do numeatomic=1,118
            ELEMmolecfix=adjustl(trim(ELEMmolec2(j,i)))

          if  (ELEMmolecfix.eq.numeatomicread(numeatomic)) then
             numeatomico=numeatomic
             !exit
          endif
         enddo
!-----------------------------------------------------------------------------

        do nj=1,nresi
            currentseg=adjustl(trim(nkindseg(nj)))
            currentsegmolec=adjustl(trim(resNmolec2(i)))
            write(*,*)currentseg,'  ',currentsegmolec, nj
            !pause
         if (currentsegmolec.eq. currentseg) then

            IMOLEC=nJ
               DO IJ=1,NATOM(IMOLEC)
		write(*,*)imolec, ' imolec', natom(imolec)
            !pause
               atommolec(ij,nj)=adjustl(trim(atommolec(ij,nj)))
               ATOMNAMEMOLEC=adjustl(trim(ATOMNAMEmolec2(j,i)))
               !write(*,*)atommolec(ij,nj),' ',   ATOMNAMEMOLEC
               !pause
               
               IF(ATOMNAMEMOLEC.EQ.atommolec(ij,nj)) THEN
               ATOMINMOLEC=IJ
               do inj=1,natomkind
               atomtype(ij,nj)=adjustl(trim(atomtype(ij,nj)))
               atomkind(inj)=adjustl(trim(atomkind(inj)))
               if (ATOMKIND(inj).eq.atomtype(ij,nj)) then
               ntipo=inj




      write(45,*)xmolec2(j,i),ymolec2(j,i),zmolec2(j,i),inj,numeatomico
               endif
               enddo
               endif
               enddo
               endif
      
      enddo
       enddo
      close(45)
       enddo

      close(imolinfile2)

      return
      end

