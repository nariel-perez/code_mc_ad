      IMPLICIT NONE
      INTEGER ISOT
      INTEGER NMOLEC
      INTEGER INMOLEC
      INTEGER NATOM(1000)
      INTEGER IISOT
      INTEGER I
      INTEGER ICNF,JCNF,KCNF,NCELLMAT
      INTEGER NATOMKINDI
      REAL CNF(-25:25,-25:25,-25:25,50,10)
      CHARACTER CONFIG*3
      CHARACTER CONFAT*3
      CHARACTER CONFNAT*3
      INTEGER ICNFi,JCNFi,KCNFi
      
      INTEGER IMIN,IMAX
      INTEGER JMIN,JMAX
      INTEGER KMIN,KMAX
      
      REAL XY(-25:25,-25:25,50,10)
      REAL XZ(-25:25,-25:25,50,10)
      REAL YZ(-25:25,-25:25,50,10)
      
      REAL X(-25:25,50,10)
      REAL Y(-25:25,50,10)
      REAL Z(-25:25,50,10)
      write(*,*)'-----------------------------------'
      WRITE(*,*)'NUMERO DE PUNTOS DE LA ISOTERMA?'
      READ(*,*)ISOT
      WRITE(*,*)'NUMERO DE MOLECULAS A CONSIDERAR?'
      READ(*,*)NMOLEC
      
      WRITE(*,*)'Z MIN?'
      READ(*,*)KMIN
      WRITE(*,*)'Z MAX?'
      READ(*,*)KMAX
      !
      WRITE(*,*)'X MIN?'
      READ(*,*)IMIN
      WRITE(*,*)'X MAX?'
      READ(*,*)IMAX
      
      WRITE(*,*)'Y MIN?'
      READ(*,*)JMIN
      WRITE(*,*)'Y MAX?'
      READ(*,*)JMAX
      
      !imin=-25
      !imax=25
      !jmin=-25
      !jmax=25
      !kmin=-25
      !kmax=25
      !nmolec=1
      !isot=3
      NCELLMAT=50

      DO INMOLEC=1,NMOLEC
          WRITE(*,*) 'NUMERO DE ATOMOS MOLECULA', INMOLEC
          READ(*,*) NATOM(INMOLEC)
      ENDDO
      
      

      DO IISOT=1,ISOT
          !------------------------------------------------------------------------------
      !RESET CONTADORES
      !------------------------------------------------------------------------
            Do I=1,nmolec
              DO NATOMKINDI=1,NATOM(I)
                  DO ICNF=-NCELLMAT/2,NCELLMAT/2
                      DO JCNF=-NCELLMAT/2,NCELLMAT/2
                          DO KCNF=-NCELLMAT/2,NCELLMAT/2
                              XY(ICNF,JCNF,I,NATOMKINDI)=0
                              XZ(ICNF,KCNF,I,NATOMKINDI)=0
                              YZ(JCNF,KCNF,I,NATOMKINDI)=0
                              
                              X(ICNF,I,NATOMKINDI)=0
                              Z(KCNF,I,NATOMKINDI)=0
                              Y(JCNF,I,NATOMKINDI)=0
                          ENDDO
                      ENDDO
                      
                  enddo
                  CLOSE(101)
              ENDDO
      ENDDO
!------------------------------------------------------------------------------------
      CONFIG='CONFIG'
	write(CONFIG,'(i3)') IISOT
      
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
	write(CONFNAT,'(i1)') NATOMKINDI
      open(101,file='CNF'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
      READ(101,*)NCELLMAT
                  DO ICNF=-NCELLMAT/2,NCELLMAT/2
                      DO JCNF=-NCELLMAT/2,NCELLMAT/2
                          DO KCNF=-NCELLMAT/2,NCELLMAT/2
          READ(101,*)icnfi,jcnfi,kcnfi,CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
!          write(*,*)icnfi,jcnfi,kcnfi,CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
!          if(CNF(ICNF,JCNF,KCNF,I,NATOMKINDI).ne.0) pause
          IF (KCNF.GE.KMIN.AND.KCNF.LE.KMAX) THEN
          XY(ICNF,JCNF,I,NATOMKINDI)=
     +XY(ICNF,JCNF,I,NATOMKINDI)+CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
          ENDIF
          IF (JCNF.GE.JMIN.AND.JCNF.LE.JMAX) THEN
                    XZ(ICNF,KCNF,I,NATOMKINDI)=
     +XZ(ICNF,KCNF,I,NATOMKINDI)+CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
          ENDIF
          IF (ICNF.GE.IMIN.AND.ICNF.LE.IMAX) THEN
                    YZ(JCNF,KCNF,I,NATOMKINDI)=
     +YZ(JCNF,KCNF,I,NATOMKINDI)+CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
          ENDIF
          
      IF (KCNF.GE.KMIN.AND.KCNF.LE.KMAX) THEN
          IF (JCNF.GE.JMIN.AND.JCNF.LE.JMAX) THEN
          X(ICNF,I,NATOMKINDI)=
     +X(ICNF,I,NATOMKINDI)+CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
          ENDIF
      ENDIF
      IF (ICNF.GE.IMIN.AND.ICNF.LE.IMAX) THEN
          IF (JCNF.GE.JMIN.AND.JCNF.LE.JMAX) THEN
      
                    Z(KCNF,I,NATOMKINDI)=
     +Z(KCNF,I,NATOMKINDI)+CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
                    
          ENDIF
      ENDIF
      IF (ICNF.GE.IMIN.AND.ICNF.LE.IMAX) THEN
          IF (KCNF.GE.KMIN.AND.KCNF.LE.JMAX) THEN      
                    Y(JCNF,I,NATOMKINDI)=
     +Y(JCNF,I,NATOMKINDI)+CNF(ICNF,JCNF,KCNF,I,NATOMKINDI)
          ENDIF
      ENDIF
      

                    
                          ENDDO
                      ENDDO
                      
                  enddo
                  CLOSE(101)
              ENDDO
      ENDDO
      !---------------------------------------------------------------------------
      !ESCRITURA DE ARCHIVOS
      !Archivo en Z
!---------------------------------------------------------------------------
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
                              write(CONFNAT,'(i1)') NATOMKINDI
          open(102,file='Z'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
                              DO KCNF=KMIN,KMAX
                                  write(102,*)KCNF,Z(KCNF,I,NATOMKINDI)
                              enddo
                              CLOSE(102)
              ENDDO
      ENDDO
      
!--------------------------------------------------------------------------------
      !Archivos X
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
                              write(CONFNAT,'(i1)') NATOMKINDI
          open(102,file='X'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
                              DO iCNF=iMIN,iMAX
                                  write(102,*)iCNF,x(iCNF,I,NATOMKINDI)
                              enddo
                              CLOSE(102)
              ENDDO
      ENDDO
      
!--------------------------------------------------------------------------------
      !Archivos Y
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
                              write(CONFNAT,'(i1)') NATOMKINDI
          open(102,file='Y'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
                              DO JCNF=JMIN,JMAX
                                  write(102,*)JCNF,Y(iCNF,I,NATOMKINDI)
                              enddo
                              CLOSE(102)
              ENDDO
      ENDDO
      
!--------------------------------------------------------------------------------
      !Archivos XY
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
	write(CONFNAT,'(i1)') NATOMKINDI
      open(102,file='XY'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
                  DO ICNF=IMIN,IMAX
                      DO JCNF=JMIN,JMAX
          WRITE(102,*)icnf,jcnf,XY(ICNF,JCNF,I,NATOMKINDI)
                      ENDDO
                      
                  enddo
                  CLOSE(102)
              ENDDO
      ENDDO

!------------------------------------------------------------------------------
      !ARCHIVO XZ
!------------------------------------------------------------------------
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
	write(CONFNAT,'(i1)') NATOMKINDI
      open(102,file='XZ'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
                  DO ICNF=IMIN,IMAX
                      DO KCNF=KMIN,KMAX
          WRITE(102,*)icnf,Kcnf,XZ(ICNF,KCNF,I,NATOMKINDI)
                      ENDDO
                      
                  enddo
                  CLOSE(102)
              ENDDO
      ENDDO
!------------------------------------------------------------------------------
      !ARCHIVO YZ
!------------------------------------------------------------------------
      Do I=1,nmolec
          	CONFAT='CONFIG'
	write(CONFAT,'(i1)') i
              DO NATOMKINDI=1,NATOM(I)
                            	CONFNAT='CONFIG'
	write(CONFNAT,'(i1)') NATOMKINDI
      open(102,file='YZ'//CONFIG//'-'//CONFAT//'-'//CONFNAT//'.TXT')
                  DO JCNF=IMIN,IMAX
                      DO KCNF=KMIN,KMAX
          WRITE(102,*)Jcnf,Kcnf,YZ(JCNF,KCNF,I,NATOMKINDI)
                      ENDDO
                      
                  enddo
                  CLOSE(102)
              ENDDO
      ENDDO
      ENDDO

      

      
      
      END
      