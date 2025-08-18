MODULE RotationModule
  USE PBC_mod, only: rk

  
  IMPLICIT NONE
  

    ! Constants
    INTEGER, PARAMETER :: TABLA_SIZE = 1000  ! Number of precomputed values
    REAL(rk), PARAMETER :: PI = 3.14159_rk
    REAL(rk), PARAMETER :: PASO = (2.0_rk*PI )/ (TABLA_SIZE - 1)  ! Step size

    ! Lookup tables for sine and cosine values
    REAL(rk), DIMENSION(TABLA_SIZE) :: TABLE_SEN, TABLE_COS

CONTAINS

    !------------------------------------------------------
    ! Initialize the precomputed trigonometric tables
    !------------------------------------------------------
    SUBROUTINE InitRotationTables()
        IMPLICIT NONE
        INTEGER :: I
        REAL(rk) :: ANGULO

        DO I = 1, TABLA_SIZE
            ANGULO = -PI + (REAL(I-1) * PASO)  ! Angle in [-pi , pi]
            TABLE_SEN(I) = SIN(ANGULO)
            TABLE_COS(I) = COS(ANGULO)
        END DO
    END SUBROUTINE InitRotationTables

    !------------------------------------------------------
    ! Compute a rotation matrix using the lookup table
    !------------------------------------------------------
    SUBROUTINE GetRotationMatrix(DX, DY, DZ, R)
        IMPLICIT NONE
        REAL(rk), INTENT(IN) :: DX, DY, DZ
        REAL(rk), DIMENSION(3,3), INTENT(OUT) :: R
        INTEGER :: INDICE_DX, INDICE_DY, INDICE_DZ
        REAL(rk) :: SEN_DX, COS_DX, SEN_DY, COS_DY, SEN_DZ, COS_DZ


        
        ! Compute indices for lookup
        INDICE_DX = INT((DX +PI) / PASO) + 1
        INDICE_DY = INT((DY +PI) / PASO) + 1
        INDICE_DZ = INT((DZ +PI) / PASO) + 1

        ! Get precomputed sine and cosine values
        SEN_DX = TABLE_SEN(INDICE_DX)
        COS_DX = TABLE_COS(INDICE_DX)
        SEN_DY = TABLE_SEN(INDICE_DY)
        COS_DY = TABLE_COS(INDICE_DY)
        SEN_DZ = TABLE_SEN(INDICE_DZ)
        COS_DZ = TABLE_COS(INDICE_DZ)

        ! Compute rotation matrix
        R(1,1) = COS_DZ * COS_DX - SEN_DZ * COS_DY * SEN_DX
        R(1,2) = SEN_DZ * COS_DX + COS_DZ * COS_DY * SEN_DX
        R(1,3) = SEN_DY * SEN_DX

        R(2,1) = -COS_DZ * SEN_DX - SEN_DZ * COS_DY * COS_DX
        R(2,2) = -SEN_DZ * SEN_DX + COS_DZ * COS_DY * COS_DX
        R(2,3) = SEN_DY * COS_DX

        R(3,1) = SEN_DZ * SEN_DY
        R(3,2) = -COS_DZ * SEN_DY
        R(3,3) = COS_DY
        
    END SUBROUTINE GetRotationMatrix

END MODULE RotationModule
