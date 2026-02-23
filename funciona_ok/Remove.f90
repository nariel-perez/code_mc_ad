SUBROUTINE REMOVE(NLOC, IPULL, MOLKIND)

    USE SimulationData, ONLY: LOCATE, N
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER :: NLOC, IPULL, MOLKIND

    ! Variables locales
    INTEGER :: K

    ! Verificar si el índice NLOC es menor que el número de moléculas de tipo MOLKIND
    IF (NLOC < N(MOLKIND)) THEN

        ! Cerrar el array LOCATE después de la eliminación
        DO K = NLOC + 1, N(MOLKIND)
            LOCATE(K - 1, MOLKIND) = LOCATE(K, MOLKIND)
        END DO

        ! Colocar el átomo fantasma IPULL justo fuera del rango activo del array LOCATE
        LOCATE(N(MOLKIND), MOLKIND) = IPULL
    END IF

    RETURN
END SUBROUTINE REMOVE
