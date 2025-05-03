SUBROUTINE gasdev(harvest)
	IMPLICIT NONE

	REAL,INTENT(INOUT) :: harvest
	REAL :: rsq, v1,v2
	
	CALL RANDOM_SEED()	
	
	DO
		CALL RANDOM_NUMBER(v1)
		CALL RANDOM_NUMBER(V2)
		v1 = 2.0*v1-1.0
		v2 = 2.0*v2-1.0
		rsq=v1**2 + v2**2
		IF (rsq>0.0 .AND. rsq<1.0) EXIT
	END DO
	
	rsq=SQRT(-2.0*LOG(rsq)/rsq)
	harvest=v1*rsq
	
END SUBROUTINE gasdev