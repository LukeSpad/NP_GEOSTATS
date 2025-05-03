SUBROUTINE mann(db,n,controlpts,cvest,weights,est)
	USE global
	USE interfaces
	IMPLICIT NONE
	
	!Dummy argument declarations
	INTEGER,INTENT(in) :: db,n
	REAL,INTENT(out) :: est
	REAL,DIMENSION(:),INTENT(in) :: cvest
	REAL,DIMENSION(:),INTENT(in) :: weights
	TYPE(point),DIMENSION(:),INTENT(in) :: controlpts
	
	!Local variable declarations
	INTEGER :: i,j,est_oid
	REAL,DIMENSION(n,n) :: manns
	TYPE(point),DIMENSION(n) :: working 
	REAL :: kendalls,q,denom,rem
	CHARACTER(4) :: nn,ii
	CHARACTER(n+34) :: fmt
	CHARACTER(24+(n*9)) :: fmt2
	CHARACTER(50+(n*9)) :: fmt3,fmt4
	
	!Manns test for spatial structure
	working=controlpts(1:n)
	working%h = sill_gl - cvest
	CALL sort(working,working,0)
	manns=0.
	denom=0.
	DO i=1,(n-1)
		DO j=(i+1),n
			IF (working(i)%value>working(j)%value) THEN
				manns(i,j) = 1 * (weights(i)*weights(j))
			ELSE IF (working(i)%value<working(j)%value) THEN
				manns(i,j) = -1 * (weights(i)*weights(j))
			ELSE
				manns(i,j) = 0
			END IF
			denom = denom + (weights(i)*weights(j))
		END DO
	END DO
	
	kendalls = SUM(manns) / denom
	kendalls = MIN(1.,kendalls)
	kendalls = MAX(-1.,kendalls)
	q = 0.5 + (0.5*kendalls)
		
	!Retrieve ith quantile from the data as the estimate
	CALL sort(working,working,1)
	rem = (REAL(n)*q)-FLOOR(REAL(n)*q)
	est = ABS((rem*working(CEILING(REAL(n)*q))%value) + ((1-rem)*working(FLOOR(REAL(n)*q))%value))
	est_oid = working(NINT(REAL(n)*q))%oid
			 
	!Debug info
	CALL sort(working,working,0)
	WRITE(nn,'(I4)') n-1
	fmt2 = "(4X,'Value',8X,'OID',5X,"//nn//"(I6,3X))"
	fmt3 = "('------------------------',"//nn//"('---------'))"
	fmt4 = "('========================',"//nn//"('========='))"
	IF(db==0) THEN
		PRINT*,""
		PRINT '("Spatial Manns Test:")'
		WRITE (*,fmt3)
		WRITE (*,fmt2) working(2:n)%oid
		WRITE (*,fmt4)
		DO i=1,n
			WRITE(ii,'(I4)') 3+((i-1)*9)
			IF(i < n .AND. working(i)%oid.NE.est_oid) THEN
				fmt  = "(F10.1,4X,I6,':',"//ii//"X,"//nn//"(F7.3,2X))"
				WRITE (*,fmt) working(i)%value,working(i)%oid,manns(i,i+1:n)
				
			ELSE IF(i < n .AND. working(i)%oid==est_oid) THEN
				fmt  = "(F10.1,'***',1X,I6,':',"//ii//"X,"//nn//"(F7.3,2X))"
				WRITE (*,fmt) working(i)%value,working(i)%oid,manns(i,i+1:n)
				
			ELSE IF(i == n .AND. working(i)%oid.NE.est_oid) THEN
				PRINT "(F10.1,4X,I6,':')",working(i)%value,working(i)%oid
				
			ELSE
				PRINT "(F10.1,'***',1X,I6,':')",working(i)%value,working(i)%oid
			END IF	
		END DO
		WRITE (*,fmt3)
		PRINT '("Kendalls Tau: ",F5.2,/"Quantile: ",F5.2)',kendalls,q
	END IF

END SUBROUTINE mann