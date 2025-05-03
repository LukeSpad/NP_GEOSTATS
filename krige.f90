SUBROUTINE krige(cvobs,cvest,weights,controlpts,est,var,ck,n,localmu,alpha,beta,db)
	USE global
	USE interfaces
	IMPLICIT NONE

	!Dummy argument declarations
	INTEGER,INTENT(in) :: n,db
	REAL,INTENT(out) :: localmu,est,var,ck,alpha,beta
	REAL,DIMENSION(:),INTENT(in) :: cvest
	REAL,DIMENSION(:),INTENT(out) :: weights
	REAL,DIMENSION(:,:),INTENT(in) :: cvobs
	TYPE(point),DIMENSION(:),INTENT(in) :: controlpts

	!Local variable declarations
	INTEGER :: i,n_negs
	REAL :: mu,lambda_hat,c_hat,r,rem,q_l,q_u
	REAL,DIMENSION(n) :: old
	REAL,DIMENSION(n+method) :: musv
	TYPE(point),DIMENSION(n) :: jk_var
	
	!A subroutine to predict the data value of an unsampled location, using geostatistical
	!methods.

	ef=0

	!Derive Weights
	weights=MATMUL(cvobs,cvest)
		
	!Check weights sum to 1 (the last 2 elements are the lagrange multipliers)
	ck=SUM(weights(1:n))
	IF (NINT(ck) .NE. 1) THEN
		PRINT '(/"ERROR: Kriging weights do not sum to 1!")'
		PRINT '("Check validity of semivariogram model."/)'
		STOP
	END IF

	!Correction for negative weights: Only use for indicator Kriging or when <0 
	!estimates are physically impossible. Correction keeps estimate within data range
	!(i.e. a convex estimator). The Method described in:
	!Deutsch, C.V. (1996). Computers & Geosciences 22(7):765-773
	IF (correction==.TRUE.) THEN
		!Count negative elwments of the weights vector
		n_negs = COUNT(weights(1:n)<0.)
		IF(db==0) old=weights(1:n)
		
		!Calculate mean absolute negative weight and mean covariance of negative weighted points
		lambda_hat = SUM(ABS(PACK(weights(1:n),weights(1:n)<0.)))/n_negs
		c_hat = SUM(PACK(cvest(1:n),weights(1:n)<0.))/n_negs
		
		!Run the correction
		WHERE (weights(1:n) < 0.)
			weights(1:n) = 0.
		ELSEWHERE (cvest(1:n)<c_hat .AND. weights(1:n)<lambda_hat)
			weights(1:n) = 0.
		END WHERE
		
		!Recale weights to sum to 1
		weights(1:n) = weights(1:n)/SUM(weights(1:n))
		IF(db==0) THEN
			PRINT '(/"---------------------------------------------")'
			PRINT '("           *** Weights Correction ***        ")'
			PRINT '("  n       Covariance   Original   Corrected")'
			PRINT '("=============================================")'
			DO i=1,n
				PRINT '(I3,7X,F10.4,5X,F6.4,6X,F6.4)',i,cvest(i),old(i),weights(i)
			END DO
			PRINT '("----------------------------------------------")'
			PRINT '("w/C hat:",2X,F10.4,5X,F6.4)',c_hat,lambda_hat
		END IF
	END IF

	!Estimate the local mean: derive new set of weights from the observtions covariance
	!matrix and an estimation covariance vector which is equal to zero. This filters
	!the residual component of the estimate, returning the local mean.
	IF (method > 0) THEN
		musv = (/ (0.0,i=1,n),cvest(n+1:n+method) /)
		musv=MATMUL(cvobs,musv)

		IF (correction==.TRUE.) THEN
			!Run the correction for negative weights
			WHERE (musv(1:n) < 0.)
				musv(1:n) = 0.
			ELSEWHERE (cvest(1:n)<c_hat .AND. musv(1:n)<lambda_hat)
				musv(1:n) = 0.
			END WHERE
		
			!Recale weights to sum to 1
			musv(1:n) = musv(1:n)/SUM(musv(1:n))
		END IF
		
		localmu=DOT_PRODUCT(musv(1:n),controlpts%value)
	END IF

	!Estimate the regression parameters by modifying the estimateion covariance vector,
	!see Wackernagel's Multivariate geostatistics, pp 290.	
	IF(method==2) THEN
		musv = (/ (0.0,i=1,n),1.0,0.0 /)
		musv=MATMUL(cvobs,musv)
		alpha = DOT_PRODUCT(musv(1:n),controlpts%value)
	
		musv = (/ (0.0,i=1,n),0.0,1.0 /)
		musv=MATMUL(cvobs,musv)
		beta  = DOT_PRODUCT(musv(1:n),controlpts%value)
	END IF

	!Spatial manns to find variance via the jack knife
	DO i=1,n
		CALL mann(1,n-1,CSHIFT(controlpts,i),CSHIFT(cvest,i),CSHIFT(weights,i),est)
		jk_var(i)%value = est
	END DO
	var = SUM((jk_var%value-est)**2)/REAL(n-1)
	
	!Spatial manns test to find estimate
	CALL mann(db,n,controlpts,cvest,weights,est)
	
	!Sample from jacknife sample distribution
	CALL gasdev(r)
	r = (r+2)/4
	q_l = MAX(1.,(REAL(n)*r))
	q_l = MIN(q_l,REAL(n))
	rem = q_l-FLOOR(q_l)
	q_l = FLOOR(q_l)
	q_u = FLOOR(q_u)
	CALL sort(jk_var,jk_var,1)
	est = (rem*jk_var(q_u)%value)+((1-rem)*jk_var(q_l)%value)

	
END SUBROUTINE krige