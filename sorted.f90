SUBROUTINE sorted(input,ns,nt,p,controlpts,hobs,last)
	USE global
	USE interfaces
	IMPLICIT NONE
	
	!Dummy argument declarations
	INTEGER,INTENT(in) :: ns,nt,p
	INTEGER,DIMENSION(:),INTENT(in) :: last
	TYPE(point),DIMENSION(:),INTENT(inout) :: input
	TYPE(point),DIMENSION(:),INTENT(out) :: controlpts
	TYPE(pairs),DIMENSION(:,:),INTENT(out) :: hobs

	!Local variable declarations
	INTEGER :: i,from,to=0,candidates,nodat=0
	TYPE(point),DIMENSION(:),POINTER :: temp

	!A subroutine to retrieve the nearest neighbors of the estimation location
	!(referred to as control points). The observations separation matrix is
	!returned (as a precursor to producing the observations covariance matrix)
	!along with the vector of nearest neighbbors.

	from=0

	!Loop through temporal separations and retrieve ns spatial nearest neighbors from each 
	DO i=1,nt
		
		!Count number of candidate data points for this temporal separation
		candidates = COUNT(ABS(input%dt)>=i-1 .AND. ABS(input%dt)<i)
	
		!Check enough data exists to fulfil search strategy
		IF(candidates < ns) THEN
			IF (msg==0) THEN
				PRINT '(/"Insufficient data at plus/minus ",I3," time step(s) to fulfil search strategy:")',i-1
				PRINT '("Attempting to fill data deficit from subsequent time steps..."/)'
				msg=1
			END IF
			nodat=nodat+1
			CYCLE
		END IF
		
		!Store data for current temporal separation in temp pointer
		ALLOCATE(temp(candidates))
		temp = PACK(input,ABS(input%dt)>=i-1 .AND. ABS(input%dt)<i)

		!Sort temp by h to find the nearest spatial neighbors
		CALL sort(temp,temp,0)
			
		!Select control points for kriging: the ns nearest spatial neighbors
		to=(ns*i)+(nodat*ns)
		controlpts(from+1:ns*i) = temp(1:ns+(ns*nodat))
		from=ns*i
		
		IF(nodat>0) nodat=0
	
		DEALLOCATE(temp)
	
	END DO
	
	IF(nodat>0) THEN
		PRINT*,"Insufficient data to fulfil search strategy - please check data file."
		STOP
		PAUSE
	END IF

	!If controlpoints the same as last prediction datum, then skip building Hobs
	IF (ANY(last .NE. controlpts%oid)) THEN
	
		!Construct observations distance matrix (hobs) from the control points
		CALL distance(hobs,controlpts,p)
	
	END IF

END SUBROUTINE sorted