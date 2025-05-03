PROGRAM geostats
	USE IFPORT
	USE global
	USE interfaces
	IMPLICIT NONE

	!****************************************************************************************
	!A program to generate sequential Gaussian simulations from a random function
	!conditioned on a set of data points on a specified set of locations by
	!Kriging. The program requires a set of data points and a semivariogram 
	!model. A grid file must also be provided, containing the space-time coordinates of the
	!required estimation locations, and the values of the external drift variable.
	!A list of coordinates with estimate and estimation variance is returned. The 
	!method returns the best linear unbiased estimator of the process. Estimation proceeds 
	!using the product or product-sum covariance model. Simple Kriging, Ordinary Kriging, 
	!and Kriging with an external drift methods are supported.
	!
	!If compiling on anything other than Compaq Visual Fortran, comment out timing functions
	!and the USE dfport statements (lines 2,78-79 and 193-197).
	!
	!Written by Luke Spadavecchia, November 2006. Last Modified, May 2007.
	!
	!****************************************************************************************

	!Parameter declarations
	INTEGER,PARAMETER :: in=30

	!Variable declarations
	INTEGER :: ios,sim,i,j,ns_h,ns_t,neighbors,window,p,db,cv_flag,candidates,n,ct=0,c1,c2
	INTEGER,DIMENSION(:),ALLOCATABLE :: last
	REAL :: est,var,ck,timer,localmu,alpha,beta
	REAL,DIMENSION(:),ALLOCATABLE :: cvest,weights
	REAL,DIMENSION(:,:),ALLOCATABLE :: cvobs
	CHARACTER(100) :: datafile,outfile
	TYPE(point) :: prediction
	TYPE(point),DIMENSION(:),POINTER :: controlpts,temp
	TYPE(pairs),DIMENSION(:,:),ALLOCATABLE :: hobs
	
	!Copyright notice
	PRINT '("Gaussim (version 1.2): Spatio-Temporal Sequential Gaussian Simulation."//&
    "Copyright (C) 2006  Luke Spadavecchia"//&
    "This program is free software; you can redistribute it and/or modify"/&
    "it under the terms of the GNU General Public License as published by"/&
    "the Free Software Foundation; either version 2 of the License."//&
    "This program is distributed in the hope that it will be useful,"/&
    "but WITHOUT ANY WARRANTY; without even the implied warranty of"/&
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"/&
    "GNU General Public License for more details."//&
    "You should have received a copy of the GNU General Public License along"/&
    "with this program; if not, write to the Free Software Foundation, Inc.,"/&
	"51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA."//)'

	!Get parameters and prepare semivariogram arrays
	CALL getparams(datafile,outfile,db,cv_flag,ns_h,ns_t,neighbors,window,n)

	!Total number of control points to use for kriging:
	window=window+1
	p=neighbors*window
	IF (n<p) THEN
		PRINT '(/"Number of nearest neighbors exceeds number of data points!")'
		PRINT '("Check parameter file!"/)'
		PAUSE
		STOP
	END IF

	!Allocate arrays for kriging
	ALLOCATE(weights(p+method))
	ALLOCATE(cvest(p+method))
	ALLOCATE(controlpts(p))
	ALLOCATE(last(p))
	ALLOCATE(hobs(p,p))
	ALLOCATE(cvobs(p+method,p+method))

	!Timing functions
	PRINT '(/"Running predictions on",1X,I3,1X,"nearest neighbors in space,")',neighbors
	PRINT '("each from a temporal window of +/-",1X,I2,1X,"timesteps.")',window-1
	PRINT '("This is a total of",1X,I3,1X,"observations."/)',p
	PRINT '("Time started:",1X,A8//)',CLOCK()
	timer=SECNDS(0.0)	
		
	!Loop over sim simulations
	DO sim=1,nsim
	
	!Reset grid values and randomise visit schedule: new seed for each simulation
	CALL RANDOM_SEED()
	CALL resetgrid(outfile,sim)

	!Loop through prediction locations
	DO ct=1,ngrid
		!Get information about the prediction datum
		IF(db .NE. 2) THEN
			prediction = grid(ct)
			
			!Select obsevations/Simulations within temporal window, and store in temp array
			c1 = COUNT(observations%t >= prediction%t-window .AND. &
						observations%t <  prediction%t+window)
			c2 = COUNT(grid%t >= prediction%t-window .AND. &
						grid%t <  prediction%t+window .AND. grid%value .NE. -999)
			candidates=c1+c2
			
			ALLOCATE(temp(candidates))
			
			temp(1:c1) = PACK(observations,observations%t >= prediction%t-window .AND. &
								observations%t <  prediction%t+window)
			temp(c1+1:candidates) = PACK(grid,grid%t >= prediction%t-window .AND. &
										 grid%t <  prediction%t+window .AND. & 
										 grid%value .NE. -999)
		ELSE
				PRINT*,"Jack knife mode not supported under simulation"
				PAUSE
				STOP
		END IF

		!Build Estimation distance Vector
		CALL separation(temp,prediction,temp,db)

		!Sort observations to find the nearest neighbors of estimation location, and return
		!the observed/unobserved h vector (controlpts) and observed/observed h matrix (hobs).
		CALL sorted(temp,neighbors,window,p,controlpts,hobs,last)
		
		!Clean up
		DEALLOCATE(temp)

		!Calculate est covariance from semivariance model, using the estimate to data
		!distance vector
		CALL covariance(controlpts%h,controlpts%theta,ABS(controlpts%dt),cvest(1:p),&
						ns_h,ns_t,cv_flag)

		!Subtract (local) mean from control point values if using Simple Kriging: Here
		!the residual component is interpolated, and the known mean added to the result.
		IF (method==0) controlpts%value = controlpts%value - prediction%secondary
		
		!Add extra elements for lagrangian minimisation if method requires
		IF (method>0) cvest(p+1)=1
		
		!Trend component
		IF (method>1) cvest(p+2)=prediction%secondary

		!If control points for this datum are the same as for the last, the next steps
		!can be skipped over: Since distances between the observations are unchanged
		!the observation covariace matrix and its inverse are also unchanged. 
		IF (ANY(last .NE. controlpts%oid)) THEN
			
			!Calculate observation covariance matrix from specified semivariance model, 
			!using the observation distance matrix (hobs).
			CALL covariance2d(hobs%h,hobs%theta,ABS(hobs%dt),cvobs(1:p,1:p),ns_h,ns_t,cv_flag)

			!Add extra rows and columns for lagrangian minimisation
			IF(method>0) THEN
				cvobs(p+1,:)=1
				cvobs(:,p+1)=1
			
				!Trend component
				IF(method>1) THEN
					cvobs(p+2,1:p)=controlpts%secondary
					cvobs(1:p,p+2)=controlpts%secondary
				END IF
			
				!Fill bottom right of matrix zeros 
				cvobs(p+1:p+method,p+1:p+method)=0
			END IF

			!Invert Observation covariance Matrix
			CALL matinv(cvobs,cvobs)
		ELSE
			IF (db == 0) PRINT*,"Control points unchanged: "&
								"Skipping observation covariance matrix inversion!"
		END IF

		!Produce estimate for location
		CALL krige(cvobs,cvest,weights,controlpts,est,var,ck,p,localmu,alpha,beta,db)
		IF (ef==1) PRINT '("Problem coordinate: East =",F11.2", North =",F11.2,", Time ="F7.2)',prediction%x,prediction%y,prediction%t

		!Add (local) mean to estimate if using simple kriging
		IF (method == 0) THEN
			est = est + prediction%secondary
			localmu = prediction%secondary
		END IF

		!Return calculation information if debug mode is selected
		IF (db == 0) &
		CALL debug(weights,controlpts,prediction,est,var,p,cvest,localmu,alpha,beta)

		!Write results to output file
		CALL writeout(prediction,localmu,est,var,alpha,beta,controlpts,p,db)

		!store last set of control points
		last=controlpts%oid

		!Pass values onto grid
		grid(ct)%value = est
	END DO

	END DO

	!Return some timing stats if not in debug mode
	IF (db .NE. 0) THEN
		timer=SECNDS(timer)
		PRINT '(/"Processing took",1X,I3,1X,"Hours",1X,I2,1X,"minutes,",1X,I2,1X,"seconds."/)',&
		FLOOR(timer/3600),FLOOR(timer/60)-(FLOOR(timer/3600)*60),NINT(MOD(timer,60.0))
	END IF

	PRINT '(/"Finished!"/)'
	PAUSE

	!Format definitions
	!Inputs
	102 FORMAT (I1,1X,F6.2,1X,F6.2)

END PROGRAM geostats