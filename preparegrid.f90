SUBROUTINE preparegrid(gridfile,datafile,db)
USE GLOBAL
USE Interfaces
	IMPLICIT NONE

	!Dummy argument definitions
	CHARACTER(100),INTENT(IN) :: gridfile,datafile
	INTEGER,INTENT(IN) :: db
	
	!Local variable declarations
	INTEGER,PARAMETER :: in=30
	INTEGER :: i,ios
	
	!Open grid file
	OPEN(UNIT=in,FILE=gridfile,STATUS="OLD",IOSTAT=ios)
	IF(ios .NE. 0) THEN
		PRINT '("Invalid prediction locations file!"/)'
		PAUSE
		STOP
	END IF

	!Skip header
	READ(in,*)

	!Prepare the grid vector
	ngrid=0

	!Loop over grid file to count observations
	DO
		READ(in,*,IOSTAT=ios)
		IF(ios .NE. 0) EXIT
		ngrid=ngrid+1
	END DO
	REWIND(in)

	ALLOCATE(grid(ngrid))

	!Skip header
	READ(in,*)

	!Read in grid locations
	DO i=1,ngrid
		READ(in,*) grid(i)%x,grid(i)%y,grid(i)%t,grid(i)%secondary
	END DO

	CLOSE (in)

END SUBROUTINE preparegrid