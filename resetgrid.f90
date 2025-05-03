SUBROUTINE resetgrid(outfile,sim)
	USE global
	USE interfaces
	IMPLICIT NONE

	!Dummy variable definition
	INTEGER,INTENT(in) :: sim
	CHARACTER(100),INTENT(in) :: outfile

	!Local variable declarations
	INTEGER :: ios,i
	CHARACTER(4) :: n
	CHARACTER(100) :: out

	!Convert sim integer value to character, and concatenate with outfile
	WRITE(UNIT=n,FMT=101) sim
	i = 4-LEN_TRIM(ADJUSTL(n))
	n(1:i) = REPEAT('0',i)
	out = TRIM(outfile)//n//'.txt'

	!Initialise nodata values
	grid%value = -999
	grid%sid   = -999
	grid%oid   = -999

	!Randomise grid by sorting on random number vector stored in h
	CALL RANDOM_NUMBER(grid%h)
	CALL sort(grid,grid,0)

	!Close last output file and create new file
	CLOSE(40)
	OPEN(UNIT=40,FILE=out,IOSTAT=ios)
	IF(ios .NE. 0) THEN
		PRINT '("Invalid output file name!"/)'
		PAUSE
		STOP
	END IF

	!Write header on output file
	WRITE(40,FMT=201) TRIM(names(7))
	PRINT '(/"Simulation #",A4)',n

	!Format definitions
	101 FORMAT (I4)
	201 FORMAT (3X,"East",3X,"North",3X,A20,3X,"Time",3X,"Local_Mu",3X,"Estimate",3X,&
				"Variance",3X,"St.Dev."3X,"Alpha",3X,"Beta",3X,"Nearest",3X,"Mu_Dist",3X,&
				"Mu_Drift",3X,"Nearest.ID")

END SUBROUTINE resetgrid