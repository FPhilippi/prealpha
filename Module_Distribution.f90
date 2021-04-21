
!This Module calculates special combined distribution functions - such as cylindrical or polar.
MODULE DISTRIBUTION ! Copyright (C) !RELEASEYEAR! Frederik Philippi
	USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	INTEGER,PARAMETER :: sampling_interval_default=10
	INTEGER,PARAMETER :: bin_count_default=100
	LOGICAL,PARAMETER :: subtract_uniform_default=.FALSE.
	LOGICAL,PARAMETER :: weigh_charge_default=.FALSE.
	LOGICAL,PARAMETER :: use_COC_default=.FALSE.
	LOGICAL,PARAMETER :: normalise_CLM_default=.FALSE.
	REAL,PARAMETER :: maxdist_default=10.0
	!variables
	CHARACTER (LEN=10) :: operation_mode="NONE"!operation mode of the distribution module.
	INTEGER :: number_of_references !number of different references. '0 0 1' uses the z-axis, etc.
	INTEGER :: bin_count_b=bin_count_default!number of bins in the final 2D combined distribution functions (in "a" direction!)
	INTEGER :: bin_count_a=bin_count_default!number of bins in the final 2D combined distribution functions (in "b" direction!)
	INTEGER :: sampling_interval=sampling_interval_default!every this many steps are used for the distribution functions
	REAL :: maxdist=maxdist_default!only molecules with a distance < maxdist are considered. Should be half the box size.
	REAL :: foolsproof_ratio!for the cdf, consider rare sampling of corners.
	LOGICAL :: subtract_uniform=subtract_uniform_default!subtract the uniform density
	LOGICAL :: weigh_charge=weigh_charge_default!weigh the distribution functions by charges.
	LOGICAL :: use_COC=use_COC_default!use centre of charge
	LOGICAL :: normalise_CLM=normalise_CLM_default!divide charge arm by product of mass and radius of gyration squared.
	INTEGER,DIMENSION(:,:),ALLOCATABLE :: references ! (x y z molecule_type_index_ref molecule_type_index_obs)=1st dim, (nreference)=2nd dim
	REAL,DIMENSION(:,:),ALLOCATABLE :: distribution_function
	REAL :: step_a,step_b !step sizes in angströms or radians
	REAL :: ideal_density
	!PRIVATE/PUBLIC declarations
	PRIVATE :: operation_mode,number_of_references,sampling_interval,references,weigh_charge,subtract_uniform,maxdist
	PRIVATE :: distribution_function
	PUBLIC :: user_distribution_input,perform_distribution_analysis,write_simple_sumrules,write_simple_charge_arm

	CONTAINS

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_distribution'
		SUBROUTINE write_simple_sumrules()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: n,ios
			FILENAME_DISTRIBUTION_INPUT="prealpha_simple.inp"
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),EXIST=file_exists)
			IF (file_exists) CALL report_error(114)
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),IOSTAT=ios)!input path is added for the MSD file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'("pdf ",I0)') give_number_of_molecule_types()
			DO n=1,give_number_of_molecule_types(),1
				WRITE(8,'("+0 +0 +1 ",I0," -1")') n
			ENDDO
			WRITE(8,'("maxdist_optimize")')
			WRITE(8,'("subtract_uniform F")')
			WRITE(8,'("weigh_charge T")')
			WRITE(8,'("sampling_interval 1")')
			WRITE(8,'("bin_count_a 300")')
			WRITE(8,'("bin_count_b 1")')
			WRITE(8,'("quit")')
			CLOSE(UNIT=8)
			IF (VERBOSE_OUTPUT) PRINT *,"File 'prealpha_simple.inp' written."
		END SUBROUTINE write_simple_sumrules

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_distribution'
		SUBROUTINE write_simple_charge_arm(normalise)
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		LOGICAL,INTENT(IN),OPTIONAL :: normalise
		INTEGER :: n,ios
			FILENAME_DISTRIBUTION_INPUT="prealpha_simple.inp"
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),EXIST=file_exists)
			IF (file_exists) CALL report_error(114)
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),IOSTAT=ios)!input path is added for the MSD file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'("charge_arm ",I0)') give_number_of_molecule_types()
			DO n=1,give_number_of_molecule_types(),1
				WRITE(8,'("+0 +0 +1 ",I0)') n
			ENDDO
			WRITE(8,'("maxdist_optimize")')
			WRITE(8,'("subtract_uniform F")')
			WRITE(8,'("sampling_interval 1")')
			WRITE(8,'("bin_count_a 100")')
			WRITE(8,'("bin_count_b 100")')
			IF (PRESENT(normalise)) THEN
				IF (normalise) WRITE(8,'("normalise_clm T")')
			ENDIF
			WRITE(8,'("quit")')
			CLOSE(UNIT=8)
			IF (VERBOSE_OUTPUT) PRINT *,"File 'prealpha_simple.inp' written."
		END SUBROUTINE write_simple_charge_arm

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_distribution'
		SUBROUTINE user_distribution_input(parallelisation_possible,parallelisation_requested,number_of_molecules,filename_distribution)
		IMPLICIT NONE
		LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
		CHARACTER (LEN=*) :: filename_distribution
		LOGICAL :: connected,optimize_distance,separate_bins
		INTEGER,INTENT(IN) :: number_of_molecules
		INTEGER :: allocstatus,deallocstatus,maxmol,n,ios
			parallelisation_possible=.TRUE.
			PRINT *,"Generating distribution function input."
			IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
				PRINT *,"First of all, the calculation of distribution functions benefits from parallelisation."
				PRINT *,"Would you like to turn on parallelisation? (y/n)"
				IF (user_input_logical()) parallelisation_requested=.TRUE.
			ENDIF
			PRINT *,"This feature calculates two-dimensional distribution functions in cylindrical or polar coordinates."
			PRINT *,"The main purpose is to investigate anisotropy, such as channel formation."
			PRINT *,"In addition, it is possible to print the polar distribution function of the charge arm / dipole vector."
			PRINT *,"Thus, the distribution functions are calculated with respect to a specified reference vector."
			PRINT *,"In cylindrical coordinates, the independent variables are"
			PRINT *,"  a) the distance perpendicular to the reference vector and"
			PRINT *,"  b) the distance collinear to the reference vector."
			PRINT *,"In polar coordinates, the independent variables are"
			PRINT *,"  a) the distance between the two molecules in the considered pair and"
			PRINT *,"  b) the angle between reference vector and the vector connecting the molecules."
			PRINT *,"Would you like to request a distribution function in polar (y) or cylindrical (n) coordinates?"
			IF (user_input_logical()) THEN
				PRINT *,"Would you like to calculate the polar distribution function only for the charge arm? (y/n)"
				IF (user_input_logical()) THEN
					operation_mode="charge_arm"
				ELSE
					operation_mode="pdf"
				ENDIF
			ELSE
				operation_mode="cdf"
			ENDIF
			IF (TRIM(operation_mode)=="charge_arm") THEN
				PRINT *,"You now need to specify the combinations of reference vector and the molecule type for the distribution."
			ELSE
				PRINT *,"You now need to specify the combinations of reference vector and the two molecule types for the distribution."
			ENDIF
			PRINT *,"How many of these combinations would you like to specify?"
			number_of_references=user_input_integer(1,100000)
			!Allocate memory for intermediate storage...
			ALLOCATE(references(5,number_of_references),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			maxmol=number_of_molecules
			IF (maxmol<1) maxmol=10000!unknown molecule number... expect the worst.
			!THEN, read references from the standard input.
			DO n=1,number_of_references,1
				WRITE(*,'(" Reading reference number ",I0,"...")') n
				PRINT *,"Please give the molecule type index / number of the reference molecule:"
				PRINT *,"(This molecule will be the origin for the distribution analysis)"
				references(4,n)=user_input_integer(1,maxmol)
				IF (.NOT.(TRIM(operation_mode)=="charge_arm")) THEN
					PRINT *,"It is possible to consider *all* molecule types around this reference molecules."
					PRINT *,"Would you like to use all (y) molecules or only those of a specific type (n)?"
					IF (user_input_logical()) THEN
						references(5,n)=-1
					ELSE
						PRINT *,"Which molecule type would you like to observe / use for the distribution function?"
						references(5,n)=user_input_integer(1,maxmol)
					ENDIF
				ENDIF
				PRINT *,"You now have to choose the reference vector for this molecule type."
				PRINT *,"For example, '0 0 1' is the z direction."
				PRINT *,"If you specify '0 0 0', then the reference vector will be randomised for *every pair*."
				PRINT *,"Please provide the x-component as integer:"
				references(1,n)=user_input_integer(-10,10)
				PRINT *,"Please provide the y-component as integer:"
				references(2,n)=user_input_integer(-10,10)
				PRINT *,"Please provide the z-component as integer:"
				references(3,n)=user_input_integer(-10,10)
			ENDDO
			!At this point, references should be initialised.
			!Thus, continue with reading in the switches:
			PRINT *,"Would you like to change advanced settings?"
			IF (TRIM(operation_mode)=="pdf") THEN
				PRINT *,"These are bin counts, sampling interval, maximum distance, and subtraction of uniform density."
			ELSE
				PRINT *,"These are bin counts, sampling interval, and maximum distance."
			ENDIF
			IF (user_input_logical()) THEN
				PRINT *,"would you like to set the bin counts for a) and b) separately?"
				WRITE(*,'(" (default for bin counts is currently set to ",I0,")")') bin_count_default
				separate_bins=user_input_logical()
				IF (separate_bins) THEN
					PRINT *,"Please enter bin count for independent variable a):"
					bin_count_a=user_input_integer(10,1000)
					PRINT *,"Please enter bin count for independent variable b):"
					bin_count_b=user_input_integer(10,1000)
				ELSE
					PRINT *,"Please enter bin count (to be used for both independent variables):"
					bin_count_a=user_input_integer(10,1000)
				ENDIF
				PRINT *,"Every how many steps would you like to use to construct the histogram?"
				WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
				sampling_interval=user_input_integer(1,1000)
				IF (TRIM(operation_mode)=="charge_arm") THEN
					PRINT *,"You need to choose a maximum distance - vectors whose length exceed that threshold will be disregarded."
					PRINT *,"It is possible to use the maximum vector length from the first time step. Would you like to use that shortcut? (y/n)"
				ELSE
					PRINT *,"You need to choose a maximum distance - molecule pairs beyond that threshold will be disregarded."
					PRINT *,"It is recommended to use half the box size where available. Would you like to use that shortcut? (y/n)"
				ENDIF
				optimize_distance=user_input_logical()
				IF (.NOT.(optimize_distance)) THEN
					PRINT *,"Please enter the maximum distance."
					maxdist=user_input_real(5.0,1000.0)
				ENDIF
				IF (TRIM(operation_mode)=="pdf") THEN
					PRINT *,"Would you like to subtract the uniform density? (y/n)"
					PRINT *,"(i.e. only show differences to the radial distribution function)"
					subtract_uniform=user_input_logical()
				ENDIF
			ELSE
				optimize_distance=.TRUE.
				CALL set_defaults()
			ENDIF
			IF (TRIM(operation_mode)=="charge_arm") THEN
				PRINT *,"Would you like to normalise the charge arm, using the charge lever moment instead?"
				PRINT *,"(i.e. dividing the charge arm by M*Rgy**2 - see  J. Chem. Phys., 2008, 129, 124507)"
				PRINT *,"(Here, M is the molar mass of the molecule, and Rgy is the radius of gyration)"
				normalise_CLM=user_input_logical()
				weigh_charge=.FALSE.
			ELSE
				PRINT *,"Should the entries for the histogram be weighed by their charge?"
				PRINT *,"(Useful only in combination with 'all surrounding molecules')"
				weigh_charge=user_input_logical()
				PRINT *,"Would you like to use the centres of charge instead of centres of mass?"
				PRINT *,"(Makes only sense if you have ions AND if you have specified atomic charges)"
				use_COC=user_input_logical()
			ENDIF
			WRITE(*,FMT='(A)',ADVANCE="NO") " writing distribution input file..."
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_distribution),IOSTAT=ios)!input path is added for the MSD file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			SELECT CASE (TRIM(operation_mode))
			CASE ("cdf")
				WRITE(8,'(" cdf ",I0," ### cylindrical distribution function.")') number_of_references
			CASE ("pdf")
				WRITE(8,'(" pdf ",I0," ### polar distribution function.")') number_of_references
			CASE ("charge_arm")
				WRITE(8,'(" charge_arm ",I0," ### charge arm distribution.")') number_of_references
			CASE DEFAULT
				CALL report_error(0)
			END SELECT
			DO n=1,number_of_references,1
				IF (TRIM(operation_mode)=="charge_arm") THEN
					WRITE(8,ADVANCE="NO",FMT='(" ",SP,I0," ",I0," ",I0,SS," ",I0," ### ")') references(1:4,n)
					WRITE(8,ADVANCE="NO",FMT='(" molecule type ",I0)') references(4,n)
				ELSE
					WRITE(8,ADVANCE="NO",FMT='(" ",SP,I0," ",I0," ",I0,SS," ",I0," ",I0," ### ")') references(:,n)
					IF (references(5,n)<1) THEN
						WRITE(8,ADVANCE="NO",FMT='(" *all* molecules around ")')
					ELSE
						WRITE(8,ADVANCE="NO",FMT='(" molecules of type ",I0," around ")') references(5,n)
					ENDIF
					WRITE(8,ADVANCE="NO",FMT='(" those of type ",I0)') references(4,n)
				ENDIF
				IF (ALL(references(1:3,n)==0)) THEN
					WRITE(8,'(", randomised reference vector.")')
				ELSE
					WRITE(8,'(".")')
				ENDIF
			ENDDO
			IF (optimize_distance) THEN
				WRITE(8,'(" maxdist_optimize ### optimise maximum distance.")')
			ELSE
				WRITE(8,'(" maxdist ",F0.3," ### set maximum distance")')maxdist
			ENDIF
			IF (subtract_uniform) THEN
				WRITE(8,'(" subtract_uniform T ### subtract anisotropic density - only for pdf")')
			ELSE
				WRITE(8,'(" subtract_uniform F ### do not subtract anisotropic density.")')
			ENDIF
			IF (TRIM(operation_mode)=="charge_arm") THEN
				IF (normalise_CLM) THEN
					WRITE(8,'(" normalise_CLM T ### use charge lever moment")')
				ELSE
					WRITE(8,'(" normalise_CLM F ### use charge arm, no charge lever moment correction")')
				ENDIF
			ELSE
				IF (weigh_charge) THEN
					WRITE(8,'(" weigh_charge T ### weigh with charge - only sensible when *all* molecules are used!")')
				ELSE
					WRITE(8,'(" weigh_charge F ### no charge weighing")')
				ENDIF
				IF (use_COC) THEN
					WRITE(8,'(" center_of_charge T ### use centres of charge - remember to specify atomic charges!")')
				ELSE
					WRITE(8,'(" center_of_charge F ### use centre of mass, as usual")')
				ENDIF
			ENDIF
			WRITE(8,'(" sampling_interval ",I0," ### sample every ",I0," timesteps.")') sampling_interval,sampling_interval
			IF (separate_bins) THEN
				WRITE(8,'(" bin_count_a ",I0," ### set bin counts for a) to ",I0,".")') bin_count_a,bin_count_a
				WRITE(8,'(" bin_count_b ",I0," ### set bin counts for b) to ",I0,".")') bin_count_b,bin_count_b
			ELSE
				WRITE(8,'(" bin_count ",I0," ### set bin counts for both variables to ",I0,".")') bin_count_a,bin_count_a
			ENDIF
			WRITE(8,*) "quit"
			WRITE(8,*)
			WRITE(8,*) "This is an input file for the calculation of distribution functions."
			WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
			ENDFILE 8
			CLOSE(UNIT=8)
			WRITE(*,*) "done"
			DEALLOCATE(references,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		END SUBROUTINE user_distribution_input

		SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
		IMPLICIT NONE
			sampling_interval=sampling_interval_default
			bin_count_a=bin_count_default
			bin_count_b=bin_count_default
			subtract_uniform=subtract_uniform_default
			weigh_charge=weigh_charge_default
			maxdist=maxdist_default
			use_COC=use_COC_default
		END SUBROUTINE set_defaults

		!initialises the distribution module by reading the specified input file.
		SUBROUTINE initialise_distribution()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,allocstatus
			! first, check if file exists.
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),EXIST=file_exists)
			IF (file_exists) THEN
				!setting defaults to start with.
				CALL set_defaults()
				IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),"'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTRIBUTION_INPUT),&
				&ACTION='READ',IOSTAT=ios)
				IF (ios/=0) CALL report_error(99,exit_status=ios)
				READ(3,IOSTAT=ios,FMT=*) operation_mode!read the operation mode.
				IF (ios/=0) CALL report_error(99,exit_status=ios)
				!support for synonyms
				IF (TRIM(operation_mode)=="cydf") operation_mode="cdf"
				IF (TRIM(operation_mode)=="podf") operation_mode="pdf"
				!Now read the body of the distribution input file in line with the requested operation mode:
				SELECT CASE (TRIM(operation_mode))
				CASE ("cdf")
					!The 'corners' of the cylinder are sampled more rarely than the middle.
					!Thus, only those pairs with distances < 0.71*maxdist are binned.
					foolsproof_ratio=2.0**(-0.5)
					IF (VERBOSE_OUTPUT) WRITE(*,*) "calculating cylindrical distribution function."
					IF (VERBOSE_OUTPUT) WRITE(*,*) "reading user-specified references."
					CALL read_references()!uses unit 3!!
					IF ((ERROR_CODE)==33) RETURN !UNIT 3 is closed by report_error
					!allocating memory - first, check for sensible input.
					CALL allocate_memory()
					!step sizes in angströms, depending on bin counts and boundary values:
					! a is the distance in the "xy" plane, i.e. perpendicular to the reference vector
					step_a=(maxdist*foolsproof_ratio)/FLOAT(bin_count_a)
					! b is the vertical distance (in "z", with z being collinear with the reference vector).
					! twice size, because this one has to go in both positive and negative direction.
					step_b=(2.0*maxdist*foolsproof_ratio)/FLOAT(bin_count_b)
				CASE ("pdf")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "calculating polar distribution functions."
					IF (VERBOSE_OUTPUT) WRITE(*,*) "reading user-specified references."
					CALL read_references()!uses unit 3!!
					IF ((ERROR_CODE)==33) RETURN !UNIT 3 is closed by report_error
					!allocating memory - first, check for sensible input.
					CALL allocate_memory()
					!step sizes in angströms or radians, depending on bin counts and boundary values:
					! a is the real space distance between the molecules in the considered pair
					step_a=(maxdist)/FLOAT(bin_count_a)
					! b is the polar angle in radians.
					step_b=(pi)/FLOAT(bin_count_b)
				CASE ("charge_arm")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "calculating charge arm polar distribution functions."
					IF (VERBOSE_OUTPUT) WRITE(*,*) "reading user-specified references."
					CALL read_references()!uses unit 3!!
					IF ((ERROR_CODE)==33) RETURN !UNIT 3 is closed by report_error
					!allocating memory - first, check for sensible input.
					CALL allocate_memory()
					!step sizes in angströms or radians, depending on bin counts and boundary values:
					! a is the real space distance between the molecules in the considered pair
					step_a=(maxdist)/FLOAT(bin_count_a)
					! b is the polar angle in radians.
					step_b=(pi)/FLOAT(bin_count_b)
				CASE DEFAULT
					CALL report_error(99)
				END SELECT
				CLOSE(UNIT=3)
			ELSE
				CALL report_error(132)!No input - no output. easy as that.
			ENDIF
			CONTAINS

				SUBROUTINE read_references()!This subroutine is responsible for reading the body of the distribution input file, connected as unit 3.
				IMPLICIT NONE
				INTEGER :: n,inputbin,molecule_index,m
				CHARACTER(LEN=32) :: inputstring
				REAL :: chargearm
				REAL(KIND=WORKING_PRECISION) :: rgy_sq
					IF ((TRIM(operation_mode)=="cdf").OR.(TRIM(operation_mode)=="pdf")) THEN
						!read user-specified references
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,number_of_references
						!Allocate memory to store the references and the molecule_type_index
						ALLOCATE(references(5,number_of_references),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
						!Try to read all the references from the input file.
						DO n=1,number_of_references,1
							READ(3,IOSTAT=ios,FMT=*) references(:,n)
							IF (ios/=0) CALL report_error(99,exit_status=ios)!ERROR 99: incorrect format in distribution.inp
							IF ((references(4,n)>give_number_of_molecule_types()).OR.(references(4,n)<1)) THEN
								!the specified molecule type doesn't exist. Stop execution.
								CALL report_error(33,exit_status=references(4,n))
								CLOSE(UNIT=3)
								RETURN
							ENDIF
							!the last reference can be <1, which means all molecule types have to be considered
							IF (references(5,n)>give_number_of_molecule_types()) THEN
								!the specified molecule type doesn't exist. Stop execution.
								CALL report_error(33,exit_status=references(5,n))
								CLOSE(UNIT=3)
								RETURN
							ENDIF
						ENDDO
					ELSEIF (TRIM(operation_mode)=="charge_arm") THEN
						!read user-specified references
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,number_of_references
						!Allocate memory to store the references and the molecule_type_index
						ALLOCATE(references(4,number_of_references),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
						!Try to read all the references from the input file.
						DO n=1,number_of_references,1
							READ(3,IOSTAT=ios,FMT=*) references(:,n)
							IF (ios/=0) CALL report_error(99,exit_status=ios)!ERROR 99: incorrect format in distribution.inp
							IF ((references(4,n)>give_number_of_molecule_types()).OR.(references(4,n)<1)) THEN
								!the specified molecule type doesn't exist. Stop execution.
								!<1 not allowed here. Too much work later on.
								CALL report_error(33,exit_status=references(4,n))
								CLOSE(UNIT=3)
								RETURN
							ENDIF
						ENDDO
					ELSE !unknown operation mode.
						CALL report_error(0)
					ENDIF
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  End-of-file condition in ",TRIM(FILENAME_DISTRIBUTION_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DISTRIBUTION_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("sampling_interval")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
							IF (ios/=0) THEN
								CALL report_error(100,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'sampling_interval' to default (=",sampling_interval_default,")"
								sampling_interval=sampling_interval_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'sampling_interval' to ",sampling_interval
							ENDIF
						CASE ("bin_count")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputbin
							IF (ios/=0) THEN
								CALL report_error(100,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'bin_count_a' and 'bin_count_b' to default (=",bin_count_default,")"
								bin_count_a=bin_count_default
								bin_count_b=bin_count_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'bin_count_a' and 'bin_count_b' to ",inputbin
								bin_count_a=inputbin
								bin_count_b=inputbin
							ENDIF
						CASE ("bin_count_a")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,bin_count_a
							IF (ios/=0) THEN
								CALL report_error(100,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'bin_count_a' to default (=",bin_count_default,")"
								bin_count_a=bin_count_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'bin_count_a' to ",bin_count_a
							ENDIF
						CASE ("bin_count_b")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,bin_count_b
							IF (ios/=0) THEN
								CALL report_error(100,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'bin_count_b' to default (=",bin_count_default,")"
								bin_count_b=bin_count_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'bin_count_b' to ",bin_count_b
							ENDIF
						CASE ("maxdist")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,maxdist
							IF (ios/=0) THEN
								CALL report_error(100,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'maxdist' to default (=",maxdist_default,")"
								maxdist=maxdist_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3)') "   setting 'maxdist' to ",maxdist
							ENDIF
						CASE ("maxdist_optimize")
							IF (TRIM(operation_mode)=="charge_arm") THEN
								maxdist=0.0
								DO m=1,number_of_references,1
									DO molecule_index=1,give_number_of_molecules_per_step(references(4,m)),1
										chargearm=SQRT(SUM(charge_arm(1,references(4,m),molecule_index)**2))
										IF (normalise_CLM) THEN
											CALL compute_squared_radius_of_gyration(1,references(4,m),molecule_index,rgy_sq)
											chargearm=chargearm/(give_mass_of_molecule(references(4,m))*rgy_sq)
										ENDIF
										IF (chargearm>maxdist) maxdist=chargearm
									ENDDO
								ENDDO
								IF (normalise_CLM) THEN
									maxdist=maxdist*100000.0
									maxdist=FLOAT(CEILING(maxdist))
									maxdist=maxdist/100000.0
								ELSE
									maxdist=maxdist*10.0
									maxdist=FLOAT(CEILING(maxdist))
									maxdist=maxdist/10.0
								ENDIF
								IF (maxdist<0.0d0) THEN
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3)')&
									&"  setting 'maxdist' to default (=",maxdist_default,")"
									maxdist=maxdist_default
								ELSEIF (maxdist<1.0d0) THEN
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,E10.3,A)') "   setting 'maxdist' to",maxdist," (based on maximum charge arm)"
								ELSE
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3,A)') "   setting 'maxdist' to ",maxdist," (based on maximum charge arm)"
								ENDIF
							ELSE
								maxdist=give_box_limit()/2.0
								IF (maxdist<0.0d0) THEN
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3)')&
									&"  setting 'maxdist' to default (=",maxdist_default,")"
									maxdist=maxdist_default
								ELSE
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3,A)') "   setting 'maxdist' to ",maxdist," (half the box length)"
								ENDIF
							ENDIF
						CASE ("subtract_uniform")
							IF ((TRIM(operation_mode)=="pdf").OR.(TRIM(operation_mode)=="charge_arm") ) THEN
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,subtract_uniform
								IF (ios/=0) THEN
									CALL report_error(100,exit_status=ios)
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
									&"  setting 'subtract_uniform' to default (=",subtract_uniform_default,")"
									subtract_uniform=subtract_uniform_default
								ELSE
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'subtract_uniform' to ",subtract_uniform
								ENDIF
							ELSE
								WRITE(*,*) "  subtract_uniform is only available for polar distribution functions."
							ENDIF
						CASE ("weigh_charge")
							IF (TRIM(operation_mode)=="charge_arm") THEN
								WRITE(*,*) "charge weighing not available for charge_arm analysis."
							ELSE
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,weigh_charge
								IF (ios/=0) THEN
									CALL report_error(100,exit_status=ios)
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') "   setting 'weigh_charge' to default (=",weigh_charge_default,")"
									weigh_charge=weigh_charge_default
								ELSE
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'weigh_charge' to ",weigh_charge
								ENDIF
							ENDIF
						CASE ("normalise_CLM","normalise_charge_arm","normalize_CLM","normalize_charge_arm","normalise_clm")
							IF (TRIM(operation_mode)=="charge_arm") THEN
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,normalise_CLM
								IF (ios/=0) THEN
									CALL report_error(100,exit_status=ios)
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
									&"  setting 'normalise_CLM' to default (=",normalise_CLM_default,")"
									normalise_CLM=normalise_CLM_default
								ELSE
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'normalise_CLM' to ",normalise_CLM
								ENDIF
							ELSE
								WRITE(*,*) "The charge lever moment normalisation is only available for charge_arm analysis."
							ENDIF
						CASE ("use_coc","use_COC","center_of_charge","centre_of_charge")
							IF (TRIM(operation_mode)=="charge_arm") THEN
								WRITE(*,*) "center_of_charge not available for charge_arm analysis."
							ELSE
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,use_COC
								IF (ios/=0) THEN
									CALL report_error(100,exit_status=ios)
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') "   setting 'use_COC' to default (=",use_COC_default,")"
									use_COC=use_COC_default
								ELSE
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'use_COC' to ",use_COC
								ENDIF
							ENDIF
						CASE ("quit")
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DISTRIBUTION_INPUT)
							EXIT
						CASE DEFAULT
							IF (VERBOSE_OUTPUT) WRITE(*,*) "  can't interpret line - continue streaming"
						END SELECT
						!check bin counts
						IF (bin_count_a<1) THEN
							bin_count_a=1
							CALL report_error(104,exit_status=1)
						ELSEIF (bin_count_a>1000) THEN
							bin_count_a=1000
							CALL report_error(104,exit_status=1000)
						ENDIF
						IF (bin_count_b<1) THEN
							bin_count_b=1
							CALL report_error(104,exit_status=1)
						ELSEIF (bin_count_b>1000) THEN
							bin_count_b=1000
							CALL report_error(104,exit_status=1000)
						ENDIF
					ENDDO
				END SUBROUTINE read_references

				SUBROUTINE allocate_memory()
					ALLOCATE(distribution_function(bin_count_a,bin_count_b),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				END SUBROUTINE allocate_memory

		END SUBROUTINE initialise_distribution

		!finalises the distribution module.
		SUBROUTINE finalise_distribution()
		IMPLICIT NONE
		INTEGER :: deallocstatus
			DEALLOCATE(references,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(distribution_function,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		END SUBROUTINE finalise_distribution

		!This subroutine generates the required distribution functions.
		SUBROUTINE make_distribution_functions(reference_number)
		USE ANGLES
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER,INTENT(IN) :: reference_number
		INTEGER(KIND=DP) :: distribution_histogram(bin_count_a,bin_count_b)
		INTEGER(KIND=DP) :: within_bin!how many pairs are in the bin?
		INTEGER(KIND=DP) :: out_of_bounds!how many pairs have exceeded the bounds? should be zero for pdf.
		INTEGER :: molecule_type_index_ref,molecule_type_index_obs,timestep
		LOGICAL :: all_surrounding_molecules!IF TRUE, THEN all surrounding molecules are considered, of any type.
		LOGICAL :: randomise!IF TRUE, THEN the reference vector is randomised each step
		LOGICAL :: aligned !IF TRUE, THEN reference vector is the z direction already
		LOGICAL :: use_charge_arm
		REAL :: maxdist_squared
		REAL(KIND=WORKING_PRECISION) :: reference_vector(3)
		REAL :: shift !vertical distance shift for cydf
			shift=maxdist*foolsproof_ratio
			out_of_bounds=0
			within_bin=0
			maxdist_squared=maxdist**2
			randomise=(ALL(references(1:3,reference_number)==0))
			IF (randomise) THEN
				CALL srand(86456)
				IF (VERBOSE_OUTPUT) WRITE(*,*) "  randomising reference vector."
			ELSE
				reference_vector(:)=FLOAT(references(1:3,reference_number))
				CALL normalize3D(reference_vector)
				CALL prepare_rotation(reference_vector,(/0.0d0,0.0d0,1.0d0/),aligned=aligned)
				IF (VERBOSE_OUTPUT) THEN
					IF (aligned) THEN
						WRITE(*,*) "  using z-direction as reference vector."
					ELSE
						WRITE(*,*) "  prepared rotation matrix for mapping onto reference vector."
					ENDIF
				ENDIF
			ENDIF
			!the molecule_type_index can be initialised from the references array:
			molecule_type_index_ref=references(4,reference_number)
			IF (TRIM(operation_mode)=="charge_arm") THEN
				all_surrounding_molecules=.FALSE.
				use_charge_arm=.TRUE.
			ELSE
				use_charge_arm=.FALSE.
				IF (references(5,reference_number)<1) THEN
					!invalid molecule type index, i.e. look at all available!
					all_surrounding_molecules=.TRUE.
				ELSE
					molecule_type_index_obs=references(5,reference_number)
					all_surrounding_molecules=.FALSE.
				ENDIF
			ENDIF
			CALL iterate_timesteps_parallelised(randomise)
			CALL transfer_histogram()

			CONTAINS

				SUBROUTINE iterate_timesteps_parallelised(random_debug)
				IMPLICIT NONE
				INTEGER :: observed_type
				INTEGER(KIND=DP) :: within_bin_local,out_of_bounds_local
				INTEGER(KIND=DP) :: distribution_histogram_local(bin_count_a,bin_count_b)
				LOGICAL,INTENT(IN) :: random_debug
					!initialise the histogram
					distribution_histogram(:,:)=0
					!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
					!$OMP PRIVATE(observed_type,distribution_histogram_local,within_bin_local,out_of_bounds_local)&
					!$OMP PRIVATE(reference_vector)
					!$OMP SINGLE
					!$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
					!$ 	WRITE(*,'(A,I0,A)') "   ### Parallel execution on ",OMP_get_num_threads()," threads (distribution function)"
					!$ 	CALL timing_parallel_sections(.TRUE.)
					!$ ENDIF
					IF (VERBOSE_OUTPUT) CALL print_progress&
					&(MAX((give_number_of_timesteps()-1+sampling_interval)/sampling_interval,0))
					!$OMP END SINGLE
					distribution_histogram_local(:,:)=0
					observed_type=molecule_type_index_obs
					!$OMP DO SCHEDULE(STATIC,1)
					DO timestep=1,give_number_of_timesteps(),sampling_interval
						IF (use_charge_arm) THEN
							CALL fill_histogram_charge_arm&
							&(timestep,within_bin_local,out_of_bounds_local,distribution_histogram_local,random_debug)
							!$OMP ATOMIC
							within_bin=within_bin+within_bin_local
							!$OMP ATOMIC
							out_of_bounds=out_of_bounds+out_of_bounds_local
						ELSE
							IF (all_surrounding_molecules) THEN
								DO observed_type=1,give_number_of_molecule_types(),1
									CALL fill_histogram&
									&(timestep,within_bin_local,out_of_bounds_local,distribution_histogram_local,random_debug,observed_type)
									!$OMP ATOMIC
									within_bin=within_bin+within_bin_local
									!$OMP ATOMIC
									out_of_bounds=out_of_bounds+out_of_bounds_local
								ENDDO
							ELSE
								CALL fill_histogram&
								&(timestep,within_bin_local,out_of_bounds_local,distribution_histogram_local,random_debug,observed_type)
								!$OMP ATOMIC
								within_bin=within_bin+within_bin_local
								!$OMP ATOMIC
								out_of_bounds=out_of_bounds+out_of_bounds_local
							ENDIF
						ENDIF
						!$OMP CRITICAL
						IF (VERBOSE_OUTPUT) CALL print_progress()
						!$OMP END CRITICAL
					ENDDO
					!$OMP END DO
					!$OMP CRITICAL
					distribution_histogram(:,:)=distribution_histogram(:,:)+distribution_histogram_local(:,:)
					!$OMP END CRITICAL
					!$OMP END PARALLEL
					IF (((MAX((give_number_of_timesteps()-1+sampling_interval)/sampling_interval,0))>100)&
					&.AND.(VERBOSE_OUTPUT)) WRITE(*,*)
				 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
				 !$ 	WRITE(*,ADVANCE="NO",FMT='("   ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
				 !$ ENDIF
				END SUBROUTINE iterate_timesteps_parallelised

				SUBROUTINE fill_histogram&
				&(timestep_in,within_bin_local,out_of_bounds_local,distribution_histogram_local,random_debug,observed_type)
				IMPLICIT NONE
				INTEGER :: binpos_a,binpos_b,molecule_index_ref,molecule_index_obs
				INTEGER,INTENT(IN) :: timestep_in,observed_type
				REAL :: current_distance_squared,a,b
				REAL(KIND=WORKING_PRECISION) :: link_vector(3),local_reference(3)
				INTEGER(KIND=DP),INTENT(INOUT) :: distribution_histogram_local(:,:)
				INTEGER(KIND=DP),INTENT(OUT) :: within_bin_local,out_of_bounds_local
				LOGICAL,INTENT(IN) :: random_debug
					within_bin_local=0
					out_of_bounds_local=0
					!iterate over all reference molecules
					DO molecule_index_ref=1,give_number_of_molecules_per_step(molecule_type_index_ref),1
						!iterate over all observed molecules
						DO molecule_index_obs=1,give_number_of_molecules_per_step(observed_type),1
							IF ((molecule_type_index_ref==observed_type).AND.(molecule_index_ref==molecule_index_obs)) CYCLE
							!check if the two molecules are within distance
							current_distance_squared=give_smallest_distance_squared&
							&(timestep_in,timestep_in,molecule_type_index_ref,observed_type,molecule_index_ref,molecule_index_obs,&
							&translation=link_vector(:))
							IF (current_distance_squared<maxdist_squared) THEN
								!molecule pair is close enough.
								IF (use_COC) THEN
									link_vector(:)=link_vector(:)+&
									&give_center_of_charge(timestep_in,observed_type,molecule_index_obs)-&
									&give_center_of_charge(timestep_in,molecule_type_index_ref,molecule_index_ref)
								ELSE
									link_vector(:)=link_vector(:)+&
									&give_center_of_mass(timestep_in,observed_type,molecule_index_obs)-&
									&give_center_of_mass(timestep_in,molecule_type_index_ref,molecule_index_ref)
								ENDIF
								!Calculate bin limits.
								SELECT CASE (TRIM(operation_mode))
								CASE ("cdf")
									!cylindrical distribution function: a will be distance in plane perpendicular to reference vector,
									! b will be vertical distance. needs to be corrected, because can be negative
									!IF necessary, randomise rotation matrix - expensive!
									IF (random_debug) THEN
										CALL rotate_random(link_vector(:))
									ELSE
										CALL rotate(link_vector(:))
									ENDIF
									a=SQRT(link_vector(1)**2+link_vector(2)**2)
									b=link_vector(3)
									binpos_a=(INT(a/step_a)+1)
									binpos_b=(INT((b+shift)/step_b))+1
									IF ((binpos_b==1).AND.((b+shift)<0)) binpos_b=0
									! a=(binpos_a-1)*step_a+0.5*step_a
									! b=(binpos_b-1)*step_b+0.5*step_b-shift
								CASE ("pdf")
									!If necessary, randomise uvector (debug mostly)
									IF (random_debug) THEN
										local_reference(1)=RAND()-0.5
										local_reference(2)=RAND()-0.5
										local_reference(3)=RAND()-0.5
										CALL normalize3D(local_reference)
									ELSE
										local_reference(:)=reference_vector(:)
									ENDIF
									!polar distribution function: a will be distance, b will be angle in radians
									a=SQRT(current_distance_squared)
									CALL normalize3D(link_vector(:))
									b=ACOS(DOT_PRODUCT(link_vector(:),local_reference(:)))
									binpos_a=INT(a/step_a)+1
									binpos_b=INT(b/step_b)+1
								CASE DEFAULT
									CALL report_error(0)
								END SELECT
								IF (((binpos_a>0).AND.(binpos_a<=bin_count_a)).AND.((binpos_b>0).AND.(binpos_b<=bin_count_b))) THEN
									within_bin_local=within_bin_local+1
									IF (weigh_charge) THEN
										distribution_histogram_local(binpos_a,binpos_b)=distribution_histogram_local(binpos_a,binpos_b)+&
										&give_charge_of_molecule(observed_type)
									ELSE
										distribution_histogram_local(binpos_a,binpos_b)=distribution_histogram_local(binpos_a,binpos_b)+1
									ENDIF
								ELSE
									out_of_bounds_local=out_of_bounds_local+1
								ENDIF
							ENDIF
						ENDDO
					ENDDO
				END SUBROUTINE fill_histogram

				SUBROUTINE fill_histogram_charge_arm&
				&(timestep_in,within_bin_local,out_of_bounds_local,distribution_histogram_local,random_debug)
				IMPLICIT NONE
				INTEGER :: binpos_a,binpos_b,molecule_index
				INTEGER,INTENT(IN) :: timestep_in
				REAL :: current_distance_squared,a,b
				REAL(KIND=WORKING_PRECISION) :: link_vector(3),local_reference(3),rgy_sq
				INTEGER(KIND=DP),INTENT(INOUT) :: distribution_histogram_local(:,:)
				INTEGER(KIND=DP),INTENT(OUT) :: within_bin_local,out_of_bounds_local
				LOGICAL,INTENT(IN) :: random_debug
					within_bin_local=0
					out_of_bounds_local=0
					!iterate over all reference molecules
					DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index_ref),1
						!iterate over all observed molecules
						!Calculate bin limits.
						link_vector(:)=charge_arm(timestep_in,molecule_type_index_ref,molecule_index)
						IF (normalise_CLM) THEN
							CALL compute_squared_radius_of_gyration(timestep_in,molecule_type_index_ref,molecule_index,rgy_sq)
							link_vector(:)=link_vector(:)/(give_mass_of_molecule(molecule_type_index_ref)*rgy_sq)
						ENDIF
						!If necessary, randomise uvector (debug mostly)
						IF (random_debug) THEN
							local_reference(1)=RAND()-0.5
							local_reference(2)=RAND()-0.5
							local_reference(3)=RAND()-0.5
							CALL normalize3D(local_reference)
						ELSE
							local_reference(:)=reference_vector(:)
						ENDIF
						!polar distribution function: a will be distance, b will be angle in radians
						a=SQRT(SUM(link_vector(:)**2))
						IF (a<1E-8) THEN !decreased threshold because of CLM correction...
							out_of_bounds_local=out_of_bounds_local+1
						ELSE
							!normalise
							link_vector(:)=link_vector(:)/a
							b=ACOS(DOT_PRODUCT(link_vector(:),local_reference(:)))
							binpos_a=INT(a/step_a)+1
							binpos_b=INT(b/step_b)+1
							IF (((binpos_a>0).AND.(binpos_a<=bin_count_a)).AND.((binpos_b>0).AND.(binpos_b<=bin_count_b))) THEN
								within_bin_local=within_bin_local+1
								distribution_histogram_local(binpos_a,binpos_b)=distribution_histogram_local(binpos_a,binpos_b)+1
							ELSE
								out_of_bounds_local=out_of_bounds_local+1
							ENDIF
						ENDIF
					ENDDO
				END SUBROUTINE fill_histogram_charge_arm

				SUBROUTINE transfer_histogram()
				IMPLICIT NONE
				REAL :: total_volume,chunk_area,chunk_volume
				REAL :: a,b
				REAL :: within_bin_real,share
				INTEGER :: binpos_a,binpos_b
					distribution_function(:,:)=0.0
					IF (within_bin<0) THEN
						CALL report_error(101)
						within_bin_real=SUM(FLOAT(distribution_histogram(:,:)))
					ELSE
						within_bin_real=FLOAT(within_bin)
					ENDIF
					IF (VERBOSE_OUTPUT) THEN
						WRITE(*,'("   within bin ",I0,", out of bounds ",I0,".")') within_bin,out_of_bounds
						share=100.0*(FLOAT(out_of_bounds)/FLOAT(out_of_bounds+within_bin))
						IF (share>1.0) WRITE(*,'("   ",F0.1,"% were within the sphere, but beyond the array bounds.")') share
					ENDIF
					!Compute total volume
					SELECT CASE (TRIM(operation_mode))
					CASE ("cdf")
						total_volume=twopi*(shift**3) !pi*height*radius^2, height is 2*radius
					CASE ("pdf","charge_arm")
						total_volume=(4.0/3.0)*pi*(maxdist**3)
					CASE DEFAULT
						CALL report_error(0)
					END SELECT
					ideal_density=within_bin_real/total_volume
					!correct histogram for ideal density
					distribution_function(:,:)=FLOAT(distribution_histogram(:,:))/ideal_density
					IF (VERBOSE_OUTPUT) THEN
						ideal_density=ideal_density/&
						&FLOAT((give_number_of_timesteps()/sampling_interval)*&
						&give_number_of_molecules_per_step(references(4,reference_number)))
						WRITE(*,'("   ideal density:",ES10.3," particles/Angström^3")')ideal_density
						IF (TRIM(operation_mode)=="charge_arm") &
						&WRITE(*,'("   (",F5.3," particles in the sphere)")')ideal_density*total_volume
						IF (DEVELOPERS_VERSION) WRITE(*,'("  ! total volume:  ",ES9.3," cubic Angströms")') total_volume
					ENDIF
					DO binpos_a=1,bin_count_a,1
						!Compute chunk area / partial integral
						SELECT CASE (TRIM(operation_mode))
						CASE ("cdf")
							! a is the distance in the 'xy' plane, i.e. perpendicular to the reference vector.
							a=(FLOAT(binpos_a)*step_a) !outer limit of a distance
							chunk_area=pi*((a**2)-((a-step_a)**2))!radial and azimuthal part of volume
						CASE ("pdf","charge_arm")
							! a is the real space distance between the molecules.
							a=(FLOAT(binpos_a)*step_a) !outer limit of radius
							chunk_area=((a**3)-((a-step_a)**3))!radial part of volume
							chunk_area=(twopi/3.0)*chunk_area!azimuthal part and constants
							a=(a-(step_a*0.5)) !Correct to the middle of the range
						CASE DEFAULT
							CALL report_error(0)
						END SELECT
						DO binpos_b=1,bin_count_b,1
							!Compute chunk volume
							SELECT CASE (TRIM(operation_mode))
							CASE ("cdf")
								! b is the vertical distance (in "z", with z being collinear with the reference vector).
								!The difference will always be delta(b), which is step_b
								chunk_volume=chunk_area*step_b
							CASE ("pdf")
								! b is the polar angle in radians.
								b=FLOAT(binpos_b)*step_b !upper limit
								chunk_volume=chunk_area*(COS(b-step_b)-COS(b)) !lower minus upper for polar part because of minus sign
							CASE ("charge_arm")
								! b is the polar angle in radians.
								b=FLOAT(binpos_b)*step_b !upper limit
								!do NOT correct for the radial part.
								chunk_volume=(COS(b-step_b)-COS(b)) !lower minus upper for polar part because of minus sign
							CASE DEFAULT
								CALL report_error(0)
							END SELECT
							!correct distribution function for chunk volume
							distribution_function(binpos_a,binpos_b)=distribution_function(binpos_a,binpos_b)/chunk_volume
						ENDDO
					ENDDO
					!for charge arm, normalisation is different.
					IF (TRIM(operation_mode)=="charge_arm") distribution_function(:,:)=distribution_function(:,:)/SUM(distribution_function(:,:))
				END SUBROUTINE transfer_histogram

		END SUBROUTINE make_distribution_functions

		!This subroutine is responsible for writing the output into a file.
		SUBROUTINE write_distribution_functions(reference_number)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: reference_number
		REAL :: shift_b,trace,number_integral,energy_integral
		REAL :: a,b,previous
		INTEGER :: binpos_a,binpos_b,ios,reference_charge
		CHARACTER(LEN=1024) :: filename_trace,filename_distribution_output,filename_number_integral,filename_energy_integral
		CHARACTER(LEN=1024) :: headerline,unitline
		LOGICAL :: connected,write_energy_integral
			reference_charge=give_charge_of_molecule(references(4,reference_number))
			write_energy_integral=((reference_charge/=0).AND.(weigh_charge))
			INQUIRE(UNIT=3,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=3)
			INQUIRE(UNIT=4,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=4)
			INQUIRE(UNIT=10,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=10)
			IF (ALL(references(1:3,reference_number)==0)) THEN
				WRITE(filename_distribution_output,'(A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"randomised_"
				WRITE(headerline,'(A)') "randomised vector, "
			ELSE
				WRITE(filename_distribution_output,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"reference_",reference_number,"_"
				WRITE(headerline,'(A,SP,I0," ",I0," ",I0,SS,", ")') "reference vector ",references(1:3,reference_number)
			ENDIF
			IF (TRIM(operation_mode)=="charge_arm") THEN
				IF (normalise_CLM) THEN
					WRITE(filename_distribution_output,'(A,"CLM-normalised_")')&
					&TRIM(filename_distribution_output)
				ENDIF
				!only charge arm
				IF (give_charge_of_molecule(references(4,reference_number))/=0) THEN
					WRITE(filename_distribution_output,'(A,"charge_arm_type_",I0)')&
					&TRIM(filename_distribution_output),references(4,reference_number)
					WRITE(headerline,'(A," charge arm for molecule type ",I0,", ")')&
					&TRIM(headerline),references(4,reference_number)
				ELSE
					WRITE(filename_distribution_output,'(A,"dipole_moment_type_",I0)')&
					&TRIM(filename_distribution_output),references(4,reference_number)
					WRITE(headerline,'(A," dipole moment for molecule type ",I0,", ")')&
					&TRIM(headerline),references(4,reference_number)
				ENDIF
			ELSEIF (references(5,reference_number)<1) THEN
				!all surroundings...
				WRITE(filename_distribution_output,'(A,"all_around_",I0)')&
				&TRIM(filename_distribution_output),references(4,reference_number)
				WRITE(headerline,'(A," all molecules around type ",I0,", ")')&
				&TRIM(headerline),references(4,reference_number)
			ELSE
				WRITE(filename_distribution_output,'(A,"type_",I0,"_around_",I0)')&
				&TRIM(filename_distribution_output),references(5,reference_number),references(4,reference_number)
				WRITE(headerline,'(A," all molecules of type ",I0," around type ",I0,", ")')&
				&TRIM(headerline),references(5,reference_number),references(4,reference_number)
			ENDIF
			!do specific 'pretreatment'
			SELECT CASE (TRIM(operation_mode))
			CASE ("cdf")
				shift_b=maxdist*foolsproof_ratio
				IF (weigh_charge) WRITE(filename_distribution_output,'(A,"_xcharge")') TRIM(filename_distribution_output)
				WRITE(filename_distribution_output,'(A,"_cdf")') TRIM(filename_distribution_output)
				WRITE(headerline,'(A," cylindrical distribution function")') TRIM(headerline)
				WRITE(unitline,'("a(perpendicular_distance) b(collinear_distance) g(a,b)")')
				IF (VERBOSE_OUTPUT) WRITE(*,*) "  writing distribution data (cdf) into file: ",TRIM(filename_distribution_output)
				OPEN(UNIT=3,FILE=TRIM(filename_distribution_output),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
			CASE ("pdf")
				shift_b=0.0
				!Write trace, i.e. the radial distribution function
				WRITE(filename_trace,'(A,"_rdf")') TRIM(filename_distribution_output)
				IF (weigh_charge) THEN
					WRITE(filename_trace,'(A,"_xcharge")') TRIM(filename_trace)
					WRITE(filename_distribution_output,'(A,"_xcharge")') TRIM(filename_distribution_output)
				ENDIF
				WRITE(filename_number_integral,'(A,"_numberintegral")') TRIM(filename_trace)
				WRITE(filename_energy_integral,'(A,"_energyintegral")') TRIM(filename_trace)
				IF (subtract_uniform) WRITE(filename_distribution_output,'(A,"_-trace")') TRIM(filename_distribution_output)
				WRITE(filename_distribution_output,'(A,"_pdf")') TRIM(filename_distribution_output)
				!filenames are complete now.
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,*) "  writing trace (rdf) into file: ",TRIM(filename_trace)
					WRITE(*,*) "  writing number integral into file: ",TRIM(filename_number_integral)
					IF (write_energy_integral) WRITE(*,*) "  writing (Coulomb) energy integral into file: ",&
					&TRIM(filename_energy_integral)
				ENDIF
				OPEN(UNIT=3,FILE=TRIM(filename_trace),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
			
				OPEN(UNIT=10,FILE=TRIM(filename_number_integral),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(3,'(A," radial distribution function (trace of pdf)")') TRIM(headerline)
				WRITE(3,'("distance g(r)")')
				WRITE(10,'(A," number integral INT(4*pi*r^2*rho*g(r))dR")') TRIM(headerline)
				WRITE(10,'("distance n(r)")')
				number_integral=0.0
				IF (write_energy_integral) THEN
					OPEN(UNIT=4,FILE=TRIM(filename_energy_integral),IOSTAT=ios)
					IF (ios/=0) CALL report_error(26,exit_status=ios)
					WRITE(4,'(A," (Coulomb) energy integral INT(z*e^2*r*rho*g(r)/(2*epsilon0))dR")') TRIM(headerline)
					WRITE(4,'("distance E(r)[J/mol]")')
					energy_integral=0.0
					WRITE(4,*) 0.0,0.0
				ENDIF
				WRITE(10,*) 0.0,0.0
				previous=0.0
				DO binpos_a=1,bin_count_a,1
					trace=SUM(distribution_function(binpos_a,:))/FLOAT(bin_count_b)
					a=(FLOAT(binpos_a)*step_a) !upper limit - integral runs over whole section
					!trapezoidal integration
					number_integral=number_integral+number_integral_trapezoid(previous,trace,(a-step_a),a)*ideal_density
					WRITE(10,*) a,number_integral
					IF (write_energy_integral) THEN
						energy_integral=energy_integral+energy_integral_trapezoid(previous,trace,(a-step_a),a)*ideal_density
						WRITE(4,*) a,energy_integral
					ENDIF
					a=a-0.5*step_a !using the middle to report everything else
					WRITE(3,*) a,trace
					previous=trace
					IF (subtract_uniform) THEN
						distribution_function(binpos_a,:)=distribution_function(binpos_a,:)-trace
					ENDIF
				ENDDO
				CLOSE(UNIT=3)
				IF (write_energy_integral) CLOSE(UNIT=4)
				CLOSE(UNIT=10)
				WRITE(unitline,'("a(distance) b(polar_angle) g(a,b)")')
				WRITE(headerline,'(A," polar distribution function")') TRIM(headerline)
				IF (VERBOSE_OUTPUT) WRITE(*,*) "  writing distribution data (pdf) into file: ",TRIM(filename_distribution_output)
				OPEN(UNIT=3,FILE=TRIM(filename_distribution_output),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
			CASE ("charge_arm")
				weigh_charge=.FALSE.
				shift_b=0.0
				!Write trace, i.e. the radial distribution function
				WRITE(filename_trace,'(A,"_rdf")') TRIM(filename_distribution_output)
				IF (subtract_uniform) WRITE(filename_distribution_output,'(A,"_-trace")') TRIM(filename_distribution_output)
				WRITE(filename_distribution_output,'(A,"_pdf")') TRIM(filename_distribution_output)
				!filenames are complete now.
				IF (VERBOSE_OUTPUT) WRITE(*,*) "  writing trace (rdf) into file: ",TRIM(filename_trace)
				OPEN(UNIT=3,FILE=TRIM(filename_trace),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(3,'(A," radial distribution function (trace of pdf)")') TRIM(headerline)
				WRITE(3,'("distance g(r)")')
				DO binpos_a=1,bin_count_a,1
					trace=SUM(distribution_function(binpos_a,:))/FLOAT(bin_count_b)
					a=(FLOAT(binpos_a)*step_a) !upper limit - integral runs over whole section
					a=a-0.5*step_a !using the middle to report everything else
					WRITE(3,*) a,trace
					IF (subtract_uniform) THEN
						distribution_function(binpos_a,:)=distribution_function(binpos_a,:)-trace
					ENDIF
				ENDDO
				CLOSE(UNIT=3)
				WRITE(unitline,'("a(distance) b(polar_angle) g(a,b)")')
				WRITE(headerline,'(A," polar distribution function")') TRIM(headerline)
				IF (VERBOSE_OUTPUT) WRITE(*,*) "  writing distribution data (pdf) into file: ",TRIM(filename_distribution_output)
				OPEN(UNIT=3,FILE=TRIM(filename_distribution_output),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
			CASE DEFAULT
				CALL report_error(0)
			END SELECT
			IF (weigh_charge) WRITE(headerline,'(A,", weighed by charge")') TRIM(headerline)
			WRITE(3,'(A,A,A)') "Output in this file is based on the input file '",TRIM(FILENAME_DISTRIBUTION_INPUT),"':"
			WRITE(3,'(A)') TRIM(headerline)
			WRITE(3,'(A)') TRIM(unitline)
			DO binpos_a=1,bin_count_a,1
				DO binpos_b=1,bin_count_b,1
					a=(FLOAT(binpos_a)*step_a)-0.5*step_a
					b=(FLOAT(binpos_b)*step_b)-0.5*step_b-shift_b
					IF ((TRIM(operation_mode)=="pdf").OR.(TRIM(operation_mode)=="charge_arm")) b=b*degrees
					WRITE(3,*) a,b,distribution_function(binpos_a,binpos_b)
				ENDDO
			ENDDO
			CLOSE(UNIT=3)

			CONTAINS

				REAL FUNCTION number_integral_trapezoid(g_a,g_b,radius_a,radius_b)
				IMPLICIT NONE
				REAL,INTENT(IN) :: g_a,g_b,radius_a,radius_b
				REAL :: m,n,upper,lower
					n=(radius_b*g_a-radius_a*g_b)/(radius_b-radius_a)
					m=(g_b-n)/(radius_b)
					upper=(m/4.0)*(radius_b**4)+(n/3.0)*(radius_b**3)
					lower=(m/4.0)*(radius_a**4)+(n/3.0)*(radius_a**3)
					!evaluate the integral
					number_integral_trapezoid=4.0*pi*(upper-lower)
				END FUNCTION number_integral_trapezoid

				REAL FUNCTION energy_integral_trapezoid(g_a,g_b,radius_a,radius_b)
				IMPLICIT NONE
				REAL,INTENT(IN) :: g_a,g_b,radius_a,radius_b
				REAL :: m,n,upper,lower
					n=(radius_b*g_a-radius_a*g_b)/(radius_b-radius_a)
					m=(g_b-n)/(radius_b)
					upper=(m/3.0)*(radius_b**3)+(n/2.0)*(radius_b**2)
					lower=(m/3.0)*(radius_a**3)+(n/2.0)*(radius_a**2)
					!evaluate the integral
					energy_integral_trapezoid=&
					&(FLOAT(reference_charge)*(elementary_charge**2)/(2.0*vacuum_permittivity))*(upper-lower)*&
					&avogadro*1.0e10!conversion to J/mol
				END FUNCTION energy_integral_trapezoid

		END SUBROUTINE write_distribution_functions

		SUBROUTINE perform_distribution_analysis()
		IMPLICIT NONE
		INTEGER :: n
		CHARACTER(LEN=1024)  :: name_obs
			IF (give_number_of_timesteps()<11) THEN
				!trajectory is too short. Abort analysis.
				CALL report_error(37)
			ELSE
				CALL initialise_distribution()
				IF ((ERROR_CODE/=31).AND.((ERROR_CODE)/=33)) THEN
					!do the actual analysis:
					!user-defined references.
					DO n=1,number_of_references,1
						IF (TRIM(operation_mode)=="charge_arm") THEN
							WRITE(*,'(" Calculating distribution for reference vector ",I0," out of ",I0,".")') n,number_of_references
							IF (give_charge_of_molecule(references(4,n))==0) THEN
								WRITE(*,'("   (Vector ",SP,I0," ",I0," ",I0,SS," vs. dipole moment of ",A,")")')&
								&references(1:3,n),TRIM(give_sum_formula(references(4,n)))
							ELSE
								WRITE(*,'("   (Vector ",SP,I0," ",I0," ",I0,SS," vs. charge arm of ",A,")")')&
								&references(1:3,n),TRIM(give_sum_formula(references(4,n)))
							ENDIF
							IF (.NOT.(check_charges(references(4,n)))) CALL report_error(127)
						ELSE
							IF (references(5,n)<1) THEN
								name_obs="all"
							ELSE
								name_obs=TRIM(give_sum_formula(references(5,n)))
							ENDIF
							WRITE(*,'(" Calculating distribution for reference vector ",I0," out of ",I0,".")') n,number_of_references
							WRITE(*,'("   (Vector ",SP,I0," ",I0," ",I0,SS,", ",A," for ",A," around ",A,")")')&
							&references(1:3,n),TRIM(operation_mode),TRIM(name_obs),TRIM(give_sum_formula(references(4,n)))
						ENDIF
						CALL make_distribution_functions(n)
						CALL write_distribution_functions(n)
					ENDDO
					CALL finalise_distribution()
				ELSE
					ERROR_CODE=ERROR_CODE_DEFAULT
				ENDIF
			ENDIF
		END SUBROUTINE perform_distribution_analysis

END MODULE DISTRIBUTION
!--------------------------------------------------------------------------------------------------------------------------------!