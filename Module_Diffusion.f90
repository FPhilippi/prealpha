
!This Module calculates (drift corrected) mean squared displacements.
!The diffusion implementation is shit. Use TRAVIS.
MODULE DIFFUSION ! Copyright (C) !RELEASEYEAR! Frederik Philippi
	USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	INTEGER,PARAMETER :: tmax_default=10000
	INTEGER,PARAMETER :: tstep_default=1
	LOGICAL,PARAMETER :: verbose_print_default=.FALSE.
	!variables
	CHARACTER (LEN=8) :: operation_mode="NONE"!operation mode of the diffusion module.
	INTEGER :: number_of_projections !number of different projections. '1 1 1' is the normal, three-dimensional self-diffusion, '0 0 1' would be along the z-axis, etc.
	INTEGER :: tmax=tmax_default!max number of timesteps into the future for the mean-squared displacement. Default is 10000
	INTEGER :: tstep=tstep_default!resolution of the diffusion functions x_num, x_squared and x_unsquared
	INTEGER,DIMENSION(:,:),ALLOCATABLE :: projections ! x y z molecule_type_index (2ns molecule_type_index)
    REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: x_squared !mean squared displacement, dimension is the time shift given in timesteps
	REAL(KIND=WORKING_PRECISION),DIMENSION(:,:),ALLOCATABLE :: x_unsquared !mean displacement (=drift), first dimension = time shift, remaining dimensions = xyz of drift
    INTEGER(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: x_num !number of averages taken for the positions
	INTEGER :: msd_exponent=2 !the exponent... 2 is for MSD, 4 can be used to manually construct the non-gaussian parameter
	LOGICAL :: verbose_print=verbose_print_default
	!PRIVATE/PUBLIC declarations
	PRIVATE :: operation_mode,initialise_diffusion,number_of_projections,tmax,projections,tstep
	PRIVATE :: finalise_diffusion,make_diffusion_functions,x_num,x_squared,x_unsquared
	PRIVATE :: tmax_default,tstep_default,verbose_print_default,verbose_print
	PRIVATE :: write_diffusion_functions
	PUBLIC :: perform_diffusion_analysis,print_helpful_info,user_msd_input,write_simple_diffusion,user_cross_input

	CONTAINS

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_distribution'
		SUBROUTINE write_simple_diffusion()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,tstep_local
			FILENAME_DIFFUSION_INPUT="prealpha_simple.inp"
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DIFFUSION_INPUT),EXIST=file_exists)
			IF (file_exists) CALL report_error(114)
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DIFFUSION_INPUT),IOSTAT=ios)!input path is added for the MSD file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			tstep_local=(give_number_of_timesteps()-1)/10000
			IF (tstep_local<5) tstep_local=1
			WRITE(8,'("default")')
			WRITE(8,'("tmax ",I0)') give_number_of_timesteps()-1
			WRITE(8,'("tstep ",I0)') tstep_local
			WRITE(8,'("print_verbose F")')
			WRITE(8,'("quit")')
			CLOSE(UNIT=8)
			IF (VERBOSE_OUTPUT) PRINT *,"File 'prealpha_simple.inp' written."
		END SUBROUTINE write_simple_diffusion

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_self_diffusion'
		SUBROUTINE user_msd_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_msd)
		IMPLICIT NONE
		LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
		CHARACTER (LEN=*) :: filename_msd
		LOGICAL :: use_default,printdrift,connected
		INTEGER,INTENT(IN) :: number_of_molecules,nsteps
		INTEGER :: nprojections,allocstatus,deallocstatus,maxmol,tstep,tmax,n,ios
	!	INTEGER,DIMENSION(:,:),ALLOCATABLE :: projections ! x y z molecule_type_index, just as in module DIFFUSION
			parallelisation_possible=.TRUE.
			PRINT *,"Generating MSD input."
			IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
				PRINT *,"First of all, the calculation of diffusion functions benefits from parallelisation."
				PRINT *,"Would you like to turn on parallelisation? (y/n)"
				IF (user_input_logical()) parallelisation_requested=.TRUE.
			ENDIF
			PRINT *,"This feature has the capability of calculating <x²> and <x> for different projections,"
			PRINT *,"which is necessary for systems with anisotropic diffusion."
			PRINT *,"(To name an example, if just the diffusion in z-direction is required, separated from x and y)"
			PRINT *,"It is possible to just calculate the 3D diffusion coefficients for every molecule."
			PRINT *,"Do you want to take this shortcut? (y/n)"
			use_default=user_input_logical()
			IF (use_default) THEN
				!shortcut - write default.
				PRINT *,"Program will use defaults - no need to specify projection. Rather robust, too."
			ELSE
				PRINT *,"Please enter the number of projections you want to specify (including molecule types)."
				nprojections=user_input_integer(1,100000)
				!Allocate memory for intermediate storage...
				ALLOCATE(projections(4,nprojections),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				maxmol=number_of_molecules
				IF (maxmol<1) maxmol=10000!unknown molecule number... expect the worst.
				!Then, read projections from the standard input.
				DO n=1,nprojections,1
					WRITE(*,'(" Reading projection number ",I0,"...")') n
					PRINT *,"Please give the molecule type index / number of the molecule to observe:"
					projections(4,n)=user_input_integer(1,maxmol)
					PRINT *,"You now have to choose the projection for this molecule type."
					PRINT *,"'1 1 1' is 3D, '0 0 1' is only in z-direction, '1 1 0' is in the x-y-plane, etc."
					PRINT *,"Please provide the x-component as integer:"
					projections(1,n)=user_input_integer(-10,10)
					PRINT *,"Please provide the y-component as integer:"
					projections(2,n)=user_input_integer(-10,10)
					PRINT *,"Please provide the z-component as integer:"
					projections(3,n)=user_input_integer(-10,10)
				ENDDO
			ENDIF
			PRINT *
			!At this point, projections should be initialised (if necessary!)
			!Thus, continue with reading in the switches:
			PRINT *,"How many steps do you want the shift of the displacement functions to be?"
			WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
			tmax=user_input_integer(10,(nsteps-1))
			PRINT *,"How fine do you want the functions <x²> and <x> to be computed?"
			PRINT *,"For example, when specifying '10', then the displacements are printed"
			PRINT *,"in intervals of 10 timesteps."
			WRITE(*,'(" The default is currently set to ",I0,".")') tstep_default
			tstep=user_input_integer(1,(tmax/10))
			PRINT *,"Would you like to use specify a custom exponent for the displacement? (y/n)"
			PRINT *,"(the default is '2' for the mean squared displacement)"
			IF (user_input_logical()) THEN
				PRINT *,"Please enter the exponent n you would like to use for <x**n>"
				msd_exponent=user_input_integer(0,4)
			ELSE
				msd_exponent=2
			ENDIF
			!tstep and tmax have sensible values.
			PRINT *,"Finally, would you like the detailed drift to be printed? (y/n)"
			printdrift=user_input_logical()
			WRITE(*,FMT='(A32)',ADVANCE="NO") " writing MSD/drift input file..."
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_msd),IOSTAT=ios)!input path is added for the MSD file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			IF (use_default) THEN
				!let the diffusion module take care of all that.
				WRITE(8,*) "default ### calculating the 3D functions <R²> / <R>² for every molecule type"
				!--> thus, no projections required.
			ELSE
				WRITE(8,'(" msd ",I0," ### mean-squared displacement for given number of projections")') nprojections
				DO n=1,nprojections,1
					WRITE(8,'(" ",I0," ",I0," ",I0," ",I0," ### x-y-z projection for molecule type ",I0)') projections(:,n),projections(4,n)
				ENDDO
				!Done writing - projections is no longer needed.
				DEALLOCATE(projections,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ENDIF
			WRITE(8,'(" tmax ",I0," ### maximum time shift of the displacement function")') tmax
			WRITE(8,'(" tstep ",I0)') tstep
			WRITE(8,FMT='(" print_verbose ",L1)',ADVANCE="NO") printdrift
			IF (printdrift) Then
				WRITE(8,'(" ### detailed drift will be printed.")')
			ELSE
				WRITE(8,'(" ### do not print detailed drift.")')
			ENDIF
			IF (msd_exponent/=2) WRITE(8,'(" exponent ",I0," ### custom exponent - computing <x**",I0,">")')&
			&msd_exponent,msd_exponent
			WRITE(8,*) "quit"
			WRITE(8,*)
			WRITE(8,*) "This is an input file for the calculation of mean-squared displacements."
			WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
			ENDFILE 8
			CLOSE(UNIT=8)
			WRITE(*,*) "done"
		END SUBROUTINE user_msd_input

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_cross_diffusion'
		SUBROUTINE user_cross_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_msd)
		IMPLICIT NONE
		LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
		CHARACTER (LEN=*) :: filename_msd
		LOGICAL :: connected
		INTEGER,INTENT(IN) :: number_of_molecules,nsteps
		INTEGER :: nprojections,allocstatus,deallocstatus,maxmol,tstep,tmax,n,ios
	!	INTEGER,DIMENSION(:,:),ALLOCATABLE :: projections ! x y z molecule_type_index, just as in module DIFFUSION
			parallelisation_possible=.TRUE.
			PRINT *,"Generating MSD input."
			IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
				PRINT *,"First of all, the calculation of diffusion functions benefits from parallelisation."
				PRINT *,"Would you like to turn on parallelisation? (y/n)"
				IF (user_input_logical()) parallelisation_requested=.TRUE.
			ENDIF
			PRINT *,"This feature has the capability of calculating <R²> for different projections,"
			PRINT *,"which is necessary for systems with anisotropic diffusion."
			PRINT *,"(To name an example, if just the diffusion in z-direction is required, separated from x and y)"
			PRINT *,"Please enter the number of projections you want to specify (including molecule types)."
			nprojections=user_input_integer(1,100000)
			!Allocate memory for intermediate storage...
			ALLOCATE(projections(5,nprojections),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			maxmol=number_of_molecules
			IF (maxmol<1) maxmol=10000!unknown molecule number... expect the worst.
			!Then, read projections from the standard input.
			DO n=1,nprojections,1
				WRITE(*,'(" Reading projection number ",I0,"...")') n
				PRINT *,"Please give the 1st molecule type index / number of the 1st molecule to observe:"
				projections(4,n)=user_input_integer(1,maxmol)
				PRINT *,"Please give the 2nd molecule type index / number of the 2nd molecule to observe:"
				projections(5,n)=user_input_integer(1,maxmol)
				PRINT *,"You now have to choose the projection for this molecule type."
				PRINT *,"'1 1 1' is 3D, '0 0 1' is only in z-direction, '1 1 0' is in the x-y-plane, etc."
				PRINT *,"Please provide the x-component as integer:"
				projections(1,n)=user_input_integer(-10,10)
				PRINT *,"Please provide the y-component as integer:"
				projections(2,n)=user_input_integer(-10,10)
				PRINT *,"Please provide the z-component as integer:"
				projections(3,n)=user_input_integer(-10,10)
			ENDDO
			PRINT *
			!At this point, projections should be initialised (if necessary!)
			!Thus, continue with reading in the switches:
			PRINT *,"How many steps do you want the shift of the displacement functions to be?"
			WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
			tmax=user_input_integer(10,(nsteps-1))
			PRINT *,"How fine do you want the functions <x²> to be computed?"
			PRINT *,"For example, when specifying '10', then the displacements are printed"
			PRINT *,"in intervals of 10 timesteps."
			WRITE(*,'(" The default is currently set to ",I0,".")') tstep_default
			tstep=user_input_integer(1,(tmax/10))
			PRINT *,"Would you like to use specify a custom exponent for the displacement? (y/n)"
			PRINT *,"(the default is '2' for the mean squared displacement)"
			IF (user_input_logical()) THEN
				PRINT *,"Please enter the exponent n you would like to use for <x**n>"
				msd_exponent=user_input_integer(0,4)
			ELSE
				msd_exponent=2
			ENDIF
			!tstep and tmax have sensible values.
			WRITE(*,FMT='(A32)',ADVANCE="NO") " writing cross MSD input file..."
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_msd),IOSTAT=ios)!input path is added for the MSD file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'(" cross ",I0," ### mean-squared displacement for given number of projections")') nprojections
			DO n=1,nprojections,1
				WRITE(8,'(" ",I0," ",I0," ",I0," ",I0," ",I0," ### x-y-z projection for molecule types ",I0," and ",I0)')&
				&projections(:,n),projections(4,n),projections(5,n)
			ENDDO
			!Done writing - projections is no longer needed.
			DEALLOCATE(projections,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			WRITE(8,'(" tmax ",I0," ### maximum time shift of the displacement function")') tmax
			WRITE(8,'(" tstep ",I0)') tstep
			IF (msd_exponent/=2) WRITE(8,'(" exponent ",I0," ### custom exponent - computing <x**",I0,">")')&
			&msd_exponent,msd_exponent
			WRITE(8,*) "quit"
			WRITE(8,*)
			WRITE(8,*) "R=(N_a+N_b)*x_a*x_b*((R_a(t)-R_b(t))-(R_a(0)-R_b(0)))"
			WRITE(8,*) "with R_a(t)=(1/N_a)*SUM(r_a(t)"
			WRITE(8,*)
			WRITE(8,*) "This is an input file for the calculation of mean-squared displacements."
			WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
			ENDFILE 8
			CLOSE(UNIT=8)
			WRITE(*,*) "done"
		END SUBROUTINE user_cross_input

		!initialises the diffusion module by reading the specified input file.
		SUBROUTINE initialise_diffusion()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,allocstatus
			! first, check if file exists.
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DIFFUSION_INPUT),EXIST=file_exists)
			IF (file_exists) THEN
				!setting defaults to start with.
				CALL set_defaults()
				IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(PATH_INPUT)//TRIM(FILENAME_DIFFUSION_INPUT),"'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DIFFUSION_INPUT),&
				&ACTION='READ',IOSTAT=ios)
				IF (ios/=0) CALL report_error(30,exit_status=ios)
				READ(3,IOSTAT=ios,FMT=*) operation_mode!read the operation mode.
				IF (ios/=0) CALL report_error(30,exit_status=ios)
				!support for synonyms
				IF (TRIM(operation_mode)=="mxd") operation_mode="msd"
				!Now read the body of the diffusion input file in line with the requested operation mode:
				SELECT CASE (TRIM(operation_mode))
				CASE ("cross")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "calculating mean-squared displacements for cross-diffusion."
					CALL read_input_for_cross_diffusion()!uses unit 3!!
					IF ((ERROR_CODE)==33) RETURN !UNIT 3 is closed by report_error
					!allocating memory - first, check for sensible input.
					CALL check_input_and_allocate_memory()
				CASE ("msd")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "calculating mean-squared displacements for self-diffusion."
					IF (VERBOSE_OUTPUT) WRITE(*,*) "reading user-specified projections."
					CALL read_input_for_self_diffusion()!uses unit 3!!
					IF ((ERROR_CODE)==33) RETURN !UNIT 3 is closed by report_error
					!allocating memory - first, check for sensible input.
					CALL check_input_and_allocate_memory()
				CASE ("default")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "calculating mean-squared displacements for self-diffusion."
					CALL read_input_for_self_diffusion()!uses unit 3!!
					IF ((ERROR_CODE)==33) RETURN !UNIT 3 is closed by report_error
					!allocating memory - first, check for sensible input.
					CALL check_input_and_allocate_memory()
					!the only difference between "default" and "msd" is how the projections are obtained, and how the result is printed.
				CASE DEFAULT
					CALL report_error(30)
				END SELECT
				CLOSE(UNIT=3)
			ELSE
				CALL report_error(31)!No input - no output. easy as that.
			ENDIF
			CONTAINS

				SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
				IMPLICIT NONE
					tmax=tmax_default
					tstep=tstep_default
					verbose_print=verbose_print_default
					msd_exponent=2
				END SUBROUTINE set_defaults

				SUBROUTINE check_input_and_allocate_memory()
				IMPLICIT NONE
					!check if tmax, tstep are sensible
					IF ((tmax>(give_number_of_timesteps()-1))) THEN
						tmax=(give_number_of_timesteps()-1)
						CALL report_error(28,exit_status=INT(tmax))
					ENDIF
					IF ((tmax<10)) THEN
						tmax=10
						CALL report_error(28,exit_status=INT(tmax))
					ENDIF
					IF (tstep>(tmax/10)) THEN
						tstep=(tmax/10)
						CALL report_error(35,exit_status=tstep)
					ENDIF
					!allocate memory for diffusion functions
					ALLOCATE(x_num(tmax/tstep),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					ALLOCATE(x_squared(tmax/tstep),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					ALLOCATE(x_unsquared(tmax/tstep,3),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				END SUBROUTINE check_input_and_allocate_memory

				SUBROUTINE read_input_for_cross_diffusion()!This subroutine is responsible for reading the body of the diffusion input file, connected as unit 3.
				IMPLICIT NONE
				INTEGER :: n
				CHARACTER(LEN=32) :: inputstring
					!read user-specified projections
					BACKSPACE 3
					READ(3,IOSTAT=ios,FMT=*) operation_mode,number_of_projections
					!Allocate memory to store the projections and the molecule_type_index
					ALLOCATE(projections(5,number_of_projections),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					!Try to read all the projections from the input file.
					DO n=1,number_of_projections,1
						READ(3,IOSTAT=ios,FMT=*) projections(:,n)
						IF (ios/=0) CALL report_error(34,exit_status=ios)!ERROR 14: incorrect format in diffusion.inp
						IF ((projections(4,n)>give_number_of_molecule_types()).OR.(projections(4,n)<1)) THEN
							!the specified molecule type doesn't exist. Stop execution.
							CALL report_error(33,exit_status=projections(4,n))
							CLOSE(UNIT=3)
							RETURN
						ENDIF
						IF ((projections(5,n)>give_number_of_molecule_types()).OR.(projections(5,n)<1)) THEN
							!the specified molecule type doesn't exist. Stop execution.
							CALL report_error(33,exit_status=projections(4,n))
							CLOSE(UNIT=3)
							RETURN
						ENDIF
						IF (projections(4,n)==projections(5,n)) THEN
							!same molecule type for cross...
							CALL report_error(133,exit_status=projections(4,n))
						ENDIF
					ENDDO
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_DIFFUSION_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DIFFUSION_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("tmax")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tmax' to default (=",tmax_default,")"
								tmax=tmax_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
							ENDIF
						CASE ("tstep")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,tstep
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tstep' to default (=",tstep_default,")"
								tstep=tstep_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tstep' to ",tstep
							ENDIF
						CASE ("exponent")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,msd_exponent
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'exponent' to default (=2)"
								msd_exponent=2
							ELSE
								IF (VERBOSE_OUTPUT) THEN
									WRITE(*,'(A,I0)') " setting 'exponent' to ",msd_exponent
									WRITE(*,ADVANCE="NO",FMT='(" (Will report <|R|^",I0)') msd_exponent
								ENDIF
							ENDIF
							SELECT CASE (msd_exponent)
							CASE (1)
								WRITE(*,'("> - mean displacement)")')
							CASE (2)
								WRITE(*,'("> - mean squared displacement)")')
							CASE (0)
								WRITE(*,'("> - unit projection for debug purposes)")')
							CASE (4)
								WRITE(*,'("> - use with MSD for non-gaussian parameter alpha2)")')
							CASE DEFAULT
								WRITE(*,'("> - exponent of unknown use)")')
							END SELECT
						CASE ("quit")
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DIFFUSION_INPUT)
							EXIT
						CASE DEFAULT
							IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
						END SELECT
					ENDDO
				END SUBROUTINE read_input_for_cross_diffusion


				SUBROUTINE read_input_for_self_diffusion()!This subroutine is responsible for reading the body of the diffusion input file, connected as unit 3.
				IMPLICIT NONE
				INTEGER :: n
				CHARACTER(LEN=32) :: inputstring
					IF (TRIM(operation_mode)=="msd") THEN
						!read user-specified projections
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) operation_mode,number_of_projections
						!Allocate memory to store the projections and the molecule_type_index
						ALLOCATE(projections(4,number_of_projections),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
						!Try to read all the projections from the input file.
						DO n=1,number_of_projections,1
							READ(3,IOSTAT=ios,FMT=*) projections(:,n)
							IF (ios/=0) CALL report_error(34,exit_status=ios)!ERROR 14: incorrect format in diffusion.inp
							IF ((projections(4,n)>give_number_of_molecule_types()).OR.(projections(4,n)<1)) THEN
								!the specified molecule type doesn't exist. Stop execution.
								CALL report_error(33,exit_status=projections(4,n))
								CLOSE(UNIT=3)
								RETURN
							ENDIF
						ENDDO
					ELSE !operation_mode is "default".
						number_of_projections=give_number_of_molecule_types()
						ALLOCATE(projections(4,number_of_projections),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
						!choose (1 1 1) as default projection for every molecule type
						projections(:,:)=1
						projections(4,:)=(/(n,n=1,number_of_projections,1)/)!remember, number_of_projections is equal to the number of molecule types.
					ENDIF
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_DIFFUSION_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DIFFUSION_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("tmax")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tmax' to default (=",tmax_default,")"
								tmax=tmax_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
							ENDIF
						CASE ("tstep")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,tstep
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tstep' to default (=",tstep_default,")"
								tstep=tstep_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tstep' to ",tstep
							ENDIF
						CASE ("print_verbose")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,verbose_print
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') "setting 'verbose_print' to default (=",verbose_print_default,")"
								verbose_print=.FALSE.
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'verbose_print' to ",verbose_print
							ENDIF
						CASE ("exponent")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,msd_exponent
							IF (ios/=0) THEN
								CALL report_error(32,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'exponent' to default (=2)"
								msd_exponent=2
							ELSE
								IF (VERBOSE_OUTPUT) THEN
									WRITE(*,'(A,I0)') " setting 'exponent' to ",msd_exponent
									WRITE(*,ADVANCE="NO",FMT='(" (Will report <|R|^",I0)') msd_exponent
								ENDIF
							ENDIF
							SELECT CASE (msd_exponent)
							CASE (1)
								WRITE(*,'("> - mean displacement)")')
							CASE (2)
								WRITE(*,'("> - mean squared displacement)")')
							CASE (0)
								WRITE(*,'("> - unit projection for debug purposes)")')
							CASE (4)
								WRITE(*,'("> - use with MSD for non-gaussian parameter alpha2)")')
							CASE DEFAULT
								WRITE(*,'("> - exponent of unknown use)")')
							END SELECT
						CASE ("quit")
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DIFFUSION_INPUT)
							EXIT
						CASE DEFAULT
							IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
						END SELECT
					ENDDO
				END SUBROUTINE read_input_for_self_diffusion

		END SUBROUTINE initialise_diffusion

		!finalises the diffusion module.
		SUBROUTINE finalise_diffusion()
		IMPLICIT NONE
		INTEGER :: deallocstatus
			IF (TRIM(operation_mode)=="msd") THEN
				DEALLOCATE(projections,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_num,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_squared,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_unsquared,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ENDIF
			IF (TRIM(operation_mode)=="default") THEN
				DEALLOCATE(projections,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_num,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_squared,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_unsquared,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ENDIF
			IF (TRIM(operation_mode)=="cross") THEN
				DEALLOCATE(projections,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_num,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_squared,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(x_unsquared,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ENDIF
		END SUBROUTINE finalise_diffusion

		!gives information about what is actually calculated and printed.
		SUBROUTINE print_helpful_info()
		IMPLICIT NONE
		CHARACTER(LEN=32) :: msd_exponent_str,power_terminology
			WRITE(msd_exponent_str,'(I0)') msd_exponent
			SELECT CASE (msd_exponent)
			CASE (1)
				power_terminology=""
			CASE (2)
				power_terminology=" squared"
			CASE (3)
				power_terminology=" cube"
			CASE (4)
				power_terminology=" fourth power"
			CASE DEFAULT
				WRITE(power_terminology,'(" power ",I0," ")') msd_exponent
			END SELECT
			! provide the user with information about what will be reported
			WRITE(*,*) "These quantities will be calculated and reported:"
			IF (TRIM(operation_mode)=="cross") THEN
				WRITE(*,*) "   'timeline': number of the timestep * time scaling factor"
				WRITE(*,*) "   '<|R|**"//TRIM(msd_exponent_str)//">': relative mean molecular",TRIM(power_terminology)," displacement"
				WRITE(*,*) "   here, R=(N_a+N_b)*x_a*x_b*((R_a(t)-R_b(t))-(R_a(0)-R_b(0)))"
				WRITE(*,*) "   with R_a(t)=(1/N_a)*SUM(r_a(t)"
			ELSE
				IF (verbose_print) THEN
					WRITE(*,*) "   'timeline':      number of the timestep * time scaling factor"
					WRITE(*,*) "   '<|R|**"//TRIM(msd_exponent_str)//">':  mean",&
					&TRIM(power_terminology)," displacement, not corrected"
					WRITE(*,*) "   '<R>':           average drift of the center of mass, calculated as <R>=SQRT(<x>²+<y>²+<z>²)"
					WRITE(*,*) "   '<|R|**"//TRIM(msd_exponent_str)//">-<R>**"//TRIM(msd_exponent_str)//"': drift corrected mean",&
					&TRIM(power_terminology)," displacement, equals to <(X-drift)**"//TRIM(msd_exponent_str)//">"
					WRITE(*,*) "   'drift_x(y/z)':  the x, y and z components of the drift vector (average over box)"
					WRITE(*,*) "   '#(N)':          number of averages taken to obtain this value"
				ELSE
					WRITE(*,*) "   'timeline': number of the timestep * time scaling factor"
					WRITE(*,*) "   '<|R|**"//TRIM(msd_exponent_str)//">': mean",TRIM(power_terminology)," displacement, not corrected"
					WRITE(*,*) "   '<R>':      average drift of the center of mass, calculated as <R>=SQRT(<x>²+<y>²+<z>²)"
				ENDIF
			ENDIF
			WRITE(*,'(" the time scaling factor is ",I0)') TIME_SCALING_FACTOR
			WRITE(*,'(" averages taken: max = ",I0,", min = ",I0)') (give_number_of_timesteps()-tstep),(give_number_of_timesteps()-tmax)
		END SUBROUTINE print_helpful_info

		!This subroutine generates the required diffusion functions.
		SUBROUTINE make_diffusion_functions(projection_number)
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: current_distance,starting_timestep,molecule_index,molecule_type_index,n,array_pos
		INTEGER,INTENT(IN) :: projection_number
		INTEGER :: number_of_timesteps,allocstatus,nmolecules,deallocstatus
		REAL(KIND=WORKING_PRECISION) :: projektionsvektor(3),vector_clip(3),squared_clip
		REAL(KIND=WORKING_PRECISION),DIMENSION(:,:),ALLOCATABLE :: initial_positions!first dimension: molecule_index, second dimension: the three coordinates.
		!local functions for parallelisation
		REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: x_squared_temp !mean squared displacement, dimension is the time shift given in timesteps
		REAL(KIND=WORKING_PRECISION),DIMENSION(:,:),ALLOCATABLE :: x_unsquared_temp !mean displacement (=drift), first dimension = time shift, remaining dimensions = xyz of drift
		INTEGER(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: x_num_temp !number of averages taken for the positions
		!get #timesteps, so that function doesn't have to be called every time.
			number_of_timesteps=give_number_of_timesteps()
			!the molecule_type_index can be initialised from the projections array:
			molecule_type_index=projections(4,projection_number)
			!get number of molecules
			nmolecules=give_number_of_molecules_per_step(molecule_type_index)
			!initialise the diffusion functions:
			x_num(:)=0
			x_squared(:)=0.0d0
			x_unsquared(:,:)=0.0d0
			!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
			!$OMP PRIVATE(initial_positions,x_num_temp,x_squared_temp,x_unsquared_temp,molecule_index)&
			!$OMP PRIVATE(current_distance,array_pos,vector_clip,squared_clip,projektionsvektor)
			!$OMP SINGLE
			!$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
			!$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (Diffusion / displacement)"
			!$ 	CALL timing_parallel_sections(.TRUE.)
			!$ ENDIF
			IF (VERBOSE_OUTPUT) CALL print_progress(number_of_timesteps)
			!$OMP END SINGLE
			!Allocate memory to store the initial positions
			ALLOCATE(initial_positions(nmolecules,3),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			!then, allocate and initialise the local diffusion functions:
			ALLOCATE(x_num_temp(tmax/tstep),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			ALLOCATE(x_squared_temp(tmax/tstep),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			ALLOCATE(x_unsquared_temp(tmax/tstep,3),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			x_num_temp(:)=0!later, in the loop, the array indices will be addressed with 'current_distance'
			x_squared_temp(:)=0.0d0
			x_unsquared_temp(:,:)=0.0d0
			!outer loop iterates over all the possible starting points
			!tmax can't be larger than #timesteps-1, in which case only the first timestep can be a starting point.
			!$OMP DO SCHEDULE(STATIC,1)
			DO starting_timestep=1,number_of_timesteps,1
				!get the initial positions, so that module MOLECULAR doesn't have to REWIND every time.
				DO molecule_index=1,nmolecules,1
					initial_positions(molecule_index,:)=give_center_of_mass(starting_timestep,molecule_type_index,molecule_index)
				ENDDO
				!inner loop iterates over the time differences for the steps.
				DO current_distance=tstep,tmax,tstep !zeroth entry is implied.
					!The starting_timestep+current_distance shall never exceed the total number of timesteps.
					IF ((starting_timestep+current_distance)>number_of_timesteps) EXIT
					!array has only tmax/tstep members!
					array_pos=current_distance/tstep
					!initialise the clipboard variable for x_unsquared
					vector_clip(:)=0.0d0
					!initialise the clipboard variable for x_squared
					squared_clip=0.0d0
					!increment x_num. It counts how many points were used for averaging. Not necessary, but less error prone than me calculating that.
					x_num_temp(array_pos)=x_num_temp(array_pos)+1
					!calculate difference between starting_timestep+current_distance and starting_timestep
					!iterate over all molecules of the given type
					DO molecule_index=1,nmolecules,1
						!This is the central part of the (self) diffusion analysis.
						!first, the shift of the center of mass of the current molecule is stored in 'projektionsvektor'.
						projektionsvektor=give_center_of_mass(starting_timestep+current_distance,molecule_type_index,molecule_index)!final position...
						projektionsvektor=projektionsvektor-initial_positions(molecule_index,:)!... minus initial position.
						!Then, the projection is applied - if all entries are one, then the 'normal', three-dimensional mean squared displacement is calculated.
						projektionsvektor(:)=DFLOAT(projections(1:3,projection_number))*projektionsvektor(:)
						!add the found distance to x_squared.
						squared_clip=squared_clip+SQRT(SUM(projektionsvektor(:)**2))**msd_exponent !squared_clip collects the quantity x²+y²+z² (or parts thereof)
						!if only the MSD is required, the following line will work just fine, and is faster:
						! --> squared_clip=squared_clip+SUM((projektionsvektor(:))**2)
						!unlike x_squared, x_unsquared is sign-sensitive and has to be collected in a clipboard variable.
						vector_clip=vector_clip+projektionsvektor
					ENDDO
					!Write squared_clip into x_squared, normalise for number of molecules.
					x_squared_temp(array_pos)=x_squared_temp(array_pos)&
					&+squared_clip/DFLOAT(nmolecules)
					!vector_clip now contains the sum of the individual, signed drifts.
					!It has to be normalised to become the drift of this molecule type:
					vector_clip(:)=vector_clip(:)/DFLOAT(nmolecules)
					!Accumulate the drift at this timestep:
					x_unsquared_temp(array_pos,:)=x_unsquared_temp(array_pos,:)+vector_clip(:)
				ENDDO
				!$OMP CRITICAL
				IF (VERBOSE_OUTPUT) CALL print_progress()
				!$OMP END CRITICAL
			ENDDO
			!$OMP END DO
			!update the original functions
			!$OMP CRITICAL
			x_num(:)=x_num(:)+x_num_temp(:)
			x_squared(:)=x_squared(:)+x_squared_temp(:)
			x_unsquared(:,:)=x_unsquared(:,:)+x_unsquared_temp(:,:)
			!$OMP END CRITICAL
			!deallocate temporary functions
			DEALLOCATE(x_num_temp,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(x_squared_temp,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(x_unsquared_temp,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			!$OMP END PARALLEL
			IF ((number_of_timesteps>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
		 !$ ENDIF
			!both x_squared and x_unsquared are normalised by the number of molecules. now, account for the averaging process:
			x_squared(:)=x_squared(:)/DFLOAT(x_num(:))
			DO n=1,3,1
				x_unsquared(:,n)=x_unsquared(:,n)/DFLOAT(x_num(:))
			ENDDO
		END SUBROUTINE make_diffusion_functions

		!This subroutine generates (N_a+N_b)*x_a*x_b*|(R_a(t)-R_b(t))-(R_a(0)-R_b(0))|^2
		!With R_a(t)=(1/N_a)*SUM(r_a(t))
		SUBROUTINE make_cross_diffusion_functions(projection_number)
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: current_distance,starting_timestep,molecule_index,array_pos,molecule_type_index_a,molecule_type_index_b
		INTEGER,INTENT(IN) :: projection_number
		INTEGER :: number_of_timesteps,allocstatus,nmolecules_a,nmolecules_b,deallocstatus
		REAL(KIND=WORKING_PRECISION) :: vector_clip(3),projektionsvektor(3),x_a,x_b
		REAL(KIND=WORKING_PRECISION),DIMENSION(:,:),ALLOCATABLE :: positions!first dimension: timestep, second dimension: the three coordinates.
		!local functions for parallelisation
		REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: x_squared_temp !mean squared displacement, dimension is the time shift given in timesteps
		INTEGER(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: x_num_temp !number of averages taken for the positions
			!get #timesteps, so that function doesn't have to be called every time.
			number_of_timesteps=give_number_of_timesteps()
			molecule_type_index_a=projections(4,projection_number)
			molecule_type_index_b=projections(5,projection_number)
			!get number of molecules
			nmolecules_a=give_number_of_molecules_per_step(molecule_type_index_a)
			nmolecules_b=give_number_of_molecules_per_step(molecule_type_index_b)
			x_a=FLOAT(nmolecules_a)/FLOAT(nmolecules_a+nmolecules_b)
			x_b=FLOAT(nmolecules_b)/FLOAT(nmolecules_a+nmolecules_b)
			!initialise the diffusion functions:
			x_num(:)=0
			x_squared(:)=0.0d0
			!Allocate memory to store the positions
			ALLOCATE(positions(number_of_timesteps,3),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			!initialise
			positions(:,:)=0.0d0
			!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
			!$OMP PRIVATE(molecule_index,vector_clip,x_num_temp,x_squared_temp,current_distance,array_pos)&
			!$OMP PRIVATE(projektionsvektor)
			!$OMP SINGLE
			!$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
			!$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (mean squared RMM displacement)"
			!$ 	CALL timing_parallel_sections(.TRUE.)
			!$ ENDIF
			!$OMP END SINGLE
			!$OMP DO SCHEDULE(STATIC,1)
			DO starting_timestep=1,number_of_timesteps,1
				!Here, the positions for the relative mean molecular diffusion coefficients are prepared
				vector_clip(:)=0.0d0
				DO molecule_index=1,nmolecules_a,1
					vector_clip(:)=vector_clip(:)+give_center_of_mass(starting_timestep,molecule_type_index_a,molecule_index)
				ENDDO
				positions(starting_timestep,:)=vector_clip(:)/DFLOAT(nmolecules_a)
				vector_clip(:)=0.0d0
				DO molecule_index=1,nmolecules_b,1
					vector_clip(:)=vector_clip(:)+give_center_of_mass(starting_timestep,molecule_type_index_b,molecule_index)
				ENDDO
				positions(starting_timestep,:)=positions(starting_timestep,:)-vector_clip(:)/DFLOAT(nmolecules_b) !Calculating R_a(t)-R_b(t)
			ENDDO
			!$OMP END DO
			!$OMP SINGLE
			!$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
			!$ 	WRITE(*,ADVANCE="NO",FMT='(" ### Preparation took ")')
			!$ 	CALL timing_parallel_sections(.FALSE.)
			!$ ENDIF
			!$OMP END SINGLE
			!then, allocate and initialise the local diffusion functions:
			ALLOCATE(x_num_temp(tmax/tstep),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			ALLOCATE(x_squared_temp(tmax/tstep),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			x_num_temp(:)=0!later, in the loop, the array indices will be addressed with 'current_distance'
			x_squared_temp(:)=0.0d0
			!outer loop iterates over all the possible starting points
			!tmax can't be larger than #timesteps-1, in which case only the first timestep can be a starting point.
			!$OMP DO SCHEDULE(STATIC,1)
			DO starting_timestep=1,number_of_timesteps,1
				!inner loop iterates over the time differences for the steps.
				DO current_distance=tstep,tmax,tstep !zeroth entry is implied.
					!The starting_timestep+current_distance shall never exceed the total number of timesteps.
					IF ((starting_timestep+current_distance)>number_of_timesteps) EXIT
					!array has only tmax/tstep members!
					array_pos=current_distance/tstep
					!increment x_num. It counts how many points were used for averaging. Not necessary, but less error prone than me calculating that analytically.
					x_num_temp(array_pos)=x_num_temp(array_pos)+1
					!calculate difference between starting_timestep+current_distance and starting_timestep
					!This is the central part of the (cross) diffusion analysis.
					!first, the shift of the center of mass of the current molecule is stored in 'projektionsvektor'.
					projektionsvektor=positions(starting_timestep+current_distance,:)-positions(starting_timestep,:)!final position minus initial position.
					!Then, the projection is applied - if all entries are one, then the 'normal', three-dimensional mean squared displacement is calculated.
					projektionsvektor(:)=DFLOAT(projections(1:3,projection_number))*projektionsvektor(:)
					!add the found distance to x_squared.
					x_squared_temp(array_pos)=x_squared_temp(array_pos)+SQRT(SUM(projektionsvektor(:)**2))**msd_exponent !squared_clip collects the quantity x²+y²+z² (or parts thereof)
				ENDDO
			ENDDO
			!$OMP END DO
			!update the original functions
			!$OMP CRITICAL
			x_num(:)=x_num(:)+x_num_temp(:)
			x_squared(:)=x_squared(:)+x_squared_temp(:)
			!$OMP END CRITICAL
			!deallocate temporary functions
			DEALLOCATE(x_num_temp,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(x_squared_temp,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			!$OMP END PARALLEL
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
		 !$ ENDIF
			DEALLOCATE(positions,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			!Account for the averaging, and the molar fractions:
			x_squared(:)=x_squared(:)*(FLOAT(nmolecules_a+nmolecules_b)*x_a*x_b/DFLOAT(x_num(:)))
		END SUBROUTINE make_cross_diffusion_functions

		!This subroutine is responsible for writing the output into a file.
		SUBROUTINE write_diffusion_functions(projection_number)
		IMPLICIT NONE
		!Output formats
		2	FORMAT (I15,2E15.6) !Module DIFFUSION - diffusion output file, no verbose_print
		3	FORMAT (I15,E15.6) !Module DIFFUSION - cross diffusion output file, no verbose_print
		6	FORMAT (I15,6E15.6,I15) !Module DIFFUSION - diffusion output file, verbose_print
		INTEGER,INTENT(IN) :: projection_number
		INTEGER :: array_pos,ios,i
		LOGICAL :: connected
		CHARACTER(LEN=1024) :: filename_diffusion_output
		CHARACTER(LEN=32) :: msd_exponent_str
			WRITE(msd_exponent_str,'(I0)') msd_exponent
			INQUIRE(UNIT=3,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=3)
			!generate filename, depending on operation mode
			SELECT CASE (TRIM(operation_mode))
			CASE ("msd")
				WRITE(filename_diffusion_output,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"projection_",projection_number,"_self_diffusion"
			CASE ("default")
				WRITE(filename_diffusion_output,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"type_",projection_number,"_diffusion"
			CASE ("cross")
				WRITE(filename_diffusion_output,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"projection_",projection_number,"_cross_diffusion"
			CASE DEFAULT
				CALL report_error(0)
			END SELECT
			IF (VERBOSE_OUTPUT) WRITE(*,*) "writing diffusion data (msd) into file: ",TRIM(filename_diffusion_output)
			OPEN(UNIT=3,FILE=TRIM(filename_diffusion_output),IOSTAT=ios)
			IF (ios/=0) CALL report_error(26,exit_status=ios)
			!Write header
			WRITE(3,*) "This file contains diffusion information based on the input file '"&
			&,TRIM(FILENAME_DIFFUSION_INPUT),"'"
			IF (TRIM(operation_mode)=="cross") THEN
				WRITE(3,'("   Projection ",I0," (",SP,I0," ",I0," ",I0,SS,"), Molecule types ",I0," (",A,") and ",I0," (",A,")")')&
				&projection_number,projections(1:4,projection_number),&
				&TRIM(give_sum_formula(projections(4,projection_number))),&
				&projections(5,projection_number),&
				&TRIM(give_sum_formula(projections(5,projection_number)))
				WRITE(3,'(2A15)') "timeline","<R**"//TRIM(msd_exponent_str)//">"!timestep*scaling, mean squared displacement
				WRITE(3,3) 0,0.0
			ELSE
				WRITE(3,'("   Projection ",I0," (",SP,I0," ",I0," ",I0,SS,"), Molecule type ",I0," (",A,")")')&
				&projection_number,projections(:,projection_number),TRIM(give_sum_formula(projections(4,projection_number)))
				IF (verbose_print) THEN
					!also print the full drift information
					WRITE(3,'(8A15)') "timeline","<R**"//TRIM(msd_exponent_str)//">",&
					&"<R>","<R**"//TRIM(msd_exponent_str)//">-<R>**"//TRIM(msd_exponent_str),"drift_x","drift_y","drift_z","#(N)"
					WRITE(3,6) 0,0.,0.,0.,0.,0.,0.,0
				ELSE
					!only write the basics
					WRITE(3,'(3A15)') "timeline","<R**"//TRIM(msd_exponent_str)//">","<R>"!timestep*scaling, mean squared displacement and mean displacement
					WRITE(3,2) 0,0.,0.
				ENDIF
			ENDIF
			DO array_pos=1,tmax/tstep,1 !zeroth entry is already there.
				IF (TRIM(operation_mode)=="cross") THEN
					WRITE(3,3) array_pos*TIME_SCALING_FACTOR*tstep,x_squared(array_pos)
				ELSEIF (verbose_print) THEN
					WRITE(3,6) array_pos*TIME_SCALING_FACTOR*tstep,x_squared(array_pos),&!timeline and <R**msd_exponent>
					&SQRT(SUM(x_unsquared(array_pos,:)**2)),&!<R>=SQRT(<x>²+<y>²+<z>²)
					&x_squared(array_pos)-SQRT(SUM(x_unsquared(array_pos,:)**2))**msd_exponent,&!<(R-drift)**msd_exponent>
					&(x_unsquared(array_pos,i),i=1,3),x_num(array_pos)!the drift vector as well as the number of taken averages.
				ELSE
					!only the basics...
					WRITE(3,2) array_pos*TIME_SCALING_FACTOR*tstep,x_squared(array_pos),&!timeline and <R**msd_exponent>
					&SQRT(SUM(x_unsquared(array_pos,:)**2))!The average drift, i.e. the drift of the center of mass of that particular molecule type.
				ENDIF
			ENDDO
			ENDFILE 3
			CLOSE(UNIT=3)
		END SUBROUTINE write_diffusion_functions

		SUBROUTINE perform_diffusion_analysis()
		IMPLICIT NONE
		INTEGER :: n
			IF (give_number_of_timesteps()<11) THEN
				!trajectory is too short. Abort analysis.
				CALL report_error(37)
			ELSE
				CALL initialise_diffusion()
				IF ((ERROR_CODE/=31).AND.((ERROR_CODE)/=33)) THEN
					!do the actual analysis:
					SELECT CASE (TRIM(operation_mode))
					CASE ("msd")
						CALL print_helpful_info()
						!user-defined projections.
						DO n=1,number_of_projections,1
							CALL make_diffusion_functions(n)
							CALL write_diffusion_functions(n)
						ENDDO
					CASE ("default")
						CALL print_helpful_info()
						!default projections, i.e. one per molecule type.
						DO n=1,number_of_projections,1
							CALL make_diffusion_functions(n)
							CALL write_diffusion_functions(n)
						ENDDO
					CASE ("cross")
						CALL print_helpful_info()
						!default projections, i.e. one per molecule type.
						DO n=1,number_of_projections,1
							CALL make_cross_diffusion_functions(n)
							CALL write_diffusion_functions(n)
						ENDDO
					CASE DEFAULT
						CALL report_error(0)
					END SELECT
					CALL finalise_diffusion()
				ELSE
					ERROR_CODE=ERROR_CODE_DEFAULT
				ENDIF
			ENDIF
		END SUBROUTINE perform_diffusion_analysis

END MODULE DIFFUSION
!--------------------------------------------------------------------------------------------------------------------------------!