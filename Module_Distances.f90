
!This Module calculates intra- and intermolecular contact distance estimates (averages of closest distance)
MODULE DISTANCE ! Copyright (C) !RELEASEYEAR! Frederik Philippi
	USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	REAL,PARAMETER :: maxdist_default=10.0
	INTEGER,PARAMETER :: nsteps_default=1 !how many steps to use from trajectory
	INTEGER,PARAMETER :: sampling_interval_default=1
	LOGICAL,PARAMETER :: calculate_stdev_default=.FALSE.
	LOGICAL,PARAMETER :: calculate_exponential_default=.FALSE.
	LOGICAL,PARAMETER :: calculate_ffc_default=.FALSE.
	!variables
	TYPE,PRIVATE :: distance_subjob
		CHARACTER(LEN=2) :: element1,element2
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: indices_1 !first dimension counts the indices, second dimension is molecule_type_index and atom_index
		INTEGER :: size_indices_1
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: indices_2 !first dimension counts the indices, second dimension is molecule_type_index and atom_index
		INTEGER :: size_indices_2
		INTEGER(KIND=WORKING_PRECISION) :: weight !how many occurrences of this distance per timestep?
		INTEGER(KIND=WORKING_PRECISION) :: missed_neighbours !how many times was the expected weight NOT obtained?
		REAL(KIND=WORKING_PRECISION) :: closest_average
		REAL(KIND=WORKING_PRECISION) :: closest_stdev
		REAL(KIND=WORKING_PRECISION) :: exponential_average
		REAL(KIND=WORKING_PRECISION) :: exponential_stdev
		REAL(KIND=WORKING_PRECISION) :: exponential_weight
		REAL(KIND=WORKING_PRECISION) :: ffc_average
		REAL(KIND=WORKING_PRECISION) :: ffc_weight
		REAL(KIND=WORKING_PRECISION) :: ffcexp_average
		CHARACTER(LEN=512) :: subjob_description
	END TYPE distance_subjob
	TYPE(distance_subjob),DIMENSION(:),ALLOCATABLE :: list_of_distances
	CHARACTER(LEN=14) :: operation_mode="NONE"!operation mode of the distance module.
	INTEGER :: nsteps
	INTEGER :: sampling_interval
	INTEGER :: number_of_subjobs
	LOGICAL :: calculate_stdev=calculate_stdev_default !if TRUE, then the standard deviation is also calculated.
	LOGICAL :: calculate_ffc=calculate_ffc_default !if TRUE, then the equivalents of dHH and rHH are calculated, see the ESI to 10.1021/acs.jpclett.0c00087.
	LOGICAL :: calculate_exponential=calculate_exponential_default !also report distance average weighed with exp(-kr)
	REAL :: maxdist=maxdist_default!only atom pairs with a distance < maxdist are considered. Should not exceed half the box size
	REAL :: maxdist_squared=maxdist_default**2
	REAL :: distance_exponent=1.0 !exponent for distance weighing with exp(-distance_exponent*distance)
	!PRIVATE/PUBLIC declarations
	PUBLIC :: write_simple_intramolecular_distances,write_simple_intermolecular_distances,perform_distance_analysis,user_distance_input

	CONTAINS

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'initialise_distance'
		SUBROUTINE user_distance_input(maxmol,filename_distance)
		IMPLICIT NONE
		CHARACTER (LEN=*) :: filename_distance
		CHARACTER (LEN=2) :: element1,element2
		LOGICAL :: connected
		INTEGER :: n,ios,molecule_type_index,atom_index
		INTEGER,INTENT(IN) :: maxmol
			PRINT *,"Generating input for average distances."

			PRINT *,"This feature reports an average of the distance between reference and observed atoms."
			PRINT *,"the reference and observed atoms can also be element names as wildcards, such as 'H' or 'F'."
			PRINT *,"Intra- or intermolecular distances are possible."
			PRINT *,"Would you like to request intramolecular (y) or intermolecular (n) distances?"
			IF (user_input_logical()) THEN
				operation_mode="intra"
			ELSE
				operation_mode="inter"
			ENDIF
			WRITE(*,'(" How many ",A,"molecular distances would you like to calculate?")') TRIM(operation_mode)
			PRINT *,"Please enter the number of subjobs as integer."
			number_of_subjobs=user_input_integer(1,100000)
			WRITE(*,'(" opening distance input file...")')
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_distance),IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'(" ",A,"molecular ",I0," ### average ",A,"molecular distances.")')&
			&TRIM(operation_mode),number_of_subjobs,TRIM(operation_mode)
			DO n=1,number_of_subjobs,1
				WRITE(*,'(" Reading subjob number ",I0,"...")') n
				PRINT *,"Do you want to use a wildcard for the reference atom?"
				PRINT *,"(For example 'all H atoms' / all atoms with element label 'H') "
				PRINT *,"(If no (n), you will be asked to give a specific atom_index."
				IF (user_input_logical()) THEN
					PRINT *,"Please enter the desired wildcard / element name for the reference atoms."
					element1=user_input_string(2)
					PRINT *,"Please enter the desired wildcard / element name for the observed atoms."
					element2=user_input_string(2)
					WRITE(8,'(" ",A," ",A," ### consider ",A,"-",A," distances.")')&
					&TRIM(element1),TRIM(element2),TRIM(element1),TRIM(element2)
				ELSE
					PRINT *,"Please enter the molecule_type_index of the reference atom."
					molecule_type_index=user_input_integer(1,maxmol)
					PRINT *,"Please enter the atom_index of the reference atom."
					atom_index=user_input_integer(1,100000)
					WRITE(8,FMT='(" ",I0," ",I0," ")',ADVANCE="NO") molecule_type_index,atom_index
					PRINT *,"Do you want to use a wildcard for the observed atom?"
					IF (user_input_logical()) THEN
						!use wildcard
						PRINT *,"Please enter the desired wildcard / element name for the observed atoms."
						element2=user_input_string(2)
						WRITE(8,'(A," ### atom #",I0," in molecule type #",I0," --> ",A," atoms")')&
						&TRIM(element2),atom_index,molecule_type_index,TRIM(element2)
					ELSE
						!use specific atom. only intermolecular analysis needs extra molecule type index.
						IF (TRIM(operation_mode)=="inter") THEN
							PRINT *,"Please enter the molecule_type_index of the observed atom."
							molecule_type_index=user_input_integer(1,maxmol)
							WRITE(8,FMT='(I0," ")',ADVANCE="NO") molecule_type_index
						ENDIF
						PRINT *,"Please enter the atom_index of the observed atom."
						atom_index=user_input_integer(1,100000)
						WRITE(8,'(I0," ### distance between two specific atoms")') atom_index
					ENDIF
				ENDIF
			ENDDO
			!At this point, subjobs should be written.
			!Thus, continue with reading in the switches:
			PRINT *,"Would you like to also calculate standard deviations?"
			calculate_stdev=user_input_logical()
			IF (calculate_stdev) THEN
				WRITE(8,'(" stdev T ### calculate standard deviations where reasonable.")')
			ELSE
				WRITE(8,'(" stdev F ### do not calculate standard deviations.")')
			ENDIF
			PRINT *,"Would you like to also calculate the FFC distances?"
			PRINT *,"(see ESI of J. Phys. Chem. Lett., 2020, 11, 2165–2170.)"
			calculate_ffc=user_input_logical()
			IF (calculate_ffc) THEN
				WRITE(8,'(" ffc T ### calculate distances as used in FFC.")')
			ELSE
				WRITE(8,'(" ffc F ### calculate distances as used in FFC.")')
			ENDIF
			PRINT *,"Would you like to calculate exponentially weighed average distances?"
			PRINT *,"i.e. average distance = SUM(exp(-k*r)*r)/SUM(exp(-k*r))"
			PRINT *,"(sum runs over valid neighbours within cutoff, r is the distance.)"
			calculate_exponential=user_input_logical()
			IF (calculate_exponential) THEN
				WRITE(8,'(" average_exponential T ### also report average distance weighed with exp(-k*r).")')
				PRINT *,"Please enter the exponent 'k'."
				PRINT *,"(if you have no idea, '1.0' will do)"
				distance_exponent=user_input_real(0.1,100.0)
				WRITE(8,'(" exponent ",E12.5," ### defines the exponent k in exp(-k*r).")') distance_exponent
			ELSE
				WRITE(8,'(" average_exponential F ### do not calculate exponentially weighed averages.")')
			ENDIF
			PRINT *,"Would you like to change advanced settings?"
			PRINT *,"These are sampling interval, range of timesteps and distance cutoff."
			IF (user_input_logical()) THEN
				PRINT *,"Every how many steps would you like to use?"
				WRITE(*,'(A,I0,A)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
				sampling_interval=user_input_integer(1,1000)
				WRITE(8,'(" sampling_interval ",I0," ### use every this many steps.")')sampling_interval
				PRINT *,"What is the largest timestep you would like to use?"
				nsteps=user_input_integer(1,10000)
				WRITE(8,'(" nsteps ",I0," ### defines largest timestep to consider.")')nsteps
				PRINT *,"You need to choose a maximum distance cutoff."
				PRINT *,"Pairs of atoms with larger distance will be ignored."
				PRINT *,"It is possible to use half the box length. Would you like to use that shortcut? (y/n)"
				IF (user_input_logical()) THEN
					WRITE(8,'(" maxdist_optimize ### take half the box length as cutoff.")')
				ELSE
					PRINT *,"Please enter the cutoff / maximum distance."
					maxdist=user_input_real(1.0,1000.0)
					WRITE(8,'(" maxdist ",E12.5," ### use specific cutoff distance.")') maxdist
				ENDIF
			ENDIF
			WRITE(8,*) "quit"
			WRITE(8,*)
			WRITE(8,*) "This is an input file for the calculation of average distances."
			WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
			ENDFILE 8
			CLOSE(UNIT=8)
			WRITE(*,*) "...done"
		END SUBROUTINE user_distance_input

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_distance'
		LOGICAL FUNCTION write_simple_intramolecular_distances()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected,H_available,F_available,both_available
		INTEGER :: n,ios,molecule_type_index
			FILENAME_DISTANCE_INPUT="prealpha_simple.inp"
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),EXIST=file_exists)
			IF (file_exists) CALL report_error(114)
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			!look if there are reasonable intramolecular pairs
			H_available=.FALSE.
			F_available=.FALSE.
			both_available=.FALSE.
			DO molecule_type_index=1,give_number_of_molecule_types(),1
				IF (give_number_of_specific_atoms_per_molecule(molecule_type_index,"H")>1) H_available=.TRUE.
				IF (give_number_of_specific_atoms_per_molecule(molecule_type_index,"F")>1) F_available=.TRUE.
				IF ((give_number_of_specific_atoms_per_molecule(molecule_type_index,"F")>1)&
				&.AND.(give_number_of_specific_atoms_per_molecule(molecule_type_index,"H")>1)) both_available=.TRUE.
			ENDDO
			n=0
			IF (H_available) n=n+1
			IF (F_available) n=n+1
			IF (both_available) n=n+2
			WRITE(8,'("intra ",I0)') n
			IF (n==0) THEN
				write_simple_intramolecular_distances=.FALSE.
				CLOSE(UNIT=8)
				RETURN
			ELSE
				write_simple_intramolecular_distances=.TRUE.
			ENDIF
			IF (H_available) WRITE(8,'("H H")')
			IF (F_available) WRITE(8,'("F F")')
			IF (both_available) THEN
				WRITE(8,'("H F")')
				WRITE(8,'("F H")')
			ENDIF
			WRITE(8,'("exponential_weight T")')
			WRITE(8,'("ffc T")')
			WRITE(8,'("stdev T")')
			IF (give_number_of_timesteps()>10) THEN
				WRITE(8,'("nsteps 10")')
			ELSE
				WRITE(8,'("nsteps ",I0)') give_number_of_timesteps()
			ENDIF
			WRITE(8,'("cutoff_optimize")')
			WRITE(8,'("quit")')
			CLOSE(UNIT=8)
			IF (VERBOSE_OUTPUT) PRINT *,"File 'prealpha_simple.inp' written."
		END FUNCTION write_simple_intramolecular_distances

		!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_input_for_distance'
		LOGICAL FUNCTION write_simple_intermolecular_distances()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected,H_available,F_available,both_available
		INTEGER :: n,ios,molecule_type_index
			FILENAME_DISTANCE_INPUT="prealpha_simple.inp"
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),EXIST=file_exists)
			IF (file_exists) CALL report_error(114)
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			!look if there are reasonable intramolecular pairs
			H_available=.FALSE.
			F_available=.FALSE.
			both_available=.FALSE.
			DO molecule_type_index=1,give_number_of_molecule_types(),1
				IF (give_number_of_specific_atoms("H")>1) H_available=.TRUE.
				IF (give_number_of_specific_atoms("F")>1) F_available=.TRUE.
				IF ((give_number_of_specific_atoms("F")>1)&
				&.AND.(give_number_of_specific_atoms("F")>1)) both_available=.TRUE.
			ENDDO
			n=0
			IF (H_available) n=n+1
			IF (F_available) n=n+1
			IF (both_available) n=n+2
			WRITE(8,'("inter ",I0)') n
			IF (n==0) THEN
				write_simple_intermolecular_distances=.FALSE.
				CLOSE(UNIT=8)
				RETURN
			ELSE
				write_simple_intermolecular_distances=.TRUE.
			ENDIF
			IF (H_available) WRITE(8,'("H H")')
			IF (F_available) WRITE(8,'("F F")')
			IF (both_available) THEN
				WRITE(8,'("H F")')
				WRITE(8,'("F H")')
			ENDIF
			WRITE(8,'("exponential_weight T")')
			WRITE(8,'("ffc T")')
			WRITE(8,'("stdev T")')
			IF (give_number_of_timesteps()>10) THEN
				WRITE(8,'("nsteps 10")')
			ELSE
				WRITE(8,'("nsteps ",I0)') give_number_of_timesteps()
			ENDIF
			WRITE(8,'("cutoff_optimize")')
			WRITE(8,'("quit")')
			CLOSE(UNIT=8)
			IF (VERBOSE_OUTPUT) PRINT *,"File 'prealpha_simple.inp' written."
		END FUNCTION write_simple_intermolecular_distances

		SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
		IMPLICIT NONE
			maxdist=maxdist_default
			maxdist_squared=maxdist_default**2
			sampling_interval=sampling_interval_default
			nsteps=nsteps_default
			calculate_stdev=calculate_stdev_default
			distance_exponent=1.0
		END SUBROUTINE set_defaults

		!initialises the distance module by reading the specified input file.
		SUBROUTINE initialise_distance()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,lines_to_skip,valid_lines,counter,allocstatus
			! first, check if file exists.
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),EXIST=file_exists)
			IF (file_exists) THEN
				!setting defaults to start with.
				CALL set_defaults()
				IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),"'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_DISTANCE_INPUT),&
				&ACTION='READ',IOSTAT=ios)
				IF (ios/=0) CALL report_error(135,exit_status=ios)
				READ(3,IOSTAT=ios,FMT=*) operation_mode,number_of_subjobs!read the operation mode.
				IF (ios/=0) CALL report_error(135,exit_status=ios)
				lines_to_skip=1
				valid_lines=0
				!support for synonyms
				IF (TRIM(operation_mode)=="intramolecular") operation_mode="intra"
				IF (TRIM(operation_mode)=="intermolecular") operation_mode="inter"
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,'(" reading ",I0," user-specified distances...")') number_of_subjobs
				ENDIF
				!open scratch file. Structure of the scratch file is:
				! number of wildcards (0,1,2) - molecule type index - atom index (or wildcard) - (molecule type index for inter) - atom index (or wildcard).
				INQUIRE(UNIT=10,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=10)
				OPEN(UNIT=10,STATUS="SCRATCH")
				REWIND 10
				!Now read the body of the distance input file in line with the requested operation mode:
				SELECT CASE (TRIM(operation_mode))
				CASE ("intra")
					CALL read_input_intra()!also checks for valid atom indices, if not valid prints warning 111 ("this line is ignored")
				CASE ("inter")
					CALL read_input_inter()!also checks for valid atom indices, if not valid prints warning 111 ("this line is ignored")
				CASE DEFAULT
					CLOSE(UNIT=10)
					CALL report_error(135)
				END SELECT
				IF (VERBOSE_OUTPUT) WRITE(*,'(" Successfully read ",I0,"/",I0," custom distance jobs.")')valid_lines,number_of_subjobs
				number_of_subjobs=valid_lines
				IF (number_of_subjobs==0) THEN
					CALL report_error(139)
					CLOSE(UNIT=3)
					CLOSE(UNIT=10)
					RETURN
				ENDIF
				!allocate memory for the subjobs,...
				IF (ALLOCATED(list_of_distances)) CALL report_error(0)
				ALLOCATE(list_of_distances(number_of_subjobs),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
				REWIND 10
				DO counter=1,number_of_subjobs,1
					CALL prepare_subjob(counter)!... then allocate memory for the subsubjobs (if any)
				ENDDO
				CALL goto_currentline()
				!read rest of input file.
				CALL read_body()
				CLOSE(UNIT=3)
				CLOSE(UNIT=10)
			ELSE
				CALL report_error(134)!No input - no output. easy as that.
			ENDIF
			CONTAINS

				SUBROUTINE read_body()
				IMPLICIT NONE
				INTEGER :: n
				CHARACTER(LEN=32) :: inputstring
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  End-of-file condition in ",TRIM(FILENAME_DISTANCE_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DISTANCE_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("sampling_interval")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'sampling_interval' to default (=",sampling_interval_default,")"
								sampling_interval=sampling_interval_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'sampling_interval' to ",sampling_interval
							ENDIF
						CASE ("exponent")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,distance_exponent
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A)') "setting 'exponent' to default (k=1 in exp(-kr))"
								distance_exponent=1.0
							ELSE
								IF (VERBOSE_OUTPUT) THEN
									WRITE(*,'(A,E10.3)') " setting 'exponent' to k=",distance_exponent
									IF (.NOT.(calculate_exponential)) WRITE(*,'(" (exponentially weighed averages are not turned on yet)")')
								ENDIF
							ENDIF
						CASE ("exponential_weight","weigh_exponential","calculate_exponential","average_exponential")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,calculate_exponential
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								calculate_exponential=calculate_exponential_default
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&" setting 'calculate_exponential' to default (",calculate_exponential,")"
							ELSE
								IF (calculate_exponential) THEN
									WRITE(*,'(" Turned on averaged distances with exponential weights exp(-k*r)")')
									WRITE(*,*)"(r=distance, k can be set with 'exponent')"
								ELSE
									WRITE(*,'(" Turned off averaged distances with exponential weights.")')
								ENDIF
							ENDIF
						CASE ("stdev","standard_deviation","standarddev")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,calculate_stdev
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								calculate_stdev=calculate_stdev_default
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&" setting 'calculate_stdev' to default (",calculate_stdev,")"
							ELSE
								IF (calculate_stdev) THEN
									WRITE(*,'(" Also calculate standard deviations, where applicable.")')
								ELSE
									WRITE(*,'(" Turned off calculation of standard deviation.")')
								ENDIF
							ENDIF
						CASE ("ffc_weight","weigh_ffc","calculate_ffc","average_ffc","ffc","dhh","rhh")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,calculate_ffc
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								calculate_ffc=calculate_ffc_default
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&" setting 'calculate_ffc' to default (",calculate_ffc,")"
							ELSE
								IF (calculate_ffc) THEN
									IF (TRIM(operation_mode)=="intra") THEN
										WRITE(*,'(" Will calculate intramolecular distances rHH as described in:")')
									ELSE
										WRITE(*,'(" Will calculate intermolecular distances dHH as described in:")')
									ENDIF
									WRITE(*,'(" J. Phys. Chem. Lett., 2020, 11, 2165–2170. DOI 10.1021/acs.jpclett.0c00087")')
								ELSE
									WRITE(*,*) "Turned off calculation of intra- or intermolecular distances for FFC/NMR"
								ENDIF
							ENDIF
						CASE ("nsteps")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,nsteps
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'nsteps' to default (=",nsteps_default,")"
								nsteps=nsteps_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'nsteps' to ",nsteps
							ENDIF
						CASE ("maxdist","cutoff")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,maxdist
							IF (ios/=0) THEN
								CALL report_error(142,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'maxdist' to default (=",maxdist_default,")"
								maxdist=maxdist_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3)') "   setting 'maxdist' to ",maxdist
							ENDIF
						CASE ("maxdist_optimize","cutoff_optimize")
							maxdist=give_box_limit()/2.0
							IF (maxdist<0.0d0) THEN
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3)')&
								&"  setting 'maxdist' to default (=",maxdist_default,")"
								maxdist=maxdist_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3,A)') "   setting 'maxdist' to ",maxdist," (half the box length)"
							ENDIF
						CASE ("quit")
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_DISTANCE_INPUT)
							EXIT
						CASE DEFAULT
							IF (VERBOSE_OUTPUT) WRITE(*,*) "  can't interpret line - continue streaming"
						END SELECT
					ENDDO
					maxdist_squared=maxdist**2
					IF (nsteps<1) THEN
						CALL report_error(57,exit_status=nsteps)
						nsteps=1
					ELSEIF (nsteps>give_number_of_timesteps()) THEN
						CALL report_error(57,exit_status=nsteps)
						nsteps=give_number_of_timesteps()
					ENDIF
				END SUBROUTINE read_body

				SUBROUTINE prepare_subjob(number_of_subjob)
				IMPLICIT NONE
				INTEGER,INTENT(IN) :: number_of_subjob
				INTEGER :: number_of_wildcards,ios
				INTEGER :: molecule_type_index_1,molecule_type_index_2,atom_index_1,atom_index_2
				CHARACTER (LEN=2) :: element1,element2 !wildcards for element names
					READ(10,IOSTAT=ios,FMT=*) number_of_wildcards
					IF (ios/=0) CALL report_error(0,exit_status=ios)
					BACKSPACE 10
					SELECT CASE (number_of_wildcards)
					CASE (0)
						READ(10,IOSTAT=ios,FMT=*)number_of_wildcards,&
						&molecule_type_index_1,atom_index_1,molecule_type_index_2,atom_index_2
						IF (ios/=0) CALL report_error(0,exit_status=ios)
						list_of_distances(number_of_subjob)%size_indices_1=1
						list_of_distances(number_of_subjob)%size_indices_2=1
						WRITE(list_of_distances(number_of_subjob)%subjob_description,&
						&'("Atom #",I0," in molecule type ",I0," -> atom #",I0," in molecule type ",I0)')&
						&atom_index_1,molecule_type_index_1,atom_index_2,molecule_type_index_2
						!Allocate and initialise memory for reference atom
						ALLOCATE(list_of_distances(number_of_subjob)%indices_1(1,2),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
						list_of_distances(number_of_subjob)%indices_1(1,1)=molecule_type_index_1
						list_of_distances(number_of_subjob)%indices_1(1,2)=atom_index_1
						!Allocate and initialise memory for observed atom
						ALLOCATE(list_of_distances(number_of_subjob)%indices_2(1,2),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
						list_of_distances(number_of_subjob)%indices_2(1,1)=molecule_type_index_2
						list_of_distances(number_of_subjob)%indices_2(1,2)=atom_index_2
						!looks alright - I have this many molecules, and one distance for each.
						list_of_distances(number_of_subjob)%weight=&
						&give_number_of_molecules_per_step(molecule_type_index_1)
					CASE (1)
						READ(10,IOSTAT=ios,FMT=*) number_of_wildcards,molecule_type_index_1,atom_index_1,element2
						IF (ios/=0) CALL report_error(0,exit_status=ios)
						list_of_distances(number_of_subjob)%size_indices_1=1
						IF (TRIM(operation_mode)=="intra") THEN
							list_of_distances(number_of_subjob)%size_indices_2=&
							&give_number_of_specific_atoms_per_molecule(molecule_type_index_1,TRIM(element2))
							WRITE(list_of_distances(number_of_subjob)%subjob_description,&
							&'("Atom #",I0," in molecule type ",I0," -> ",I0," ",A," atoms.")')&
							&atom_index_1,molecule_type_index_1,&
							&list_of_distances(number_of_subjob)%size_indices_2,&
							&TRIM(element2)
							!Allocate and initialise memory for reference atom
							ALLOCATE(list_of_distances(number_of_subjob)%indices_1(1,2),STAT=allocstatus)
							IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
							list_of_distances(number_of_subjob)%indices_1(1,1)=molecule_type_index_1
							list_of_distances(number_of_subjob)%indices_1(1,2)=atom_index_1
							!Allocate and initialise memory for observed atom
							ALLOCATE(list_of_distances(number_of_subjob)%indices_2(&
							&list_of_distances(number_of_subjob)%size_indices_2,2),STAT=allocstatus)
							IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
							!The following subroutine also initialised all entries of the first dimension to "molecule_type_index_1"
							CALL give_indices_of_specific_atoms_per_molecule(&
							&molecule_type_index_1,element2,list_of_distances(number_of_subjob)%indices_2)
						ELSE
							list_of_distances(number_of_subjob)%size_indices_2=give_number_of_specific_atoms(TRIM(element2))
							WRITE(list_of_distances(number_of_subjob)%subjob_description,&
							&'("Atom #",I0," in molecule type ",I0," -> ",I0," ",A," atoms.")')&
							&atom_index_1,molecule_type_index_1,&
							&list_of_distances(number_of_subjob)%size_indices_2,&
							&TRIM(element2)
							!Allocate and initialise memory for reference atom
							ALLOCATE(list_of_distances(number_of_subjob)%indices_1(1,2),STAT=allocstatus)
							IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
							list_of_distances(number_of_subjob)%indices_1(1,1)=molecule_type_index_1
							list_of_distances(number_of_subjob)%indices_1(1,2)=atom_index_1
							!Allocate and initialise memory for observed atom
							ALLOCATE(list_of_distances(number_of_subjob)%indices_2(&
							&list_of_distances(number_of_subjob)%size_indices_2,2),STAT=allocstatus)
							IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
							CALL give_indices_of_specific_atoms(element2,list_of_distances(number_of_subjob)%indices_2)
						ENDIF
						!there is still only one closest distance.
						list_of_distances(number_of_subjob)%weight=&
						&give_number_of_molecules_per_step(molecule_type_index_1)
					CASE (2)
						READ(10,IOSTAT=ios,FMT=*) number_of_wildcards,element1,element2
						IF (ios/=0) CALL report_error(0,exit_status=ios)
						list_of_distances(number_of_subjob)%size_indices_1=&
						&give_number_of_specific_atoms(TRIM(element1))
						list_of_distances(number_of_subjob)%size_indices_2=&
						&give_number_of_specific_atoms(TRIM(element2))
						WRITE(list_of_distances(number_of_subjob)%subjob_description,&
						&'("",I0," ",A," atoms -> ",I0," ",A," atoms.")')&
						&list_of_distances(number_of_subjob)%size_indices_1,TRIM(element1),&
						&list_of_distances(number_of_subjob)%size_indices_2,TRIM(element2)
						!Allocate and initialise memory for reference atom
						ALLOCATE(list_of_distances(number_of_subjob)%indices_1(&
						&list_of_distances(number_of_subjob)%size_indices_1,2),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
						CALL give_indices_of_specific_atoms(element1,list_of_distances(number_of_subjob)%indices_1)
						!Allocate and initialise memory for observed atom
						ALLOCATE(list_of_distances(number_of_subjob)%indices_2(&
						&list_of_distances(number_of_subjob)%size_indices_2,2),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(137,exit_status=allocstatus)
						CALL give_indices_of_specific_atoms(element2,list_of_distances(number_of_subjob)%indices_2)
						!Here, calculating the weights gets more challenging - especially for the intermolecular distances.
						list_of_distances(number_of_subjob)%weight=count_wildcard_weights(number_of_subjob)
					CASE DEFAULT
						CALL report_error(0,exit_status=number_of_wildcards)
					END SELECT
					list_of_distances(number_of_subjob)%element1=give_element_symbol(&
					&list_of_distances(number_of_subjob)%indices_1(1,1),&
					&list_of_distances(number_of_subjob)%indices_1(1,2))
					list_of_distances(number_of_subjob)%element2=give_element_symbol(&
					&list_of_distances(number_of_subjob)%indices_2(1,1),&
					&list_of_distances(number_of_subjob)%indices_2(1,2))
				END SUBROUTINE prepare_subjob

				SUBROUTINE goto_currentline()
				IMPLICIT NONE
				INTEGER :: i
					REWIND 3
					DO i=1,lines_to_skip,1
						READ(3,IOSTAT=ios,FMT=*)
						IF (ios/=0) CALL report_error(0,exit_status=ios) 
					ENDDO
				END SUBROUTINE goto_currentline

				SUBROUTINE read_input_inter()!This subroutine reads and checks input
				IMPLICIT NONE
				CHARACTER (LEN=2) :: element1,element2 !wildcards for element names
				INTEGER :: molecule_type_index_1,molecule_type_index_2,atom_index_1,atom_index_2
					!try and read the current line
					DO
						READ(3,IOSTAT=ios,FMT=*) molecule_type_index_1,atom_index_1,molecule_type_index_2,atom_index_2
						IF (ios<0) EXIT
						IF (ios/=0) THEN
							CALL goto_currentline()
							!test for first wildcard
							READ(3,IOSTAT=ios,FMT=*) element1,molecule_type_index_2,atom_index_2
							IF (ios<0) EXIT
							IF (ios/=0) THEN
								CALL goto_currentline()
								!test for second wildcard
								READ(3,IOSTAT=ios,FMT=*) molecule_type_index_1,atom_index_1,element2
								IF (ios<0) EXIT
								IF (ios/=0) THEN
									CALL goto_currentline()
									!test for both wildcards
									READ(3,IOSTAT=ios,FMT=*) element1,element2
									IF (ios<0) EXIT
									IF (ios/=0) THEN
										!possibility to print error message here.
										EXIT
									ELSE
										!two wildcards!
										lines_to_skip=lines_to_skip+1
										IF (give_number_of_specific_atoms(element1)==0) THEN
											CALL report_error(136)
										ELSEIF (give_number_of_specific_atoms(element2)==0) THEN
											CALL report_error(136)
										ELSE
											valid_lines=valid_lines+1
											WRITE(10,*) 2,element1,element2
										ENDIF
									ENDIF
								ELSE
									lines_to_skip=lines_to_skip+1
									IF ((atom_index_1<1).OR.(atom_index_1>give_number_of_atoms_per_molecule(molecule_type_index_1))) THEN
										CALL report_error(111,exit_status=atom_index_1)
									ELSEIF ((molecule_type_index_1<1).OR.(molecule_type_index_1>give_number_of_molecule_types())) THEN
										CALL report_error(110,exit_status=molecule_type_index_1)
									ELSEIF (give_number_of_specific_atoms(element2)==0) THEN
										CALL report_error(136)
									ELSE
										valid_lines=valid_lines+1
										WRITE(10,*) 1,molecule_type_index_1,atom_index_1,element2
									ENDIF
								ENDIF
							ELSE
								lines_to_skip=lines_to_skip+1
								IF ((atom_index_2<1).OR.(atom_index_2>give_number_of_atoms_per_molecule(molecule_type_index_2))) THEN
									CALL report_error(111,exit_status=atom_index_2)
								ELSEIF ((molecule_type_index_2<1).OR.(molecule_type_index_2>give_number_of_molecule_types())) THEN
									CALL report_error(110,exit_status=molecule_type_index_2)
								ELSEIF (give_number_of_specific_atoms(element1)==0) THEN
									CALL report_error(136)
								ELSE
									valid_lines=valid_lines+1
									CALL report_error(140)
									WRITE(10,*) 1,molecule_type_index_2,atom_index_2,element1
								ENDIF
							ENDIF
						ELSE
							!no wildcards
							lines_to_skip=lines_to_skip+1
							IF ((atom_index_1<1).OR.(atom_index_1>give_number_of_atoms_per_molecule(molecule_type_index_1))) THEN
								CALL report_error(111,exit_status=atom_index_1)
							ELSEIF ((atom_index_2<1).OR.(atom_index_2>give_number_of_atoms_per_molecule(molecule_type_index_2))) THEN
								CALL report_error(111,exit_status=atom_index_2)
							ELSEIF ((molecule_type_index_1<1).OR.(molecule_type_index_1>give_number_of_molecule_types())) THEN
								CALL report_error(110,exit_status=molecule_type_index_1)
							ELSEIF ((molecule_type_index_2<1).OR.(molecule_type_index_2>give_number_of_molecule_types())) THEN
								CALL report_error(110,exit_status=molecule_type_index_2)
							ELSEIF ((molecule_type_index_1==molecule_type_index_2).AND.&
							&(give_number_of_molecules_per_step(molecule_type_index_1)==1)) THEN
								CALL report_error(141)
							ELSE
								valid_lines=valid_lines+1
								WRITE(10,*) 0,molecule_type_index_1,atom_index_1,molecule_type_index_2,atom_index_2
							ENDIF
						ENDIF
						IF (lines_to_skip==(number_of_subjobs+1)) EXIT
					ENDDO
				END SUBROUTINE read_input_inter

				SUBROUTINE read_input_intra()!This subroutine reads and checks input
				IMPLICIT NONE
				CHARACTER (LEN=2) :: element1,element2 !wildcards for element names
				INTEGER :: molecule_type_index_1,atom_index_1,atom_index_2,i
				LOGICAL :: valid_elements
					!try and read the current line
					DO
						READ(3,IOSTAT=ios,FMT=*) molecule_type_index_1,atom_index_1,atom_index_2
						IF (ios<0) EXIT
						IF (ios/=0) THEN
							CALL goto_currentline()
							!test for first wildcard
							READ(3,IOSTAT=ios,FMT=*) molecule_type_index_1,element1,atom_index_2
							IF (ios<0) EXIT
							IF (ios/=0) THEN
								CALL goto_currentline()
								!test for second wildcard
								READ(3,IOSTAT=ios,FMT=*) molecule_type_index_1,atom_index_1,element2
								IF (ios<0) EXIT
								IF (ios/=0) THEN
									CALL goto_currentline()
									!test for both wildcards
									READ(3,IOSTAT=ios,FMT=*) element1,element2
									IF (ios<0) EXIT
									IF (ios/=0) THEN
										!possibility to print error message here.
										EXIT
									ELSE
										!both wildcards
										lines_to_skip=lines_to_skip+1
										valid_elements=.FALSE.
										DO i=1,give_number_of_molecule_types(),1
											IF ((give_number_of_specific_atoms_per_molecule(i,TRIM(element1))>0)&
											&.AND.(give_number_of_specific_atoms_per_molecule(i,TRIM(element2))>0)) THEN
												valid_elements=.TRUE.
												EXIT
											ENDIF
										ENDDO
										IF (valid_elements) THEN
											valid_lines=valid_lines+1
											WRITE(10,*) 2,element1,element2
										ELSE
											CALL report_error(136)
										ENDIF
									ENDIF
								ELSE
									!second wildcard
									lines_to_skip=lines_to_skip+1
									IF ((atom_index_1<1).OR.(atom_index_1>give_number_of_atoms_per_molecule(molecule_type_index_1))) THEN
										CALL report_error(111,exit_status=atom_index_1)
									ELSEIF ((molecule_type_index_1<1).OR.(molecule_type_index_1>give_number_of_molecule_types())) THEN
										CALL report_error(110,exit_status=molecule_type_index_1)
									ELSEIF (give_number_of_specific_atoms_per_molecule(molecule_type_index_1,element2)==0) THEN
										CALL report_error(136,molecule_type_index_1)
									ELSE
										valid_lines=valid_lines+1
										WRITE(10,*) 1,molecule_type_index_1,atom_index_1,element2
									ENDIF
								ENDIF
							ELSE
								!first wildcard
								lines_to_skip=lines_to_skip+1
								IF ((atom_index_2<1).OR.(atom_index_2>give_number_of_atoms_per_molecule(molecule_type_index_1))) THEN
									CALL report_error(111,exit_status=atom_index_2)
								ELSEIF ((molecule_type_index_1<1).OR.(molecule_type_index_1>give_number_of_molecule_types())) THEN
									CALL report_error(110,exit_status=molecule_type_index_1)
								ELSEIF (give_number_of_specific_atoms_per_molecule(molecule_type_index_1,element1)==0) THEN
									CALL report_error(136,molecule_type_index_1)
								ELSE
									CALL report_error(140)
									valid_lines=valid_lines+1
									WRITE(10,*) 1,molecule_type_index_1,atom_index_2,element1
								ENDIF
							ENDIF
						ELSE
							!no wildcards
							lines_to_skip=lines_to_skip+1
							IF ((atom_index_1<1).OR.(atom_index_1>give_number_of_atoms_per_molecule(molecule_type_index_1))) THEN
								CALL report_error(111,exit_status=atom_index_1)
							ELSEIF ((atom_index_2<1).OR.(atom_index_2>give_number_of_atoms_per_molecule(molecule_type_index_1))) THEN
								CALL report_error(111,exit_status=atom_index_2)
							ELSEIF (atom_index_2==atom_index_1) THEN
								CALL report_error(111,exit_status=atom_index_2)
							ELSEIF ((molecule_type_index_1<1).OR.(molecule_type_index_1>give_number_of_molecule_types())) THEN
								CALL report_error(110,exit_status=molecule_type_index_1)
							ELSE
								valid_lines=valid_lines+1
								WRITE(10,*) 0,molecule_type_index_1,atom_index_1,molecule_type_index_1,atom_index_2
							ENDIF
						ENDIF
						IF (lines_to_skip==(number_of_subjobs+1)) EXIT
					ENDDO
				END SUBROUTINE read_input_intra

		END SUBROUTINE initialise_distance

		INTEGER FUNCTION count_wildcard_weights(number_of_subjob2)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: number_of_subjob2
		INTEGER :: indices_counter_1,indices_counter_2
		INTEGER :: molecule_type_index_1,molecule_type_index_2,atom_index_1,atom_index_2
			count_wildcard_weights=0
			!two loops: over the reference atoms...
			DO indices_counter_1=1,list_of_distances(number_of_subjob2)%size_indices_1,1
				molecule_type_index_1=list_of_distances(number_of_subjob2)%indices_1(indices_counter_1,1)
				atom_index_1=list_of_distances(number_of_subjob2)%indices_1(indices_counter_1,2)
				!... and the observed atoms.
				!Each reference atom will contribute to the weight IF it has at least one valid neighbour.
				DO indices_counter_2=1,list_of_distances(number_of_subjob2)%size_indices_2,1
					molecule_type_index_2=list_of_distances(number_of_subjob2)%indices_2(indices_counter_2,1)
					atom_index_2=list_of_distances(number_of_subjob2)%indices_2(indices_counter_2,2)
					IF (TRIM(operation_mode)=="intra") THEN
						IF ((molecule_type_index_1==molecule_type_index_2).AND.(atom_index_1/=atom_index_2)) THEN
							count_wildcard_weights=count_wildcard_weights+&
							&give_number_of_molecules_per_step(molecule_type_index_1)
							EXIT
						ENDIF
					ELSE
						IF ((molecule_type_index_1/=molecule_type_index_2).OR.(atom_index_1/=atom_index_2)) THEN
							count_wildcard_weights=count_wildcard_weights+&
							&give_number_of_molecules_per_step(molecule_type_index_1)
							EXIT
						ELSEIF ((molecule_type_index_1==molecule_type_index_2).AND.(atom_index_1==atom_index_2)) THEN
							IF (give_number_of_molecules_per_step(molecule_type_index_1)>1) THEN
								count_wildcard_weights=count_wildcard_weights+&
								&give_number_of_molecules_per_step(molecule_type_index_1)
								EXIT
							ENDIF
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		END FUNCTION count_wildcard_weights

		!finalises the distance module.
		SUBROUTINE finalise_distance()
		IMPLICIT NONE
		INTEGER :: deallocstatus
			IF (ALLOCATED(list_of_distances)) THEN
				DEALLOCATE(list_of_distances,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
		END SUBROUTINE finalise_distance

		SUBROUTINE compute_intramolecular_distances(stdev_run)
		IMPLICIT NONE
		INTEGER :: stepcounter,subjob_counter,indices_counter_1,indices_counter_2,localweight
		INTEGER :: molecule_type_index_1,molecule_type_index_2,atom_index_1,atom_index_2,molecule_index
		REAL :: distance_squared,closest_distance_squared,distance
		REAL(KIND=WORKING_PRECISION) :: local_closest,local_exponential,local_ffc,weight_exponential,local_exponential_stdev,local_ffcexp
		REAL(KIND=WORKING_PRECISION) :: exponential_function,ffc_function
		LOGICAL :: valid_neighbour
		LOGICAL,INTENT(IN) :: stdev_run
		INTEGER(KIND=WORKING_PRECISION) :: localweight_ffc
			!average over timesteps, a few should be enough
			DO stepcounter=1,nsteps,sampling_interval
				!iterating through the subjobs
				DO subjob_counter=1,number_of_subjobs,1
					localweight=0
					local_closest=0.0d0
					local_exponential=0.0d0
					weight_exponential=0.0d0
					local_exponential_stdev=0.0d0
					local_ffc=0.0d0
					local_ffcexp=0.0d0
					localweight_ffc=0
					DO indices_counter_1=1,list_of_distances(subjob_counter)%size_indices_1,1 !two loops: over the reference atoms...
						molecule_type_index_1=list_of_distances(subjob_counter)%indices_1(indices_counter_1,1)
						atom_index_1=list_of_distances(subjob_counter)%indices_1(indices_counter_1,2)
						DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index_1),1
							!Find the closest neighbour (=observed atom) for each molecule_index
							closest_distance_squared=maxdist_squared
							valid_neighbour=.FALSE.
							!for the standard deviation, ensure values are reset before the loop over observed atoms
							IF ((stdev_run).AND.(calculate_exponential)) THEN
								local_exponential=0.0d0
								weight_exponential=0.0d0
							ENDIF
							DO indices_counter_2=1,list_of_distances(subjob_counter)%size_indices_2,1 !... and the observed atoms.
								molecule_type_index_2=list_of_distances(subjob_counter)%indices_2(indices_counter_2,1)
								IF (molecule_type_index_1/=molecule_type_index_2) CYCLE !intramolecular distances only...
								atom_index_2=list_of_distances(subjob_counter)%indices_2(indices_counter_2,2)
								IF (atom_index_1/=atom_index_2) THEN
									distance_squared=give_smallest_atom_distance_squared&
									&(stepcounter,stepcounter,molecule_type_index_1,molecule_type_index_1,&
									&molecule_index,molecule_index,atom_index_1,atom_index_2)
									IF (distance_squared<maxdist_squared) THEN
										!Atom pair is closer than cutoff
										valid_neighbour=.TRUE.
										IF (distance_squared<closest_distance_squared) THEN
											closest_distance_squared=distance_squared
										ENDIF
										!include exp function
										IF (calculate_exponential) THEN
											distance=SQRT(distance_squared)
											exponential_function=exp(-distance_exponent*distance)
											local_exponential=local_exponential+exponential_function*distance
											weight_exponential=weight_exponential+exponential_function
										ENDIF
										IF (calculate_ffc) THEN
											!eq (6) in the supporting information of J. Phys. Chem. Lett., 2020, 11, 2165–2170, DOI 10.1021/acs.jpclett.0c00087
											!summing 1/(r^6) = 1/((r^2)^3)
											ffc_function=1/((distance_squared)**3)
											local_ffc=local_ffc+ffc_function
											localweight_ffc=localweight_ffc+1
											IF (calculate_exponential) THEN
												!exponentially weighed ffc...
												local_ffcexp=local_ffcexp+exponential_function*ffc_function
											ENDIF
										ENDIF
									ENDIF
								ENDIF
							ENDDO
							IF (valid_neighbour) THEN
								localweight=localweight+1
								IF (stdev_run) THEN
									local_closest=local_closest+&
									&(SQRT(closest_distance_squared)-list_of_distances(subjob_counter)%closest_average)**2
									IF (calculate_exponential) THEN
										local_exponential_stdev=local_exponential_stdev+&
										&(local_exponential/weight_exponential-list_of_distances(subjob_counter)%exponential_average)**2
									ENDIF
								ELSE
									local_closest=local_closest+SQRT(closest_distance_squared)
								ENDIF
							ELSE
								list_of_distances(subjob_counter)%missed_neighbours=list_of_distances(subjob_counter)%missed_neighbours+1
							ENDIF
						ENDDO
					ENDDO
					!Update results for subjob
					IF (stdev_run) THEN
						list_of_distances(subjob_counter)%closest_stdev=&
						&list_of_distances(subjob_counter)%closest_stdev+local_closest/FLOAT(localweight)
						IF (calculate_exponential) list_of_distances(subjob_counter)%exponential_stdev=&
						&list_of_distances(subjob_counter)%exponential_stdev+local_exponential_stdev/FLOAT(localweight)
					ELSE
						list_of_distances(subjob_counter)%closest_average=&
						&list_of_distances(subjob_counter)%closest_average+local_closest/FLOAT(localweight)
						IF ((list_of_distances(subjob_counter)%missed_neighbours+localweight)&
						&/=(list_of_distances(subjob_counter)%weight)) THEN
							CALL report_error(0)
						ENDIF
						IF (calculate_exponential) THEN
							list_of_distances(subjob_counter)%exponential_average=&
							&list_of_distances(subjob_counter)%exponential_average+local_exponential
							list_of_distances(subjob_counter)%exponential_weight=&
							&list_of_distances(subjob_counter)%exponential_weight+weight_exponential
						ENDIF
						IF (calculate_ffc) THEN
							list_of_distances(subjob_counter)%ffc_average=&
							&list_of_distances(subjob_counter)%ffc_average+local_ffc
							IF (calculate_exponential) THEN
								list_of_distances(subjob_counter)%ffcexp_average=&
								&list_of_distances(subjob_counter)%ffcexp_average+local_ffcexp
							ENDIF
							list_of_distances(subjob_counter)%ffc_weight=&
							&list_of_distances(subjob_counter)%ffc_weight+FLOAT(localweight_ffc)
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		END SUBROUTINE compute_intramolecular_distances

		SUBROUTINE compute_intermolecular_distances(stdev_run)
		IMPLICIT NONE
		INTEGER :: stepcounter,subjob_counter,indices_counter_1,indices_counter_2,localweight
		INTEGER :: molecule_type_index_1,molecule_type_index_2,atom_index_1,atom_index_2,molecule_index_1,molecule_index_2
		REAL :: distance_squared,closest_distance_squared,distance
		REAL(KIND=WORKING_PRECISION) :: local_closest,local_exponential,local_ffc,weight_exponential,local_exponential_stdev,local_ffcexp
		REAL(KIND=WORKING_PRECISION) :: exponential_function,ffc_function
		LOGICAL :: valid_neighbour
		LOGICAL,INTENT(IN) :: stdev_run
		INTEGER(KIND=WORKING_PRECISION) :: localweight_ffc
			!average over timesteps, a few should be enough
			DO stepcounter=1,nsteps,sampling_interval
				!iterating through the subjobs
				DO subjob_counter=1,number_of_subjobs,1
					localweight=0
					local_closest=0.0d0
					local_exponential=0.0d0
					weight_exponential=0.0d0
					local_exponential_stdev=0.0d0
					local_ffc=0.0d0
					local_ffcexp=0.0d0
					localweight_ffc=0
					DO indices_counter_1=1,list_of_distances(subjob_counter)%size_indices_1,1 !two loops: over the reference atoms...
						molecule_type_index_1=list_of_distances(subjob_counter)%indices_1(indices_counter_1,1)
						atom_index_1=list_of_distances(subjob_counter)%indices_1(indices_counter_1,2)
						DO molecule_index_1=1,give_number_of_molecules_per_step(molecule_type_index_1),1
							!Find the closest neighbour (=observed atom) for each molecule_index_1
							closest_distance_squared=maxdist_squared
							valid_neighbour=.FALSE.
							!for the standard deviation, ensure values are reset before the loop over observed atoms
							IF ((stdev_run).AND.(calculate_exponential)) THEN
								local_exponential=0.0d0
								weight_exponential=0.0d0
							ENDIF
							DO indices_counter_2=1,list_of_distances(subjob_counter)%size_indices_2,1 !... and the observed atoms.
								molecule_type_index_2=list_of_distances(subjob_counter)%indices_2(indices_counter_2,1)
								atom_index_2=list_of_distances(subjob_counter)%indices_2(indices_counter_2,2)
								DO molecule_index_2=1,give_number_of_molecules_per_step(molecule_type_index_2),1
									IF ((molecule_index_1==molecule_index_2).AND.(molecule_type_index_1==molecule_type_index_2)) CYCLE
									distance_squared=give_smallest_atom_distance_squared&
									&(stepcounter,stepcounter,molecule_type_index_1,molecule_type_index_2,&
									&molecule_index_1,molecule_index_2,atom_index_1,atom_index_2)
									IF (distance_squared<maxdist_squared) THEN
										!Atom pair is closer than cutoff
										valid_neighbour=.TRUE.
										IF (distance_squared<closest_distance_squared) THEN
											closest_distance_squared=distance_squared
										ENDIF
										!include exp function
										IF (calculate_exponential) THEN
											distance=SQRT(distance_squared)
											exponential_function=exp(-distance_exponent*distance)
											local_exponential=local_exponential+exponential_function*distance
											weight_exponential=weight_exponential+exponential_function
										ENDIF
										IF (calculate_ffc) THEN
											!eq (13) in the supporting information of J. Phys. Chem. Lett., 2020, 11, 2165–2170, DOI 10.1021/acs.jpclett.0c00087
											!summing 1/(r^3)
											distance=SQRT(distance_squared)
											ffc_function=1/((distance)**3)
											local_ffc=local_ffc+ffc_function
											localweight_ffc=localweight_ffc+1
											IF (calculate_exponential) THEN
												!exponentially weighed ffc...
												local_ffcexp=local_ffcexp+exponential_function*ffc_function
											ENDIF
										ENDIF
									ENDIF
								ENDDO
							ENDDO
							IF (valid_neighbour) THEN
								localweight=localweight+1
								IF (stdev_run) THEN
									local_closest=local_closest+&
									&(SQRT(closest_distance_squared)-list_of_distances(subjob_counter)%closest_average)**2
									IF (calculate_exponential) THEN
										local_exponential_stdev=local_exponential_stdev+&
										&(local_exponential/weight_exponential-list_of_distances(subjob_counter)%exponential_average)**2
									ENDIF
								ELSE
									local_closest=local_closest+SQRT(closest_distance_squared)
								ENDIF
							ELSE
								list_of_distances(subjob_counter)%missed_neighbours=list_of_distances(subjob_counter)%missed_neighbours+1
							ENDIF
						ENDDO
					ENDDO
					!Update results for subjob
					IF (stdev_run) THEN
						list_of_distances(subjob_counter)%closest_stdev=&
						&list_of_distances(subjob_counter)%closest_stdev+local_closest/FLOAT(localweight)
						IF (calculate_exponential) list_of_distances(subjob_counter)%exponential_stdev=&
						&list_of_distances(subjob_counter)%exponential_stdev+local_exponential_stdev/FLOAT(localweight)
					ELSE
						list_of_distances(subjob_counter)%closest_average=&
						&list_of_distances(subjob_counter)%closest_average+local_closest/FLOAT(localweight)
						IF (calculate_exponential) THEN
							list_of_distances(subjob_counter)%exponential_average=&
							&list_of_distances(subjob_counter)%exponential_average+local_exponential
							list_of_distances(subjob_counter)%exponential_weight=&
							&list_of_distances(subjob_counter)%exponential_weight+weight_exponential
						ENDIF
						IF (calculate_ffc) THEN
							list_of_distances(subjob_counter)%ffc_average=&
							&list_of_distances(subjob_counter)%ffc_average+local_ffc
							IF (calculate_exponential) THEN
								list_of_distances(subjob_counter)%ffcexp_average=&
								&list_of_distances(subjob_counter)%ffcexp_average+local_ffcexp
							ENDIF
							list_of_distances(subjob_counter)%ffc_weight=&
							&list_of_distances(subjob_counter)%ffc_weight+FLOAT(localweight_ffc)
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		END SUBROUTINE compute_intermolecular_distances

		SUBROUTINE report_distances()
		IMPLICIT NONE
		INTEGER :: subjob_counter
			WRITE(*,'(" Done with ",A,"molecular distance calculation. Results:")') TRIM(operation_mode)
			DO subjob_counter=1,number_of_subjobs,1
				WRITE(*,'("   Subjob #",I0," out of ",I0,":")') subjob_counter,number_of_subjobs
				WRITE(*,'("     ",A,"")') TRIM(list_of_distances(subjob_counter)%subjob_description)
				IF (list_of_distances(subjob_counter)%missed_neighbours>0) THEN
					WRITE(*,'("     ",I0," atoms without neighbour - consider increasing cutoff.")')&
					&list_of_distances(subjob_counter)%missed_neighbours
				ELSE
					WRITE(*,'("     Used ",I0," reference atoms per step.")') list_of_distances(subjob_counter)%weight
				ENDIF
				CALL formatted_distance_output("average closest distance",list_of_distances(subjob_counter)%closest_average)
				IF (calculate_stdev) CALL formatted_distance_output("standard deviation",list_of_distances(subjob_counter)%closest_stdev)
				IF (calculate_exponential) THEN
					CALL formatted_distance_output("exponentially weighed distance",list_of_distances(subjob_counter)%exponential_average)
					IF (calculate_stdev) CALL formatted_distance_output("standard deviation",list_of_distances(subjob_counter)%exponential_stdev)
				ENDIF
				IF (calculate_ffc) THEN
					IF (TRIM(operation_mode)=="intra") THEN
						!"ffc" calculates rXX
						CALL formatted_distance_output(&
						&"FFC distance r"//TRIM(list_of_distances(subjob_counter)%element1)//&
						&TRIM(list_of_distances(subjob_counter)%element2),&
						&list_of_distances(subjob_counter)%ffc_average)
						CALL formatted_distance_output(&
						&"exp(-kr) weighed FFC distance r'"//TRIM(list_of_distances(subjob_counter)%element1)//&
						&TRIM(list_of_distances(subjob_counter)%element2),&
						&list_of_distances(subjob_counter)%ffcexp_average)
						CALL formatted_distance_output("Number of averages for FFC",&
						&list_of_distances(subjob_counter)%ffc_weight)
					ELSE
						!"ffc" calculates dXX
						CALL formatted_distance_output(&
						&"FFC distance d"//TRIM(list_of_distances(subjob_counter)%element1)//&
						&TRIM(list_of_distances(subjob_counter)%element2),&
						&list_of_distances(subjob_counter)%ffc_average)
						CALL formatted_distance_output(&
						&"exp(-kr) weighed FFC distance d'"//TRIM(list_of_distances(subjob_counter)%element1)//&
						&TRIM(list_of_distances(subjob_counter)%element2),&
						&list_of_distances(subjob_counter)%ffcexp_average)
						CALL formatted_distance_output("Number of averages for FFC",&
						&list_of_distances(subjob_counter)%ffc_weight)
					ENDIF
				ENDIF
			ENDDO
			WRITE(*,'(" ( reference atom(s) -> observed atom(s) )")')
		CONTAINS

			SUBROUTINE formatted_distance_output(description,distance)
			IMPLICIT NONE
			REAL(KIND=WORKING_PRECISION),INTENT(IN) :: distance
			CHARACTER(LEN=*),INTENT(IN) :: description
				WRITE(*,FMT='("     ",A,": ")',ADVANCE="NO") TRIM(description)
				IF (distance>1000.0) THEN
					WRITE(*,'(E9.3)') distance
				ELSEIF (distance>1.0) THEN
					WRITE(*,'(F0.3)') distance
				ELSEIF (distance>0.01) THEN
					WRITE(*,'(F5.3)') distance
				ELSE
					WRITE(*,'(E9.3)') distance
				ENDIF
			END SUBROUTINE formatted_distance_output

		END SUBROUTINE report_distances

		SUBROUTINE perform_distance_analysis()
		IMPLICIT NONE
			CALL initialise_distance()
			IF ((ERROR_CODE/=134).AND.(ERROR_CODE/=139)) THEN
				WRITE(*,'(" Using ",I0," timesteps in intervals of ",I0," for averaging.")')&
				&MAX((nsteps-1+sampling_interval)/sampling_interval,0),sampling_interval
				WRITE(*,'(" calculating average ",A,"molecular distances.")')TRIM(operation_mode)
				CALL refresh_IO()
				!initialise averages
				list_of_distances(:)%closest_average=0.0d0
				list_of_distances(:)%closest_stdev=0.0d0
				list_of_distances(:)%exponential_average=0.0d0
				list_of_distances(:)%exponential_stdev=0.0d0
				list_of_distances(:)%exponential_weight=0.0d0
				list_of_distances(:)%missed_neighbours=0
				list_of_distances(:)%ffc_average=0.0d0
				list_of_distances(:)%ffcexp_average=0.0d0
				list_of_distances(:)%ffc_weight=0.0d0
				SELECT CASE (TRIM(operation_mode))
				CASE ("intra")
					CALL compute_intramolecular_distances(.FALSE.)
				CASE ("inter")
					CALL compute_intermolecular_distances(.FALSE.)
				CASE DEFAULT
					CALL report_error(0)
				END SELECT
				!post-processing
				list_of_distances(:)%closest_average=&
				&list_of_distances(:)%closest_average/FLOAT(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
				IF (calculate_exponential) list_of_distances(:)%exponential_average=&
				&list_of_distances(:)%exponential_average/list_of_distances(:)%exponential_weight
				IF (calculate_ffc) THEN
					list_of_distances(:)%ffc_average=&
					&list_of_distances(:)%ffc_average/list_of_distances(:)%ffc_weight
					IF (calculate_exponential) THEN
						!exponentially weighed ffc...
						list_of_distances(:)%ffcexp_average=&
						&list_of_distances(:)%ffcexp_average/list_of_distances(:)%exponential_weight
					ENDIF
					IF (TRIM(operation_mode)=="intra") THEN
						list_of_distances(:)%ffc_average=&
						&list_of_distances(:)%ffc_average**(-1.0/6.0)
						IF (calculate_exponential) THEN
							!exponentially weighed ffc...
							list_of_distances(:)%ffcexp_average=&
							&list_of_distances(:)%ffcexp_average**(-1.0/6.0)
						ENDIF
					ELSE
						list_of_distances(:)%ffc_average=&
						&list_of_distances(:)%ffc_average**(-1.0/3.0)
						IF (calculate_exponential) THEN
							!exponentially weighed ffc...
							list_of_distances(:)%ffcexp_average=&
							&list_of_distances(:)%ffcexp_average**(-1.0/3.0)
						ENDIF
					ENDIF
				ENDIF
				!calculate standard deviation if necessary
				IF (calculate_stdev) THEN
					IF (VERBOSE_OUTPUT) WRITE(*,'(" calculating standard deviations.")')
					CALL refresh_IO()
					SELECT CASE (TRIM(operation_mode))
					CASE ("intra")
						CALL compute_intramolecular_distances(.TRUE.)
					CASE ("inter")
						CALL compute_intermolecular_distances(.TRUE.)
					CASE DEFAULT
						CALL report_error(0)
					END SELECT
					!post-processing
					list_of_distances(:)%closest_stdev=&
					&SQRT(list_of_distances(:)%closest_stdev/FLOAT(MAX((nsteps-1+sampling_interval)/sampling_interval,0)))
					IF (calculate_exponential) list_of_distances(:)%exponential_stdev=&
					&SQRT(list_of_distances(:)%exponential_stdev/FLOAT(MAX((nsteps-1+sampling_interval)/sampling_interval,0)))
				ENDIF
				CALL report_distances()
				CALL finalise_distance()
			ELSE
				ERROR_CODE=ERROR_CODE_DEFAULT
			ENDIF
		END SUBROUTINE perform_distance_analysis

END MODULE DISTANCE
!--------------------------------------------------------------------------------------------------------------------------------!