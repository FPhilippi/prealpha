
!This is the main program unit. the MOLECULAR, DEBUG and ANGLES modules are general purpose, everything else is invoked as specified in general.inp
RECURSIVE SUBROUTINE initialise_global()
USE MOLECULAR
USE DEBUG
USE SETTINGS
IMPLICIT NONE
!default values
!variables
!PRIVATE/PUBLIC declarations
LOGICAL :: file_exists,connected,smalltask!'smalltask' means that the analysis can be run on a login node (usually)
LOGICAL :: anytask !is there any task at all?
LOGICAL :: wrapping_is_sensible!this variable is initialised to .TRUE., and set to .FALSE. as soom as something like msd is requested.
CHARACTER(LEN=*),PARAMETER :: filename_rmmvcf="rmmvcf.inp",filename_dihedral="dihedral.inp",filename_reorient="reorient.inp" !correlation module standard filenames
CHARACTER(LEN=*),PARAMETER :: filename_vc_components="vc_components.inp",filename_cacf_components="cacf_components.inp" !correlation module standard filenames
CHARACTER(LEN=*),PARAMETER :: filename_conductivity="conductivity.inp"!correlation module standard filenames
CHARACTER(LEN=*),PARAMETER :: filename_msd="diffusion.inp" !diffusion module standard filename
CHARACTER(LEN=*),PARAMETER :: filename_distribution="distribution.inp" !distribution module standard filename
CHARACTER(LEN=*),PARAMETER :: filename_distance="distance.inp" !distance module standard filename
INTEGER :: i,number_of_molecules!number_of_molecules is required for some safety checks, and initialised in generate_molecular_input()
INTEGER :: nsteps!nsteps is required again for checks (tmax...), and is initialised in generate_molecular_input()
	CALL set_default_masses()
	CALL set_default_charges()
	TIME_OUTPUT=TIME_OUTPUT_DEFAULT
	VERBOSE_OUTPUT=VERBOSE_OUTPUT_DEFAULT
	ERROR_OUTPUT=ERROR_OUTPUT_DEFAULT
	PARALLEL_OPERATION=PARALLEL_OPERATION_DEFAULT
	BOX_VOLUME_GIVEN=BOX_VOLUME_GIVEN_DEFAULT
	READ_SEQUENTIAL=READ_SEQUENTIAL_DEFAULT
	ERROR_CODE=ERROR_CODE_DEFAULT
	TIME_SCALING_FACTOR=TIME_SCALING_FACTOR_DEFAULT
	HEADER_LINES=HEADER_LINES_DEFAULT
	TRAJECTORY_TYPE=TRAJECTORY_TYPE_DEFAULT
	INFORMATION_IN_TRAJECTORY="UNK"
	WRAP_TRAJECTORY=WRAP_TRAJECTORY_DEFAULT
	number_of_molecules=-1
	nsteps=1000000000!expect the worst.
	smalltask=.TRUE.
	anytask=.FALSE.
	wrapping_is_sensible=.TRUE.
	ALPHABET_small=(/ (i,i=IACHAR("a"),IACHAR("a")+25,1) /) !just lowercase letters
	ALPHABET=(/ (i,i=IACHAR("a"),IACHAR("a")+25,1),& !lowercase letters
	&(i,i=IACHAR("A"),IACHAR("A")+25,1),& !uppercase letters
	&IACHAR("_"),IACHAR("/"),IACHAR("."),& !file stuff
	&(i,i=IACHAR("0"),IACHAR("0")+9,1) /)!... and some numbers. Numbers are always good.
	IF (DEVELOPERS_VERSION) THEN
		PRINT *, "   #######################"
		PRINT *, "   # DEVELOPER'S VERSION #"
		PRINT *, "   #######################"
		PRINT *
	ENDIF
	PRINT *, "   Copyright (C) !RELEASEYEAR! Frederik Philippi (Tom Welton Group)"
	PRINT *, "   Please report any bugs."
	PRINT *, "   Suggestions and questions are also welcome. Thanks."
	PRINT *, "   Date of Release: !RELEASEDATE!"
	PRINT *
	IF (DEVELOPERS_VERSION) THEN!only people who actually read the code get my contacts.
		PRINT *, "   Imperial College London"
		PRINT *, "   MSRH Room 601"
		PRINT *, "   White City Campus"
		PRINT *, "   80 Wood Lane"
		PRINT *, "   W12 0BZ London"
		PRINT *, "   f.philippi18"," at ","imperial.ac.uk"
		PRINT *
	ENDIF
	! first, check if file exists. If not, switch to user input for this part.
	INQUIRE(FILE=TRIM(FILENAME_GENERAL_INPUT),EXIST=file_exists)
	IF (file_exists) THEN
		PRINT *, "Note: The program *only* switches to user input if no 'general.inp' is present"
		PRINT *
		CALL read_general_input_header()
		!READ_SEQUENTIAL should be initialised here.
		IF (READ_SEQUENTIAL) PARALLEL_OPERATION=.FALSE.
		CLOSE(UNIT=7) !Bear in mind that this close statement is necessary, even though it is part of finalise_global.
	ELSE!Switch to user input.
		!The following lines of this subroutine serve for generating input files from user input.
		WRITE(*,*) "### FILE '",TRIM(FILENAME_GENERAL_INPUT),"' NOT AVAILABLE - SWITCH TO USER INPUT ###"
		USER_INPUT=.TRUE.
		PRINT *
		!as default, assume that the analysis is NOT necessary:
		SKIP_ANALYSIS=.TRUE.
		CALL user_general_input()
		IF (.NOT.(SKIP_ANALYSIS)) THEN
			!Restarting the whole thing!
			IF (VERBOSE_OUTPUT) THEN
				PRINT *
				PRINT *,"Restarting from the beginning now."
				PRINT *
			ENDIF
			CALL initialise_global()
		ENDIF
	ENDIF

	CONTAINS

		!reading the first lines of general.inp (plus the read_sequential line)
		SUBROUTINE read_general_input_header()
		IMPLICIT NONE
		INTEGER :: ios,n
		CHARACTER(LEN=32) :: inputstring
		CHARACTER(LEN=3) :: trajectory_type_input
		LOGICAL :: input_condition,trajectory_statement_absent
			IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(FILENAME_GENERAL_INPUT),"'"
			INQUIRE(UNIT=7,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=7)
			OPEN(UNIT=7,FILE=TRIM(FILENAME_GENERAL_INPUT),ACTION='READ',IOSTAT=ios)
			IF (ios/=0) CALL report_error(5,exit_status=ios)
			READ(7,IOSTAT=ios,FMT=*) FILENAME_TRAJECTORY
			IF (ios/=0) CALL report_error(5,exit_status=ios)
			READ(7,IOSTAT=ios,FMT=*) FILENAME_MOLECULAR_INPUT
			IF (ios/=0) CALL report_error(5,exit_status=ios)
			READ(7,IOSTAT=ios,FMT=*) PATH_TRAJECTORY
			IF (ios/=0) CALL report_error(5,exit_status=ios)
			READ(7,IOSTAT=ios,FMT=*) PATH_INPUT
			IF (ios/=0) CALL report_error(5,exit_status=ios)
			READ(7,IOSTAT=ios,FMT=*) PATH_OUTPUT
			IF (ios/=0) CALL report_error(5,exit_status=ios)
			!optional line: requesting serial read.
			!first, initialise to default:
			READ_SEQUENTIAL=READ_SEQUENTIAL_DEFAULT
			inputstring=""
			!then, try to find the corresponding section in the input file - change READ_SEQUENTIAL if necessary.
			DO n=1,MAXITERATIONS,1
				READ(7,IOSTAT=ios,FMT=*) inputstring
				IF (ios<0) THEN
					!end of file encountered
					EXIT
				ENDIF
				IF (ios==0) THEN
					IF ((TRIM(inputstring)=="sequential_read").OR.(TRIM(inputstring)=="read_sequential")) THEN
						!"sequential_read found in input file."
						WRITE(*,'(A45,I0)') " found a 'sequential_read' statement in line ",n+HEADER_LINES
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,input_condition
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) PRINT *,"Can't interpret line - setting sequential_read to default."
							READ_SEQUENTIAL=READ_SEQUENTIAL_DEFAULT!setting to default
						ELSE
							IF ((TRIM(inputstring)=="sequential_read").OR.(TRIM(inputstring)=="read_sequential")) THEN
								READ_SEQUENTIAL=input_condition
							ELSE
								CALL report_error(0)
							ENDIF
						ENDIF
						EXIT
					ELSEIF (TRIM(inputstring)=="quit") THEN
						EXIT
					ENDIF
				ENDIF
			ENDDO
			REWIND 7
			!try to get the trajectory type
			TRAJECTORY_TYPE=TRAJECTORY_TYPE_DEFAULT
			inputstring=""
			trajectory_statement_absent=.TRUE.
			DO n=1,MAXITERATIONS,1
				READ(7,IOSTAT=ios,FMT=*) inputstring
				IF (ios<0) THEN
					!end of file encountered
					EXIT
				ENDIF
				IF (ios==0) THEN
					IF (TRIM(inputstring)=="trajectory_type") THEN
						trajectory_statement_absent=.FALSE.
						!"trajectory_type found in input file."
						WRITE(*,'(A45,I0)') " found a 'trajectory_type' statement in line ",n
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,trajectory_type_input
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) PRINT *,"Can't interpret line - setting trajectory_type to default."
							TRAJECTORY_TYPE=TRAJECTORY_TYPE_DEFAULT!setting to default
						ELSE
							SELECT CASE (trajectory_type_input)
							CASE ("lmp")
								TRAJECTORY_TYPE="lmp"
							CASE ("xyz")
								TRAJECTORY_TYPE="xyz"
							CASE DEFAULT
								CALL report_error(51)!unknown trajectory format.
							END SELECT
						ENDIF
						EXIT
					ELSEIF (TRIM(inputstring)=="quit") THEN
						EXIT
					ENDIF
				ENDIF
			ENDDO
			!If no statement was present, then try to get the type from the extension.
			IF (trajectory_statement_absent) THEN
				SELECT CASE (FILENAME_TRAJECTORY(LEN(TRIM(FILENAME_TRAJECTORY))-3:LEN(TRIM(FILENAME_TRAJECTORY))))
				CASE (".lmp",".LMP")
					PRINT *,"Assuming lammps trajectory based on file extension."
					TRAJECTORY_TYPE="lmp"
				CASE (".xyz",".XYZ")
					PRINT *,"Assuming xyz trajectory based on file extension."
					TRAJECTORY_TYPE="xyz"
				CASE DEFAULT
					CALL report_error(51)!unknown trajectory format.
				END SELECT
			ENDIF
			REWIND 7
			!search for a 'wrap' statement
			WRAP_TRAJECTORY=WRAP_TRAJECTORY_DEFAULT
			inputstring=""
			DO n=1,MAXITERATIONS,1
				READ(7,IOSTAT=ios,FMT=*) inputstring
				IF (ios<0) THEN
					!end of file encountered
					EXIT
				ENDIF
				IF (ios==0) THEN
					IF (TRIM(inputstring)=="wrap_trajectory") THEN
						!"trajectory_type found in input file."
						WRITE(*,'(A45,I0)') " found a 'wrap_trajectory' statement in line ",n
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,input_condition
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) PRINT *,"Can't interpret line - setting wrap_trajectory to default."
							WRAP_TRAJECTORY=WRAP_TRAJECTORY_DEFAULT!setting to default
						ELSE
							IF (TRIM(inputstring)=="wrap_trajectory") THEN
								WRAP_TRAJECTORY=input_condition
							ELSE
								CALL report_error(0)
							ENDIF
						ENDIF
						EXIT
					ELSEIF (TRIM(inputstring)=="quit") THEN
						EXIT
					ENDIF
				ENDIF
			ENDDO
			!error report 73 also sets WRAP_TRAJECTORY to .FALSE.
			IF ((INFORMATION_IN_TRAJECTORY=="VEL").AND.(WRAP_TRAJECTORY)) CALL report_error(73)
		END SUBROUTINE read_general_input_header

		!This function checks if an input file is available and whether the user wants to overwrite it.
		LOGICAL FUNCTION check_if_inputfile_necessary(filename)
		IMPLICIT NONE
		CHARACTER(LEN=*) :: filename
		LOGICAL :: inputfile_exists
		INQUIRE(FILE=TRIM(filename),EXIST=inputfile_exists)
		check_if_inputfile_necessary=.FALSE.
		IF (inputfile_exists) THEN
			WRITE(*,*) "There is a file ",filename," - do you want to keep it? (y/n)"
			!IF yes, jump over that part. otherwise, open as overwrite.
			IF (user_input_logical()) THEN !Keep file
				WRITE(*,*) "Existing file is kept. You are responsible for ensuring proper format."
			ELSE
				check_if_inputfile_necessary=.TRUE.
			ENDIF
			!overwrite was requested!
		ELSE !no file present - make a new one.
			check_if_inputfile_necessary=.TRUE.
		ENDIF
		END FUNCTION check_if_inputfile_necessary

		SUBROUTINE hackprint()
		IMPLICIT NONE
		INTEGER :: ionpairs,ios
			PRINT *,"Hack mode. Supported are 'BMIMTFSI', 'N2O231TFSI' and 'N4441TFSI'."
			PRINT *,"How many ion pairs do you want?"
			ionpairs=user_input_integer(1,16384)
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE="N2O131TFSI.inp",IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'("100000 steps.")')
			WRITE(8,'("2 different types of molecules.")')
			WRITE(8,'("-1 15 ",I0," NTf2")') ionpairs
			WRITE(8,'("+1 38 ",I0," N(2O1)31")') ionpairs
			CLOSE(UNIT=8)
			OPEN(UNIT=8,FILE="N4441TFSI.inp",IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'("100000 steps.")')
			WRITE(8,'("2 different types of molecules.")')
			WRITE(8,'("-1 15 ",I0," NTf2")') ionpairs
			WRITE(8,'("+1 44 ",I0," N4441")') ionpairs
			CLOSE(UNIT=8)
			OPEN(UNIT=8,FILE="BMIMTFSI.inp",IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'("100000 steps.")')
			WRITE(8,'("2 different types of molecules.")')
			WRITE(8,'("+1 25 ",I0," BMIM")') ionpairs
			WRITE(8,'("-1 15 ",I0," NTf2")') ionpairs
			CLOSE(UNIT=8)
			OPEN(UNIT=8,FILE="BMIMTFSI_drudes.inp",IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'("100000 steps.")')
			WRITE(8,'("2 different types of molecules.")')
			WRITE(8,'("+1 25 ",I0," BMIM")') ionpairs
			WRITE(8,'("-1 24 ",I0," NTf2")') ionpairs
			WRITE(8,'("masses 1")')
			WRITE(8,'("X  0.400")')
			WRITE(8,'("constraints 1")')
			WRITE(8,'("1 15")')
			CLOSE(UNIT=8)
			PRINT *,"Printed with 100000 steps."
		END SUBROUTINE hackprint

		!This routine creates a new file 'general.inp' in unit 8 from user input.
		SUBROUTINE user_general_input()
		USE SETTINGS
		IMPLICIT NONE
		INTEGER :: n,nlow
			DO n=1,MAXITERATIONS,1
				IF (n==1) THEN
					PRINT *,"This program is an analyser for MD trajectories."
					PRINT *,"Please choose one of the following:"
				ELSE
					PRINT *
					PRINT *,"Please choose one of the following:"
				ENDIF
				IF (DEVELOPERS_VERSION) THEN
					PRINT *,"-1 - hack mode"
					nlow=-1
				ELSE
					nlow=0
				ENDIF
				PRINT *," 0 - exit"
				PRINT *," 1 - show the features of this software    (come back here)"!Complete so far.
				PRINT *," 2 - explain the structure of input files  (come back here)"!Complete. Don't change without changing 'show_program_features', too!
				PRINT *," 3 - generate input files from user input  (and then exit)"
				PRINT *," 4 - explain program flow / analysis       (come back here)"!Don't change without changing 'show_program_features', too!
				PRINT *," 5 - which format does the output take?    (come back here)"
				PRINT *," 6 - how to format the input trajectory?   (come back here)"
				SELECT CASE (user_input_integer(nlow,6))
				CASE (-1)
					CALL hackprint()
				CASE (0)!exit.
					EXIT
				CASE (1)!show the features of this software"
					CALL show_program_features()
				CASE (2)!explain the structure of input files
					CALL explain_input_files()
				CASE (3)!generate input files from user input
					CALL generate_all_input_files()
					EXIT
				CASE (4)!explain program flow / analysis
					CALL explain_program_flow()
				CASE (5)!which format does the output take? 
					CALL explain_output_format()
				CASE (6)!how to format the input trajectory?
					CALL explain_trajectory_format()
				CASE DEFAULT
					CALL report_error(0)
				END SELECT
			ENDDO
			!At this point, the user might have provided sufficient information.
			!If this is the case, SKIP_ANALYSIS is set to .FALSE.
			!Here, the user has the chance to manually skip the analysis anyway:
			IF (SKIP_ANALYSIS) THEN
				!no big job specified, but maybe a small one?
				IF ((smalltask).AND.(anytask)) THEN
					PRINT *,"You have requested relatively easy tasks."
					PRINT *,"Do you want to start these *now*? (y/n)"
					SKIP_ANALYSIS=(.NOT.(user_input_logical()))
				ENDIF
			ELSE
				PRINT *,"There should be sufficient input now."
				PRINT *,"Bear in mind that some tasks might be very involved."
				PRINT *,"Do you want to start the analysis? (y/n)"
				SKIP_ANALYSIS=(.NOT.(user_input_logical()))
				!Ask user if analysis should actually be started.
			ENDIF
			IF (n==MAXITERATIONS) WRITE(*,*) "Please take this seriously."
		END SUBROUTINE user_general_input

		!" 1 - show the features of this software"
		SUBROUTINE show_program_features()
		IMPLICIT NONE
			PRINT *,"The program has the following features:"
			PRINT *,
			PRINT *,"   Dihedral Conditions:"
			PRINT *,"   Allows the user to specify a set of dihedral conditions to be fulfilled."
			PRINT *,"   These could e.g. be the two dihedrals in NTf2 ('cisoid' vs. 'transoid'),"
			PRINT *,"   or some dihedrals along a side chain (= a certain conformation like 'all-trans')."
			PRINT *,"   It is possible to 'fold' the specified dihedrals (convenient for cisoid/transoid),"
			PRINT *,"   then on top of the range 'a' to 'b', also check for (360-b) to (360-a). Can be turned off."
			PRINT *,"   For these conditions, the following analyses are available:"
			PRINT *,"    - (Independent) incidences (or 'counts') for each specified dihedral"
			PRINT *,"    - Dependent incidences, i.e. the 2D PES subset (only for 2 dihedrals)"
			PRINT *,"    - For each timestep the share of fulfilled conditions (like, '42.0% transoid')"
			PRINT *,"    - The intermittent binary autocorrelation function of the specified condition"
			PRINT *,"     (from which e.g. the lifetime of 'cisoid' can be obtained)(PARALLELISED)"
			PRINT *,"   The encountered values of the specified dihedrals can also be exported in a separate file."
			PRINT *
			PRINT *,"   Orientation correlation functions: (PARALLELISED)"
			PRINT *,"   Computes the reorientational time correlation function for a given vector."
			PRINT *,"   The base and tip point of this vector is defined as fragment of a molecule (including single atoms)"
			PRINT *,"   Alternatively, when atomic charges are available, the dipole moment / charge arm vector can be used."
			PRINT *,"   Different legendre polynomials are available, and the computed quantity is Cl(t)=<Pl(u(t)*u(t=0))>"
			PRINT *,"   (u unit vector of fragment, P legendre polynomial of order l, t time shift)"
			PRINT *,"   see also equation (11.11.1) in 'THEORY OF SIMPLE LIQUIDS' (Hansen / McDonald), 4th edition."
			PRINT *,"   Or, for ionic liquids, for example equation (2) in DOI: 10.1016/j.cplett.2007.03.084"
			PRINT *
			PRINT *,"   Relative Mean Molecular Velocity Correlation Coefficients:"
			PRINT *,"   Unlike most other modules, this one needs atomic VELOCITIES instead of coordinates."
			PRINT *,"   Computes relative mean molecular velocity correlation coefficients based on:"
			PRINT *,"   Phys. Rev. E, 1994, 50, 1162–1170. DOI: 10.1103/PhysRevE.50.1162"
			PRINT *,"   No reference frame dependent properties are calculated (who needs these?)."
			PRINT *,"   Only two molecule types at once are supported currently."
			PRINT *,"   The following quantities from this reference have been implemented:"
			PRINT *,"    - RMM-VCFs, equation (4)"
			PRINT *,"    - The integral and the normalised function C12(t)"
			PRINT *,"   Optionally, these self-contributions can also be computed:"
			PRINT *,"    - lambdas(t), equation (8), for both specified particles (PARALLELISED)"
			PRINT *,"    - The integral and the normalised functions C1(t) and C2(t)"
			PRINT *,"    - Self-velocity correlations, eq (7), as well as C0(t)"
			PRINT *,"    - All corresponding diffusion quantities based on eq (17)"
			PRINT *,"    - The delta function as in eq (19)"
			PRINT *,"    - The time-dependent delta function, calculated as delta(t)=C12(t)-C0(t)"
			PRINT *,"    - A reference frame independent combination of distinct contributions Dd12"
			PRINT *,"   By that, everything in Table II. of the above reference is available."
			PRINT *,"   Additionally, conductivities are printed:"
			PRINT *,"    - Self, distinct and total contributions to the specific electrolytical conductivity"
			PRINT *,"    - The same for the molar conductivity (based on total particle number - *2 for ILs)"
			PRINT *,"    - Based on that, the predicted Haven Ratio in this framework of theory."
			PRINT *,"   The equations for electrolytical conductivity can be found in:"
			PRINT *,"   J. Chem. Phys., 1993, 99, 3983–3989. DOI: 10.1063/1.466191"
			PRINT *,"   Note that quite a large number of averages has to be taken to obtain sensible values."
			PRINT *
			PRINT *,"   Custom components of velocity correlation functions:"
			PRINT *,"   (and electric current autocorrelation function)"
			PRINT *,"   Unlike most other modules, this one needs atomic VELOCITIES instead of coordinates."
			PRINT *,"   You are left to choose the two molecule types to correlate freely."
			PRINT *,"   Additionally, you can choose to compute self- or distinct contributions."
			PRINT *,"   A good reference is, for example, J. Phys. Chem. B, 2011, 115, 13212–13221."
			PRINT *,"   Requesting velocity correlation functions will compute the quantities in <...> brackets"
			PRINT *,"   from equations (A6) to (A10) in above reference - i.e. not including the N, but Ncat, Nan,..."
			PRINT *,"   Requesting CACFs, however, will compute the quantities in <...> brackets"
			PRINT *,"   from equation (3). Essentially they are the same, but weighed with charges rather than 1/numbers"
			PRINT *,"   Important note: the distinct contributions are extremely expensive to calculate."
			PRINT *,"   If possible at all, it is advisable to use the RMM-VCFs, or even the 'conductivity_simple' keyword."
			PRINT *
			PRINT *,"   Mean (Squared) Displacement 1: (PARALLELISED)"
			PRINT *,"   Calculates the mean squared displacement including a drift correction."
			PRINT *,"   (other exponents can be chosen as well, e.g. for the mean fourth power displacement)"
			PRINT *,"   The diffusion coefficients obtained by that can be used for comparison with the VACFs."
			PRINT *,"   Different projections can by chosen by which the displacement vector is to be multiplied."
			PRINT *,"   This could be e.g. '1 1 1' (giving the 'standard' 3D diffusion coefficient),"
			PRINT *,"   or something like '0 0 1' (which would give only the component in z-direction)."
			PRINT *,"   Two print levels are available. Default is to only print:"
			PRINT *,"    - The mean squared displacement <R²>"
			PRINT *,"    - The mean displacement <R>"
			PRINT *,"   When the verbose print is requested, then the output additionally contains:"
			PRINT *,"    - Drift corrected mean squared displacement <R²>-<R>²"
			PRINT *,"    - All three components of the drift vector, <x>, <y>, and <z>"
			PRINT *,"    - The number of averages taken to obtain these values."
			PRINT *
			PRINT *,"   Mean (Squared) Displacement 2: (PARALLELISED)"
			PRINT *,"   It is also possible to calculate the relative mean molecular diffusion coefficients,"
			PRINT *,"   using the corresponding Einstein relation."
			PRINT *
			PRINT *,"   Distribution Functions: (PARALLELISED)"
			PRINT *,"   Calculates distribution functions referenced to a fixed vector in space."
			PRINT *,"   cylindrical and polar coordinate systems are supported."
			PRINT *,"   Different references vectors can be chosen, including a randomised reference vector (debug purposes)."
			PRINT *,"   This could be e.g. '0 0 1' (z-direction is the reference)."
			PRINT *,"   In both coordinate systems, values are averaged over all azimuthal angles."
			PRINT *,"   This feature is also capable of calculating Sum Rules and Coulomb interaction energies."
			PRINT *,"   Number integrals and radial distribution functions are also printed."
			PRINT *
			PRINT *,"   Average distances:"
			PRINT *,"   Calculates the average closest distance between a reference atom and an observed atom."
			PRINT *,"   Reference and observed atom can either be specific (i.e. a certain molecule_type_index and atom_index),"
			PRINT *,"   or they can be of a certain type (such as 'H' or 'F')."
			PRINT *,"   Two modi are available: intermolecular or intramolecular distances. Also available are:"
			PRINT *,"    - exponentially weighed averaged distances ( weighed with exp(-k*r) )"
			PRINT *,"    - distances used in FFC theory, using eq (6) and (13) in the ESI of J. Phys. Chem. Lett., 2020, 11, 2165–2170."
			PRINT *,"    - standard deviations of closest distances and exponentially weighed distances."
			PRINT *
			PRINT *,"Each of these blocks is treated as a distinct feature with its own input file."
			PRINT *,"Some of the more demanding routines exist also in a parallelised version."
			PRINT *,"(see option '4' in the main menu for more information about parallelisation)"
			PRINT *
			PRINT *,"For most of these features, a number of switches are available."
			PRINT *,"These switches influence e.g. the print level, bin counts or steps to analyse."
			PRINT *,"The example input files (option '5' in the main menu) contain the most common switches."
			PRINT *,"Information about all possible switches is provided by option '2' in the main menu."
		END SUBROUTINE show_program_features

		!" 2 - explain the structure of input files"
		SUBROUTINE explain_input_files()
		IMPLICIT NONE
		CHARACTER(LEN=1),PARAMETER :: doublequote='"'
			PRINT *,"The input files have to follow the following format:"
			PRINT *
			PRINT *,"General input file:"
			PRINT *,"The 'general.inp' file is the main input file, located in the same folder as the executable."
			PRINT *,"It is possible to specify other names for the general input file as command line arguments."
			PRINT *,"When multiple general input files are specified, they will be invoked subsequently."
			PRINT *,"If the file / one of these files isn't found, then the program switches to user input."
			WRITE(*,'(" It is read line wise, with the first ",I0," lines being reserved (and strictly fixed).")') HEADER_LINES
			PRINT *,"The content of these lines is:"
			PRINT *," line 1 - the filename of the trajectory, e.g. 'trajectory.lmp'"
			PRINT *," line 2 - the name of the molecular input file, e.g. 'mymolecule.inp'"
			PRINT *," line 3 - Path to the trajectory"
			PRINT *," line 4 - Path to the input files other than the general and molecular input files."
			PRINT *," line 5 - Output folder path"
			WRITE(*,'(" Path names have to be enclosed in quotes ",A1)') doublequote
			PRINT *,"The body of 'general.inp' is read line-wise, and the"
			PRINT *,"program finishes when either 'quit' or the end of file is encountered."
			PRINT *,"Each line contains a switch or keyword, followed by an argument (if required)"
			PRINT *,"Only the necessary information is read from any line, with the rest being ignored."
			PRINT *,"Be aware that keywords affect only the lines below them."
			PRINT *,"This is with the exception of sequential_read, trajectory_type and wrap_trajectory."
			PRINT *,"These latter three act on the whole analysis, no matter where specified."
			PRINT *,"Only their first occurence matters - everything afterwards is ignored."
			PRINT *,"An incorrectly formatted 'general.inp' is not tolerated (read the error messages)."
			PRINT *,"For many switches, a simple mode that doesn't require additional input is available."
			PRINT *,"The simple mode is requested by appending '_simple' to the switch, e.g. 'gyradius_simple'."
			PRINT *,"Available switches are: (case-sensitive, everything is lowercase)"
			PRINT *," - 'sequential_read':"
			PRINT *,"    If true 'T', then the trajectory is read line by line."
			PRINT *,"    This is slow, but requires only the minimum amount of RAM."
			PRINT *,"    Not recommended for mean-squared displacement and VACFs."
			PRINT *,"    If false 'F', then the whole trajectory is read into RAM."
			PRINT *,"    This is the first switch that affects every line, not just the ones after it."
			PRINT *," - 'trajectory_type':"
			PRINT *,"    expects either 'xyz' or 'lmp' as string input."
			PRINT *,"    This is the second switch that affects every line, not just the ones after it."
			PRINT *," - 'wrap_trajectory':"
			PRINT *,"    Expects one logical. If 'T', then molecules are wrapped into the box."
			PRINT *,"    (based on their centre of mass. Might not be sensible for some analyses.)"
			PRINT *,"    This is the third and last switch affecting every line."
			PRINT *," - 'parallel_operation':"
			PRINT *,"    Turns parallelisation on (T) or off (F)."
			PRINT *,"    Parallelisation is only available with 'sequential_read F'"
			PRINT *," - 'set_threads':"
			PRINT *,"    Sets the number of threads to use. 'set_threads 0' uses all available threads."
			PRINT *," - 'error_output':"
			PRINT *,"    Turns error output on (T) or off (F)."
			PRINT *," - 'time_scaling':"
			PRINT *,"    Takes an integer value, by which the timestep is multiplied. For example,"
			PRINT *,"    specify 'time_scaling 1000' if your trajectory is dumped every 1000fs."
			PRINT *," - 'set_prefix':"
			PRINT *,"    The specified prefix is prepended to the output files."
			PRINT *,"    Useful if, for example, the dihedral analysis is specified multiple times."
			PRINT *," - 'dump_example':"
			PRINT *,"    Writes an xyz file of every specified molecule type into the output folder."
			PRINT *,"    Can be used to extract the atom numbers for the dihedral analysis."
			PRINT *," - 'dump_snapshot': (simple mode available)"
			PRINT *,"    Expects an integer and a logical and dumps the specified timestep as .xyz file."
			PRINT *,"    If the logical is 'T', then every molecule is written into a separate file."
			PRINT *," - 'dump_split': (simple mode available)"
			PRINT *,"    Splits the trajectory into separate files for every molecule type (centred to centre of mass!)."
			PRINT *,"    Expects two integers: the first timestep and the last timestep."
			PRINT *," - 'dump_single':"
			PRINT *,"    Writes a trajectory containing just one single molecule."
			PRINT *,"    This keyword expects a logical, followed by four integers in the same line:"
			PRINT *,"    If the logical is 'T', then the molecule is centred to its centre-of-mass."
			PRINT *,"    The first and second integers specify the first and last timestep to write."
			PRINT *,"    The third and fourth integers are the molecule type index and the molecule index, respectively."
			PRINT *," - 'contact_distance': (simple mode available)"
			PRINT *,"    Reports the smallest intra- and intermolecular distances and the largest intramolecular distance."
			PRINT *,"    This keyword expects two integers as input: the timestep to analyse and the molecule type index."
			PRINT *,"    If a molecule type index of 0 is specified, then all molecule types are considered."
			PRINT *," - 'dump_cut':"
			PRINT *,"    like dump_single - but the surrounding molecules are also written."
			PRINT *,"    This keyword expects a logical, followed by four integers and one real in the same line:"
			PRINT *,"    If the logical is 'T', then the molecule is centred to its centre-of-mass in every step."
			PRINT *,"    The first and second integers specify the first and last timestep to write."
			PRINT *,"    The third and fourth integers are the molecule type index and the molecule index, respectively."
			PRINT *,"    The real number defines the cutoff for centre-of-mass distance for exporting molecules"
			PRINT *,"    Note that the properly wrapped mirror images of the closest encounters are given."
			PRINT *," - 'dump_dimers'"
			PRINT *,"    Dumps the closest molecule of type X around all molecules of type Y for a certain timestep."
			PRINT *,"    Expects a logical, followed by three integers."
			PRINT *,"    If the logical is (T), then the output is combined in a single 'trajectory'-like xyz file."
			PRINT *,"    Otherwise (F), one output file is written per dimer."
			PRINT *,"    The first integer is the timestep. The following two integers are the types Y and X, respectively"
			PRINT *," - 'dump_neighbour_traj': (simple mode available)"
			PRINT *,"    Dumps the N closest molecules of type X around a certain molecule M of type Y as trajectory."
			PRINT *,"    This keyword expects six integers:"
			PRINT *,"    The first and second integers specify the first and last timestep to write."
			PRINT *,"    The third and fourth integers are the molecule type index Y and the molecule index M, respectively."
			PRINT *,"    The fith integer is the integer of the neighbour molecule to consider"
			PRINT *,"    The last integer is the number of neighbours to write."
			PRINT *," - 'cubic_box_edge':"
			PRINT *,"    this keyword expects two real values, the lower and upper bounds of the simulation box."
			PRINT *,"    i.e. cubic_box_edge 0.0 100.0 corresponds to a cubic box with side length 100.0 Angströms"
			PRINT *,"    useful if e.g. dump_cut is used with a xyz trajectory."
			PRINT *," - 'convert': (simple mode available)"
			PRINT *,"    converts the given trajectory to a centre-of-mass trajectory (per specified molecule type)."
			PRINT *,"    i.e. only the centres of mass for the molecules are printed instead of the atoms."
			PRINT *,"    This keyword expects a logical. If (T), then a new, modified molecular input file is written as well."
			PRINT *," - 'temperature': (simple mode available)"
			PRINT *,"    Computes the instantaneous temperature of a particular molecule type."
			PRINT *,"    This keyword expects exactly three integers:"
			PRINT *,"    The molecule type index, and the range of analysis, given as first step and last step."
			PRINT *,"    If a molecule type index of 0 is specified, then all molecule types are considered."
			PRINT *," - 'drude_temp': (simple mode available)"
			PRINT *,"    Computes drude, centre of mass, and total temperature of the whole box."
			PRINT *,"    This keyword computes equation (13), (14) and (15) in:"
			PRINT *,"    J. Phys. Chem. Lett., 2019, 10, 7523–7530. DOI 10.1021/acs.jpclett.9b02983"
			PRINT *,"    Support of drude particles requires them to be read in manually,"
			PRINT *,"    since the automatic drude particle assignment is only available for position trajectories."
			PRINT *,"    This keyword expects exactly two integers:"
			PRINT *,"    The range of analysis, given as first step and last step."
			PRINT *," - 'remove_drudes': (simple mode available)"
			PRINT *,"    writes a new trajectory, with drude particles merged into their respective cores."
			PRINT *,"    (requires assigned drude particles, either manually or automatically)"
			PRINT *,"    This keyword expects exactly two integers:"
			PRINT *,"    The range of analysis, given as first step and last step."
			PRINT *," - 'remove_cores': (simple mode available)"
			PRINT *,"    writes a new trajectory only with drude particles minus the positions of their respective cores."
			PRINT *,"    (requires assigned drude particles, either manually or automatically)"
			PRINT *,"    This keyword expects exactly two integers:"
			PRINT *,"    The range of analysis, given as first step and last step."
			PRINT *," - 'gyradius': (simple mode available)"
			PRINT *,"    Computes the ensemble averages and standard deviations of radius of gyration,"
			PRINT *,"    radius of gyration squared, and maximum distance of any atom in a molecule from its centre of mass."
			PRINT *,"    This keyword expects exactly three integers:"
			PRINT *,"    The molecule type index, and the range of analysis, given as first step and last step."
			PRINT *,"    If a molecule type index of 0 is specified, then all molecule types are considered."
			PRINT *," - 'jump_velocity': (simple mode available)"
			PRINT *,"    Computes a histogram / probability distribution of jump velocities."
			PRINT *,"    This keyword expects exactly three integers:"
			PRINT *,"    The molecule type index, the range of analysis, given as first step and last step,"
			PRINT *,"    and the maximum jump length / shift given as number of timesteps."
			PRINT *,"    If a molecule type index of 0 is specified, then all molecule types are considered."
			PRINT *," - 'show_settings':"
			PRINT *,"    Writes settings and useful information to the standard output"
			PRINT *," - 'print_atomic_masses':"
			PRINT *,"    Writes atomic masses to the standard output, using the molecular input file format."
			PRINT *," - 'print_atomic_charges':"
			PRINT *,"    Writes atomic charges to the standard output, using the molecular input file format."
			PRINT *," - 'print_dipole_statistics':"
			PRINT *,"    outputs average/minimum/maximum/standarddev of dipole moment for the first timestep."
			PRINT *,"    for charged molecules, the vector from center of mass to center of charge is used."
			PRINT *," - 'show_drude':"
			PRINT *,"    Writes detailed current information about drude particles."
			PRINT *," - 'switch_to_com'"
			PRINT *,"    Irreversibly switches to the barycentric reference frame."
			PRINT *,"    i.e. for all analyses below this keyword, the box' centre-of-mass is removed in every step."
			PRINT *,"    (cannot be turned off again, until the program switches to the next general input file)"
			PRINT *," - 'conductivity_simple'"
			PRINT *,"    This is a special keyword that computes the overall electrical conductivity of the whole system."
			PRINT *,"    It uses the autocorrelation module, but is much faster than the CACFs."
			PRINT *," - 'verbose_output':"
			PRINT *,"    Turned on (T) by default. If (F), then only very limited output is obtained."
			PRINT *," - 'time_output':"
			PRINT *,"    Turns the timing on (T) or off (F)."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *,"These keywords require separate input files (explained below):"
			PRINT *," - 'conductivity'  (requests feature 'cacf components' or 'conductivity_simple')"
			PRINT *," - 'velocity'      (requests feature 'vcf components')"
			PRINT *," - 'rmm-vcf'       (requests feature 'Relative Mean Molecular Velocity Correlation Coefficients')"
			PRINT *," - 'diffusion'     (requests feature 'Mean Squared Displacement') (simple mode available)"
			PRINT *," - 'dihedral'      (requests feature 'Dihedral Conditions')"
			PRINT *," - 'reorientation' (requests feature 'reorientational time correlation')"
			PRINT *," - 'distribution'  (requests feature 'polar/cylindrical distribution function) (simple mode available)"
			PRINT *," - 'distance'      (requests feature 'Average distances') (simple mode available)"
			PRINT *
			PRINT *,"Molecular input file:"
			PRINT *,"This file contains information about the system, located in the same folder as the executable"
			PRINT *,"The first line is the number of timesteps,"
			PRINT *,"followed by the number of molecule types in the second line."
			PRINT *,"For every molecule type, the following information is read:"
			PRINT *,"Charge - Number of atoms per molecule - number of molecules."
			PRINT *,"The program expects as many lines as there are molecule types."
			PRINT *,"Following this fixed section, the rest of the input file is read."
			PRINT *,"(Until either a 'quit' statement or the end of file are encountered)"
			PRINT *,"In that free-format section, the following optional subsections can be placed:"
			PRINT *," - 'default_masses':"
			PRINT *,"    this keyword triggers the specification of custom default masses."
			PRINT *,"    it expects an integer, which is the number of subsequent lines to read."
			PRINT *,"    This is available for single lowercase letters (a,b,c,...,z) and element names."
			PRINT *,"    (Including 'X' and 'D', which are treated as drude particles)"
			PRINT *,"    If e.g. the trajectory contains an anion of mass 123.4, abbreviated as 'a',"
			PRINT *,"    and a cation of mass 432.1, abbreviated as 'c', then this section should be added:"
			PRINT *,"      masses 2"
			PRINT *,"      a 123.4"
			PRINT *,"      c 432.1"
			PRINT *,"    Furthermore, the support of drude particles can be turned on by adding:"
			PRINT *,"      masses 1"
			PRINT *,"      X 0.4"
			PRINT *,"    Note that drude particles are added to N,O,C,S,P,Li,F - but not Hydrogen."
			PRINT *,"    This keyword changes the defaults, i.e. ALL atoms of type a,X,P,.... see also:"
			PRINT *," - 'atomic_masses':"
			PRINT *,"    this keyword triggers the specification of custom atomic masses."
			PRINT *,"    it expects an integer, which is the number of subsequent lines to read."
			PRINT *,"    each line must have two integers and one real:"
			PRINT *,"    molecule type index, atom index, and atomic mass."
			PRINT *,"    This allows the user to specify atomic weights, which can be used two introduce charges."
			PRINT *,"    The order is important:"
			PRINT *,"      1) The program starts with its own defaults, such as 12.011 for carbon."
			PRINT *,"      2) If necessary, these values are changed by 'default_masses'"
			PRINT *,"      3) Any positive drude mass 'X' or 'D', if present, is subtracted from N,O,C,S,P,Li,F."
			PRINT *,"      4) After that, masses of particular atoms are overwritten by 'atomic_masses'."
			PRINT *," - 'default_charges' and 'atomic_charges'"
			PRINT *,"    These two keywords can be used to specify atomic charges."
			PRINT *,"    Their syntax follows that of the keywords 'default_masses' and 'atomic_masses', respectively."
			PRINT *,"    (With the exception that lowercase letters are not accepted)"
			PRINT *,"    The order is as follows:"
			PRINT *,"      1) The program starts with its own default charge - 0.0 for every atom."
			PRINT *,"      2) If necessary, these values are changed by 'default_charges'"
			PRINT *,"      3) Charges of particular atoms are overwritten by 'atomic_charges'."
			PRINT *," - 'constraints':"
			PRINT *,"    this keyword triggers the specification of custom constraints."
			PRINT *,"    it expects an integer, which is the number of subsequent lines to read."
			PRINT *,"    each of these subsequent lines has to contain two integers:"
			PRINT *,"    First, the molecule type index, and second, the number of constraints."
			PRINT *," - 'drudes':"
			PRINT *,"    This keyword is used to manually assign drude particles to their respective core."
			PRINT *,"    it expects an integer, which is the number of subsequent lines to read."
			PRINT *,"    each drude particle is assigned by giving three integers (per line):"
			PRINT *,"    the molecule type index - atom index core - atom index drude."
			PRINT *
			PRINT *,"relative mean molecular velocity correlation input file:"
			PRINT *,"The two molecules to correlate have to be given in the first line."
			PRINT *,"'rmm-vcf' is given in the second line, indicating the type of analysis."
			PRINT *,"Switches are read from the following lines. Available are:"
			PRINT *," - 'tmax':"
			PRINT *,"    expects an integer, which is then taken as the maximum number of steps"
			PRINT *,"    into the future for the autocorrelation function (the shift, so to say)."
			PRINT *," - 'skip_autocorrelation':"
			PRINT *,"    If yes (T), then only the cross-contributions are calculated."
			PRINT *," - 'sampling_interval':"
			PRINT *,"    Expects an integer. Every so many steps will be used as origin to compute self-contributions."
			PRINT *,"    These are usually computationally more expensive, but need less averaging."
			PRINT *,"    Note that the printed tcf will always have the same time resolution as the trajectory.."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *,
			PRINT *,"custom velocity correlation AND cacf components input file:"
			PRINT *,"The first line gives the operation mode, followed by the number of custom components."
			PRINT *,"The operation mode is either 'vcf' or 'cacf'."
			PRINT *,"the custom compontents are specified below, one per line. Each line must contain, in this order:"
			PRINT *," - the molecule type index of the first, reference molecule"
			PRINT *," - the molecule type index of the second, observed molecule"
			PRINT *," - 'T' if self contributions are to calculate, and 'F' for distinct contributions."
			PRINT *,"Thus, the following 7 lines specify every unique vcf in a system with two constituents:"
			PRINT *,"  vcf 6"
			PRINT *,"  1 1 T"
			PRINT *,"  1 1 F"
			PRINT *,"  1 2 F"
			PRINT *,"  2 1 F"
			PRINT *,"  2 2 T"
			PRINT *,"  2 2 F"
			PRINT *,"(note that '1 2 F' and '2 1 F' are redundant and will/should give the same value)"
			PRINT *,"'rmm-vcf' is given in the second line, indicating the type of analysis."
			PRINT *,"Switches are then read from the following lines. Available are:"
			PRINT *," - 'tmax':"
			PRINT *,"    expects an integer, which is then taken as the maximum number of steps"
			PRINT *,"    into the future for the autocorrelation function (the shift, so to say)."
			PRINT *," - 'sampling_interval':"
			PRINT *,"    Expects an integer. Every so many steps will be used as origin of the correlation functions."
			PRINT *,"    Note that the printed tcf will always have the same time resolution as the trajectory.."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *,"A rule of thumb as final remark: everything with an 'F' will be very, very, very slow."
			PRINT *,"This is because the double sum over every distinct particle is evaluated."
			PRINT *,"It is possible (and a lot faster) to just calculate the overall electrical conductivity."
			PRINT *,"To request this, just put 'conductivity' in the first line."
			PRINT *
			PRINT *,"diffusion input file for self diffusion:"
			PRINT *,"The first line contains the expression 'msd', followed by the number of projections N."
			PRINT *,"The latter are read from the following N lines. The format of each line is:"
			PRINT *,"x - y - z - number of the molecule type."
			PRINT *,"For the 'standard' 3D diffusion of molecule type 2, the line would thus be '1 1 1 2'."
			PRINT *,"After the projections have been specified, switches can be specified in an arbitrary order."
			PRINT *,"Available are:"
			PRINT *," - 'tmax':"
			PRINT *,"    Expects an integer, which is then taken as the maximum number of steps"
			PRINT *,"    into the future for the mean squared displacement."
			PRINT *," - 'tstep':"
			PRINT *,"    The given integer is taken as the step size. i.e. if 'tstep 10' is specified,"
			PRINT *,"    then only shifts by 1,10,20,...,tmax are computed."
			PRINT *," -  'print_verbose':"
			PRINT *,"    If yes (T), then the detailed drift is printed, too."
			PRINT *," - 'exponent':"
			PRINT *,"    Expects an integer, which is used as exponent in the displacement."
			PRINT *,"    Thus, 'exponent 2' is the normal MSD, but others can be used, too."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *
			PRINT *,"diffusion input file for cross diffusion:"
			PRINT *,"The first line contains the expression 'cross', followed by the number of projections N."
			PRINT *,"The latter are read from the following N lines. The format of each line is:"
			PRINT *,"x - y - z - number of the first molecule type - number of the second molecule type."
			PRINT *,"To get the cross diffusion between the first two molecule types, the line would thus be '1 1 1 1 2'."
			PRINT *,"After the projections have been specified, switches can be specified in an arbitrary order."
			PRINT *,"Available are:"
			PRINT *," - 'tmax':"
			PRINT *,"    Expects an integer, which is then taken as the maximum number of steps"
			PRINT *,"    into the future for the mean squared displacement."
			PRINT *," - 'tstep':"
			PRINT *,"    The given integer is taken as the step size. i.e. if 'tstep 10' is specified,"
			PRINT *,"    then only shifts by 1,10,20,...,tmax are computed."
			PRINT *," - 'exponent':"
			PRINT *,"    Expects an integer, which is used as exponent in the displacement."
			PRINT *,"    Thus, 'exponent 2' is the normal MSD, but others can be used, too."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *
			PRINT *,"distribution input file:"
			PRINT *,"The first line contains the expression 'cdf', 'pdf' or 'charge_arm', followed by the number of references N."
			PRINT *,"'cdf' requests the cylindrical distribution function, 'pdf' requests the polar distribution function."
			PRINT *,"'charge_arm' requests a polar distribution function of the charge arm (ions) or dipole moment (neutral molecules)."
			PRINT *,"The references are read from the N lines following the first line. The format of each line is:"
			PRINT *,"x - y - z - number of reference molecule type - number of observed molecule type."
			PRINT *,"For molecule type 1 around 2 relative to z-direction, the line would thus contain '0 0 1 2 1'."
			PRINT *,"a 'zero' reference vector, i.e. '0 0 0', triggers the randomisation of the reference vector."
			PRINT *,"If the molecule type of the observed molecule is -1, then *all* molecule types are used."
			PRINT *,"Note that for 'charge_arm', only one molecule type is required - no observed molecule type is required,"
			PRINT *,"only the reference type. The charge arm pdf is NOT corrected for azimuthal or radial parts - only polar."
			PRINT *,"After the projections have been specified, switches can be specified in an arbitrary order."
			PRINT *,"Available are:"
			PRINT *," - 'bin_count':"
			PRINT *,"    Expects an integer, which is then used as bin count for both independent variables."
			PRINT *," - 'bin_count_a':"
			PRINT *,"    Expects an integer, which is then used as bin count for variable a)."
			PRINT *," - 'bin_count_b':"
			PRINT *,"    Expects an integer, which is then used as bin count for variable b)."
			PRINT *," - 'maxdist':"
			PRINT *,"    Expects a real value, which is taken as the cutoff distance of molecule pairs to be considered."
			PRINT *," -  'maxdist_optimize':"
			PRINT *,"    sets 'maxdist' to half the box size where available. doesn't need additional input values."
			PRINT *,"    for charge arm and dipole moment analyses, 'maxdist' will be set to the maximum value in the first timestep."
			PRINT *,"    (considering all molecule types specified as references, rounded to 1 digit)"
			PRINT *," - 'subtract_uniform'"
			PRINT *,"    If yes (T), then the uniform density / radial distribution function is subtracted."
			PRINT *,"    Only available for the polar distribution function."
			PRINT *," - 'weigh_charge'"
			PRINT *,"    The charges of observed molecules are added to the distribution histogram, rather than unity."
			PRINT *," - 'normalise_CLM'"
			PRINT *,"    The charge arm is divided by (M*Rgy**2) (charge lever moment correction)."
			PRINT *,"    Here, M is the mass of the molecule, and Rgy is the radius of gyration."
			PRINT *,"    Only available for the charge arm distribution function."
			PRINT *," - 'sampling_interval'"
			PRINT *,"    Expects an integer. Every so many steps are used for the analysis."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *
			PRINT *,"dihedral input file:"
			PRINT *,"The molecule type index (= the molecule to observe) is given in the first line."
			PRINT *,"The expession 'dihedrals', followed by the number of dihedral conditions, in 2nd line"
			PRINT *,"For every dihedral condition follows one line, giving the atoms (in that molecule)"
			PRINT *,"which are part of the dihedral, as well as the lower and upper bound."
			PRINT *,"'1 2 3 4 0.0 90.0' thus means that dihedral 1-2-3-4 has to be between 0 and 90 degrees."
			PRINT *,"Important: dihedrals are defined from 0.0° to 360.0°"
			PRINT *,"ALL specified conditions have to be fulfilled simultaneously for the h operator to become true."
			PRINT *,"After the condition section, the following switches may follow:"
			PRINT *," - 'tmax':"
			PRINT *,"    Expects an integer, which is then taken as the maximum number of steps"
			PRINT *,"    into the future for the intermittent binary autocorrelation function."
			PRINT *," - 'export':"
			PRINT *,"    Requires one integer (the index of the molecule) as input."
			PRINT *,"    All specified dihedrals for this particular molecule will be exported in an output file."
			PRINT *,"    Note that 'export' can be specified more than once!"
			PRINT *," - 'fold':"
			PRINT *,"    If true (T), then apart from the dihedrals being in the range a to b,"
			PRINT *,"    also check for the range (360-b) to (360-a)."
			PRINT *," - 'dump_verbose':"
			PRINT *,"    If true (T), then also report PES subset population and <h> as a function of the timestep."
			PRINT *,"    also check for the range (360-b) to (360-a)."
			PRINT *," - 'skip_autocorrelation':"
			PRINT *,"    If true (T), then the actual autocorrelation analysis is skipped."
			PRINT *,"    This is useful if only the PES is required."
			PRINT *," - 'bin_count':"
			PRINT *,"    Sets the bin count to the specified integer."
			PRINT *,"    e.g. 'bin_count 36' equals to binning in steps of 10°."
			PRINT *," - 'jump_analysis'"
			PRINT *,"    Performs an analysis of average jump velocities as a function of:"
			PRINT *,"      a) The number of changes between fulfilment/non-fulfilment of the dihedral condition."
			PRINT *,"      b) The share of fulfilled dihedral conditions over the specified jump time."
			PRINT *,"    Requires one integer, i.e. the desired jump time (in number of timesteps)."
			PRINT *,"    It is possible to request multiple jump analyses in one dihedral input file."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *
			PRINT *,"reorientation input file:"
			PRINT *,"The molecule type index (= the molecule to observe) is given in the first line."
			PRINT *,"The expession 'reorientation' in the second line to request the appropriate analysis"
			PRINT *,"For the reorientation analysis, the user has the choice between two vectors:"
			PRINT *,"  a) the vector from a base fragment to a tip fragment. Both *must* be defined as outlined below."
			PRINT *,"  b) the charge arm or dipole moment vector. This required atomic charges (cf. molecular input file)."
			PRINT *,"Regarding method a):"
			PRINT *,"A fragment is defined by the expression 'base' or 'tip', followed by the number of atoms in this fragment."
			PRINT *,"The immediately following line must contain a list of the atom indices in this fragment."
			PRINT *,"For example, these lines define atom 16 as the base fragment and atoms 1, 3 and 4 as the tip fragment:"
			PRINT *,"  base 1"
			PRINT *,"  16"
			PRINT *,"  tip 3"
			PRINT *,"  3 1 4"
			PRINT *,"The two fragments must appear before the quit statement (if applicable)."
			PRINT *,"Method b) only requires 'charge_arm' or 'dipole', but no further input (apart from the atomic charges)."
			PRINT *,"The following switches may be used as well:"
			PRINT *," - 'tmax':"
			PRINT *,"    Expects an integer, which is then taken as the maximum number of steps"
			PRINT *,"    into the future for the time correlation function."
			PRINT *," - 'legendre':"
			PRINT *,"    Expects an integer, which defines the order of the legendre polynomial to use."
			PRINT *," - 'sampling_interval':"
			PRINT *,"    Expects an integer. Every so many steps will be used as starting point for the tcf."
			PRINT *,"    Note that the printed tcf will always have the same time resolution as the trajectory."
			PRINT *," - 'quit'"
			PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
			PRINT *
			PRINT *,"Distance Input File:"
			PRINT *,"The first line must state the type of analysis, i.e. either 'intramolecular' or 'intermolecular',"
			PRINT *,"followed by the number of subjobs as an integer."
			PRINT *,"The following lines, one for each subjob, describe which atoms are reference and observed atoms."
			PRINT *,"For intramolecular analysis, you need to specify one of these for every subjob:"
			PRINT *," - Three integers: the molecule type index and the two atom indices for reference and observed atom, for example "
			PRINT *,"   '3 1 2' calculates closest intramolecular distances for the atoms with indices 1 and 2 in molecule type 3."
			PRINT *," - two integers and an element name, acting as wildcard."
			PRINT *,"   The program will automatically consider all atoms of that type."
			PRINT *,"   for example '3 1 H' uses all H atoms in molecule type 3 as observed atoms."
			PRINT *," - two element names, such as 'H H'."
			PRINT *,"   Note that this type of analysis averages over all possible intramolecular combinations,"
			PRINT *,"   if you have more than one molecule type containing hydrogen atoms,"
			PRINT *,"   then you will obtain an average over all of these."
			PRINT *,"The same principles apply for intermolecular distance analyses,"
			PRINT *,"only that the second molecule type index needs to be given:"
			PRINT *," - Four integers: the molecule type index and atom index for the reference atom"
			PRINT *,"   and the molecule type index and atom index for the observed atom."
			PRINT *,"   Note that even if the molecule type indices are the same, the analysis is still intermolecular."
			PRINT *,"   (and thus makes only sense if there are at least two molecules of that type)."
			PRINT *," - two integers and an element name, acting as wildcard."
			PRINT *,"   The program will automatically consider all atoms of that type."
			PRINT *,"   For example '3 1 H' uses all H atoms in all molecule types as observed atoms."
			PRINT *," - two element names, such as 'H H'."
			PRINT *,"   Note that this type of analysis averages over all possible intermolecular combinations."
			PRINT *,"After the subjob section, the following switches may be used:"
			PRINT *," - 'quit' Terminates the analysis. Lines after this switch are ignored."
			PRINT *," - 'maxdist': Expects one real number, which is the maximum distance / cutoff to consider."
			PRINT *," - 'maxdist_optimize': sets the cutoff to half the box size."
			PRINT *," - 'nsteps': Followed by one integer, which is interpreted as the highest timestep number to consider."
			PRINT *," - 'sampling_interval': Expects one integer - the sampling interval,"
			PRINT *,"    i.e. every this many steps of the trajectory will be used."
			PRINT *," - 'ffc': expects a logical ('T' or 'F'). If 'T', then the average distances used for FFC are calculated"
			PRINT *,"    (i.e. dHH and rHH depending on the operation mode.)"
			PRINT *,"    Note that this average does not converge, which is not my fault."
			PRINT *," - 'calculate_exponential': expects a logical ('T' or 'F'). If 'T',"
			PRINT *,"    then the weighted average distances are also calculated."
			PRINT *,"    The weights are exp(-kr), where r is the distance and k is a constant defaulting to 1."
			PRINT *," - 'exponent': sets the exponent 'k' for the exponentially weighed average."
			PRINT *," - 'standard_deviation': expects a logical ('T' or 'F')."
			PRINT *,"    If 'T', then the standard deviation is - where reasonable - calculated."
			PRINT *
		END SUBROUTINE explain_input_files

		!" 3 - generate input files from user input"
		SUBROUTINE generate_all_input_files()
		IMPLICIT NONE
			PRINT *,"Generating input files from user input. Please answer the following questions."
			PRINT *,"First, the molecular input file will be generated. Do you want to give it a name? (y/n)"
			IF (user_input_logical()) THEN
				PRINT *,"Please type the name you want to give your molecular input file, such as 'mysystem.inp'"
				FILENAME_MOLECULAR_INPUT=TRIM(user_input_string(128))
			ENDIF
			PRINT *,"The file will be named '",TRIM(FILENAME_MOLECULAR_INPUT),"'"
			!number_of_molecules is '-1' now, will be initialised after the next line
			IF (check_if_inputfile_necessary(TRIM(FILENAME_MOLECULAR_INPUT))) CALL generate_molecular_input()
			PRINT *
			PRINT *,"Now, the general input file will be generated as '",TRIM(FILENAME_GENERAL_INPUT),"'"
			CALL generate_general_input()
		END SUBROUTINE generate_all_input_files

		!This subroutine generates the molecular input file from user input (in unit 8)
		SUBROUTINE generate_molecular_input()
		IMPLICIT NONE
		INTEGER :: n,deallocstatus,allocstatus,ios,totalcharge,number_of_constraints
		INTEGER :: total_number_of_drudes,number_of_assigned_drudes,local_drude,m,position_in_list
		INTEGER :: natoms!steps and natoms are required for the parallelisation memory estimate in case of suitable user input.
		INTEGER,ALLOCATABLE :: molecule_list(:,:),constraints_list(:,:),drude_list(:,:)
		LOGICAL :: turn_on_drudes,turn_on_constraints,manual_drude_assignment
		REAL :: drude_mass
			PRINT *,"How many timesteps are in your trajectory?"
			nsteps=user_input_integer(1,1000000000)
			PRINT *,"How many different molecule types do you have? (2 for a pure IL)"
			number_of_molecules=user_input_integer(1,10000)
			!allocate memory to store user input.
			ALLOCATE(molecule_list(number_of_molecules,3),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			totalcharge=0
			natoms=0
			DO n=1,number_of_molecules,1
				WRITE(*,'(A,I0,A)') " Please enter information for molecule number '",n,"':"
				PRINT *,"Which charge does this molecule have?"
				molecule_list(n,1)=user_input_integer(-99,99)
				PRINT *,"How many atoms are in one molecule (of this type)?"
				molecule_list(n,2)=user_input_integer(1,10000)
				PRINT *,"How many molecules (of this type) are in one timestep?"
				molecule_list(n,3)=user_input_integer(1,10000)
				natoms=natoms+molecule_list(n,3)*molecule_list(n,2)
				totalcharge=totalcharge+(molecule_list(n,1)*molecule_list(n,3))
			ENDDO
			!check for total charge neutrality
			IF (totalcharge/=0) CALL report_error(52,exit_status=totalcharge)
			WRITE(*,ADVANCE="NO",FMT='(A30)') " You would need approximately "
			CALL print_memory_requirement(DFLOAT(natoms)*DFLOAT(nsteps)*(12.0/1024.0d0))
			WRITE(*,FMT='(A38)') " of RAM to store the whole trajectory."
			PRINT *
			PRINT *,"Would you like to skip over the advanced settings (y/n)?"
			PRINT *,"These are currently: support for drude particles and constraints."
			IF (user_input_logical()) THEN
				PRINT *,"Not using advanced settings."
				turn_on_drudes=.FALSE.
				turn_on_constraints=.FALSE.
				manual_drude_assignment=.FALSE.
			ELSE
				PRINT *,"Would you like to turn on support for drude particles? (y/n)"
				turn_on_drudes=user_input_logical()
				IF (turn_on_drudes) THEN
					PRINT *,"The drude particles have to be represented by the capital letter 'X'."
					PRINT *,"Please enter the mass of the drude particles."
					drude_mass=user_input_real(0.0e0,1.0e0)
					PRINT *,"The temperature calculation of the drude particles requires drude-core assignment."
					PRINT *,"It is possible to use automatic assignment with a trajectory that contains cartesian coordinates."
					PRINT *,"The corresponding molecular input file section is then printed and can be used for the velocity trajectory."
					PRINT *,"Would you like to manually assign the drude particles instead?"
					PRINT *,"(Please note that the drudes have to be in the same molecule type as their cores)"
					manual_drude_assignment=user_input_logical()
					IF (manual_drude_assignment) THEN
						WRITE(*,'(A,I0,A)') " How many drude particles would you like to manually assign in total?"
						total_number_of_drudes=user_input_integer(1,natoms)
						number_of_assigned_drudes=0
						position_in_list=1
						ALLOCATE(drude_list(total_number_of_drudes,3),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
						DO n=1,number_of_molecules-1,1
							WRITE(*,'(A,I0,A)') " How many drude particles would you like to assign for molecule number '",n,"'?"
							WRITE(*,'(" (There are ",I0," drude particles left to distribute in total)")')&
							&(total_number_of_drudes-number_of_assigned_drudes)
							local_drude=user_input_integer(0,total_number_of_drudes-number_of_assigned_drudes)
							number_of_assigned_drudes=number_of_assigned_drudes+local_drude
							DO m=1,local_drude,1
								drude_list(position_in_list,1)=n !molecule type index
								WRITE(*,'(" Please enter the atom index of the drude particle (number ",I0," of ",I0,"):")')&
								&m,local_drude
								drude_list(position_in_list,3)=user_input_integer(1,molecule_list(n,2)) !atom index drude
								PRINT *,"Please enter the atom index of the core to which this drude particle is attached to:"
								drude_list(position_in_list,2)=user_input_integer(1,molecule_list(n,2)) !atom index core
								position_in_list=position_in_list+1
							ENDDO
							IF (total_number_of_drudes==number_of_assigned_drudes) EXIT
						ENDDO
						IF (total_number_of_drudes/=number_of_assigned_drudes) THEN
							!still some drudes left to distribute!
							PRINT *,"The remaining drudes have to be assigned for the last molecule type."
							local_drude=(total_number_of_drudes-number_of_assigned_drudes)
							WRITE(*,'(" (There are ",I0," drude particles left to distribute)")') local_drude
							!molecule type index of last molecule...
							n=number_of_molecules
							DO m=1,local_drude,1
								drude_list(position_in_list,1)=n !molecule type index
								WRITE(*,'(" Please enter the atom index of the drude particle (number ",I0," of ",I0,"):")')&
								&m,local_drude
								drude_list(position_in_list,3)=user_input_integer(1,molecule_list(n,2)) !atom index drude
								PRINT *,"Please enter the atom index of the core to which this drude particle is attached to:"
								drude_list(position_in_list,2)=user_input_integer(1,molecule_list(n,2)) !atom index core
								position_in_list=position_in_list+1
							ENDDO
							number_of_assigned_drudes=number_of_assigned_drudes+local_drude
						ENDIF
						IF (total_number_of_drudes==number_of_assigned_drudes) PRINT *,"All drudes have been assigned."
					ENDIF
				ENDIF
				PRINT *,"Would you like to add constraints to some molecules? (y/n)"
				PRINT *,"This influences the keyword 'temperature' by decreasing the degrees of freedom"
				turn_on_constraints=user_input_logical()
				IF (turn_on_constraints) THEN
					PRINT *,"The constraints have to specified per molecule."
					PRINT *,"Thus, if you specify 10 constraints for a molecule of which there are 64 in a box,"
					PRINT *,"a total of 10*64=640 constraints will be subtracted from the degrees of freedom."
					PRINT *,"Please enter the number of constraints you would like to specify."
					number_of_constraints=user_input_integer(0,number_of_molecules)
					IF (number_of_constraints==0) THEN
						PRINT *,"no constraints specified - turning constraints off again."
						turn_on_constraints=.FALSE.
					ELSE
						ALLOCATE(constraints_list(number_of_constraints,2),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
						DO n=1,number_of_constraints,1
							WRITE(*,'("Constraint ",I0," out of ",I0,":")') n,number_of_constraints
							PRINT *,"Please enter the molecule type index of this constraint."
							constraints_list(n,1)=user_input_integer(1,number_of_molecules)
							PRINT *,"How many constraints are there in *one* molecule of this type?"
							!In the limiting case of an entirely rigid molecule, this will be all of the 3*N-6 internal coordinates.
							constraints_list(n,2)=user_input_integer(1,molecule_list(constraints_list(n,1),2)*3-6)
						ENDDO
					ENDIF
				ENDIF
			ENDIF
			!Now, write the molecular input file.
			WRITE(*,FMT='(A32)',ADVANCE="NO") " writing molecular input file..."
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(FILENAME_MOLECULAR_INPUT),IOSTAT=ios)!no input path is added for the molecular file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,'(" ",I0," ### number of timesteps in trajectory")') nsteps
			WRITE(8,'(" ",I0," ### number of different types of molecules. Followed by list of molecules.")') number_of_molecules
			DO n=1,number_of_molecules,1
				WRITE(8,'(SP,I3,SS," ",I0," ",I0," ### There are ",I0," molecules per step with charge ",SP,I2,SS," and ",I0," atoms each.")')&
				& molecule_list(n,:),molecule_list(n,3),molecule_list(n,1:2)
			ENDDO
			IF (turn_on_drudes) THEN
				WRITE(8,'(" default_masses 1 ### The following line specifies a custom mass.")')
				WRITE(8,'(" X",F7.3," ### By that, the support for drude particles is turned on.")') drude_mass
			ENDIF
			IF (manual_drude_assignment) THEN
				WRITE(8,'(" drudes ",I0," ### manual assignment of drude particles to their respective cores")') total_number_of_drudes
				DO position_in_list=1,total_number_of_drudes,1
					WRITE(8,'(" ",I0," ",I0," ",I0)') drude_list(position_in_list,:)
				ENDDO
			ENDIF
			IF (turn_on_constraints) THEN
				WRITE(8,'(" constraints ",I0)') number_of_constraints
				DO n=1,number_of_constraints,1
					WRITE(8,'(" ",I0," ",I0," ### Putting ",I0," constraints on every molecule of type ",I0,".")')&
					&constraints_list(n,:),constraints_list(n,2),constraints_list(n,1)
				ENDDO
			ENDIF
			CLOSE(UNIT=8)
			WRITE(*,*) "done"
			DEALLOCATE(molecule_list,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			IF (turn_on_constraints) THEN
				DEALLOCATE(constraints_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ENDIF
			IF (manual_drude_assignment) THEN
				DEALLOCATE(drude_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ENDIF
		END SUBROUTINE generate_molecular_input

		SUBROUTINE trajectory_filename_input()
		IMPLICIT NONE
		LOGICAL :: valid_trajectory_type_extension
			!first, ask for trajectory name
			PRINT *,"Please type the name of your trajectory file, such as 'trajectory.lmp'"
			FILENAME_TRAJECTORY=TRIM(user_input_string(128))
			!Try to get type from extension.
			valid_trajectory_type_extension=.TRUE.
			BOX_VOLUME_GIVEN=.FALSE.
			SELECT CASE (FILENAME_TRAJECTORY(LEN(TRIM(FILENAME_TRAJECTORY))-3:LEN(TRIM(FILENAME_TRAJECTORY))))
			CASE (".lmp",".LMP")
				PRINT *,"assuming lammps trajectory based on file extension."
				TRAJECTORY_TYPE="lmp"
				BOX_VOLUME_GIVEN=.TRUE.
			CASE (".xyz",".XYZ")
				PRINT *,"assuming trajectory in xyz format based on file extension."
				TRAJECTORY_TYPE="xyz"
				BOX_VOLUME_GIVEN=.FALSE.
			CASE DEFAULT
				valid_trajectory_type_extension=.FALSE.
			END SELECT
			! ask the user for the trajectory type if that didn't work.
			IF (.NOT.(valid_trajectory_type_extension)) THEN
				PRINT *,"Please enter which type of trajectory you have. Currently supported are 'lmp' and 'xyz'."
				SELECT CASE (TRIM(user_input_string(4)))!'4' to suppress the error.
				CASE ("lmp")
					TRAJECTORY_TYPE="lmp"
					BOX_VOLUME_GIVEN=.TRUE.
				CASE ("xyz")
					TRAJECTORY_TYPE="xyz"
					BOX_VOLUME_GIVEN=.FALSE.
				CASE DEFAULT
					CALL report_error(51)
				END SELECT
			ENDIF
		END SUBROUTINE trajectory_filename_input

		!This subroutine generates the general input file from user input (in unit 8)
		SUBROUTINE generate_general_input()
		IMPLICIT NONE
		INTEGER :: ios
		REAL :: lower,upper
		LOGICAL :: manual_box_edge
			!get the trajectory filename
			CALL trajectory_filename_input()
			manual_box_edge=.FALSE.
			!take care of the custom box volume
			IF (.NOT.(BOX_VOLUME_GIVEN)) THEN
				PRINT *,"Box boundaries are not available from this type of trajectory file."
				PRINT *,"Would you like to specify them manually? (y/n)"
				PRINT *,"(only cubic boxes supported here)"
				manual_box_edge=user_input_logical()
				IF (manual_box_edge) THEN
					PRINT *,"Please enter lower boundary."
					lower=user_input_real(-2000.0,+1000.0)
					PRINT *,"Please enter upper boundary."
					upper=user_input_real(lower,+2000.0)
				ENDIF
			ENDIF
			!Get the paths
			PRINT *,"Do you want to specify the paths were the files are located? (y/n)"
			IF (user_input_logical()) THEN
				PRINT *,"Please enter the following Paths:"
				PRINT *," - Path to the trajectory:"
				PATH_TRAJECTORY=TRIM(user_input_string(128))
				PRINT *," - Path to the input files (other than '",TRIM(FILENAME_GENERAL_INPUT),"' and '",TRIM(FILENAME_MOLECULAR_INPUT),"')."
				PATH_INPUT=TRIM(user_input_string(128))
				PRINT *," - Output path:"
				PATH_OUTPUT=TRIM(user_input_string(128))
				PRINT *,"you are responsible for these to be the correct paths of existing folders."
			ELSE
				PRINT *,"Standard location for all files is the executable's folder (",'"./"',")."
				PATH_TRAJECTORY='./'
				PATH_INPUT='./'
				PATH_OUTPUT='./'	
			ENDIF
			PRINT *
			PRINT *,"Do you want to define a time scaling factor? (y/n)"
			IF (user_input_logical()) THEN
				PRINT *,"Please specify the time scaling factor (as an integer)."
				PRINT *,"(e.g. when dumping every 100fs, type '100' to get the right time unit)"
				TIME_SCALING_FACTOR=user_input_integer(1,1000000000)
			ELSE
				TIME_SCALING_FACTOR=-1
			ENDIF
			!Now, write the general input file.
			PRINT *,"The file '",TRIM(FILENAME_GENERAL_INPUT),"' will be opened now. Please do not change it during write."
			WRITE(*,FMT='(A40)',ADVANCE="YES") " Starting to write general input file..."
			PRINT *
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(FILENAME_GENERAL_INPUT),IOSTAT=ios)!no input path is added for the general input file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			WRITE(8,*) '"',TRIM(FILENAME_TRAJECTORY),'"'," ### trajectory filename"
			WRITE(8,*) '"',TRIM(FILENAME_MOLECULAR_INPUT),'"'," ### inputfile for module MOLECULAR"
			WRITE(8,*) '"',TRIM(PATH_TRAJECTORY),'"'," ### path to trajectory"
			WRITE(8,*) '"',TRIM(PATH_INPUT),'"'," ### path to other input files"
			WRITE(8,*) '"',TRIM(PATH_OUTPUT),'"'," ### output folder"
			IF (TIME_SCALING_FACTOR>0) WRITE(8,'(A,I0,A)') " time_scaling ",TIME_SCALING_FACTOR,&
			&" ### factor to scale the timestep with to arrive at a useful time unit."
			IF (manual_box_edge) WRITE(8,'(A,2E16.8,A)') " cubic_box_edge ",lower,upper," ### lower and upper cubic box boundaries."
			CLOSE(UNIT=8)
			PRINT *,"Does your trajectory contain velocities (y) or cartesian coordinates (n)?"
			IF (user_input_logical()) THEN
				!Velocities in trajectory - read input for VACFs.
				wrapping_is_sensible=.FALSE.!well, duh.
				CALL velocity_user_input()
			ELSE
				!Cartesian coordinates in trajectory - read other input.
				CALL coordinates_user_input()
			ENDIF
			anytask=.TRUE.
			!append the trajectory type.
			CALL append_string("trajectory_type "//TRAJECTORY_TYPE//" ### which format the (input) trajectory has")
			!finally, append "quit" statement
			CALL append_string("quit")
			!...and some comments.
			CALL append_string("")
			CALL append_string("This is the general input file.")
			CALL append_string("It controls the behaviour of the trajectory analyser.")
			!change SKIP_ANALYSIS to .FALSE. if possible
			PRINT *,"...done writing general input file."
			PRINT *
		END SUBROUTINE generate_general_input

		!This routine writes into unit 8 from user input - assuming that velocities are in the trajectory.
		SUBROUTINE velocity_user_input()
		USE AUTOCORRELATION
		IMPLICIT NONE
		LOGICAL :: parallelisation_possible,parallelisation_requested,own_prefix
		INTEGER :: nthreads,analysis_number,n,molecule_type_index,maxmol,startstep,endstep
		CHARACTER(LEN=1024) :: fstring
			parallelisation_possible=.FALSE.!
			parallelisation_requested=.FALSE.!only changes to .TRUE. from here, never back.
			analysis_number=1
			nthreads=0
			own_prefix=.FALSE.
			DO n=1,MAXITERATIONS,1
				IF (.NOT.(own_prefix)) WRITE(OUTPUT_PREFIX,'("out_",I0,"_")') analysis_number
				IF (n==1) THEN
					PRINT *,"Please choose an action you want the program to take later:"
				ELSE
					PRINT *
					PRINT *,"Please choose another action you want the program to take later:"
				ENDIF
				PRINT *," 0  - No more actions needed."
				PRINT *," 1  - Velocity (auto)correlation functions / components thereof."
				PRINT *," 2  - Change prefix. (Currently '",TRIM(OUTPUT_PREFIX),"')"
				PRINT *," 3  - Print settings."
				IF (TIME_OUTPUT) THEN
					PRINT *," 4  - turn off time output."
				ELSE
					PRINT *," 4  - turn on time output."
				ENDIF
				IF (VERBOSE_OUTPUT) THEN
					PRINT *," 5  - turn off verbose output."
				ELSE
					PRINT *," 5  - turn on verbose output."
				ENDIF
				IF (ERROR_OUTPUT) THEN
					PRINT *," 6  - turn off error output."
				ELSE
					PRINT *," 6  - turn on error output."
				ENDIF
				PRINT *," 7  - Specify the number of threads to use."
				PRINT *," 8  - Reduce the trajectory to centres of mass."
				PRINT *," 9  - Compute instantaneous temperature."
				PRINT *," 10 - Print information about drude particles."
				PRINT *," 11 - Compute temperature for drude particles."
				PRINT *," 12 - Write trajectory with drude particles merged into cores."
				PRINT *," 13 - Print atomic masses (in format suitable for a molecular input file)"
				PRINT *," 14 - Electrical conductivity via CACF / components thereof."
				PRINT *," 15 - Print atomic charges (in format suitable for a molecular input file)"
				PRINT *," 16 - Write trajectory only with drude particles (minus velocity of their cores)"
				SELECT CASE (user_input_integer(0,16))
				CASE (0)!done here.
					EXIT
				CASE (1)!compute VACFs...
					PRINT *,"There are two possible ways to approach this:"
					PRINT *,"  a) manually specify the components to be calculated."
					PRINT *,"  b) calculate RMM-VCFs for two components - includes conductivity and Haven Ratio."
					PRINT *,"Would you like to use the manual specification (y), or rather the RMM-VCFs (n)?"
					smalltask=.FALSE.
					IF (user_input_logical()) THEN
						CALL user_vcf_components_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,&
						&filename_vc_components,"vcf")
						IF (parallelisation_requested) THEN
							CALL append_string("parallel_operation T ### turn on parallel operation")
							WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
							CALL append_string(fstring)
						ENDIF
						CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
						CALL append_string('correlation "'//TRIM(OUTPUT_PREFIX)//TRIM(filename_vc_components)//&
						&'" ### compute custom velocity cross correlations')
						IF (own_prefix) THEN
							own_prefix=.FALSE.
						ELSE
							analysis_number=analysis_number+1
						ENDIF
						!enough information for the analysis.
						SKIP_ANALYSIS=.FALSE.
					ELSE
						IF (number_of_molecules==1) THEN
							PRINT *,"This module needs at least two types of molecules to be present."
							PRINT *,"Only one is specified. If you want to force it:"
							PRINT *,"Specify the molecule twice with half the number."
						ELSE!number_of_molecules is unknown (-1) or large enough.
							CALL user_vacf_input&
							&(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_rmmvcf)
							IF (parallelisation_requested) THEN
								CALL append_string("parallel_operation T ### turn on parallel operation")
								WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
								CALL append_string(fstring)
							ENDIF
							CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
							CALL append_string('correlation "'//TRIM(OUTPUT_PREFIX)//TRIM(filename_rmmvcf)//'" ### compute velocity cross correlations')
							IF (own_prefix) THEN
								own_prefix=.FALSE.
							ELSE
								analysis_number=analysis_number+1
							ENDIF
							!enough information for the analysis.
							SKIP_ANALYSIS=.FALSE.
						ENDIF
					ENDIF
				CASE (2)!set own prefix
					own_prefix=.TRUE.
					PRINT *,"Please enter the leading string for output files:"
					OUTPUT_PREFIX=TRIM(user_input_string(64))
				CASE (3)!show the settings at this point
					CALL append_string("show_settings ### show the values of important variables/parameters at this point")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (4)!switch TIME_OUTPUT
					TIME_OUTPUT=(.NOT.(TIME_OUTPUT))
					IF (TIME_OUTPUT) THEN
						CALL append_string("time_output T ### give time output")
					ELSE
						CALL append_string("time_output F ### switch off time output")
					ENDIF
				CASE (5)!switch VERBOSE_OUTPUT
					VERBOSE_OUTPUT=(.NOT.(VERBOSE_OUTPUT))
					IF (VERBOSE_OUTPUT) THEN
						CALL append_string("verbose_output T ### verbose output")
					ELSE
						CALL append_string("verbose_output F ### terse output")
					ENDIF
				CASE (6)!switch ERROR_OUTPUT
					ERROR_OUTPUT=(.NOT.(ERROR_OUTPUT))
					IF (ERROR_OUTPUT) THEN
						CALL append_string("error_output T ### print error messages")
					ELSE
						CALL append_string("error_output F ### turn off error messages")
					ENDIF
				CASE (7)!change number of threads
					PRINT *,"Please give the number of threads you want to use."
					PRINT *,"If you type '0' the program will (later) try to use the (permitted) maximum."
					nthreads=user_input_integer(0,64)
				CASE (8)!convert to centre of mass
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("convert_simple ### centre of mass trajectory with default values.")
						CYCLE
					ENDIF
					smalltask=.FALSE.
					PRINT *,"Would you also like to write the appropriately adjusted molecular input file? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("convert T ### reduce trajectory to centre of mass, write new molecular.inp")
					ELSE
						CALL append_string("convert F ### reduce trajectory to centre of mass, don't write new molecular.inp")
					ENDIF
				CASE (9)!compute instantaneous temperature.
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("temperature_simple ### calculate temperature with default values.")
						CYCLE
					ENDIF
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					IF (number_of_molecules==1) THEN
						PRINT *,"Only one molecule type available, which will be observed."
						molecule_type_index=1
					ELSE
						PRINT *,"Would you compute the temperature of just one molecule type (y), or for every type (n)?"
						IF (user_input_logical()) THEN
							PRINT *,"Please enter the index of the molecule type you wish to observe."
							molecule_type_index=user_input_integer(1,maxmol)
						ELSE
							molecule_type_index=-1
						ENDIF
					ENDIF
					PRINT *,"It is necessary to provide a range of timesteps which are to be analysed."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>10) smalltask=.FALSE.
					WRITE(fstring,'("temperature ",I0," ",I0," ",I0)') molecule_type_index,startstep,endstep
					IF (molecule_type_index==-1) THEN
						WRITE(fstring,'(A," ### calculate temperature of every molecule type for timesteps ",I0,"-",I0)')&
						&TRIM(fstring),startstep,endstep
					ELSE
						WRITE(fstring,'(A," ### calculate temperature of molecule type ",I0," for timesteps ",I0,"-",I0)')&
						&TRIM(fstring),molecule_type_index,startstep,endstep
					ENDIF
					CALL append_string(fstring)
				CASE (10)!show the drude settings at this point
					CALL append_string("show_drude ### print detailed information about drude assignments")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (11) !drude temperature
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("drude_temp_simple ### calculate drude temperature with default values.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide a range of timesteps which are to be analysed."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>50) smalltask=.FALSE.
					WRITE(fstring,'("drude_temp ",I0," ",I0)') startstep,endstep
					WRITE(fstring,'(A," ### calculate drude temperature for timesteps ",I0,"-",I0)')&
					&TRIM(fstring),startstep,endstep
					CALL append_string(fstring)
				CASE (12) !remove drudes
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("remove_drudes_simple ### merge drude particles into cores with default values.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide a range of timesteps which are to be analysed."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>50) smalltask=.FALSE.
					WRITE(fstring,'("remove_drudes ",I0," ",I0)') startstep,endstep
					WRITE(fstring,'(A," ### write trajectory without drude particles for timesteps ",I0,"-",I0)')&
					&TRIM(fstring),startstep,endstep
					CALL append_string(fstring)
				CASE (13)!show the atomic masses at this point
					CALL append_string("print_atomic_masses ### print atomic masses")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (14)
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There are two possible ways to approach this:"
					PRINT *,"  a) manually specify the CACF components to be calculated."
					PRINT *,"  b) Only calculate the overall conductivity (very fast)."
					PRINT *,"Would you like to use the manual specification (y), or rather the overall conductivity (n)?"
					PRINT *,"(If you need the distinct contributions / the Haven ratio,"
					PRINT *," and have a 2 component sytem, consider using RMM-VCFs.)"
					smalltask=.FALSE.
					IF (user_input_logical()) THEN
						CALL user_vcf_components_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,&
						&filename_cacf_components,"cacf")
						IF (parallelisation_requested) THEN
							CALL append_string("parallel_operation T ### turn on parallel operation")
							WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
							CALL append_string(fstring)
						ENDIF
						CALL append_string('correlation "'//TRIM(OUTPUT_PREFIX)//TRIM(filename_cacf_components)//&
						&'" ### compute custom velocity cross correlations')
						!enough information for the analysis.
						SKIP_ANALYSIS=.FALSE.
					ELSE
						PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
						PRINT *,"Would you like to take this shortcut? (y/n)"
						IF (user_input_logical()) THEN
							CALL append_string("conductivity_simple ### calculate electrical conductivity with default values.")
							CYCLE
							
						ENDIF
						CALL user_conductivity_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,&
						&filename_conductivity)
						IF (parallelisation_requested) THEN
							CALL append_string("parallel_operation T ### turn on parallel operation")
							WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
							CALL append_string(fstring)
						ENDIF
						CALL append_string('correlation "'//TRIM(OUTPUT_PREFIX)//TRIM(filename_conductivity)//&
						&'" ### compute custom velocity cross correlations')
						!enough information for the analysis.
						SKIP_ANALYSIS=.FALSE.
					ENDIF
				CASE (15)!show the atomic masses at this point
					CALL append_string("print_atomic_charges ### print atomic charges")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (16) !remove cores
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("remove_cores_simple ### subtract cores, write drude particles.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide a range of timesteps which are to be analysed."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>50) smalltask=.FALSE.
					WRITE(fstring,'("remove_cores ",I0," ",I0)') startstep,endstep
					WRITE(fstring,'(A," ### subtract cores, write drude particles for timesteps ",I0,"-",I0)')&
					&TRIM(fstring),startstep,endstep
					CALL append_string(fstring)
				CASE DEFAULT
					CALL report_error(0)
				END SELECT
			ENDDO
			CALL toggle_sequential_read(parallelisation_possible,parallelisation_requested)
		END SUBROUTINE velocity_user_input

		!This routine writes into unit 8 from user input - assuming that cartesian coordinates are in the trajectory.
		SUBROUTINE coordinates_user_input()
		USE DIFFUSION
		USE AUTOCORRELATION
		USE DISTRIBUTION
		USE DISTANCE
		IMPLICIT NONE
		LOGICAL :: parallelisation_possible,parallelisation_requested,own_prefix
		INTEGER :: nthreads,analysis_number,n,snap,startstep,endstep,molecule_type_index,molecule_index
		INTEGER :: maxmol,deltasteps,molecule_type_index_2,inputnumber
		REAL :: cutoff
		CHARACTER(LEN=1024) :: fstring
			parallelisation_possible=.FALSE.!
			parallelisation_requested=.FALSE.!only changes to .TRUE. from here, never back.
			analysis_number=1
			nthreads=0
			own_prefix=.FALSE.
			DO n=1,MAXITERATIONS,1
				IF (.NOT.(own_prefix)) WRITE(OUTPUT_PREFIX,'("out_",I0,"_")') analysis_number
				IF (n==1) THEN
					PRINT *,"Please choose an action you want the program to take later:"
				ELSE
					PRINT *
					PRINT *,"Please choose another action you want the program to take later:"
				ENDIF
				PRINT *," 0 - No more actions needed."
				PRINT *," 1 - Dihedral condition analysis"
				PRINT *," 2 - Calculate mean squared displacements (self or cross)"
				PRINT *," 3 - Change prefix. (Currently '",TRIM(OUTPUT_PREFIX),"')"
				PRINT *," 4 - Print settings."
				IF (TIME_OUTPUT) THEN
					PRINT *," 5 - turn off time output."
				ELSE
					PRINT *," 5 - turn on time output."
				ENDIF
				IF (VERBOSE_OUTPUT) THEN
					PRINT *," 6 - turn off verbose output."
				ELSE
					PRINT *," 6 - turn on verbose output."
				ENDIF
				IF (ERROR_OUTPUT) THEN
					PRINT *," 7 - turn off error output."
				ELSE
					PRINT *," 7 - turn on error output."
				ENDIF
				PRINT *," 8 - Specify the number of threads to use."
				PRINT *," 9 - Write an xyz file for every molecule type."
				PRINT *," 10 - Write a snapshot of the whole box."
				PRINT *," 11 - Split the trajectory according to molecule type"
				PRINT *," 12 - Convert the trajectory to centre of mass."
				PRINT *," 13 - Write a trajectory subset for just one specific molecule."
				PRINT *," 14 - Write a trajectory for a specific molecule and its neighbours"
				PRINT *," 15 - vector reorientation dynamics (time correlation function)"
				PRINT *," 16 - Print information about drude particles."
				PRINT *," 17 - Compute ensemble average of radius of gyration."
				PRINT *," 18 - Write trajectory with drude particles merged into cores."
				PRINT *," 19 - Calculate close contact distances (inter- and intramolecular). See also #28"
				PRINT *," 20 - Calculate polar or cylindrical distribution function."
				PRINT *," 21 - Jump analysis / jump velocity distribution"
				PRINT *," 22 - Print atomic masses (in format suitable for a molecular input file)"
				PRINT *," 23 - Write trajectory with neighbouring molecules around one selected molecule"
				PRINT *," 24 - Dump close contact dimers for one step"
				PRINT *," 25 - Print atomic charges (in format suitable for a molecular input file)"
				PRINT *," 26 - Calculate dipole moment statistics from first timestep."
				PRINT *," 27 - Write trajectory only with drude particles (minus position of their cores)"
				PRINT *," 28 - Calculate average distances (closest or weighed, intra- or intermolecular)."
				SELECT CASE (user_input_integer(0,28))
				CASE (0)!done here.
					EXIT
				CASE (1)!dihedral condition analysis
					smalltask=.FALSE.
					CALL user_dihedral_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_dihedral)
					IF (parallelisation_requested) THEN
						CALL append_string("parallel_operation T ### turn on parallel operation")
						WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
						CALL append_string(fstring)
					ENDIF
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					CALL append_string('dihedral "'//TRIM(OUTPUT_PREFIX)//TRIM(filename_dihedral)//'" ### invoke dihedral condition analysis')
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					!enough information for the analysis.
					SKIP_ANALYSIS=.FALSE.
				CASE (2)!mean squared displacement section.
					wrapping_is_sensible=.FALSE.
					smalltask=.FALSE.
					IF (nsteps<11) THEN
						PRINT *,"Your trajectory is really too short for that. Please use more timesteps."
						PRINT *,"The minimum of steps - even if only for debugging purposes - should be 11."
					ELSE
						CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
						IF (own_prefix) THEN
							own_prefix=.FALSE.
						ELSE
							analysis_number=analysis_number+1
						ENDIF
						PRINT *,"There is a default mode (self-diffusion) available for this analysis that doesn't require additional input."
						PRINT *,"Would you like to take this shortcut? (y/n)"
						IF (user_input_logical()) THEN
							CALL append_string("diffusion_simple ### mean squared displacement with default values.")
							CYCLE
						ENDIF
						PRINT *,"Would you like to calculate self components (y) or cross components (n)?"
						PRINT *,"('cross' components are calculated using the Einstein Relation equivalent of the RMM-VCFs)"
						IF (user_input_logical()) THEN
							CALL user_msd_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_msd)
						ELSE
							CALL user_cross_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_msd)
						ENDIF
						IF (parallelisation_requested) THEN
							CALL append_string("parallel_operation T ### turn on parallel operation")
							WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
							CALL append_string(fstring)
						ENDIF
						CALL append_string('diffusion "'//TRIM(OUTPUT_PREFIX)//&
						&TRIM(filename_msd)//'" ### compute (drift-corrected) mean squared displacements')
						!enough information for the analysis.
						SKIP_ANALYSIS=.FALSE.
					ENDIF
				CASE (3)!set own prefix
					own_prefix=.TRUE.
					PRINT *,"Please enter the leading string for output files:"
					OUTPUT_PREFIX=TRIM(user_input_string(64))
				CASE (4)!show the settings at this point
					CALL append_string("show_settings ### show the values of important variables/parameters at this point")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (5)!switch TIME_OUTPUT
					TIME_OUTPUT=(.NOT.(TIME_OUTPUT))
					IF (TIME_OUTPUT) THEN
						CALL append_string("time_output T ### give time output")
					ELSE
						CALL append_string("time_output F ### switch off time output")
					ENDIF
				CASE (6)!switch VERBOSE_OUTPUT
					VERBOSE_OUTPUT=(.NOT.(VERBOSE_OUTPUT))
					IF (VERBOSE_OUTPUT) THEN
						CALL append_string("verbose_output T ### verbose output")
					ELSE
						CALL append_string("verbose_output F ### terse output")
					ENDIF
				CASE (7)!switch ERROR_OUTPUT
					ERROR_OUTPUT=(.NOT.(ERROR_OUTPUT))
					IF (ERROR_OUTPUT) THEN
						CALL append_string("error_output T ### print error messages")
					ELSE
						CALL append_string("error_output F ### turn off error messages")
					ENDIF
				CASE (8)!specify number of threads
					PRINT *,"Please give the number of threads you want to use."
					PRINT *,"If you type '0' the program will (later) try to use the (permitted) maximum."
					nthreads=user_input_integer(0,64)
				CASE (9)!dump example files.
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					CALL append_string("dump_example ### write xyz files for all the molecule types")
					PRINT *,"The corresponding section has been added to the input file."
					smalltask=.TRUE.
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
				CASE (10)!dump a snapshot of a given timestep
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("dump_snapshot_simple ### dump snapshot with default values.")
						CYCLE
					ENDIF
					PRINT *,"Please give the timestep of which you want the snapshot to be:"
					snap=user_input_integer(1,nsteps)
					PRINT *,"Do you want the molecules to be exported as separate files? (y/n)"
					WRITE(fstring,'("dump_snapshot ",I0," ",L1," ### write snapshot of timestep ",I0," into xyz file(s)")')&
					&snap,user_input_logical(),snap
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
					smalltask=.TRUE.
				CASE (11)!splits the trajectory
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("dump_split_simple ### split trajectory per molecule type.")
						CYCLE
					ENDIF
					smalltask=.FALSE.
					PRINT *,"Please give the first timestep in the trajectory:"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Please type the final timestep you wish to export:"
					endstep=user_input_integer(startstep,nsteps)
					WRITE(fstring,'("dump_split ",I0," ",I0," ### split trajectory (per molecule type)")') startstep,endstep
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
				CASE (12)!converts to centre of mass
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("convert_simple ### centre of mass trajectory with default values.")
						CYCLE
					ENDIF
					smalltask=.FALSE.
					PRINT *,"Would you also like to produce an adjusted molecular input file? (y/n)"
					IF (user_input_logical()) THEN
						WRITE(fstring,'("convert T ### produce centre of mass trajectory and molecular input file")')
					ELSE
						WRITE(fstring,'("convert F ### produce centre of mass trajectory, but no molecular input file")')
					ENDIF
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
				CASE (13)!writes a single molecule trajectory.
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					!Write the corresponding section to 'fstring'.
					PRINT *,"Please give the first timestep in the trajectory:"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Please type the final timestep you wish to export:"
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>10) smalltask=.FALSE.
					PRINT *,"Now, you have to specify the molecule type and the number of the molecule."
					PRINT *,"Please enter the molecule type index:"
					molecule_type_index=user_input_integer(1,maxmol)
					PRINT *,"Please enter the molecule index: ('Which one is it?')"
					molecule_index=user_input_integer(1,10000)
					PRINT *,"Should the output trajectory be referenced to its centre of mass? (y/n)"
					IF (user_input_logical()) THEN
						WRITE(fstring,'(A,I0," ",I0," ",I0," ",I0," ",A,I0,A,I0,A,I0,A,I0)') &
						&"dump_single T ",startstep,endstep,molecule_type_index,molecule_index,&
						&" ### produce centre-of-mass trajectory for steps ",startstep,"-",endstep,&
						&" and molecule ",molecule_index," of type ",molecule_type_index
					ELSE
						WRITE(fstring,'(A,I0," ",I0," ",I0," ",I0," ",A,I0,A,I0,A,I0,A,I0)') &
						&"dump_single F ",startstep,endstep,molecule_type_index,molecule_index,&
						&" ### dump subset of original trajectory for steps ",startstep,"-",endstep,&
						&" and molecule ",molecule_index," of type ",molecule_type_index
					ENDIF
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					CALL append_string(fstring)
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"The corresponding section has been added to the input file."
				CASE (14)!writes a trajectory for a specified molecule and its neighbours.
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					!Write the corresponding section to 'fstring'.
					PRINT *,"Please give the first timestep in the trajectory:"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Please type the final timestep you wish to export:"
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>10) smalltask=.FALSE.
					PRINT *,"Now, you have to specify the molecule type and the number of the molecule."
					PRINT *,"Please enter the molecule type index:"
					molecule_type_index=user_input_integer(1,maxmol)
					PRINT *,"Please enter the molecule index: ('Which one is it?')"
					molecule_index=user_input_integer(1,10000)
					PRINT *,"Please enter the cutoff you want to use to recognise neighbours:"
					PRINT *,"(This number should be smaller than half the box size)"
					cutoff=user_input_real(0.1e0,100.0e0)
					PRINT *,"Should the output trajectory be referenced to its centre of mass? (y/n)"
					IF (user_input_logical()) THEN
						WRITE(fstring,'(A,I0," ",I0," ",I0," ",I0," ",F0.2,A,I0,A,I0,A,I0,A,I0)') &
						&"dump_cut T ",startstep,endstep,molecule_type_index,molecule_index,cutoff,&
						&" ### produce centred trajectory for neighbours for steps ",startstep,"-",endstep,&
						&" and molecule ",molecule_index," of type ",molecule_type_index
					ELSE
						WRITE(fstring,'(A,I0," ",I0," ",I0," ",I0," ",F0.2,A,I0,A,I0,A,I0,A,I0)') &
						&"dump_cut F ",startstep,endstep,molecule_type_index,molecule_index,cutoff,&
						&" ### dump neighbours from original trajectory for steps ",startstep,"-",endstep,&
						&" and molecule ",molecule_index," of type ",molecule_type_index
					ENDIF
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					CALL append_string(fstring)
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"The corresponding section has been added to the input file."
				CASE (15)!reorientational time correlation function
					smalltask=.FALSE.
					CALL user_reorientation_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_reorient)
					IF (parallelisation_requested) THEN
						CALL append_string("parallel_operation T ### turn on parallel operation")
						WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
						CALL append_string(fstring)
					ENDIF
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					CALL append_string('correlation "'//TRIM(OUTPUT_PREFIX)//TRIM(filename_reorient)//'" ### vector reorientation tcf')
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					!enough information for the analysis.
					SKIP_ANALYSIS=.FALSE.
				CASE (16)!show the drude settings at this point
					CALL append_string("show_drude ### print detailed information about drude assignments")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (17)!radius of gyration
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("gyradius_simple ### radius of gyration / maxdist analysis with default values.")
						CYCLE
					ENDIF
					PRINT *,"Please give the first timestep in the trajectory:"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Please type the final timestep you wish to export:"
					endstep=user_input_integer(startstep,nsteps)
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					IF (number_of_molecules==1) THEN
						PRINT *,"Only one molecule type available, which will be observed."
						molecule_type_index=1
					ELSE
						PRINT *,"Would you compute the radius of gyration of just one molecule type (y), or for every type (n)?"
						IF (user_input_logical()) THEN
							PRINT *,"Please enter the index of the molecule type you wish to observe."
							molecule_type_index=user_input_integer(1,maxmol)
						ELSE
							molecule_type_index=-1
						ENDIF
					ENDIF
					WRITE(fstring,'(A,I0," ",I0," ",I0,A)')&
					&"gyradius ",molecule_type_index,startstep,endstep,&
					&" ### compute radius of gyration and maximum distance of any atom from centre of mass."
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
				CASE (18)!merge drude particles into cores
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("remove_drudes_simple ### merge drude particles into cores with default values.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide a range of timesteps which are to be analysed."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>50) smalltask=.FALSE.
					WRITE(fstring,'("remove_drudes ",I0," ",I0)') startstep,endstep
					WRITE(fstring,'(A," ### write trajectory without drude particles for timesteps ",I0,"-",I0)')&
					&TRIM(fstring),startstep,endstep
					CALL append_string(fstring)
				CASE (19)!close contact distance, inter- and intramolecular
					parallelisation_possible=.TRUE.
					IF (.NOT.(parallelisation_requested)) THEN! ask for parallelisation, if not yet requested.
						PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
						IF (user_input_logical()) parallelisation_requested=.TRUE.
					ENDIF
					IF (parallelisation_requested) THEN
						CALL append_string("parallel_operation T ### turn on parallel operation")
						WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
						CALL append_string(fstring)
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("contact_distance_simple ### inter/intramolecular distance analysis with default values.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide the timestep to be analysed."
					PRINT *,"Please enter this timestep as integer:"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Would you like to compute the contact distances for just one molecule type (y), or for every type (n)?"
					IF (user_input_logical()) THEN
						PRINT *,"Please enter the index of the molecule type you wish to observe."
						molecule_type_index=user_input_integer(1,maxmol)
					ELSE
						molecule_type_index=-1
					ENDIF
					WRITE(fstring,'(" contact_distance ",I0," ",I0)') startstep,molecule_type_index
					WRITE(fstring,'(A," ### compute largest/smallest intramolecular and smallest intermolecular distances")') TRIM(fstring)
					CALL append_string(fstring)
				CASE (20)!cylindrical and polar distribution function
					smalltask=.FALSE.
					IF (nsteps<11) THEN
						PRINT *,"Your trajectory is really too short for that. Please use more timesteps."
						PRINT *,"The minimum of steps - even if only for debugging purposes - should be 11."
					ELSE
						CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
						IF (own_prefix) THEN
							own_prefix=.FALSE.
						ELSE
							analysis_number=analysis_number+1
						ENDIF
						PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
						PRINT *,"Would you like to take this shortcut? (y/n)"
						IF (user_input_logical()) THEN
							CALL append_string("distribution_simple ### sum rules and coulomb energy integral.")
							CYCLE
						ENDIF
						CALL user_distribution_input(parallelisation_possible,parallelisation_requested,number_of_molecules,filename_distribution)
						IF (parallelisation_requested) THEN
							CALL append_string("parallel_operation T ### turn on parallel operation")
							WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
							CALL append_string(fstring)
						ENDIF
						CALL append_string('distribution "'//TRIM(OUTPUT_PREFIX)//&
						&TRIM(filename_distribution)//'" ### compute distribution function')
						!enough information for the analysis.
						SKIP_ANALYSIS=.FALSE.
					ENDIF
				CASE (21)!jump analysis
					IF (nsteps==1) THEN
						PRINT *,"Not enough timesteps available."
						CYCLE
					ENDIF
					IF (.NOT.(parallelisation_requested)) THEN
						PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
						IF (user_input_logical()) parallelisation_requested=.TRUE.
					ENDIF
					IF (parallelisation_requested) THEN
						CALL append_string("parallel_operation T ### turn on parallel operation")
						WRITE(fstring,'("set_threads ",I0," ### set the number of threads to use to ",I0)') nthreads,nthreads
						CALL append_string(fstring)
					ENDIF
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("jump_velocity_simple ### jump analysis with default values.")
						CYCLE
					ENDIF
					smalltask=.FALSE.
					PRINT *,"Please give the first timestep to consider for the analysis:"
					startstep=user_input_integer(1,nsteps-1)
					PRINT *,"Please type the final timestep:"
					endstep=user_input_integer(startstep+1,nsteps)
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					IF (number_of_molecules==1) THEN
						PRINT *,"Only one molecule type available, which will be observed."
						molecule_type_index=1
					ELSE
						PRINT *,"Would you perform the jump analysis for just one molecule type (y), or for every type (n)?"
						IF (user_input_logical()) THEN
							PRINT *,"Please enter the index of the molecule type you wish to observe."
							molecule_type_index=user_input_integer(1,maxmol)
						ELSE
							molecule_type_index=-1
						ENDIF
					ENDIF
					deltasteps=nsteps-startstep
					PRINT *,"Please specify the maximum jump length (in timesteps):"
					deltasteps=user_input_integer(1,deltasteps)
					WRITE(fstring,'(A,I0," ",I0," ",I0," ",I0,A)')&
					&"jump_velocity ",molecule_type_index,startstep,endstep,deltasteps,&
					&" ### compute histogram / probability distribution of jump velocities."
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
				CASE (22)!show the atomic masses at this point
					CALL append_string("print_atomic_masses ### print atomic masses")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (23)!dump neighbour trajectory
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("dump_neighbour_traj_simple ### dump trajectory of nearest neighbours with default values.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide the range of timesteps for which to write the trajectory."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>50) smalltask=.FALSE.
					PRINT *,"Now, you have to specify molecule type and the number of the reference molecule."
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					PRINT *,"Please enter the molecule type index:"
					molecule_type_index=user_input_integer(1,maxmol)
					PRINT *,"Please enter the molecule index: ('Which one is it?')"
					molecule_index=user_input_integer(1,10000)
					PRINT *,"Please enter the index of the molecule type you wish to observe (as neighbours)."
					molecule_type_index_2=user_input_integer(1,maxmol)
					PRINT *,"How many neighbours would you like to observe?"
					PRINT *,"('10' will print the 10 closest neighbours around the reference molecule)"
					inputnumber=user_input_integer(1,100)
					WRITE(fstring,'(A,I0," ",I0," ",I0," ",I0," ",I0," ",I0,A)')&
					&"dump_neighbour_traj ",startstep,endstep,molecule_type_index,molecule_index,molecule_type_index_2,inputnumber,&
					&" ### dump neighbour trajectory of nearest neighbours."
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
				CASE (24)! dump dimers / closest neighbours
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"Writes all dimers (closest neighbours) in a certain timestep."
					PRINT *,"It is necessary to provide the timestep to be analysed."
					PRINT *,"Please enter this timestep as integer:"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"You now need to specify two molecule types as integers:"
					PRINT *,"First, the reference molecule type (the one that will be the centre)."
					IF (number_of_molecules<1) THEN
						maxmol=10000!unknown molecule number... expect the worst.
					ELSE
						maxmol=number_of_molecules
					ENDIF
					molecule_type_index=user_input_integer(1,maxmol)
					PRINT *,"Second, the molecule type of the neighbour molecule."
					molecule_type_index_2=user_input_integer(1,maxmol)
					PRINT *,"Would you like to dump all those dimer in one trajectory (y) or as separate files (n)?"
					IF (user_input_logical()) THEN
						!all in one trajectory file - append T
						WRITE(fstring,'("dump_dimers T ",I0," ",I0," ",I0,A,I0," around ",I0," for step ",I0,".")')&
						&startstep,molecule_type_index,molecule_type_index_2," ### dimers as trajectory file: type ",&
						&molecule_type_index_2,molecule_type_index,startstep
					ELSE
						!separate files - append F
						WRITE(fstring,'("dump_dimers F ",I0," ",I0," ",I0,A,I0," around ",I0," for step ",I0,".")')&
						&startstep,molecule_type_index,molecule_type_index_2," ### dimers as separate files: type ",&
						&molecule_type_index_2,molecule_type_index,startstep
					ENDIF
					CALL append_string(fstring)
					PRINT *,"The corresponding section has been added to the input file."
				CASE (25)!show the atomic charges at this point
					CALL append_string("print_atomic_charges ### print atomic charges")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (26)!show the atomic charges at this point
					CALL append_string("print_dipole_statistics ### print atomic charges")
					PRINT *,"The corresponding section has been added to the input file."
				CASE (27) !remove cores
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("remove_cores_simple ### subtract cores, write drude particles.")
						CYCLE
					ENDIF
					PRINT *,"It is necessary to provide a range of timesteps which are to be analysed."
					PRINT *,"To this end, please enter the first timestep to analyse"
					startstep=user_input_integer(1,nsteps)
					PRINT *,"Now, enter the last timestep of the range."
					endstep=user_input_integer(startstep,nsteps)
					IF ((endstep-startstep)>50) smalltask=.FALSE.
					WRITE(fstring,'("remove_cores ",I0," ",I0)') startstep,endstep
					WRITE(fstring,'(A," ### subtract cores, write drude particles for timesteps ",I0,"-",I0)')&
					&TRIM(fstring),startstep,endstep
					CALL append_string(fstring)
				CASE (28)!distances, usually rHH, dHH, rFF, dFF
					smalltask=.FALSE.
					CALL append_string("set_prefix "//TRIM(OUTPUT_PREFIX)//" ### This prefix will be used subsequently.")
					IF (own_prefix) THEN
						own_prefix=.FALSE.
					ELSE
						analysis_number=analysis_number+1
					ENDIF
					PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
					PRINT *,"Would you like to take this shortcut? (y/n)"
					IF (user_input_logical()) THEN
						CALL append_string("distances_simple ### intra- and intermolecular H-H and F-F distances, where available.")
						CYCLE
					ELSE
						IF (number_of_molecules<1) THEN
							maxmol=10000!unknown molecule number... expect the worst.
						ELSE
							maxmol=number_of_molecules
						ENDIF
						CALL user_distance_input(maxmol,filename_distance)
						CALL append_string('distance "'//TRIM(OUTPUT_PREFIX)//&
						&TRIM(filename_distance)//'" ### compute distances')
						!enough information for the analysis.
						SKIP_ANALYSIS=.FALSE.
					ENDIF
				CASE DEFAULT
					CALL report_error(0)
				END SELECT
			ENDDO
			IF (WRAP_TRAJECTORY) THEN!parallelisation is available...
				parallelisation_possible=.TRUE.
				IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
					PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
					IF (user_input_logical()) parallelisation_requested=.TRUE.
				ENDIF
			ENDIF
			CALL toggle_wrapping(parallelisation_possible,parallelisation_requested)!not necessary for the velocities, anyway. I'll nevertheless keep the 'wrapping_is_sensible=.FALSE.' there just in case.
			CALL toggle_sequential_read(parallelisation_possible,parallelisation_requested)
		END SUBROUTINE coordinates_user_input

		!This subroutine takes care of the sequential_read statement, and nothing else.
		SUBROUTINE toggle_sequential_read(parallelisation_possible,parallelisation_requested)
		IMPLICIT NONE
		LOGICAL,INTENT(IN) :: parallelisation_possible,parallelisation_requested
			IF (parallelisation_requested) THEN
				PRINT *,"Will load whole trajectory into RAM for parallel access."
				CALL append_string("sequential_read F ### load trajectory into RAM.")
				!parallel_operation has already been appended!
			ELSE
				PRINT *,"Do you want to read the trajectory file step by step? (y/n)"
				IF (parallelisation_possible) THEN
					PRINT *,"This uses very little RAM. Note that parallelisation is not turned on!."
					IF (user_input_logical()) THEN
						CALL append_string("sequential_read T ### read one timestep after the other")
					ELSE
						CALL append_string("sequential_read F ### load trajectory into RAM.")
						PRINT *,"Parts of the code have been parallelised for a reason."
						PRINT *,"You appear to believe that you have enough RAM - please consider using parallelisation."
					ENDIF
				ELSE
					PRINT *,"This uses very little RAM. Recommended, as parallelisation is not available anyways."
					IF (user_input_logical()) THEN
						CALL append_string("sequential_read T ### read one timestep after the other")
					ELSE
						CALL append_string("sequential_read F ### load trajectory into RAM.")
					ENDIF
				ENDIF
			ENDIF
		END SUBROUTINE toggle_sequential_read

		SUBROUTINE toggle_wrapping(parallelisation_possible,parallelisation_requested)
		IMPLICIT NONE
		LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
			IF ((wrapping_is_sensible).AND.(BOX_VOLUME_GIVEN)) THEN
				PRINT *,"Do you want to wrap the molecules you specified back into the box? (y/n)"
				PRINT *,"This could be useful if you want to export a snapshot."
				WRAP_TRAJECTORY=user_input_logical()
				IF (WRAP_TRAJECTORY) THEN
					parallelisation_possible=.TRUE.
					IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
						PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
						IF (user_input_logical()) THEN
							parallelisation_requested=.TRUE.
							CALL append_string("parallel_operation T ### turn on parallel operation")
						ENDIF
					ENDIF
					CALL append_string("wrap_trajectory T ### wrapping molecules (their centre of mass) into the box.")
				ELSE
					CALL append_string("wrap_trajectory F ### using unwrapped coordinates.")
				ENDIF
			ENDIF
		END SUBROUTINE toggle_wrapping

		!This subroutine appends the given string to the general input unit.
		SUBROUTINE append_string(inputstring)
		IMPLICIT NONE
		CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: inputstring
		INTEGER :: ios
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(FILENAME_GENERAL_INPUT),IOSTAT=ios,STATUS="OLD",POSITION="APPEND")!no input path is added for the general input file!
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			IF (PRESENT(inputstring)) THEN
				WRITE(8,FMT=*,IOSTAT=ios) TRIM(inputstring)
				IF (ios/=0) CALL report_error(46,exit_status=ios)
			ELSE
				ENDFILE 8
			ENDIF
			CLOSE(UNIT=8)
		END SUBROUTINE append_string

		!" 4 - explain program flow / analysis"
		SUBROUTINE explain_program_flow()
		IMPLICIT NONE
			PRINT *,"When invoked, any command line arguments will be treated as general input files."
			PRINT *,"For every valid one of these, the main program is invoked once."
			PRINT *,"The main program itself consists of two largely independent parts:"
			PRINT *,"The first part (this one) is the user interface."
			PRINT *,"It starts only if no general input file could be found."
			PRINT *,"(This can also be one of those specified in the command line)"
			PRINT *,"If there is a general input file, then the actual analysis starts,"
			PRINT *,"and the general input file is read line by line."
			PRINT *,"During the analysis, no further input is required."
			PRINT *,"The result is that errors in the general input file are not tolerated."
			PRINT *,"If a line cannot be read, then the execution stops and an error is printed."
			PRINT *,"When a keyword linking to a separate input file is found,"
			PRINT *,"then the corresponding module is invoked to read and run the separate input."
			PRINT *,"Some of the calculations are quite involved and might run for hours or days."
			PRINT *,"The required real time can be greatly reduced by using parallelisation."
			PRINT *,"Parallelisation uses the OpenMP library to share the workload."
			PRINT *,"To use the parallelisation, it is necessary to load the whole trajectory into RAM."
			PRINT *,"The reason for this is that the bottleneck is usually fileIO."
		END SUBROUTINE explain_program_flow

		!" 5 - which format does the output take?"
		SUBROUTINE explain_output_format()
		IMPLICIT NONE
			PRINT *,"All output file, including structures, are written to the output folder."
			PRINT *,"These files contain a header with variable names and the reference if necessary."
			PRINT *,"The names of the output files from a certain type of calculation are fixed."
			PRINT *,"This means that if you want to perform, say, two dihedral condition analyses,"
			PRINT *,"the files will be overwritten. To avoid this, request a prefix such as 'cation_'."
			PRINT *,"The 'timeline' printed in some files is obtained from timestep * time scaling factor."
		END SUBROUTINE explain_output_format

		!" 6 - how to format the input trajectory?"
		SUBROUTINE explain_trajectory_format()
		IMPLICIT NONE
			PRINT *,"The trajectory is expected to be in lammps format."
			PRINT *,"For each timestep, there is a header and a body."
			PRINT *,"The header should look like this:"
			PRINT *
			PRINT *,"ITEM: TIMESTEP"
			PRINT *,"0"
			PRINT *,"ITEM: NUMBER OF ATOMS"
			PRINT *,"20480"
			PRINT *,"ITEM: BOX BOUNDS pp pp pp"
			PRINT *,"0 63.9223"
			PRINT *,"0 63.9223"
			PRINT *,"0 63.9223"
			PRINT *,"ITEM: ATOMS element xu yu zu"
			PRINT *
			PRINT *,"After this follows the body, consisting of one line per atom."
			PRINT *,"each line begins with the element label (e.g. 'C'),"
			PRINT *,"followed by three floating point (=real) number."
			PRINT *,"Depending on the type of analysis you need, these have"
			PRINT *,"to be either cartesian coordinates or velocities."
			PRINT *,"For coordinates, your lammps input file should include something like:"
			PRINT *
			PRINT *,"dump TRAJECTORY all custom 1000 trajectory.lmp element xu yu zu"
			PRINT *,"dump_modify TRAJECTORY element C F N O S C C C C H H N sort id"
			PRINT *
			PRINT *,"whereas for velocities, 'xu yu zu' has to be changed to 'vx vy vz'."
			PRINT *,"to obtain sensible results, consistent ordering is imperative."
			PRINT *,"This is the purpose of the second line given in the example above."
			PRINT *,"Important final note:"
			PRINT *,"For performance issues, the format is not checked during read!"
		END SUBROUTINE explain_trajectory_format

END SUBROUTINE initialise_global

SUBROUTINE finalise_global()
USE MOLECULAR
USE SETTINGS
IMPLICIT NONE
42	FORMAT ("    ########   #####")
44	FORMAT ("        ########")
43	FORMAT ("        #      #")
28	FORMAT ("    #   ##  #      #")
45	FORMAT ("    #      #       #")
13	FORMAT ("    ##### #### #####")
	CALL finalise_molecular()!also closes unit 9, if still open.
	CLOSE(UNIT=7)
	IF (DEVELOPERS_VERSION) THEN
		PRINT *
		WRITE(*,44)
		WRITE(*,43)
		WRITE(*,42)
		WRITE(*,45)
		WRITE(*,45)
		WRITE(*,28)
		WRITE(*,13)
		WRITE(*,43)
		WRITE(*,44)
		PRINT *
	ENDIF
END SUBROUTINE finalise_global

!The analysis tasks are read line wise from general.inp, and executed subsequently.
SUBROUTINE run_analysis()
USE MOLECULAR
USE DEBUG
USE AUTOCORRELATION
USE DIFFUSION
USE SETTINGS
USE DISTRIBUTION
USE DISTANCE
IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ 	FUNCTION OMP_get_max_threads()
	 !$ 	INTEGER :: OMP_get_max_threads
	 !$ 	END FUNCTION OMP_get_max_threads
	 !$ 	FUNCTION OMP_get_num_procs()
	 !$ 	INTEGER :: OMP_get_num_procs
	 !$ 	END FUNCTION OMP_get_num_procs
	 !$ 	SUBROUTINE OMP_set_num_threads(number_of_threads)
	 !$ 	INTEGER,INTENT(IN) :: number_of_threads
	 !$ 	END SUBROUTINE OMP_set_num_threads
	 !$ END INTERFACE
CHARACTER(LEN=1024) :: inputstring,dummy
INTEGER :: ios,n
	!open file, read head
	CALL read_head()
	!read body of general.inp
	CALL read_body()
	!unit 7 will be closed by finalise_global. This is so it is not kept open in error cases.
	CONTAINS
	
		SUBROUTINE read_head()
		IMPLICIT NONE
		LOGICAL :: connected
			WRITE(*,*) "streaming input from general input file"
			IF (INFORMATION_IN_TRAJECTORY=="UNK") CALL report_error(55)
			WRITE(*,*)
			INQUIRE(UNIT=7,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=7)
			OPEN(UNIT=7,FILE=TRIM(FILENAME_GENERAL_INPUT),ACTION='READ',IOSTAT=ios)
			IF (ios/=0) CALL report_error(18,exit_status=ios)
			REWIND 7
			!Skip over head of general.inp
			DO n=1,HEADER_LINES,1
				READ(7,*)
			ENDDO
		END SUBROUTINE read_head

		SUBROUTINE read_body()
		IMPLICIT NONE
		INTEGER :: inputinteger,startstep,endstep,inputinteger2,inputinteger3,inputinteger4,counter
		LOGICAL :: inputlogical
		REAL(KIND=WORKING_PRECISION) :: inputreal
		REAL :: lower,upper
			DO n=1,MAXITERATIONS,1
				READ(7,IOSTAT=ios,FMT=*) inputstring
				IF (ios<0) THEN
					IF (VERBOSE_OUTPUT) THEN 
						WRITE(*,*) "End of file condition in general input file."
						WRITE(*,*) "--> terminating analysis peacefully."
					ENDIF
					EXIT
				ENDIF
				IF (ios/=0) THEN
					CALL report_error(19,exit_status=ios)
					EXIT
				ENDIF
				IF (TRIM(inputstring)=="autocorrelation") inputstring="correlation"!support for older nomenclature
				IF (TRIM(inputstring)=="dihedral") inputstring="correlation"!support for synonyms
				IF (TRIM(inputstring)=="reorientation") inputstring="correlation"!support for synonyms
				IF (TRIM(inputstring)=="rmm-vcf") inputstring="correlation"!support for synonyms
				IF (TRIM(inputstring)=="conductivity") inputstring="correlation"!support for synonyms
				IF (TRIM(inputstring)=="velocity") inputstring="correlation"!support for synonyms
				IF (TRIM(inputstring)=="show_drude_settings") inputstring="show_drude"!support for synonyms
				IF (TRIM(inputstring)=="show_drudes") inputstring="show_drude"!support for synonyms
				IF (TRIM(inputstring)=="drude_temperature") inputstring="drude_temp"!support for synonyms
				IF (TRIM(inputstring)=="remove_drude") inputstring="remove_drudes"!support for synonyms
				IF (TRIM(inputstring)=="remove_core") inputstring="remove_cores"!support for synonyms
				IF (TRIM(inputstring)=="exit") inputstring="quit"!support for synonyms
				IF (TRIM(inputstring)=="read_sequential") inputstring="sequential_read"!support for synonyms
				IF (TRIM(inputstring)=="time_scaling_factor") inputstring="time_scaling"!support for synonyms
				IF (TRIM(inputstring)=="parallel_execution") inputstring="parallel_operation"!support for synonyms
				IF (TRIM(inputstring)=="print_atom_masses") inputstring="print_atomic_masses"!support for synonyms
				IF (TRIM(inputstring)=="print_atomic_weights") inputstring="print_atomic_masses"!support for synonyms
				IF (TRIM(inputstring)=="print_atom_weights") inputstring="print_atomic_masses"!support for synonyms
				IF (TRIM(inputstring)=="print_atom_charges") inputstring="print_atomic_charges"!support for synonyms
				IF (TRIM(inputstring)=="print_dipole_properties") inputstring="print_dipole_statistics"!support for synonyms
				!so far, only error handling has occurred. Now, check what the corresponding task was, re-read with the appropriate formatting, and start analysis.
				IF (VERBOSE_OUTPUT) WRITE(*,'(" Line ",I0,": ",A)')  n+HEADER_LINES,"'"//TRIM(inputstring)//"'"
				SELECT CASE (TRIM(inputstring))
				CASE ("quit")
					WRITE(*,*) "exiting analysis."
					EXIT
				CASE ("show_settings")
					CALL show_settings()
				CASE ("print_atomic_masses")
					CALL print_atomic_masses()
				CASE ("print_atomic_charges")
					CALL print_atomic_charges()
				CASE ("print_dipole_statistics")
					CALL print_dipole_statistics()
				CASE ("show_drude")
					CALL show_drude_settings()
				CASE ("switch_to_com","barycentric","barycenter","barycentre") !Module MOLECULAR
					CALL switch_to_barycenter()
				CASE ("verbose_output")
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ELSE
						VERBOSE_OUTPUT=inputlogical
						WRITE(*,*) "setting 'VERBOSE_OUTPUT' to ",VERBOSE_OUTPUT
					ENDIF
				CASE ("time_output")
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ELSE
						TIME_OUTPUT =inputlogical
						WRITE(*,*) "setting 'TIME_OUTPUT' to ",TIME_OUTPUT
					ENDIF
				CASE ("sequential_read")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "skip line (sequential_read)"
				CASE ("wrap_trajectory")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "skip line (wrap_trajectory)"
				CASE ("trajectory_type")
					IF (VERBOSE_OUTPUT) WRITE(*,*) "skip line (trajectory_type)"
				CASE ("dump_snapshot") !Module DEBUG
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger,inputlogical
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					CALL check_timestep(inputinteger)
					CALL dump_snapshot(inputinteger,inputlogical)
					WRITE(*,'(A,I0,A)') " Snapshot of step ",inputinteger," written to output folder."
				CASE ("dump_snapshot_simple") !Module DEBUG
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					CALL dump_snapshot(1,.FALSE.)
					WRITE(*,'(A,I0,A)') " Snapshot of first step written to output folder."
				CASE ("dump_split") !Module DEBUG
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,startstep,endstep
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					CALL check_timesteps(startstep,endstep)
					WRITE(*,'(A,I0,A)') " Trajectory will be split into ",&
					&give_number_of_molecule_types()," molecule types."
					WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
					CALL dump_split(startstep,endstep,"xyz")
				CASE ("dump_split_simple") !Module DEBUG
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					WRITE(*,'(A,I0,A)') " Trajectory will be split into ",&
					&give_number_of_molecule_types()," molecule types - simple mode."
					startstep=1
					endstep=give_number_of_timesteps()
					WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
					CALL dump_split(startstep,endstep,"xyz")
				CASE ("remove_drudes") !Module DEBUG
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,startstep,endstep
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					IF (are_drudes_assigned()) THEN
						CALL check_timesteps(startstep,endstep)
						WRITE(*,*) "Writing trajectory with drude particles merged into cores."
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
						IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory type will be '",TRAJECTORY_TYPE,"'"
						CALL remove_drudes(startstep,endstep,TRAJECTORY_TYPE,.TRUE.)
					ELSE
						CALL report_error(91)
					ENDIF
				CASE ("remove_drudes_simple") !Module DEBUG
					IF (are_drudes_assigned()) THEN
						WRITE(*,*) "Writing trajectory with drude particles merged into cores - simple mode."
						startstep=1
						endstep=give_number_of_timesteps()
						CALL check_timesteps(startstep,endstep)
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
						IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory type will be '",TRAJECTORY_TYPE,"'"
						CALL remove_drudes(startstep,endstep,TRAJECTORY_TYPE,.TRUE.)
					ELSE
						CALL report_error(91)
					ENDIF
				CASE ("remove_cores") !Module DEBUG
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,startstep,endstep
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					IF (are_drudes_assigned()) THEN
						CALL check_timesteps(startstep,endstep)
						WRITE(*,*) "Writing trajectory with only drude particles (minus cores)."
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
						IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory type will be '",TRAJECTORY_TYPE,"'"
						CALL remove_cores(startstep,endstep,TRAJECTORY_TYPE)
					ELSE
						CALL report_error(91)
					ENDIF
				CASE ("remove_cores_simple") !Module DEBUG
					IF (are_drudes_assigned()) THEN
						WRITE(*,*) "Writing trajectory with only drude particles (minus cores) - simple mode."
						startstep=1
						endstep=give_number_of_timesteps()
						CALL check_timesteps(startstep,endstep)
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
						IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory type will be '",TRAJECTORY_TYPE,"'"
						CALL remove_cores(startstep,endstep,TRAJECTORY_TYPE)
					ELSE
						CALL report_error(91)
					ENDIF
				CASE ("contact_distance") !Module DEBUG
					IF (BOX_VOLUME_GIVEN) THEN
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger,inputinteger2
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						CALL check_timestep(inputinteger)
						CALL contact_distance(inputinteger,inputinteger2)
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("contact_distance_simple")
					IF (BOX_VOLUME_GIVEN) THEN
						inputinteger=1
						inputinteger2=-1
						CALL check_timestep(inputinteger)
						CALL contact_distance(inputinteger,inputinteger2)
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("dump_single") !Module DEBUG
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical,startstep,endstep,inputinteger,inputinteger2
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					CALL check_timesteps(startstep,endstep)
					WRITE(*,'(A,I0,A,I0,A)') " Trajectory for molecule ",inputinteger2," of type ",inputinteger," will be written."
					IF (inputlogical) THEN
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,", using centre of mass)"
					ELSE
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,", NOT using centre of mass)"
					ENDIF
					CALL dump_single(inputlogical,startstep,endstep,inputinteger,inputinteger2)
				CASE ("dump_single_simple")
					CALL report_error(113)
				CASE ("dump_cut") !Module DEBUG
					IF (BOX_VOLUME_GIVEN) THEN
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical,startstep,endstep,inputinteger,inputinteger2,inputreal
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						CALL check_timesteps(startstep,endstep)
						WRITE(*,'(A,I0,A,I0,A)')&
						&" Trajectory for molecule ",inputinteger2," of type ",inputinteger," will be written, including neighbours."
						IF (inputlogical) THEN
							WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,&
							&", referencing to centre of mass of reference molecule in every step)"
						ELSE
							WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,&
							&", referencing to centre of mass of reference molecule in first step)"
						ENDIF
						IF (inputreal>=1.0d0) THEN
							WRITE(*,'(" The threshold for neighbourhood is currently set to ",F0.2," Angström.")') inputreal
						ELSE
							WRITE(*,*) "Small threshold used (less than 1 Angström)"
						ENDIF
						CALL dump_cut(inputlogical,startstep,endstep,inputinteger,inputinteger2,inputreal)
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("dump_cut_simple")
					CALL report_error(113)
				CASE ("dump_dimers") !Module DEBUG
					IF (BOX_VOLUME_GIVEN) THEN
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical,startstep,inputinteger,inputinteger2
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						CALL check_timestep(startstep)
						WRITE(*,'(A,I0,A,I0,A,I0,A)') " Dumping dimers, i.e. closest molecules of type ",&
						&inputinteger2," around ",inputinteger," for timestep ",startstep,"."
						IF (inputlogical) THEN
							WRITE(*,*) "(Combined in a single 'trajectory' file)"
						ELSE
							WRITE(*,*) "(As separate files for each dimer)"
						ENDIF
						CALL dump_dimers(startstep,inputlogical,inputinteger,inputinteger2)
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("dump_dimers_simple")
					CALL report_error(113)
				CASE ("dump_neighbour_traj")
					IF (BOX_VOLUME_GIVEN) THEN
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,startstep,endstep,&
						&inputinteger,inputinteger2,inputinteger3,inputinteger4
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						CALL check_timesteps(startstep,endstep)
						WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') " Dumping neighbour trajectory, i.e. closest molecules of type ",&
						&inputinteger3," around ",inputinteger," for timestep ",startstep," to ",endstep,"."
						WRITE(*,'(" The molecule with index ",I0," (and type ",I0,") is the reference.")') inputinteger2,inputinteger
						WRITE(*,'(" ",I0," neighbours (of type ",I0,") will be written for each step.")') inputinteger4,inputinteger3
						CALL dump_neighbour_traj(.FALSE.,startstep,endstep,&
						&inputinteger,inputinteger2,inputinteger3,inputinteger4)
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("dump_neighbour_traj_simple")
					IF (BOX_VOLUME_GIVEN) THEN
						startstep=1 !first timestep
						endstep=give_number_of_timesteps() !second timestep
						inputinteger2=1 !molecule index (reference)
						inputinteger4=2 !number of neighbours
						IF (give_number_of_molecule_types()==2) THEN
							WRITE(*,'(A,I0,A,I0,A)') " Dumping neighbour trajectories, i.e. closest molecules of type ",&
							&inputinteger3," around ",inputinteger," (and vice versa) for all timesteps."
							WRITE(*,'(" The first molecule is always the reference, and 2 neighbours will be written.")')
							inputinteger=1 !molecule type index 1(reference)
							inputinteger3=2 !molecule type index 2 (the observed one)
							CALL dump_neighbour_traj(.FALSE.,startstep,endstep,&
							&inputinteger,inputinteger2,inputinteger3,inputinteger4)
							inputinteger=2 !molecule type index 1(reference)
							inputinteger3=1 !molecule type index 2 (the observed one)
							CALL dump_neighbour_traj(.FALSE.,startstep,endstep,&
							&inputinteger,inputinteger2,inputinteger3,inputinteger4)
						ELSE
							WRITE(*,*) "Dumping neighbour trajectories, i.e. closest molecules of all types"&
							&//" around themselves for all timesteps."
							WRITE(*,'(" The first molecule is always the reference, and 2 neighbours will be written.")')
							DO counter=1,give_number_of_molecule_types(),1
								inputinteger=counter !molecule type index 1(reference)
								inputinteger3=counter !molecule type index 2 (the observed one)
								CALL dump_neighbour_traj(.FALSE.,startstep,endstep,&
								&inputinteger,inputinteger2,inputinteger3,inputinteger4)
							ENDDO
						ENDIF
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("cubic_box_edge","cubic_box")
					IF (BOX_VOLUME_GIVEN) CALL report_error(92)
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,lower,upper
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					IF (lower>upper) THEN
						CALL report_error(93)
					ELSE
						CALL set_cubic_box(lower,upper)
						WRITE(*,'(" Box boundaries set to:")') 
						WRITE(*,'("   lower bound: ",E12.6)') lower
						WRITE(*,'("   upper bound: ",E12.6)') upper
					ENDIF
				CASE ("cubic_box_edge_simple")
					CALL report_error(113)
				CASE ("convert") !Module DEBUG
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					WRITE(*,*) "Reduce Trajectory to centre of mass for each molecule type."
					IF (inputlogical) WRITE(*,*) "An adjusted molecular input file will be written, too."
					IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory type will be '",TRAJECTORY_TYPE,"'"
					CALL convert(inputlogical,TRAJECTORY_TYPE)
				CASE ("convert_simple") !Module DEBUG
					WRITE(*,*) "Reduce Trajectory to centre of mass for each molecule type."
					WRITE(*,*) "An adjusted molecular input file will be written, too."
					IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory type will be '",TRAJECTORY_TYPE,"'"
					CALL convert(.TRUE.,TRAJECTORY_TYPE)
				CASE ("temperature") !Module DEBUG
					IF (WRAP_TRAJECTORY) THEN
						CALL report_error(72)
					ELSE
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger,startstep,endstep
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						WRITE(*,*) "Calculating the kinetic temperature, based on NkT=mv²."
						CALL check_timesteps(startstep,endstep)
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
						CALL report_temperature(inputinteger,startstep,endstep)
					ENDIF
				CASE ("temperature_simple") !Module DEBUG
					IF (WRAP_TRAJECTORY) THEN
						CALL report_error(72)
					ELSE
						WRITE(*,*) "Calculating the kinetic temperature for first timestep, based on NkT=mv² (simple mode)."
						inputinteger=-1
						startstep=1
						endstep=1
						CALL check_timesteps(startstep,endstep)
						CALL report_temperature(inputinteger,startstep,endstep)
					ENDIF
				CASE ("gyradius") !Module DEBUG
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger,startstep,endstep
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					WRITE(*,*) "Calculating ensemble average of radius of gyration."
					CALL check_timesteps(startstep,endstep)
					WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
					CALL report_gyradius(inputinteger,startstep,endstep)
				CASE ("gyradius_simple") !Module DEBUG
					WRITE(*,*) "Calculating ensemble average of radius of gyration - simple mode."
					WRITE(*,*) "(First 100 timesteps - all molecule types.)"
					startstep=1
					endstep=100
					CALL check_timesteps(startstep,endstep)
					CALL report_gyradius(-1,startstep,endstep)
				CASE ("drude_temp") !Module DEBUG
					IF (WRAP_TRAJECTORY) THEN
						CALL report_error(72)
					ELSE
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,startstep,endstep
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						WRITE(*,*) "Calculating the kinetic temperature, based on NkT=mv²."
						WRITE(*,*) "(Extended support of drude particles requires manual drude assignment)"
						CALL check_timesteps(startstep,endstep)
						WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
						CALL report_drude_temperature(startstep,endstep)
					ENDIF
				CASE ("drude_temp_simple") !Module DEBUG
					IF (WRAP_TRAJECTORY) THEN
						CALL report_error(72)
					ELSE
						WRITE(*,*) "Calculating the kinetic temperature of first timestep, based on NkT=mv² (simple mode)."
						WRITE(*,*) "(Extended support of drude particles requires manual drude assignment)"
						CALL report_drude_temperature(1,1)
					ENDIF
				CASE ("set_threads") !Module DEBUG
				 !$ IF (.FALSE.) THEN
						WRITE(*,*) "keyword 'set_threads' has no effect (Compiler not OpenMP compliant)"
				 !$ ENDIF
				 !$ BACKSPACE 7
				 !$ READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger
				 !$ IF (ios/=0) THEN
				 !$ 	CALL report_error(19,exit_status=ios)
				 !$ 	EXIT
				 !$ ENDIF
				 !$ IF (inputinteger<1) THEN
				 !$ 	inputinteger=OMP_get_max_threads()
				 !$ ELSEIF (inputinteger>OMP_get_num_procs()) THEN
				 !$ 	CALL report_error(42,exit_status=inputinteger)
				 !$ 	inputinteger=OMP_get_num_procs()
				 !$ ENDIF
				 !$ CALL OMP_set_num_threads(inputinteger)
				 !$ WRITE(*,'(" number of threads set to ",I0)') inputinteger
				CASE ("set_threads_simple")
				 !$ IF (.FALSE.) THEN
						WRITE(*,*) "keyword 'set_threads_simple' has no effect (Compiler not OpenMP compliant)"
				 !$ ENDIF
				 !$ inputinteger=OMP_get_max_threads()
				 !$ CALL OMP_set_num_threads(inputinteger)
				 !$ WRITE(*,'(" number of threads set to ",I0)') inputinteger
				CASE ("time_scaling") !Module SETTINGS
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ELSE
						IF (inputinteger<1) THEN
							CALL report_error(106,exit_status=inputinteger)
							inputinteger=TIME_SCALING_FACTOR_DEFAULT
						ENDIF
						WRITE(*,'(" scaling timelines with ",I0)') inputinteger
						TIME_SCALING_FACTOR=inputinteger
					ENDIF
				CASE ("set_prefix") !Module SETTINGS
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,OUTPUT_PREFIX
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					WRITE(*,*) "prefix set to '",TRIM(ADJUSTL(OUTPUT_PREFIX)),"'"
				CASE ("dump_example","dump_example_simple") !Module DEBUG
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					CALL dump_example()
					WRITE(*,*) "Example molecules written to output folder."
				CASE ("correlation") !Module AUTOCORRELATION
					!the (INFORMATION_IN_TRAJECTORY=="VEL") test is done in perform_autocorrelation()!
					!same for the wrap test.
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,dummy
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					FILENAME_AUTOCORRELATION_INPUT=dummy
					WRITE(*,*) "(Auto)correlation module invoked."
					CALL perform_autocorrelation()
				CASE ("correlation_simple")
					CALL report_error(113)
				CASE ("distance","distances") !Module DISTANCE
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,dummy
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					FILENAME_DISTANCE_INPUT=dummy
					WRITE(*,*) "distance module invoked."
					CALL perform_distance_analysis()
				CASE ("distance_simple","distances_simple") !Module DISTANCE
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					WRITE(*,*) "Calculating, where possible, intra- and intermolecular distances for H and F."
					IF (write_simple_intramolecular_distances()) THEN
						WRITE(*,*) "Intramolecular distance calculation, simple mode:"
						CALL perform_distance_analysis()
					ELSE
						WRITE(*,*) "Can't do intramolecular distances. Do manually if necessary."
					ENDIF
					IF (write_simple_intermolecular_distances()) THEN
						WRITE(*,*) "Intermolecular distance calculation, simple mode:"
						CALL perform_distance_analysis()
					ELSE
						WRITE(*,*) "Can't do intermolecular distances. Do manually if necessary."
					ENDIF
				CASE ("dihedral") !Module AUTOCORRELATION
					!the (INFORMATION_IN_TRAJECTORY=="VEL") test is done in perform_autocorrelation()!
					!same for the wrap test.
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,dummy
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					FILENAME_AUTOCORRELATION_INPUT=dummy
					WRITE(*,*) "(Auto)correlation module invoked."
					CALL perform_autocorrelation()
				CASE ("dihedral_simple")
					CALL report_error(113)
				CASE ("diffusion") !Module DIFFUSION
					IF (WRAP_TRAJECTORY) THEN
						CALL report_error(72)
					ELSE
						IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,dummy
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
						FILENAME_DIFFUSION_INPUT=dummy
						WRITE(*,*) "Diffusion module invoked."
						CALL perform_diffusion_analysis()
					ENDIF
				CASE ("diffusion_simple")
					IF (WRAP_TRAJECTORY) THEN
						CALL report_error(72)
					ELSE
						IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
						CALL write_simple_diffusion()
						WRITE(*,*) "Diffusion module invoked."
						CALL perform_diffusion_analysis()
					ENDIF
				CASE ("distribution") !Module DISTRIBUTION
					IF (BOX_VOLUME_GIVEN) THEN
						IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
						BACKSPACE 7
						READ(7,IOSTAT=ios,FMT=*) inputstring,dummy
						IF (ios/=0) THEN
							CALL report_error(19,exit_status=ios)
							EXIT
						ENDIF
							FILENAME_DISTRIBUTION_INPUT=dummy
						WRITE(*,*) "Distribution module invoked."
						CALL perform_distribution_analysis()
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("distribution_simple")
					IF (BOX_VOLUME_GIVEN) THEN
						IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
						WRITE(*,*) "Distribution module invoked - simple mode:"
						WRITE(*,*) "Only sum rules and Coulomb energy integral."
						CALL write_simple_sumrules()
						CALL perform_distribution_analysis()
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("charge_arm_simple")
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					WRITE(*,*) "Distribution module invoked - simple mode:"
					WRITE(*,*) "Charge Arm distribution. Requires charges to be initialised."
					CALL write_simple_charge_arm()
					CALL perform_distribution_analysis()
				CASE ("conductivity_simple")
					IF (BOX_VOLUME_GIVEN) THEN
						WRITE(*,*) "Autocorrelation module invoked - simple mode:"
						WRITE(*,*) "Calculating overall conductivity."
						CALL write_simple_conductivity()
						CALL perform_autocorrelation()
					ELSE
						CALL report_error(41)
					ENDIF
				CASE ("error_output")
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ELSE
						ERROR_OUTPUT=inputlogical
						WRITE(*,*) "setting 'ERROR_OUTPUT' to ",ERROR_OUTPUT
					ENDIF
				CASE ("parallel_operation")
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputlogical
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ELSE
						PARALLEL_OPERATION=inputlogical
						IF (READ_SEQUENTIAL) THEN
							WRITE(*,*) "parallel operation is not available with sequential read."
							PARALLEL_OPERATION=.FALSE.
						ENDIF
					 !$ IF (.FALSE.) THEN
							WRITE(*,*) "keyword 'parallel_operation' has no effect (Compiler not OpenMP compliant)"
					 !$ ENDIF
					 !$ WRITE(*,*) "setting 'PARALLEL_OPERATION' to ",PARALLEL_OPERATION
					ENDIF
				CASE ("jump_velocity_simple")
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					WRITE(*,*) "Constructing a histogram of jump duration vs jump velocity."
					inputinteger=give_number_of_timesteps()/10
					IF (inputinteger<10) inputinteger=give_number_of_timesteps()
					CALL jump_analysis(inputinteger,100,-1,1,give_number_of_timesteps(),.TRUE.)
				CASE ("jump_velocity")
					IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
					BACKSPACE 7
					READ(7,IOSTAT=ios,FMT=*) inputstring,inputinteger,startstep,endstep,inputinteger2
					IF (ios/=0) THEN
						CALL report_error(19,exit_status=ios)
						EXIT
					ENDIF
					WRITE(*,*) "Constructing a histogram of jump duration vs jump velocity."
					CALL check_timesteps(startstep,endstep)
					WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
					CALL jump_analysis(inputinteger2,100,inputinteger,startstep,endstep,.TRUE.)
				CASE ("DEBUG")
					!Here is some space for testing stuff
					WRITE(*,*) "################################"
					CALL perform_distance_analysis()
					!FILENAME_DISTANCE_INPUT="intra.inp"
					!CALL perform_distance_analysis()
					WRITE(*,*) "################################"
				CASE DEFAULT
					IF ((inputstring(1:1)=="#").OR.(inputstring(1:1)=="!")) THEN
						IF (VERBOSE_OUTPUT) WRITE(*,'(" (is commented out)")')
					ELSE
						CALL report_error(20,n+HEADER_LINES)!HEADER_LINES = number of fixed lines in general input file
					ENDIF
				END SELECT
				CALL timing()
			ENDDO
			!note that an EXIT condition in this loop effectively equals the soft stop in some error reports.

		END SUBROUTINE read_body

		!the following subroutine prints the settings and how to influence them.
		SUBROUTINE show_settings()
		IMPLICIT NONE
	15	FORMAT ("    ",A," ",A)
	16	FORMAT ("    ",A," ",I0)
			WRITE(*,*) "Printing current global settings."
			WRITE(*,15) "VERBOSE_OUTPUT      ",TRIM(logical_to_yesno(VERBOSE_OUTPUT))
			WRITE(*,15) "TIME_OUTPUT         ",TRIM(logical_to_yesno(TIME_OUTPUT))
			WRITE(*,15) "DEVELOPERS_VERSION  ",TRIM(logical_to_yesno(DEVELOPERS_VERSION))
			WRITE(*,15) "ERROR_OUTPUT        ",TRIM(logical_to_yesno(ERROR_OUTPUT))
			WRITE(*,15) "READ_SEQUENTIAL     ",TRIM(logical_to_yesno(READ_SEQUENTIAL))
			WRITE(*,15) "BOX_VOLUME_GIVEN    ",TRIM(logical_to_yesno(BOX_VOLUME_GIVEN))
			WRITE(*,15) "WRAP_TRAJECTORY     ",TRIM(logical_to_yesno(WRAP_TRAJECTORY))
			WRITE(*,15) "DISCONNECTED        ",TRIM(logical_to_yesno(DISCONNECTED))
			WRITE(*,16) "GLOBAL_ITERATIONS   ",GLOBAL_ITERATIONS
			WRITE(*,16) "TIME_SCALING_FACTOR ",TIME_SCALING_FACTOR
			WRITE(*,16) "HEADER_LINES_GINPUT ",HEADER_LINES
			WRITE(*,16) "CURRENT_ERROR_CODE  ",ERROR_CODE
			WRITE(*,16) "ERROR_COUNT         ",give_error_count()
			WRITE(*,*) '   OUTPUT_PREFIX        "',TRIM(OUTPUT_PREFIX),'"'
			WRITE(*,*) '   TRAJECTORY_TYPE      "',TRIM(TRAJECTORY_TYPE),'"'
			WRITE(*,*) '   INFO_IN_TRAJECTORY   "',TRIM(INFORMATION_IN_TRAJECTORY),'"'
			WRITE(*,15) "PARALLEL_OPERATION  ",TRIM(logical_to_yesno(PARALLEL_OPERATION))
		 !$ IF(.FALSE.) THEN
			WRITE(*,*) "-fopenmp flag not set! PARALLEL_OPERATION has no effect."
		 !$ ENDIF
		 !$ WRITE(*,*) "OMP flag set!"
		 !$ WRITE(*,16) "NUMBER_OF_THREADS     ",OMP_get_num_threads()
		 !$ WRITE(*,16) "NUMBER_OF_PROCS       ",OMP_get_num_procs()
		 !$ WRITE(*,16) "MAX_NUMBER_OF_THREADS ",OMP_get_max_threads()
			CALL show_molecular_settings()
			WRITE(*,*) "current default (input) Filenames:"
			WRITE(*,*) '   TRAJECTORY      "',TRIM(FILENAME_TRAJECTORY),'"'
			WRITE(*,*) '   GENERAL         "',TRIM(FILENAME_GENERAL_INPUT),'"'
			WRITE(*,*) '   MOLECULAR       "',TRIM(FILENAME_MOLECULAR_INPUT),'"'
			WRITE(*,*) '   AUTOCORRELATION "',TRIM(FILENAME_AUTOCORRELATION_INPUT),'"'
			WRITE(*,*) '   DIFFUSION       "',TRIM(FILENAME_DIFFUSION_INPUT),'"'
			WRITE(*,*) '   DISTRIBUTION    "',TRIM(FILENAME_DISTRIBUTION_INPUT),'"'
			WRITE(*,*) '   DISTANCE        "',TRIM(FILENAME_DISTANCE_INPUT),'"'
			WRITE(*,*) "Paths:"
			WRITE(*,*) '   TRAJECTORY "',TRIM(PATH_TRAJECTORY),'"'
			WRITE(*,*) '   INPUT      "',TRIM(PATH_INPUT),'"'
			WRITE(*,*) '   OUTPUT     "',TRIM(PATH_OUTPUT),'"'
			WRITE(*,'(" currently assuming ",I0," as time difference between steps")') TIME_SCALING_FACTOR
			WRITE(*,*) "assuming Angström as distance unit and femtoseconds as time unit."
			WRITE(*,*) "Using these constants:"
			WRITE(*,'("    boltzmann constant:  ",E15.9E2," J/K")') boltzmann
			WRITE(*,'("    avogadro constant:   ",E15.9E2," 1/mol")') avogadro
			WRITE(*,'("    elementary charge:   ",E15.9E2," C")') elementary_charge
			WRITE(*,'("    archimedes constant: ",E15.9E2)') pi
			WRITE(*,'("    vacuum permittivity: ",E15.9E2," F/m")') vacuum_permittivity
		END SUBROUTINE show_settings

END SUBROUTINE run_analysis

SUBROUTINE initialise_command_line_arguments()
USE SETTINGS
USE RECOGNITION
IMPLICIT NONE
CHARACTER(LEN=1024) :: inputstring,trajstring
INTEGER :: allocstatus,file_counter,i,charnum
LOGICAL :: valid_filename,file_exists,command_line_used,skip_next,missing_argument
	IF (COMMAND_ARGUMENT_COUNT()>0) THEN
		GLOBAL_ITERATIONS=0
		!The maximum of sensible members in the GENERAL_INPUT_FILENAMES list is the number of passed arguments.
		ALLOCATE(GENERAL_INPUT_FILENAMES(COMMAND_ARGUMENT_COUNT()),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(66,exit_status=allocstatus)
		!Iterate over the input arguments
		command_line_used=.FALSE.
		skip_next=.FALSE.
		DO file_counter=1,COMMAND_ARGUMENT_COUNT(),1
			IF (skip_next) THEN
				skip_next=.FALSE.
				CYCLE
			ENDIF
			CALL GET_COMMAND_ARGUMENT(file_counter,inputstring)
			inputstring=ADJUSTL(inputstring)
			!Check for command line switches.
			SELECT CASE (inputstring(1:2))
			CASE ("-d")
				IF (LEN(TRIM(inputstring))>2) THEN
					REDIRECTED_OUTPUT=TRIM(inputstring(3:))
				ELSE
					skip_next=.TRUE.
					!use following command line argument as input.
					missing_argument=(file_counter==COMMAND_ARGUMENT_COUNT())
					IF (missing_argument) THEN
						CALL report_error(102)
					ELSE
						CALL GET_COMMAND_ARGUMENT(file_counter+1,REDIRECTED_OUTPUT)
						REDIRECTED_OUTPUT=TRIM(ADJUSTL(REDIRECTED_OUTPUT))
					ENDIF
				ENDIF
				valid_filename=.FALSE.
				IF (.NOT.((DISCONNECTED).OR.(missing_argument))) THEN
					DISCONNECTED=.TRUE.
					WRITE(*,*) " # REDIRECTING UNIT 6 TO '",TRIM(REDIRECTED_OUTPUT),"'"
					WRITE(*,*)
					OPEN(UNIT=6,FILE=TRIM(REDIRECTED_OUTPUT))
					WRITE(*,*) " # UNIT 6 REDIRECTED TO '",TRIM(REDIRECTED_OUTPUT),"'"
				ENDIF
				IF (missing_argument) command_line_used=.TRUE.
			CASE ("-r")
				IF (LEN(TRIM(inputstring))>2) THEN
					trajstring=TRIM(inputstring(3:))
				ELSE
					skip_next=.TRUE.
					!use following command line argument as input.
					missing_argument=(file_counter==COMMAND_ARGUMENT_COUNT())
					IF (missing_argument) THEN
						CALL report_error(102)
					ELSE
						CALL GET_COMMAND_ARGUMENT(file_counter+1,trajstring)
						trajstring=TRIM(ADJUSTL(trajstring))
					ENDIF
				ENDIF
				valid_filename=.FALSE.
				IF (.NOT.(missing_argument)) THEN
					WRITE(*,*) " # CALLING MOLECULE RECOGNITION MODULE"
					WRITE(*,*)
					CALL molecule_recognition(trajstring)
				ENDIF
				command_line_used=.TRUE.
			CASE DEFAULT
				!Check if there are some weird characters
				valid_filename=.TRUE.
				DO i=1,LEN(TRIM(inputstring)),1
					charnum=IACHAR(inputstring(i:i))
					IF (.NOT.(ANY(ALPHABET==charnum))) THEN
						CALL report_error(67,charnum)
						valid_filename=.FALSE.
						EXIT
					ENDIF
				ENDDO
			END SELECT
			IF (valid_filename) THEN
				GLOBAL_ITERATIONS=GLOBAL_ITERATIONS+1
				GENERAL_INPUT_FILENAMES(GLOBAL_ITERATIONS)=TRIM(inputstring)
			ENDIF
		ENDDO
		IF ((GLOBAL_ITERATIONS==0).AND.(.NOT.(command_line_used))) THEN
			CALL report_error(68)
			WRITE(*,*)
		ELSEIF (GLOBAL_ITERATIONS>1) THEN
			PRINT *," ** The following general input files will be used successively:"
			DO file_counter=1,GLOBAL_ITERATIONS,1
				INQUIRE(FILE=TRIM(GENERAL_INPUT_FILENAMES(file_counter)),EXIST=file_exists)
				IF (file_exists) THEN
					WRITE(inputstring,'("existing file")')
				ELSE
					WRITE(inputstring,'("nonexistent")')
				ENDIF
				PRINT *," **  - '",TRIM(GENERAL_INPUT_FILENAMES(file_counter)),"' (",TRIM(inputstring),")"
			ENDDO
		ENDIF
	ELSE
		GLOBAL_ITERATIONS=GLOBAL_ITERATIONS_DEFAULT
	ENDIF
END SUBROUTINE initialise_command_line_arguments

SUBROUTINE finalise_command_line_arguments()
USE SETTINGS
IMPLICIT NONE
INTEGER :: deallocstatus
	IF (COMMAND_ARGUMENT_COUNT()>0) THEN
		DEALLOCATE(GENERAL_INPUT_FILENAMES,STAT=deallocstatus)
		IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
	ENDIF
	IF (DISCONNECTED) CLOSE(UNIT=6)
END SUBROUTINE finalise_command_line_arguments

PROGRAM PREALPHA ! Copyright (C) !RELEASEYEAR! Frederik Philippi
USE SETTINGS
USE MOLECULAR
USE DEBUG
IMPLICIT NONE
INTEGER :: global_iteration_counter

INTERFACE

	SUBROUTINE run_analysis()
	IMPLICIT NONE
	END SUBROUTINE run_analysis

	SUBROUTINE initialise_global()
	IMPLICIT NONE
	END SUBROUTINE initialise_global

	SUBROUTINE finalise_global()
	IMPLICIT NONE
	END SUBROUTINE finalise_global

	SUBROUTINE initialise_command_line_arguments()
	IMPLICIT NONE
	END SUBROUTINE initialise_command_line_arguments

	SUBROUTINE finalise_command_line_arguments()
	IMPLICIT NONE
	END SUBROUTINE finalise_command_line_arguments
	
END INTERFACE
	!begin timing here
	CALL timing()
	CALL initialise_command_line_arguments()
	DO global_iteration_counter=1,GLOBAL_ITERATIONS,1
		IF (GLOBAL_ITERATIONS>1) WRITE(*,*) " ** current general input file is '",&
		&TRIM(GENERAL_INPUT_FILENAMES(global_iteration_counter)),"'"
		WRITE(*,*)
		IF (COMMAND_ARGUMENT_COUNT()>0) FILENAME_GENERAL_INPUT=TRIM(GENERAL_INPUT_FILENAMES(global_iteration_counter))
		!first, load all the necessary information and allocate memory for the trajectory.
		CALL initialise_global()
		!is it necessary to actually start the analysis?
		IF (SKIP_ANALYSIS) THEN
			WRITE(*,*) "Analysis is skipped."
		ELSE
			CALL timing()
			!Call the initialisation of the module MOLECULAR.
			CALL initialise_molecular()
			!then, read the trajectory from the lammps file.
			CALL load_trajectory()
			CALL timing()
			!Perform all the analyses requested in the body of general.inp
			CALL run_analysis()
			!deallocate, close open units
			CALL finalise_global()
		ENDIF
		CALL timing(total=.TRUE.)
		IF (GLOBAL_ITERATIONS>1) WRITE(*,*) " ** done with general input file '",&
		&TRIM(GENERAL_INPUT_FILENAMES(global_iteration_counter)),"'"
		IF ((GLOBAL_ITERATIONS>1).AND.(global_iteration_counter<GLOBAL_ITERATIONS))&
		&WRITE(*,*) " ** changing to next general input file"
		IF (GLOBAL_ITERATIONS>1) WRITE(*,*)
	ENDDO
	IF (give_error_count()==0) THEN
		Print *,"No warnings or errors encountered."
	ELSE
		WRITE(*,'(" Encountered ",I0," errors/warnings globally.")') give_error_count()
	ENDIF
	WRITE(*,*)
	CALL finalise_command_line_arguments()
END