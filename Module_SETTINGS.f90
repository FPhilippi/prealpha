
    ! prealpha - a tool to extract information from molecular dynamics trajectories.
    ! Copyright (C) !RELEASEYEAR! Frederik Philippi
    ! This work is funded by the Imperial President's PhD Scholarship.
	! This work is funded by the Japan Society for the Promotion of Science

    ! This program is free software: you can redistribute it and/or modify
    ! it under the terms of the GNU General Public License as published by
    ! the Free Software Foundation, either version 3 of the License, or
    ! (at your option) any later version.

    ! This program is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! GNU General Public License for more details.

    ! You should have received a copy of the GNU General Public License
    ! along with this program.  If not, see <https://www.gnu.org/licenses/>

! TO USE:
! Please use a compiler that supports OpenMP.
! For the GNU fortran compiler, set the corresponding flag:
! > gfortran pre-alpha.f03 -fopenmp
! The code is designed to also run without that library, but is not as nice.

!Length unit: Angström
!mass unit: Dalton
!Can only handle one-character element names (Extended to two characters now)
!consistent ordering assumed, in the format Elementname x y z
!Path names have to be enclosed in quotes. Probably even double quotes.
MODULE SETTINGS !This module contains important globals and subprograms.
	IMPLICIT NONE
	PUBLIC
	!default values - basic parameters
	INTEGER,PARAMETER :: SP=KIND(1.e0)
	INTEGER,PARAMETER :: DP=KIND(1.d0)
	INTEGER,PARAMETER :: WORKING_PRECISION=DP!'Standard' precision is double precision; change to single precision here if required.
	INTEGER,PARAMETER :: GENERAL_PRECISION=DP!'Standard' precision for positions and dihedrals.
	INTEGER,PARAMETER :: STORAGE_PRECISION=SP!'Standard' precision for the storage of positions in the trajectory.
	REAL(KIND=WORKING_PRECISION),PARAMETER :: degrees=57.295779513082320876798154814105d0 !constant: 360/2*Pi
	REAL(KIND=GENERAL_PRECISION),PARAMETER :: avogadro=6.02214076d23!avogadro's constant
	REAL(KIND=GENERAL_PRECISION),PARAMETER :: elementary_charge=1.602176634E-19!elementary_charge in Coulomb
	REAL(KIND=GENERAL_PRECISION),PARAMETER :: vacuum_permittivity=8.8541878128E-12 !vacuum permittivity in farads/metre
	REAL(KIND=GENERAL_PRECISION),PARAMETER :: boltzmann=1.380649E-23!boltzmann constant in joules per kelvin
	REAL(KIND=WORKING_PRECISION),PARAMETER :: pi=3.141592653589793238462643383279502884197169
	REAL(KIND=WORKING_PRECISION),PARAMETER :: twopi=2.0*pi
	!First, set default values
	LOGICAL,PARAMETER :: VERBOSE_OUTPUT_DEFAULT=.TRUE.
	LOGICAL,PARAMETER :: TIME_OUTPUT_DEFAULT=.TRUE.
	LOGICAL,PARAMETER :: DEVELOPERS_VERSION_DEFAULT=.TRUE.
	LOGICAL,PARAMETER :: ERROR_OUTPUT_DEFAULT=.TRUE.
	LOGICAL,PARAMETER :: PARALLEL_OPERATION_DEFAULT=.TRUE.
	LOGICAL,PARAMETER :: READ_SEQUENTIAL_DEFAULT=.FALSE.
	LOGICAL,PARAMETER :: BOX_VOLUME_GIVEN_DEFAULT=.FALSE.!If (T), then xlo,xhi,ylo, etc should be correct, too.
	INTEGER,PARAMETER :: ERROR_CODE_DEFAULT=-1
	INTEGER,PARAMETER :: TIME_SCALING_FACTOR_DEFAULT=1
	LOGICAL,PARAMETER :: WRAP_TRAJECTORY_DEFAULT=.FALSE.
	LOGICAL,PARAMETER :: UNWRAP_TRAJECTORY_DEFAULT=.FALSE.
	LOGICAL,PARAMETER :: EXTRA_VELOCITY_DEFAULT=.FALSE.
	INTEGER,PARAMETER :: HEADER_LINES_DEFAULT=5
	INTEGER,PARAMETER :: MAXITERATIONS=500
	INTEGER,PARAMETER :: GLOBAL_ITERATIONS_DEFAULT=1
	REAL,PARAMETER :: CUTOFF_INTERMOLECULAR_DEFAULT=2.0
	REAL,PARAMETER :: VDW_RATIO_INTERMOLECULAR_DEFAULT=2.8 !good values usually 2.2-3.4
	INTEGER :: q !nobody uses 'q' anywhere else, hopefully.
	INTEGER :: ALPHABET(26*2+3+10)=(/ (q,q=IACHAR("a"),IACHAR("a")+25,1),(q,q=IACHAR("A"),IACHAR("A")+25,1),IACHAR("_"),&
	&IACHAR("/"),IACHAR("."),(q,q=IACHAR("0"),IACHAR("0")+9,1) /)!The 'benign' alphabet, re-initialised in 'initialise_global'
	INTEGER :: ALPHABET_small(26)!only lowercase letters, i.e. a,b,c...
	CHARACTER(LEN=3),PARAMETER  :: TRAJECTORY_TYPE_DEFAULT="lmp"
	CHARACTER(LEN=*),PARAMETER :: FILENAME_GENERAL_INPUT_DEFAULT="general.inp"
	!variables
	LOGICAL :: VERBOSE_OUTPUT=VERBOSE_OUTPUT_DEFAULT!Controls how detailed the output is - useful for debugging.
	LOGICAL :: TIME_OUTPUT=TIME_OUTPUT_DEFAULT!Give timing information
	LOGICAL :: DEVELOPERS_VERSION=DEVELOPERS_VERSION_DEFAULT!Turns on suspicious stuff
	LOGICAL :: ERROR_OUTPUT=ERROR_OUTPUT_DEFAULT!Report any encountered Errors
	LOGICAL :: PARALLEL_OPERATION=PARALLEL_OPERATION_DEFAULT!Operation mode set to parallel
	LOGICAL :: READ_SEQUENTIAL=READ_SEQUENTIAL_DEFAULT!reading the trajectory in a serial way rather than everything at once.
	LOGICAL :: BOX_VOLUME_GIVEN=BOX_VOLUME_GIVEN_DEFAULT!is there a box volume available?
	LOGICAL :: WRAP_TRAJECTORY=WRAP_TRAJECTORY_DEFAULT!Wrap the trajectory?
	LOGICAL :: UNWRAP_TRAJECTORY=UNWRAP_TRAJECTORY_DEFAULT!manually UNWrap the trajectory?
	LOGICAL :: EXTRA_VELOCITY=EXTRA_VELOCITY_DEFAULT!are there extra 3 columns to read? (i.e. "E vx vy vz xu yu zu", or vice versa)
	LOGICAL :: SKIP_ANALYSIS!don't do the actual analysis...
	LOGICAL :: USER_INPUT=.FALSE.!Turns on as soon as user input started...
	LOGICAL :: DISCONNECTED=.FALSE. !If true, then the standard output is redirected into 'output.dat' (or REDIRECTED_OUTPUT)
	INTEGER :: ERROR_CODE=ERROR_CODE_DEFAULT!latest error is stored in here.
	INTEGER :: TIME_SCALING_FACTOR=TIME_SCALING_FACTOR_DEFAULT !integer value to scale the timelines.
	INTEGER :: HEADER_LINES=HEADER_LINES_DEFAULT!number of fixed lines in general.inp, without the optional 'sequential read'.
	INTEGER :: GLOBAL_ITERATIONS=GLOBAL_ITERATIONS_DEFAULT!number of time the whole program (including global initialisation and finalisation) is repeated.
	INTEGER :: error_count=0 !number of warnings and errors encountered
	REAL :: CUTOFF_INTERMOLECULAR=CUTOFF_INTERMOLECULAR_DEFAULT ! Everything bigger than that is considered intermolecular
	REAL :: VDW_RATIO_INTERMOLECULAR=VDW_RATIO_INTERMOLECULAR_DEFAULT ! same as the cutoff, just defined using the covalence radii
	CHARACTER(LEN=1024),DIMENSION(:),ALLOCATABLE :: GENERAL_INPUT_FILENAMES !the general input filenames passed from the command line.
	CHARACTER(LEN=3) :: TRAJECTORY_TYPE=TRAJECTORY_TYPE_DEFAULT!type of the trajectory, e.g. lmp or xyz
	CHARACTER(LEN=1024) :: FILENAME_TRAJECTORY,PATH_TRAJECTORY,PATH_INPUT,PATH_OUTPUT
	CHARACTER(LEN=1024) :: FILENAME_GENERAL_INPUT=FILENAME_GENERAL_INPUT_DEFAULT
	CHARACTER(LEN=1024) :: FILENAME_MOLECULAR_INPUT="molecular.inp"
	CHARACTER(LEN=1024) :: FILENAME_AUTOCORRELATION_INPUT="autocorrelation.inp"
	CHARACTER(LEN=1024) :: FILENAME_DIFFUSION_INPUT="diffusion.inp"
	CHARACTER(LEN=1024) :: FILENAME_DISTRIBUTION_INPUT="distribution.inp"
	CHARACTER(LEN=1024) :: FILENAME_DISTANCE_INPUT="distance.inp"
	CHARACTER(LEN=1024) :: FILENAME_SPECIATION_INPUT="speciation.inp"
	CHARACTER(LEN=1024) :: FILENAME_DISPERSION_INPUT="dispersion.inp"
	CHARACTER(LEN=1024) :: REDIRECTED_OUTPUT="output.dat"
	CHARACTER(LEN=1024) :: OUTPUT_PREFIX="" !prefix added to output, for example to distinguish between different autocorrelation analyses
	CHARACTER(LEN=3) :: INFORMATION_IN_TRAJECTORY="UNK"!does the trajectory contain velocities (VEL) or coordinates (POS)????
	!dealing with references
	TYPE,PRIVATE :: reference_entry
        CHARACTER(LEN=1024) :: reference_name
		CHARACTER(LEN=1024) :: reference_DOI
		INTEGER :: reference_number
		LOGICAL :: cited !TRUE if cited
    END TYPE reference_entry
	TYPE(reference_entry) :: LIST_OF_REFERENCES(6)
	!LIST OF ERRORS HANDLED BY THE ROUTINES:
	!0 unspecified error. These errors should (in theory) never be encountered.
	!1 divided by zero in normalize3D
	!2 tried to calculated arcuscosinus of a number <-1 or >1
	!3 divided by zero in normalize2D
	!4 Unknown Element
	!5 general.inp is not correctly formatted
	!6 couldn't allocate memory in initialise_molecular
	!7 molecular.inp is not correctly formatted
	!8 couldn't deallocate memory in finalise_molecular
	!9 couldn't find trajectory file
	!10 mismatch in atom number
	!11 couldn't allocate memory in initialise_dihedrals or initialise_autocorrelation
	!12 User requested overwriting dihedral_member_indices by repeatedly calling initialise_dihedrals.
	!13 mismatch in array size
	!14 autocorrelation.inp is not correctly formatted
	!15 couldn't deallocate memory in finalise_autocorrelation
	!16 couldn't allocate memory in initialise_autocorrelation
	!17 molecular.inp doesn't exist
	!18 streaming general.inp not successful.
	!19 positive exit status while reading unit 7 (general.inp)
	!20 unknown input keyword in general.inp
	!21 autocorrelation.inp is not available
	!22 problem with memory allocation, unspecified.
	!23 probably with memory deallocation, unspecified.
	!24 problem when streaming autocorrelation input.
	!25 array bounds have been exceeded during binning procedure.
	!26 problem writing output, unspecified.
	!27 unit already connected
	!28 tmax is larger than the number of timesteps available. (or smaller than 100)
	!29 Could not center molecules
	!30 diffusion.inp is not correctly formatted
	!31 diffusion.inp is not available
	!32 problem when streaming diffusion input.
	!33 invalid molecule_type_index
	!34 diffusion.inp is not correctly formatted
	!35 tstep is too large, has been reset to maximum permissible value
	!36 tmax is very small
	!37 trajectory is too short for the diffusion or distribution analysis.
	!38 couldn't open trajectory file
	!39 cannot invoke rmm-vcf with only one molecule type.
	!40 remind user to provide velocities.
	!41 Box volume not available
	!42 Invalid number of threads (parallelisation issue)
	!43 Bad integer in user input
	!44 Bad Logical in user input
	!45 incompatible analyses requested!
	!46 Couldn't write input file for some reason.
	!47 nonstandard character
	!48 filename contains blanks
	!49 input string too long? might have been truncated
	!50 invalid user input, unspecified
	!51 unknown trajectory type
	!52 box is not charge neutral
	!53 check format of dump statement: couldn't find expression 'element'
	!54 could neither find 'xu yu zu' nor 'vx vy vz' in the trajectory header.
	!55 Be careful - unknown trajectory format!
	!56 mismatch between information stored in trajectory and requested analysis.
	!57 requested step out of bounds (less than 1, more than number of steps)
	!58 overwriting mass while reading molecular input file
	!59 negative mass has been specified
	!60 compromised format of masses section in molecular input file.
	!61 Failed reading the custom default mass input - masses set to zero.
	!62 particle with vanishing mass detected.
	!63 Failed reading the custom constraints input - constraints removed.
	!64 compromised format of constraints section in molecular input file.
	!65 overwriting constraint while reading molecular input file
	!66 could not allocate memory to store filenames. Abort immediately!
	!67 nonstandard character not accepted.
	!68 No useful filenames provided. reset to default
	!69 invalid molecule_index
	!70 element name mismatch while checking the first step.
	!71 couldn't read trajectory lines...
	!72 wrapping requested where it wouldn't make sense.
	!73 some genious requested wrapping velocities.
	!74 constraints not physical
	!75 more molecules specified to export than available.
	!76 molecule index has been specified twice / more than once
	!77 legendre polynomial not available
	!78 User requested overwriting fragment_list by repeatedly calling initialise_fragments.
	!79 couldn't allocate memory for fragment_list
	!80 atom_index out of bounds...
	!81 couldn't read base or tip atom input for fragments
	!82 invalid number of fragments
	!83 Cannot initialise fragments - not enough information.
	!84 Cannot read head of trajectory. (recoverable error)
	!85 Cannot read head of trajectory. (occurred while reading a snapshot - not recoverable)
	!86 lonely drude particle found
	!87 Failed reading the drude particles from molecular input file.
	!88 not all drude particles have been assigned.
	!89 using sequential read
	!90 t-value not available
	!91 this feature requires assigned drude particles
	!92 box volume will be overwritten
	!93 boundaries not sensible.
	!94 Error count exceeds maximum
	!95 Cannot perform molecule recognition
	!96 Atoms for one of the molecular units were separated.
	!97 need at least two different molecules for dimers
	!98 distribution.inp is not available
	!99 distribution.inp is not correctly formatted
	!100 problem when streaming distribution input.
	!101 integer overflow while binning
	!102 Not enough arguments in command line
	!103 Not enough steps in trajectory for jump analysis
	!104 bin count set to nearest sensible value
	!105 endstep<startstep, switch the two
	!106 invalid TIME_SCALING_FACTOR
	!107 Cannot read ITEM: TIMESTEP
	!108 inconsistency in ITEM: TIMESTEP detected
	!109 Failed reading the custom atomic masses input. will try to use element names where available.
	!110 invalid molecule_type_index, ignore line
	!111 invalid atom_index, ignore line
	!112 atomic mass already specified, will be overwritten.
	!113 "simple" keyword not available
	!114 "prealpha_simple.inp" will be overwritten
	!115 couldn't allocate memory for vc_components
	!116 vc_components contains invalid molecule type
	!117 vc_components contains invalid "self/distinct" flag
	!118 component number < 1
	!119 molecule_index out of bounds...
	!120 no charged particles - cannot calculate conductivity.
	!121 number of dihedral conditions less than 1
	!122 default charges only for defined elements, not lowercase letters.
	!123 Could not read default charges
	!124 Could not read atomic charges
	!125 atomic charge already specified, will be overwritten.
	!126 mismatch between sum of atomic charges and charge of whole molecule.
	!127 something is wrong with the atomic charges of the specified molecule type index
	!128 custom_atom_charges and custom_default_charges are both FALSE
	!129 jump time in dihedral autocorrelation has been specified twice / more than once
	!130 requested number of jump analyses exceeded permissible maximum.
	!131 jump length is a tiny bit too short
	!132 distribution.inp is not available.
	!133 Same molecule type for cross interactions
	!134 distance.inp is not available
	!135 distance.inp is not correctly formatted
	!136 the molecule (cf molecule_type_index) does not have an element with this name
	!137 Couldn't allocate memory for distance list
	!138 mismatch in array size (atomic indices)
	!139 no valid subjobs in distance module
	!140 distance module: swapped "wildcard -> specific" to "specific -> wildcard"
	!141 not enough molecules (molecule_index)
	!142 problem when streaming distance input.
	!143 no charged particles - centre of charge trajectory will be empty.
	!144 the specified trajectory output or input format is not supported (.gro, .xyz, .lmp)
	!145 print notice about nm convention in gromacs
	! DISPERSION MODULE CURRENTLY NOT FURTHER DEVELOPED
	!146 dispersion.inp is not available
	!147 dispersion.inp is not correctly formatted
	!148 no valid subjobs in dispersion module
	!149 problem when streaming dispersion input.
	!150 Couldn't allocate memory for dispersion list
	!151 requested both wrapped and unwrapped trajectory
	!152 requested unwrapping with sequential read
	!153 requested change of exponent with alpha2. But we won't let them.
	!154 non-unit vector requested as projection for alpha2.
	!155 speciation.inp is not correctly formatted
	!156 Couldn't allocate memory for speciation list
	!157 speciation.inp is not available
	!158 problem when streaming speciation input.
	!159 no valid donors or acceptors left
	!160 neighbour overflow in speciation module
	!161 species overflow in speciation module
	!162 overwriting existing atom groups in speciation module
	!163 existing element groups are kept / overwritten
	!164 no group defined for this atom
	!165 atom assigned twice to a group
	!166 need time_series for autocorrelation
	!167 could not read the time_series file - return.
	!PRIVATE/PUBLIC declarations
	PUBLIC :: normalize2D,normalize3D,crossproduct,report_error,timing_parallel_sections,legendre_polynomial
	PUBLIC :: FILENAME_TRAJECTORY,PATH_TRAJECTORY,PATH_INPUT,PATH_OUTPUT,user_friendly_time_output
	PUBLIC :: user_input_string,user_input_integer,user_input_logical,user_input_real
	PUBLIC :: student_t_value,covalence_radius,give_error_count,DEVELOPERS_VERSION,initialise_references,add_reference,print_references
	PRIVATE :: q !that's not a typo!
	PRIVATE :: error_count
	CONTAINS

		SUBROUTINE initialise_references()
		IMPLICIT NONE
			LIST_OF_REFERENCES(:)%reference_number=0
			LIST_OF_REFERENCES(:)%cited=.FALSE.
			!Here, add the references. need to manually increase the max number for LIST_OF_REFERENCES in the declaration.
			!Daniel's curled cation paper
			LIST_OF_REFERENCES(1)%reference_name="Phys. Chem. Chem. Phys., 2021, 23, 21042–21064."
			LIST_OF_REFERENCES(1)%reference_DOI="10.1039/D1CP02889H"
			!The CTPOL paper
			LIST_OF_REFERENCES(2)%reference_name="Phys. Chem. Chem. Phys., 2022, 24, 3144–3162."
			LIST_OF_REFERENCES(2)%reference_DOI="10.1039/D1CP04592J"
			!yongji's paper
			LIST_OF_REFERENCES(3)%reference_name="J. Chem. Phys., 2022, 156, 204312."
			LIST_OF_REFERENCES(3)%reference_DOI="10.1063/5.0091322"
			!The fluorination/flexibility paper
			LIST_OF_REFERENCES(4)%reference_name="Chem. Sci., 2022, 13, 9176–9190."
			LIST_OF_REFERENCES(4)%reference_DOI="10.1039/d2sc03074h"
			!Julian's FFC paper
			LIST_OF_REFERENCES(5)%reference_name="J. Phys. Chem. B., 2022, 126, 7143–7158."
			LIST_OF_REFERENCES(5)%reference_DOI="10.1021/acs.jpcb.2c01372"
			!The battery paper
			LIST_OF_REFERENCES(6)%reference_name="Chem. Sci., 2024, 15, 7342–7358."
			LIST_OF_REFERENCES(6)%reference_DOI="10.1039/D4SC01492H"
		END SUBROUTINE initialise_references

		SUBROUTINE add_reference(reference_number_to_add)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: reference_number_to_add
			IF (reference_number_to_add>SIZE(LIST_OF_REFERENCES)) THEN
				CALL report_error(0)
				RETURN
			ENDIF
			IF (LIST_OF_REFERENCES(reference_number_to_add)%cited) THEN
				WRITE(*,'(" See also reference [",I0,"].")')&
				&LIST_OF_REFERENCES(reference_number_to_add)%reference_number
			ELSE
				LIST_OF_REFERENCES(reference_number_to_add)%reference_number=&
				&MAXVAL(LIST_OF_REFERENCES(:)%reference_number)+1
				LIST_OF_REFERENCES(reference_number_to_add)%cited=.TRUE.
				WRITE(*,'(" Please cite reference [",I0,"].")')&
				&LIST_OF_REFERENCES(reference_number_to_add)%reference_number
			ENDIF
		END SUBROUTINE add_reference

		SUBROUTINE print_references()
		IMPLICIT NONE
		INTEGER :: reference_counter,current_reference
			IF (ALL(LIST_OF_REFERENCES(:)%cited.eqv..FALSE.)) THEN
				WRITE(*,'(" Please consider citing our first paper with prealpha:")')
				WRITE(*,'("   [1] ",A)') &
				&TRIM(LIST_OF_REFERENCES(1)%reference_name)
				WRITE(*,'("       dx.doi.org/",A)') &
				&TRIM(LIST_OF_REFERENCES(1)%reference_DOI)
			ELSE
				WRITE(*,'(" References:")')
				current_reference=0
				DO
					current_reference=current_reference+1
					DO reference_counter=1,SIZE(LIST_OF_REFERENCES),1
						IF (LIST_OF_REFERENCES(reference_counter)%reference_number==current_reference) THEN
							WRITE(*,'("   [",I0,"] ",A)') current_reference,&
							&TRIM(LIST_OF_REFERENCES(reference_counter)%reference_name)
							WRITE(*,'("       dx.doi.org/",A)') &
							&TRIM(LIST_OF_REFERENCES(reference_counter)%reference_DOI)
							LIST_OF_REFERENCES(reference_counter)%cited=.FALSE.
						ENDIF
					ENDDO
					IF (ALL(LIST_OF_REFERENCES(:)%cited.eqv..FALSE.)) EXIT
				ENDDO
			ENDIF
		END SUBROUTINE print_references

		SUBROUTINE report_error(local_error_code,exit_status)!Routine for error handling. Severe Errors cause the program to stop.
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: local_error_code
		INTEGER,INTENT(IN),OPTIONAL :: exit_status
			error_count=error_count+1
			ERROR_CODE=local_error_code
			IF ((VERBOSE_OUTPUT).AND.(PRESENT(exit_status))) WRITE(*,'("  #  EXIT STATUS ",I0," was passed to report_error")') exit_status
			IF (ERROR_OUTPUT) THEN
				SELECT CASE (local_error_code)
				CASE (0)
					WRITE(*,*) " #  ERROR 0. Something internal went wrong that really shouldn't."
					WRITE(*,*) " #  Results might be biased! Please report this error."
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (1)
					WRITE(*,*) " #  ERROR 1 in normalize2D: dividing by zero"
					WRITE(*,*) " #  Results might be biased!"
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (2)
					WRITE(*,*) " #  ERROR 2 in dihedral_angle: arccos(x) with x<-1 or x>1 beyond tolerance"
					WRITE(*,*) "--> set to nearest sensible value"
				CASE (3)
					WRITE(*,*) " #  ERROR 3 in normalize3D: dividing by zero"
					WRITE(*,*) " #  Results might be biased!"
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (4)
					WRITE(*,*) " #  SEVERE ERROR 4 in atomic_weight/atomic_charge: unknown element"
					WRITE(*,*) " #  If necessary, add element to function atomic_weight and atomic_charge in module MOLECULAR"
					CALL finalise_global()
					STOP
				CASE (5)
					WRITE(*,*) " #  SEVERE ERROR 5: couldn't read '",TRIM(FILENAME_GENERAL_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CALL finalise_global()
					STOP
				CASE (6)
					WRITE(*,*) " #  SEVERE ERROR 6: couldn't allocate memory during initialisation"
					WRITE(*,*) " #  give more RAM!"
					CALL finalise_global()
					STOP
				CASE (7)
					WRITE(*,*) " #  SEVERE ERROR 7: couldn't read '",TRIM(FILENAME_MOLECULAR_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CALL finalise_global()
					STOP
				CASE (8)
					WRITE(*,*) " #  ERROR 8: couldn't deallocate memory during finalisation (module MOLECULAR)"
					WRITE(*,*) "--> Program will continue"
				CASE (9)
					WRITE(*,*) " #  SEVERE ERROR 9: couldn't find trajectory file"
					WRITE(*,*) " #  check path and filename!"
					CALL finalise_global()
					STOP
				CASE (10)
					WRITE(*,*) " #  ERROR 10: mismatch in atom number between input and lammps trajectory"
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (11)
					WRITE(*,*) " #  SEVERE ERROR 11: couldn't allocate memory for list of dihedrals."
					WRITE(*,*) " #  give more RAM!"
					CALL finalise_global()
					STOP
				CASE (12)
					IF (PRESENT(exit_status)) THEN !If an exit_status is present, then the error has been called after the deallocation, which is bad.
						WRITE(*,*) "SEVERE ERROR 12: deallocation during reinitialisation unsuccessful."
						CALL finalise_global()
						STOP
					ENDIF
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 12: dihedral_member_indices is already initialised."
					WRITE(*,*) "--> Program will try to reinitialise."
				CASE (13)
					WRITE(*,*) " #  ERROR 13: array of wrong size was passed to give_dihedrals."
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (14)
					WRITE(*,*) " #  SEVERE ERROR 14: couldn't read '",TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the autocorrelation input file
					CALL finalise_global()
					STOP
				CASE (15)
					WRITE(*,*) " #  ERROR 15: couldn't deallocate memory during finalisation (module AUTOCORRELATION)"
					WRITE(*,*) "--> Program will continue"
				CASE (16)
					WRITE(*,*) " #  SEVERE ERROR 16: couldn't allocate memory for autocorrelation array."
					WRITE(*,*) " #  give more RAM!"
					CALL finalise_global()
					STOP
				CASE (17)
					WRITE(*,*) " #  SEVERE ERROR 17: The specified molecular input file doesn't exist."
					IF (.NOT.(PRESENT(exit_status))) STOP
				CASE (18)
					WRITE(*,*) " #  SEVERE ERROR 18: problem streaming '",TRIM(FILENAME_GENERAL_INPUT),"'"
					CALL finalise_global()
					STOP
				CASE (19)
					WRITE(*,*) " #  ERROR 19: problem streaming '",TRIM(FILENAME_GENERAL_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_GENERAL_INPUT),"'!"
				CASE (20)
					WRITE(*,*) " #  WARNING 20: Couldn't interpret the specified line of ",TRIM(FILENAME_GENERAL_INPUT)
				CASE (21)
					WRITE(*,*) " #  ERROR 21: couldn't find '",TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (22)
					WRITE(*,*) " #  SEVERE ERROR 22: couldn't allocate memory."
					WRITE(*,*) " #  give more RAM? Could also be a genuine issue - please report."
					CALL finalise_global()
					STOP
				CASE (23)
					WRITE(*,*) " #  ERROR 23: couldn't deallocate memory."
					WRITE(*,*) "--> Program will continue. Could also be a genuine issue - please report."
				CASE (24)
					WRITE(*,*) " #  ERROR 24: problem streaming '",TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_AUTOCORRELATION_INPUT),"'!"
				CASE (25)
					WRITE(*,*) " #  ERROR 25: array bounds have been exceeded (see 'EXIT STATUS') during binning."
					WRITE(*,*) "--> Program will continue, and most likely crash. Please report this error."
				CASE (26)
					WRITE(*,*) " #  ERROR 26: unspecified problem while writing data."
					WRITE(*,*) "--> Program will continue."
				CASE (27)
					WRITE(*,*) " #  ERROR 27: The unit which is specified as 'EXIT STATUS' is already open."
					WRITE(*,*) "--> Program will continue, and most likely crash. Please report this error."
				CASE (28)
					WRITE(*,*) " #  WARNING 28: tmax is too large (or too small)."
					WRITE(*,*) "--> Program will continue, tmax is set to its maximum (see 'EXIT STATUS')."
				CASE (29)
					WRITE(*,*) " #  WARNING 29: couldn't center molecules."
					WRITE(*,*) "--> Program will continue."
				CASE (30)
					WRITE(*,*) " #  SEVERE ERROR 30: couldn't read '",TRIM(FILENAME_DIFFUSION_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the diffusion input file
					CALL finalise_global()
					STOP
				CASE (31)
					WRITE(*,*) " #  ERROR 31: couldn't find '",TRIM(FILENAME_DIFFUSION_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (32)
					WRITE(*,*) " #  ERROR 32: problem streaming '",TRIM(FILENAME_DIFFUSION_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_DIFFUSION_INPUT),"'!"
				CASE (33)
					WRITE(*,*) " #  ERROR 33: the specified molecule type index doesn't exist."
					WRITE(*,*) "--> Main program will continue, but this analysis is aborted."
				CASE (34)
					WRITE(*,*) " #  SEVERE ERROR 34: couldn't read '",TRIM(FILENAME_DIFFUSION_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the diffusion input file
					CALL finalise_global()
					STOP
				CASE (35)
					WRITE(*,*) " #  WARNING 35: tstep is too large."
					WRITE(*,*) "--> Program will continue, but tstep is set to its maximum (see 'EXIT STATUS')."
				CASE (36)
					WRITE(*,*) " #  WARNING 36: tmax is smaller than 10."
					WRITE(*,*) "--> Program will continue, tmax is set to its minimum (see 'EXIT STATUS')."
				CASE (37)
					WRITE(*,*) " #  ERROR 37: trajectory is too short for this module (needs > 10 steps)."
					WRITE(*,*) "--> Main program will continue, diffusion/distribution analysis is aborted."
				CASE (38)
					WRITE(*,*) "SEVERE ERROR 38: Trajectory file couldn't be opened."
					CALL finalise_global()
					STOP
				CASE (39)
					WRITE(*,*) " #  ERROR 39: need at least two molecule types for *CROSS*-correlation."
					WRITE(*,*) "--> Main program will continue, VCF analysis is aborted."
				CASE (40)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 40: This module requires *velocities* as input instead of cartesian coordinates."
				CASE (41)
					WRITE(*,*) " #  ERROR 41: Box volume required, but not available."
					WRITE(*,*) "--> Try setting the box volume using the keyword 'cubic_box_edge'."
				CASE (42)
					WRITE(*,*) " #  WARNING 42: Invalid number of threads requested (see 'EXIT STATUS')."
					WRITE(*,*) "--> Program will continue, number_of_threads is set to its maximum."
				CASE (43)
					WRITE(*,*) " #  WARNING 43: Invalid user input. Please enter a valid number."
					WRITE(*,*) "--> Program will continue."
				CASE (44)
					WRITE(*,*) " #  WARNING 44: Invalid user input. Please type 'yes'/'y' or 'no'/'n'."
					WRITE(*,*) "--> Program will continue."
				CASE (45)
					WRITE(*,*) " #  REALLY SERIOUS WARNING 45: Two incompatible analyses were requested."
					WRITE(*,*) " #  AT LEAST one of the results will not only be biased, but just wrong."
					WRITE(*,*) " #  Carefully think about your analyses: do you need coordinates or velocities?"
					WRITE(*,*) "--> Program will continue anyway."
				CASE (46)
					WRITE(*,*) " #  ERROR 46: Couldn't write input file."
					WRITE(*,*) "--> Check for sensible filename, blanks, etc."
				CASE (47)
					IF (PRESENT(exit_status)) THEN !If an exit_status is present, then the character has been passed.
						WRITE(*,*) " #  WARNING 47: String contains nonstandard character '",CHAR(exit_status),"'."
					ELSE
						WRITE(*,*) " #  WARNING 47: String contains nonstandard character."
					ENDIF
					WRITE(*,*) "--> Program will try to continue, might crash later."
				CASE (48)
					WRITE(*,*) " #  WARNING 48: String contains blanks."
					WRITE(*,*) "--> Program will try to continue, might crash later."
				CASE (49)
					WRITE(*,*) " #  WARNING 49: Input string length matches requested size."
					WRITE(*,*) " #  Check output: string might have been truncated."
				CASE (50)
					WRITE(*,*) " #  WARNING 50: Invalid user input. Please type again."
					WRITE(*,*) "--> Program will continue."
				CASE (51)
					WRITE(*,*) " #  ERROR 51: Trajectory type unknown or not supported."
					WRITE(*,*) "--> TRAJECTORY_TYPE set to default value (",TRAJECTORY_TYPE_DEFAULT,")"
					TRAJECTORY_TYPE=TRAJECTORY_TYPE_DEFAULT
				CASE (52)
					WRITE(*,*) " #  WARNING 52: The specified system is NOT charge neutral."
					IF (PRESENT(exit_status)) WRITE(*,*) " #  check total charge (given in EXIT STATUS)"!If an exit_status is present, then the total charge has been passed.
				CASE (53)
					WRITE(*,*) " #  WARNING 53: Could not find 'element' in the trajectory."
					WRITE(*,*) " #  Please assure that your trajectory has the correct format, e.g. in lammps:"
					WRITE(*,*) "--> use 'dump TRAJ all custom 1 trajectory.lmp element xu yu zu'"
					WRITE(*,*) " #  (for velocities / VACF analysis, use '...element vx vy vz' instead)"
				CASE (54)
					WRITE(*,*) " #  WARNING 54: Your trajectory might have the wrong format."
					WRITE(*,*) " #  Please assure that your trajectory has the correct format, e.g. in lammps:"
					WRITE(*,*) "--> use 'dump TRAJ all custom 1 trajectory.lmp element xu yu zu'"
					WRITE(*,*) " #  (for velocities / VACF analysis, use '...element vx vy vz' instead)"
				CASE (55)
					WRITE(*,*) " #  WARNING 55: Unknown trajectory format! Please carefully check results and input."
					WRITE(*,*) " #  (Velocity or Cartesian coordinates)"
				CASE (56)
					WRITE(*,*) " #  SERIOUS WARNING 56: This analysis/feature is not meaningful with the information in the trajectory (",&
					&INFORMATION_IN_TRAJECTORY,")."
					WRITE(*,*) "--> Program will continue. Results might be useless."
				CASE (57)
					WRITE(*,*) " #  WARNING 57: The requested timestep (see EXIT STATUS) is out of range."
					WRITE(*,*) "--> Timestep is set to nearest sensible value."
				CASE (58)
					WRITE(*,*) " #  WARNING 58: custom mass for '",CHAR(exit_status),"' has been overwritten."
					WRITE(*,*) " #  If masses are specified twice, then only the last one will be used!"
				CASE (59)
					WRITE(*,*) " #  WARNING 59: negative mass specified for '",CHAR(exit_status),"'."
				CASE (60)
					WRITE(*,*) " #  WARNING 60: The character '",CHAR(exit_status),"' has been ignored."
					WRITE(*,*) " #  Only lowercase letters (a,b,c,...,z) and element names including 'X' (for drudes) are allowed."
				CASE (61)
					WRITE(*,*) " #  ERROR 61: Couldn't read 'default_masses' section from molecular input file."
					WRITE(*,*) "--> Check format of molecular input file! COM masses have been set to zero."
				CASE (62)
					WRITE(*,*) " #  WARNING 62: A massless particle has been specified."
					WRITE(*,*) " #  Quantities like the kinetic temperature will be biased!"
				CASE (63)
					WRITE(*,*) " #  ERROR 63: Couldn't read 'constraints' section from molecular input file."
					WRITE(*,*) "--> Check format of molecular input file! All constraints have been removed."
				CASE (64)
					WRITE(*,*) " #  WARNING 64: ignored incorrectly formatted constraints section line (see EXIT STATUS)"
				CASE (65)
					WRITE(*,*) " #  WARNING 65: constraints have been overwritten! (see molecule type in EXIT STATUS)."
					WRITE(*,*) " #  If constraints are specified twice, then only the last one will be used!"
				CASE (66)
					WRITE(*,*) " #  SEVERE ERROR 66: couldn't allocate memory for filenames from command line!"
					WRITE(*,*) " #  Program will stop immediately. Please report this issue."
					STOP
				CASE (67)
					IF (PRESENT(exit_status)) THEN !If an exit_status is present, then the character has been passed.
						WRITE(*,*) " #  ERROR 67: String with nonstandard character '",CHAR(exit_status),"' not accepted."
					ELSE
						WRITE(*,*) " #  WARNING 67: String contains nonstandard character."
					ENDIF
				CASE (68)
					WRITE(*,*) " #  ERROR 68: No valid inputfiles. Allowed characters are: "
					WRITE(*,*) " #  ",CHAR(ALPHABET(:))
					FILENAME_GENERAL_INPUT=FILENAME_GENERAL_INPUT_DEFAULT
					WRITE(*,*) "--> analysis skipped."
				CASE (69)
					WRITE(*,*) " #  ERROR 69: the specified molecule index doesn't exist."
					WRITE(*,*) "--> Main program will continue, but this analysis is aborted."
				CASE (70)
					WRITE(*,*) " #  SERIOUS WARNING 70: A mismatch of element names in the trajectory has been detected!"
					WRITE(*,*) " #  (molecule_type_index given in EXIT STATUS)"
					WRITE(*,*) "--> Carefully check your trajectory and your molecular input file!"
				CASE (71)
					WRITE(*,*) " #  SEVERE ERROR 71: FileIO error while reading the trajectory."
					WRITE(*,*) "--> Check format, path, filename etc of your trajectory."
					CALL finalise_global()
					STOP
				CASE (72)
					WRITE(*,*) " #  ERROR 72: This feature is not available with a wrapped trajectory."
					WRITE(*,*) "--> Main program will continue, but this analysis is aborted."
				CASE (73)
					WRITE(*,*) " #  ERROR 73: Trajectory contains velocities - (UN)wrapping not available."
					WRAP_TRAJECTORY=.FALSE.
					UNWRAP_TRAJECTORY=.FALSE.
				CASE (74)
					WRITE(*,*) " #  ERROR 74: unphysical number of constraints."
					WRITE(*,*) " #  Temperature values don't include the constraints correction!"
				CASE (75)
					WRITE(*,*) " #  ERROR 75: specified molecules for export exceeds total number."
					WRITE(*,*) " #  Ignoring keyword 'export' for this molecule!"
				CASE (76)
					WRITE(*,*) " #  WARNING 76: redundant molecule indices specified by 'export'."
				CASE (77)
					WRITE(*,*) " #  ERROR 77: Legendre polynomials of this order are not available."
				CASE (78)
					IF (PRESENT(exit_status)) THEN !If an exit_status is present, then the error has been called after the deallocation, which is bad.
						WRITE(*,*) "SEVERE ERROR 78: deallocation during reinitialisation unsuccessful."
						CALL finalise_global()
						STOP
					ENDIF
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 78: fragment_list is already initialised."
					WRITE(*,*) "--> Program will try to reinitialise."
				CASE (79)
					WRITE(*,*) " #  SEVERE ERROR 79: couldn't allocate memory for fragment_list."
					WRITE(*,*) " #  give more RAM!"
					CALL finalise_global()
					STOP
				CASE (80)
					WRITE(*,*) " #  ERROR 80: specified atom index (see EXIT STATUS) is out of bounds."
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (81)
					WRITE(*,*) " #  ERROR 81: cannot read fragment atoms (base/tip) from input file."
					WRITE(*,*) "--> check format, including fragment record"
					WRITE(*,*) "--> check line (given as EXIT STATUS)"
				CASE (82)
					WRITE(*,*) " #  ERROR 82: invalid number of fragments (see EXIT STATUS)."
				CASE (83)
					WRITE(*,*) " #  ERROR 83: Cannot initialise fragments - not enough information."
					WRITE(*,*) "--> assure that both tip and base fragments have been specified in the correct format."
					WRITE(*,*) "--> Main program will continue, reorientational analysis is aborted."
				CASE (84)
					WRITE(*,*) " #  ERROR 84: Cannot read trajectory head. Check format of trajectory and molecular input file."
					IF (PRESENT(exit_status)) THEN
						IF (exit_status<0) WRITE(*,*) " #  end of file encountered!"
					ENDIF
				CASE (85)
					WRITE(*,*) " #  SEVERE ERROR 85: Cannot read trajectory head. Check format of trajectory and molecular input file."
					IF (PRESENT(exit_status)) THEN
						IF (exit_status<0) WRITE(*,*) " #  end of file encountered!"
					ENDIF
					WRITE(*,*) "--> use correct number of steps or switch to reading the whole trajectory."
					CALL finalise_global()
					STOP
				CASE (86)
					WRITE(*,*) " #  WARNING 86: Couldn't find core for drude particle (see EXIT STATUS)."
					WRITE(*,*) " #  This drude particle will be ignored in some analyses (e.g. temperature)."
					WRITE(*,*) "--> Program will continue. Check your trajectory and (molecular) input files."
				CASE (87)
					WRITE(*,*) " #  ERROR 87: Couldn't read 'drude' section from molecular input file."
					WRITE(*,*) "--> Check format of molecular input file! Drudes are not yet assigned."
					IF (INFORMATION_IN_TRAJECTORY=="POS") THEN
						WRITE(*,*) "--> Will try to use first step of the trajectory later."
					ELSE
						WRITE(*,*) "--> Drude particles can be automatically assigned, if positions are given."
					ENDIF
				CASE (88)
					WRITE(*,*) " #  SERIOUS WARNING 88: Manual assignment didn't include all drude particles."
					WRITE(*,*) "--> carefully check your molecular input file and your trajectory."
				CASE (89)
					WRITE(*,*) " # WARNING 89: Reading the trajectory sequentially is extremely slow for some analyses."
					WRITE(*,*) "--> if possible, load into RAM."
				CASE (90)
					WRITE(*,*) " # WARNING 90: Value for student's t distribution for given N not available."
					WRITE(*,*) "--> Using 1.960 instead."
				CASE (91)
					WRITE(*,*) " #  ERROR 91: Drude particles are not assigned to their respective cores."
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (92)
					WRITE(*,*) " #  WARNING 92: Box volume will be overwritten."
				CASE (93)
					WRITE(*,*) " #  WARNING 93: lower bound is not smaller than higher bound, box volume is not changed."
				CASE (95)
					WRITE(*,*) " #  ERROR 95: Molecule recognition failed."
					WRITE(*,*) " #  check path, format, and filename of given trajectory!"
					WRITE(*,*)
				CASE (96)
					WRITE(*,*) " #  ERROR 96: Atoms of one of the molecular units (see EXIT STATUS) were separated in the trajectory."
					WRITE(*,*) " #  ensure that trajectory is sorted, and that the used cutoff is meaningful for your system!"
				CASE (97)
					WRITE(*,*) " #  ERROR 91: Cannot dump neighbours - need enough different molecules!"
					WRITE(*,*) " #  (either different molecule types, or different molecules of one type)"
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (98)
					WRITE(*,*) " #  ERROR 98: couldn't find '",TRIM(FILENAME_DISTRIBUTION_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (99)
					WRITE(*,*) " #  SEVERE ERROR 99: couldn't read '",TRIM(FILENAME_DISTRIBUTION_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the distribution input file
					CALL finalise_global()
					STOP
				CASE (100)
					WRITE(*,*) " #  ERROR 100: problem streaming '",TRIM(FILENAME_DISTRIBUTION_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_DISTRIBUTION_INPUT),"'!"
				CASE (101)
					WRITE(*,*) " #  WARNING 101: integer overflow - switch to floating point operation?"
					WRITE(*,*) "--> Program will continue. Check Results."
				CASE (102)
					WRITE(*,*) " #  ERROR 102: command line argument missing."
				CASE (103)
					WRITE(*,*) " #  ERROR 103: Not enough timesteps in trajectory for the requested jump analysis."
					WRITE(*,*) " #  (see EXIT STATUS for required minimum number of steps)"
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (104)
					WRITE(*,*) " #  ERROR 104: bin count out of sensible range."
					WRITE(*,*) "--> bin count is set to nearest sensible value. (see EXIT STATUS)"
				CASE (105)
					WRITE(*,*) " #  WARNING 105: startstep is bigger than endstep."
					WRITE(*,*) "--> startstep and endstep are swapped."
				CASE (106)
					WRITE(*,*) " #  ERROR 106: Invalid TIME_SCALING_FACTOR (see EXIT STATUS)."
					WRITE(*,*) "--> reset to default."
				CASE (107)
					WRITE(*,*) " #  WARNING 107: Cannot read timestep item from trajectory head."
				CASE (108)
					WRITE(*,*) " #  WARNING 108: Inconsistency in timestep item detected."
					WRITE(*,*) "--> Check for double timesteps!"
				CASE (109)
					WRITE(*,*) " #  ERROR 109: Couldn't read 'atomic_masses' section from molecular input file."
					WRITE(*,*) "--> Check format of molecular input file! Will use default masses where possible."
				CASE (110)
					WRITE(*,*) " #  WARNING 110: Invalid molecule type. This line is ignored."
				CASE (111)
					WRITE(*,*) " #  WARNING 111: Invalid atom index. This line is ignored."
				CASE (112)
					WRITE(*,*) " #  WARNING 112: Atomic mass already specified - will be overwritten."
				CASE (113)
					WRITE(*,*) " #  WARNING 113: Simple mode not (yet) available for this keyword."
				CASE (114)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 114: Overwriting existing 'prealpha_simple.inp'."
				CASE (115)
					WRITE(*,*) " #  SEVERE ERROR 115: Couldn't allocate memory for vc_components."
					WRITE(*,*) " #  give more RAM!"
					CALL finalise_global()
					STOP
				CASE (116)
					WRITE(*,*) " #  ERROR 116: Invalid molecule type in velocity correlation component (see EXIT STATUS)."
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (117)
					WRITE(*,*) " #  ERROR 117: Invalid self/distinct flag in velocity correlation component."
					WRITE(*,*) " #  self contributions require the two molecule type indices to be the same."
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (118)
					WRITE(*,*) " #  ERROR 118: Number of components < 1 (see EXIT STATUS)."
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (119)
					WRITE(*,*) " #  ERROR 119: specified molecule index (see EXIT STATUS) is out of bounds."
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (120)
					WRITE(*,*) " #  ERROR 120: No charged particles - cannot calculate electrical conductivity."
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (121)
					WRITE(*,*) " #  ERROR 121: Number of dihedral conditions < 1 (see EXIT STATUS). "
					WRITE(*,*) "--> Main program will continue, this analysis is aborted."
				CASE (122)
					WRITE(*,*) " #  ERROR 122: Custom default charges only for defined elements, not lowercase letters."
				CASE (123)
					WRITE(*,*) " #  ERROR 123: Couldn't read 'default_charges' section from molecular input file."
					WRITE(*,*) "--> Check format of molecular input file!"
				CASE (124)
					WRITE(*,*) " #  ERROR 124: Couldn't read 'atomic_charges' section from molecular input file."
					WRITE(*,*) "--> Check format of molecular input file! Will use default charges where possible."
				CASE (125)
					WRITE(*,*) " #  WARNING 125: Atomic charge already specified - will be overwritten."
				CASE (126)
					WRITE(*,*) " #  SERIOUS WARNING 126: Sum of atomic charges NOT EQUAL to specified molecular charge!"
					WRITE(*,*) "--> Check molecular input file. Don't use atomic charges if not required."
				CASE (127)
					WRITE(*,*) " #  WARNING 127: Something is wrong with the atomic charges of this molecule type index."
					WRITE(*,*) "--> Main program will continue anyway - check molecular input file."
					WRITE(*,*) "--> specify charges in the molecular input file?"
				CASE (128)
					WRITE(*,*) " #  WARNING 128: Neither default nor atomic charges have been defined by the user."
					WRITE(*,*) "--> Main program will continue anyway - check molecular input file."
				CASE (129)
					WRITE(*,*) " #  WARNING 129: redundant jump length specified by 'jump_analysis'."
				CASE (130)
					WRITE(*,*) " #  ERROR 130: Requested more jump analyses than current maximum (see EXIT STATUS)."
				CASE (131)
					WRITE(*,*) " #  WARNING 131: jump length should be large enough for accurate results - at least 100 time steps."
					WRITE(*,*) "--> This is currently not the case - consider choosing a trajectory with shorter stride"
				CASE (132)
					WRITE(*,*) " #  ERROR 132: couldn't find '",TRIM(FILENAME_DISTRIBUTION_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (133)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 133: Two times the same molecule_type_index - cross interactions will be zero."
				CASE (134)
					WRITE(*,*) " #  ERROR 134: couldn't find '",TRIM(FILENAME_DISTANCE_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (135)
					WRITE(*,*) " #  SEVERE ERROR 135: couldn't read '",TRIM(FILENAME_DISTANCE_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the distance input file
					CALL finalise_global()
					STOP
				CASE (136)
					WRITE(*,*) " #  WARNING 136: Invalid element name (for this molecule type). This line is ignored."
				CASE (137)
					WRITE(*,*) " #  SEVERE ERROR 137: couldn't allocate memory for distance list (or atom_indices)."
					WRITE(*,*) " #  Program will stop immediately. Please report this issue."
					STOP
				CASE (138)
					WRITE(*,*) " #  ERROR 138: indices_array of wrong size was passed to module MOLECULAR (atomic indices)."
					WRITE(*,*) "--> Program will try to continue anyway, probably crashes."
				CASE (139)
					WRITE(*,*) " #  ERROR 139: No (remaining) valid subjobs in module DISTANCE."
					WRITE(*,*) "--> Check your input files."
				CASE (140)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 140: Swapped reference and observed atoms (wildcard -> specific not allowed)."
				CASE (141)
					WRITE(*,*) " #  WARNING 141: Not enough molecules available (as defined by molecule_index). This line is ignored."
				CASE (142)
					WRITE(*,*) " #  ERROR 142: problem streaming '",TRIM(FILENAME_DISTANCE_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_DISTANCE_INPUT),"'!"
				CASE (143)
					WRITE(*,*) " #  WARNING 143: No charged particles - centre of charge trajectory will be empty."
					WRITE(*,*) "--> Please specify both atomic (real) and molecular (integer) charges."
				CASE (144)
					WRITE(*,*) " #  ERROR 144: The requested trajectory format is not supported."
					WRITE(*,*) "--> Use 'gro' for GROMACS, 'lmp' for LAMMPS, 'xyz' for xyz format."
					WRITE(*,*) "--> Main program will continue, but this analysis is aborted."
				CASE (145)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 145: GROMACS assumes nm, prealpha assumes Angström."
					WRITE(*,*) "--> Coordinates will be divided by 10 for '.gro' output."
				CASE (146)
					WRITE(*,*) " #  ERROR 146: couldn't find '",TRIM(FILENAME_DISPERSION_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (147)
					WRITE(*,*) " #  SEVERE ERROR 147: couldn't read '",TRIM(FILENAME_DISPERSION_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the distance input file
					CALL finalise_global()
					STOP
				CASE (148)
					WRITE(*,*) " #  ERROR 148: No (remaining) valid subjobs in module DISPERSION."
					WRITE(*,*) "--> Check your input files."
				CASE (149)
					WRITE(*,*) " #  ERROR 149: problem streaming '",TRIM(FILENAME_DISPERSION_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_DISPERSION_INPUT),"'!"
				CASE (150)
					WRITE(*,*) " #  SEVERE ERROR 150: couldn't allocate memory for dispersion list (or atom_indices)."
					WRITE(*,*) " #  Program will stop immediately. Please report this issue."
					STOP
				CASE (151)
					WRITE(*,*) " #  ERROR 151: Requested both wrapping and unwrapping."
					WRITE(*,*) " #  Trajectory will be left as is - no (un)wrapping action."
				CASE (152)
					WRITE(*,*) " #  ERROR 152: Requested unwrapping with sequential read."
					WRITE(*,*) " #  Trajectory will not be unwrapped! Do do so, please load in RAM."
				CASE (153)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 153: 'exponent' has no meaning for alpha2 calculations."
				CASE (154)
					WRITE(*,*) " #  WARNING 154: alpha2 only makes sense with the unit vector."
					WRITE(*,*) " #  (i.e. the projection 1 1 1, analysis continues but results are likely meaningless)"
				CASE (155)
					WRITE(*,*) " #  SEVERE ERROR 155: couldn't read '",TRIM(FILENAME_SPECIATION_INPUT),"'"
					WRITE(*,*) " #  Check format of input file!"
					CLOSE(UNIT=3)!unit 3 is the speciation input file
					CALL finalise_global()
					STOP
				CASE (156)
					WRITE(*,*) " #  SEVERE ERROR 156: couldn't allocate memory for species list (or clipboards)."
					WRITE(*,*) " #  Program will stop immediately. Please report this issue."
					STOP
				CASE (157)
					WRITE(*,*) " #  ERROR 157: couldn't find '",TRIM(FILENAME_SPECIATION_INPUT),"'"
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (158)
					WRITE(*,*) " #  ERROR 158: problem streaming '",TRIM(FILENAME_SPECIATION_INPUT),"'"
					WRITE(*,*) " #  check format of '",TRIM(FILENAME_SPECIATION_INPUT),"'!"
				CASE (159)
					WRITE(*,*) " #  ERROR 159: not enough valid acceptor and donor atoms."
					WRITE(*,*) "--> redelivering control to main unit"
				CASE (160)
					WRITE(*,*) " #  SERIOUS WARNING 160: Neighbour overflow occurred in Module SPECIATION."
					WRITE(*,*) "--> consider increasing N_neighbours"
				CASE (161)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 161: Species overflow occurred in Module SPECIATION."
					WRITE(*,*) "--> consider increasing N_species"
				CASE (162)
					WRITE(*,*) " #  WARNING 162: Overwriting existing atom groups in Module SPECIATION."
					WRITE(*,*) "--> carefully check that your atom groups are what you think they are."
				CASE (163)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 163: Atoms already grouped by elements."
					WRITE(*,*) "--> carefully check that your atom groups are what you think they are."
				CASE (164)
					WRITE(*,*) " #  ERROR 164: Could not add atom to group - no group defined."
					WRITE(*,*) "--> please first make a new group with new_acceptor_group or new_donor_group"
				CASE (165)
					error_count=error_count-1
					WRITE(*,*) " #  NOTICE 165: An atom_index was assigned twice."
					WRITE(*,*) "--> carefully check that your atom groups are what you think they are."
				CASE (166)
					WRITE(*,*) " #  ERROR 166: 'time_series' is required for 'autocorrelation'."
					WRITE(*,*) "--> please add 'time_series T' to the body of your speciation input file."
				CASE (167)
					WRITE(*,*) " #  ERROR 167: 'xxx_time_series.dat' could not be read."
					WRITE(*,*) "--> did the previous analysis terminate normally? check if file is compromised."
				CASE DEFAULT
					WRITE(*,*) " #  ERROR: Unspecified error"
				END SELECT
			ENDIF
			IF (error_count>MAXITERATIONS) THEN
				WRITE(*,*) " #  SEVERE ERROR 94: Error count exceeds maximum."
				error_count=0
				CALL finalise_global()
				STOP
			ENDIF
		END SUBROUTINE report_error

		INTEGER FUNCTION give_error_count()
			give_error_count=error_count
		END FUNCTION give_error_count

		! THE FOLLOWING PART IS REQUIRED BY THE MODULE "ANGLES"
		FUNCTION crossproduct(a,b) !This function returns the crossproduct of the vectors a and b.
			REAL(KIND=WORKING_PRECISION) :: crossproduct(3) !higher precision, because intermediate result.
			REAL(KIND=WORKING_PRECISION),INTENT(IN) :: a(3),b(3)
			crossproduct(1)=a(2)*b(3)-a(3)*b(2)
			crossproduct(2)=a(3)*b(1)-a(1)*b(3)
			crossproduct(3)=a(1)*b(2)-a(2)*b(1)
		END FUNCTION crossproduct

		SUBROUTINE normalize3D(vector)!normalizes a 3D vector
		IMPLICIT NONE
		REAL(KIND=WORKING_PRECISION),INTENT(INOUT) :: vector(3)
		REAL(KIND=WORKING_PRECISION) :: length
			length=SQRT((vector(1))**2+(vector(2))**2+(vector(3))**2)
			IF (length==0.0_WORKING_PRECISION) CALL report_error(1)
			vector=(vector/length)
		END SUBROUTINE normalize3D

		!legendre_polynomial computes the legendre polynomial (of the given order) of x_value
		REAL(KIND=WORKING_PRECISION) FUNCTION legendre_polynomial(x_value,order)
		IMPLICIT NONE
		REAL(KIND=WORKING_PRECISION),INTENT(IN) :: x_value
		INTEGER,INTENT(IN) :: order
			SELECT CASE (order)
			CASE (0)
				legendre_polynomial=0.0d0
			CASE (1)
				legendre_polynomial=x_value
			CASE (2)
				legendre_polynomial=1.5d0*x_value**2-0.5d0
			CASE (3)
				legendre_polynomial=2.5d0*x_value**3-1.5d0*x_value
			CASE (4)
				legendre_polynomial=4.375d0*x_value**4-3.75d0*x_value**2+0.375
			CASE DEFAULT
				CALL report_error(0,exit_status=order)
				legendre_polynomial=0.0d0
			END SELECT
		END FUNCTION legendre_polynomial

		SUBROUTINE normalize2D(vector)!normalizes a 2D vector
		IMPLICIT NONE
		REAL(KIND=WORKING_PRECISION),INTENT(INOUT) :: vector(2)
		REAL(KIND=WORKING_PRECISION) :: length
			length=SQRT((vector(1))**2+(vector(2))**2)
			IF (length==0.0_WORKING_PRECISION) CALL report_error(3)
			vector=(vector/length)
		END SUBROUTINE normalize2D
		! END OF THE PART THAT BELONGS TO THE MODULE ANGLES

		!asks the user to input an integer in the specified boundaries.
		INTEGER FUNCTION user_input_integer(low,high)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: low,high
		INTEGER :: inputinteger,ios		
			DO
				WRITE(*,ADVANCE="NO",FMT='("  > ")')
				READ(*,IOSTAT=ios,FMT=*) inputinteger
				IF (ios/=0) THEN
					CALL report_error(43)
					WRITE(*,'(" Please enter an integer.")')
				ELSE
					IF ((inputinteger<low).OR.(inputinteger>high)) THEN
						WRITE(*,'(" Please enter an integer between ",I0," and ",I0)') low,high
					ELSE
						EXIT
					ENDIF
				ENDIF
			ENDDO
			user_input_integer=inputinteger
			WRITE(*,*)
		END FUNCTION user_input_integer

		!asks the user to input a real number in the specified boundaries.
		REAL FUNCTION user_input_real(low,high)
		IMPLICIT NONE
		REAL,INTENT(IN) :: low,high
		INTEGER :: ios
		REAL :: inputreal
			DO
				WRITE(*,ADVANCE="NO",FMT='("  > ")')
				READ(*,IOSTAT=ios,FMT=*) inputreal
				IF (ios/=0) THEN
					CALL report_error(53)
					WRITE(*,'(" Please enter a real number.")')
				ELSE
					IF ((inputreal<low).OR.(inputreal>high)) THEN
						IF ((low<0.001).OR.(high>999.9)) THEN
							WRITE(*,'(" Please enter a real number between ",E16.8," and ",E16.8)') low,high
						ELSE
							WRITE(*,'(" Please enter a real number between ",F5.1," and ",F5.1)') low,high
						ENDIF
					ELSE
						EXIT
					ENDIF
				ENDIF
			ENDDO
			user_input_real=inputreal
			WRITE(*,*)
		END FUNCTION user_input_real


		!asks the user to input a logical - "yes" = .TRUE. / "no" = .FALSE.
		LOGICAL FUNCTION user_input_logical()
		IMPLICIT NONE
		CHARACTER(LEN=1) :: inputstring
		INTEGER :: ios
			DO
				WRITE(*,ADVANCE="NO",FMT='("  > ")')
				READ(*,IOSTAT=ios,FMT=*) inputstring
				IF (ios/=0) THEN
					CALL report_error(44)
				ELSE
					IF ((inputstring=="y").OR.(inputstring=="Y")) THEN
						user_input_logical=.TRUE.
						EXIT
					ELSEIF ((inputstring=="n").OR.(inputstring=="N")) THEN
						user_input_logical=.FALSE.
						EXIT
					ELSE
						CALL report_error(44)
					ENDIF
				ENDIF
			ENDDO
			WRITE(*,*)
		END FUNCTION user_input_logical

		!asks the user to input a string of given length
		FUNCTION user_input_string(length)
		IMPLICIT NONE
		INTEGER :: ios,length,i,charnum,input_size
		CHARACTER(LEN=length) :: user_input_string
		CHARACTER(LEN=length) :: inputstring
			DO
				WRITE(*,ADVANCE="NO",FMT='("  > ")')
				READ(*,IOSTAT=ios,FMT='(A)') inputstring
				IF (ios/=0) THEN
					CALL report_error(50)!how do you even enter a string incorrectly???
				ELSE
					input_size=LEN(TRIM(inputstring))
					IF (LEN(TRIM(inputstring))==length) CALL report_error(49,exit_status=input_size)
					IF (input_size/=LEN(ADJUSTL(TRIM(inputstring))).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "leading spaces removed."
					!remove blanks from the beginning...
					inputstring=ADJUSTL(inputstring)
					!get length of string
					input_size=LEN(TRIM(inputstring))
					!check for blanks
					DO i=1,input_size,1
						IF (" "==(inputstring(i:i))) THEN
							CALL report_error(48)
							EXIT
						ENDIF
					ENDDO
					!check for non-standard characters
					DO i=1,input_size,1
						charnum=IACHAR(inputstring(i:i))
						IF (.NOT.(ANY(ALPHABET==charnum))) THEN
							IF (" "/=(inputstring(i:i))) THEN
								CALL report_error(47,charnum)
								EXIT
							ENDIF
						ENDIF
					ENDDO
					WRITE(*,*)
					user_input_string=inputstring
					EXIT
				ENDIF
			ENDDO
		END FUNCTION user_input_string

		SUBROUTINE user_friendly_time_output(seconds)
		IMPLICIT NONE
		REAL(8) :: seconds
		IF (seconds<(999.0d-6)) THEN
			WRITE(*,'(F5.1,A)') seconds*(1.0d6)," microseconds"
		ELSEIF (seconds<(999.0d-3)) THEN
			WRITE(*,'(F5.1,A)') seconds*(1.0d3)," milliseconds"
		ELSEIF (seconds>(86400.0d0)) THEN
			WRITE(*,'(F5.1,A)') seconds/(86400.0d0)," days"
		ELSEIF (seconds>(3600.0d0)) THEN
			WRITE(*,'(F5.1,A)') seconds/(3600.0d0)," hours"
		ELSEIF (seconds>(60.0d0)) THEN
			WRITE(*,'(F5.1,A)') seconds/(60.0d0)," minutes"
		ELSE
			WRITE(*,'(F5.1,A)') seconds," seconds"
		ENDIF
		END SUBROUTINE user_friendly_time_output

		!This function returns the values for student's t distribution at 95% confidence level for a given N
		REAL FUNCTION student_t_value(N)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: N
			SELECT CASE (N)
			CASE (1)
				student_t_value=12.71
			CASE (2)
				student_t_value=4.303
			CASE (3)
				student_t_value=3.182
			CASE (4)
				student_t_value=2.776
			CASE (5)
				student_t_value=2.571
			CASE (6)
				student_t_value=2.447
			CASE (7)
				student_t_value=2.365
			CASE (8)
				student_t_value=2.306
			CASE (9)
				student_t_value=2.262
			CASE (10)
				student_t_value=2.228
			CASE (11)
				student_t_value=2.201
			CASE (12)
				student_t_value=2.179
			CASE (13)
				student_t_value=2.160
			CASE (14)
				student_t_value=2.145
			CASE (15)
				student_t_value=2.131
			CASE (16)
				student_t_value=2.120
			CASE (17)
				student_t_value=2.110
			CASE (18)
				student_t_value=2.101
			CASE (19)
				student_t_value=2.093
			CASE (20)
				student_t_value=2.086
			CASE (21)
				student_t_value=2.080
			CASE (22)
				student_t_value=2.074
			CASE (23)
				student_t_value=2.069
			CASE (24)
				student_t_value=2.064
			CASE (25)
				student_t_value=2.060
			CASE (26)
				student_t_value=2.056
			CASE (27)
				student_t_value=2.052
			CASE (28)
				student_t_value=2.048
			CASE (29)
				student_t_value=2.045
			CASE (30)
				student_t_value=2.042
			CASE (31)
				student_t_value=2.040
			CASE (32)
				student_t_value=2.037
			CASE (48)
				student_t_value=2.011
			CASE DEFAULT
				CALL report_error(90,exit_status=N)
				student_t_value=1.960
			END SELECT
		END FUNCTION student_t_value

		REAL FUNCTION covalence_radius(element_name)
		IMPLICIT NONE
		CHARACTER(LEN=*),INTENT(IN) :: element_name
			!IF you change this part, THEN change Module_Molecular, too!
			SELECT CASE (TRIM(element_name))
			CASE ("H")
				covalence_radius=0.33
			CASE ("F")
				covalence_radius=0.71
			CASE ("B")
				covalence_radius=0.82
			CASE ("Cl")
				covalence_radius=0.99
			CASE ("Br")
				covalence_radius=1.14
			CASE ("I")
				covalence_radius=1.33
			CASE ("N")
				covalence_radius=0.71
			CASE ("O")
				covalence_radius=0.73
			CASE ("C")
				covalence_radius=0.77
			CASE ("S")
				covalence_radius=1.02
			CASE ("P")
				covalence_radius=1.06
			CASE ("Li")
				covalence_radius=1.34
			CASE ("Na")
				covalence_radius=1.54
			CASE DEFAULT
				covalence_radius=1.00
			END SELECT
		END FUNCTION covalence_radius

		SUBROUTINE timing_parallel_sections(start)
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_wtime()
	 !$ 	REAL(8) :: OMP_get_wtime
	 !$ 	END FUNCTION OMP_get_wtime
	 !$ END INTERFACE
		LOGICAL :: start
	 !$ REAL(8) :: clipboard_real
	 !$ REAL(8),SAVE :: timeline_real=0.0d0
		 !$ clipboard_real=OMP_get_wtime()
		 !$ IF (start) THEN
		 !$ 	timeline_real=clipboard_real
		 !$ ELSE
		 !$ 	IF ((timeline_real>0.0d0).AND.(TIME_OUTPUT)) THEN
		 !$ 		CALL user_friendly_time_output(clipboard_real-timeline_real)
		 !$ 	ENDIF
		 !$ ENDIF
			!Flush I/O to ease identification of bottlenecks
			CALL refresh_IO()
		END SUBROUTINE timing_parallel_sections

		CHARACTER(LEN=3) FUNCTION logical_to_yesno(input)
		IMPLICIT NONE
		LOGICAL,INTENT(IN) :: input
			IF (input) THEN
				logical_to_yesno="yes"
			ELSE
				logical_to_yesno="no"
			ENDIF
		END FUNCTION logical_to_yesno

		!prints the progress, is initialised by passing the number of total iterations.
		!Thus, the calling structure should be the following:
		!	IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-start+sampling)/sampling,0))
		!	DO something=start,nsteps,sampling
		!		IF (VERBOSE_OUTPUT) CALL print_progress()
		!	ENDDO
		!	IF (((MAX((nsteps-1+sampling)/sampling,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
		SUBROUTINE print_progress(total_iterations_in)
		IMPLICIT NONE
		INTEGER,INTENT(IN),OPTIONAL :: total_iterations_in
		INTEGER,SAVE :: progress_counter,total_iterations,iteration_counter
		LOGICAL,SAVE :: printsteps
			IF (.NOT.(VERBOSE_OUTPUT)) RETURN
			IF (PRESENT(total_iterations_in)) THEN
				!initialise
				IF (total_iterations_in>100) THEN
					WRITE(*,*) "0% - - - - - 50% - - - - - 100%"
					CALL refresh_IO()
					WRITE(*,ADVANCE="NO",FMT='(" ")')
					progress_counter=0
					total_iterations=total_iterations_in
					printsteps=.TRUE.
					iteration_counter=0
				ELSE
					progress_counter=40
					printsteps=.FALSE.
				ENDIF
			ELSE
				iteration_counter=iteration_counter+1
				IF (printsteps) THEN
					IF ((MOD(iteration_counter,total_iterations/31)==0).AND.(progress_counter<31)) THEN
						WRITE(*,ADVANCE="NO",FMT='("#")')
						CALL refresh_IO()
						progress_counter=progress_counter+1
					ENDIF
				ENDIF
			ENDIF
		END SUBROUTINE print_progress

		SUBROUTINE refresh_IO()
		IMPLICIT NONE
			CALL FLUSH()
		END SUBROUTINE refresh_IO

END MODULE SETTINGS
!--------------------------------------------------------------------------------------------------------------------------------!