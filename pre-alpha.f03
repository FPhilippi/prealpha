! RELEASED ON 14_May_2020 AT 17:44

    ! prealpha - a tool to extract information from molecular dynamics trajectories.
    ! Copyright (C) 2020 Frederik Philippi
    ! This work is funded by the Imperial President's PhD Scholarship.

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
 LOGICAL,PARAMETER :: DEVELOPERS_VERSION_DEFAULT=.FALSE.
 LOGICAL,PARAMETER :: ERROR_OUTPUT_DEFAULT=.TRUE.
 LOGICAL,PARAMETER :: PARALLEL_OPERATION_DEFAULT=.TRUE.
 LOGICAL,PARAMETER :: READ_SEQUENTIAL_DEFAULT=.FALSE.
 LOGICAL,PARAMETER :: BOX_VOLUME_GIVEN_DEFAULT=.FALSE.!If (T), then xlo,xhi,ylo, etc should be correct, too.
 INTEGER,PARAMETER :: ERROR_CODE_DEFAULT=-1
 INTEGER,PARAMETER :: TIME_SCALING_FACTOR_DEFAULT=1
 LOGICAL,PARAMETER :: WRAP_TRAJECTORY_DEFAULT=.FALSE.
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
 CHARACTER(LEN=1024) :: REDIRECTED_OUTPUT="output.dat"
 CHARACTER(LEN=1024) :: OUTPUT_PREFIX="" !prefix added to output, for example to distinguish between different autocorrelation analyses
 CHARACTER(LEN=3) :: INFORMATION_IN_TRAJECTORY="UNK"!does the trajectory contain velocities (VEL) or coordinates (POS)????
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

 !PRIVATE/PUBLIC declarations
 PUBLIC :: normalize2D,normalize3D,crossproduct,report_error,timing_parallel_sections,legendre_polynomial
 PUBLIC :: FILENAME_TRAJECTORY,PATH_TRAJECTORY,PATH_INPUT,PATH_OUTPUT,user_friendly_time_output
 PUBLIC :: user_input_string,user_input_integer,user_input_logical,user_input_real
 PUBLIC :: student_t_value,covalence_radius,give_error_count,DEVELOPERS_VERSION
 PRIVATE :: q !that's not a typo!
 PRIVATE :: error_count
 CONTAINS

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
     WRITE(*,*) " #  SEVERE ERROR 4 in atomic_weight: unknown element"
     WRITE(*,*) " #  If necessary, add element to function atomic_weight in module MOLECULAR"
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
     WRITE(*,*) " #   ERROR 69: the specified molecule index doesn't exist."
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
     WRITE(*,*) " #  ERROR 73: Trajectory contains velocities - wrapping not available."
     WRAP_TRAJECTORY=.FALSE.
    CASE (74)
     WRITE(*,*) " #  ERROR 74: unphysical number of constraints."
     WRITE(*,*) " #  Temperature values don't include the constraints correction!"
    CASE (75)
     WRITE(*,*) " #  ERROR 75: specified molecules to export dihedrals exceed total number."
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
   CASE DEFAULT
    covalence_radius=1.00
   END SELECT
  END FUNCTION covalence_radius

  SUBROUTINE timing_parallel_sections(start)
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_wtime()
  !$  REAL(8) :: OMP_get_wtime
  !$  END FUNCTION OMP_get_wtime
  !$ END INTERFACE
  LOGICAL :: start
  !$ REAL(8) :: clipboard_real
  !$ REAL(8),SAVE :: timeline_real=0.0d0
   !$ clipboard_real=OMP_get_wtime()
   !$ IF (start) THEN
   !$  timeline_real=clipboard_real
   !$ ELSE
   !$  IF ((timeline_real>0.0d0).AND.(TIME_OUTPUT)) THEN
   !$   CALL user_friendly_time_output(clipboard_real-timeline_real)
   !$  ENDIF
   !$ ENDIF
   !Flush I/O to ease identification of bottlenecks
   CALL refresh_IO()
  END SUBROUTINE timing_parallel_sections
  
  !prints the progress, is initialised by passing the number of total iterations.
  SUBROUTINE print_progress(total_iterations_in)
  IMPLICIT NONE
  INTEGER,INTENT(IN),OPTIONAL :: total_iterations_in
  INTEGER,SAVE :: progress_counter,total_iterations,iteration_counter
  LOGICAL,SAVE :: printsteps
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
    IF ((printsteps).AND.(MOD(iteration_counter,total_iterations/31)==0).AND.(progress_counter<31)) THEN
     WRITE(*,ADVANCE="NO",FMT='("#")')
     CALL refresh_IO()
     progress_counter=progress_counter+1
    ENDIF
   ENDIF
  END SUBROUTINE print_progress

  SUBROUTINE refresh_IO()
  IMPLICIT NONE
   CALL FLUSH()
  END SUBROUTINE refresh_IO

END MODULE SETTINGS
!--------------------------------------------------------------------------------------------------------------------------------!
!This Module can be used to perform rotations and also to turn a set of coordinates into a dihedral angle.
MODULE ANGLES ! Copyright (C) 2020 Frederik Philippi
    USE SETTINGS
 IMPLICIT NONE
 PRIVATE
 !default values
 REAL(KIND=WORKING_PRECISION),PARAMETER :: tolerance=0.0001d0!tolerance in degrees used to define identity of angles. Also used in some other case as distance criterion.
 !variables
    REAL(KIND=WORKING_PRECISION) :: uvector(3) !uvector is the unit vector around which the rotation will occur
 REAL(KIND=WORKING_PRECISION) :: sina,cosa !sina and cosa are actually sinus(alpha) and cosinus(alpha), with alpha being the angle of the rotation.
    REAL(KIND=WORKING_PRECISION) :: rotation(3,3) !the rotation matrix
 !PRIVATE/PUBLIC declarations
 PUBLIC :: prepare_rotation,rotate,dihedral_angle,rotate_random!dihedral_angle changes the rotation matrix!
 PRIVATE :: make_rotation_matrix,uvector,sina,cosa,tolerance,rotation!Thou shalt not temper with the matrix.
 CONTAINS

  !brute-force generation of a uniformly distributed random rotation matrix. needs srand to be called before
  SUBROUTINE rotate_random(vector)
  IMPLICIT NONE
  REAL(KIND=GENERAL_PRECISION),INTENT(INOUT) ::  vector(3)
  REAL(KIND=WORKING_PRECISION) :: x,y,z
  REAL(KIND=WORKING_PRECISION) :: cosg,sing,cosb,sinb,alpha,beta,gamma,sina2,cosa2
  REAL(KIND=WORKING_PRECISION) :: rotation_local(3,3) !the rotation matrix - needs to be a local variable for parallelisation
   alpha=RAND()*twopi
   beta=RAND()*twopi
   gamma=RAND()*twopi
   sina2=SIN(alpha)
   cosa2=COS(alpha)
   sinb=SIN(beta)
   cosb=COS(beta)
   sing=SIN(gamma)
   cosg=COS(gamma)
   rotation_local(1,1)=cosb*cosg
   rotation_local(1,2)=-sing*cosb
   rotation_local(1,3)=sinb
   rotation_local(2,1)=sina2*sinb*cosg+sing*cosa2
   rotation_local(2,2)=-sina2*sinb*sing+cosa2*cosg
   rotation_local(2,3)=-sina2*cosb
   rotation_local(3,1)=-sinb*cosa2*cosg+sina2*sing
   rotation_local(3,2)=sinb*sing*cosa2+sina2*cosg
   rotation_local(3,3)=cosa2*cosb
   x=vector(1)
   y=vector(2)
   z=vector(3)
   vector(1)=x*rotation_local(1,1)+y*rotation_local(1,2)+z*rotation_local(1,3)
   vector(2)=x*rotation_local(2,1)+y*rotation_local(2,2)+z*rotation_local(2,3)
   vector(3)=x*rotation_local(3,1)+y*rotation_local(3,2)+z*rotation_local(3,3)
  END SUBROUTINE rotate_random

  SUBROUTINE make_rotation_matrix()! This subroutine builds the rotation matrix. Not accessible globally.
  IMPLICIT NONE
   rotation(1,1)=cosa+(uvector(1)*uvector(1))*(1.0d0-cosa)
   rotation(1,2)=uvector(1)*uvector(2)*(1.0d0-cosa)-uvector(3)*sina
   rotation(1,3)=uvector(1)*uvector(3)*(1.0d0-cosa)+uvector(2)*sina
   rotation(2,1)=uvector(1)*uvector(2)*(1.0d0-cosa)+uvector(3)*sina
   rotation(2,2)=cosa+(uvector(2)*uvector(2))*(1.0d0-cosa)
   rotation(2,3)=uvector(2)*uvector(3)*(1.0d0-cosa)-uvector(1)*sina
   rotation(3,1)=uvector(1)*uvector(3)*(1.0d0-cosa)-uvector(2)*sina
   rotation(3,2)=uvector(3)*uvector(2)*(1.0d0-cosa)+uvector(1)*sina
   rotation(3,3)=cosa+(uvector(3)*uvector(3))*(1.0d0-cosa)
  END SUBROUTINE make_rotation_matrix

  SUBROUTINE prepare_rotation(startvector,targetvector,aligned)!prepares the rotation matrix that maps 'startvector' onto 'targetvector' by rotation around an axis perpendicular to both these vectors.
  IMPLICIT NONE
  LOGICAL,INTENT(OUT),OPTIONAL :: aligned
  REAL(KIND=GENERAL_PRECISION),INTENT(IN) :: startvector(3),targetvector(3)
  REAL(KIND=WORKING_PRECISION) :: a(3),b(3),angle,dummy_axis(3)
   a=startvector
   b=targetvector
   CALL normalize3D(a)
   CALL normalize3D(b)
   cosa=DOT_PRODUCT(a,b)!get the angle between a and b
   sina=SQRT(1-cosa*cosa)!also sinus(angle), needed for the rotation matrix later.
   angle=(ACOS(cosa)*degrees)!calculate angle in degrees, mainly for debugging purposes
   IF (PRESENT(aligned)) aligned=.FALSE. !Some external procedures need to know whether the vectors are aligned, hence the optional variable.
   IF (angle<=tolerance) THEN !angle is zero - vectors are aligned!
    IF (DEVELOPERS_VERSION) WRITE(*,'(A,E7.1)') "  ! vectors are aligned, angle = ",angle
    IF (PRESENT(aligned)) aligned=.TRUE.
    !The rotation matrix is the identity matrix, because nothing should change.
    rotation(:,:)=0.0d0
    rotation(1,1)=1.0d0
    rotation(2,2)=1.0d0
    rotation(3,3)=1.0d0
   ELSE
    IF ((180.0-angle)<=tolerance) THEN
     IF (DEVELOPERS_VERSION) WRITE(*,'(A,F7.3)') "  ! vectors are antiparallel, angle = ",angle
     IF ((a(1)<=(1.0d0-tolerance)).AND.(a(1)>=(-1.0d0+tolerance))) THEN !This part is just to handle the rare case of antiparallel a and b along the x axis...
      dummy_axis(:)=0.0d0
      dummy_axis(1)=1.0d0
     ELSE
      dummy_axis(:)=0.0d0
      dummy_axis(3)=1.0d0
     ENDIF
     uvector=crossproduct(a,dummy_axis)
    ELSE !The two vectors are neither aligned nor antiparallel. This should be the normal case.
     uvector=crossproduct(a,b)
    ENDIF
    CALL normalize3D(uvector)
    CALL make_rotation_matrix()
   ENDIF
  END SUBROUTINE prepare_rotation

  SUBROUTINE rotate(vector)!subroutine that rotates the given vector, using the rotation matrix ("Rotation")
  IMPLICIT NONE
  REAL(KIND=GENERAL_PRECISION),INTENT(INOUT) ::  vector(3)
  REAL(KIND=WORKING_PRECISION) :: x,y,z
    x=vector(1)
    y=vector(2)
    z=vector(3)
    !I don't like MATMUL.
    vector(1)=x*Rotation(1,1)+y*Rotation(1,2)+z*Rotation(1,3)
    vector(2)=x*Rotation(2,1)+y*Rotation(2,2)+z*Rotation(2,3)
    vector(3)=x*Rotation(3,1)+y*Rotation(3,2)+z*Rotation(3,3)
  END SUBROUTINE rotate

  !dihedral_angle returns the dihedral from 0 to 360° if minusplus=FALSE (default), or from -180 to +180 if minusplus=TRUE. minusplus is an optional variable.
  REAL(KIND=GENERAL_PRECISION) FUNCTION dihedral_angle(dihedral_members,minusplus) ! dihedral_members is an array containing the four positions of the atoms.
  IMPLICIT NONE
  LOGICAL,INTENT(IN),OPTIONAL :: minusplus
  REAL(KIND=WORKING_PRECISION) :: dihedral_members(4,3),connection_vector(3),uvector1(3),uvector2(3)
  REAL(KIND=WORKING_PRECISION) :: projection1(2),projection2(2),ccw(2),clip,length
   connection_vector=(dihedral_members(2,:)-dihedral_members(3,:)) !vector connecting the second and third atom
   uvector1=(dihedral_members(2,:)-dihedral_members(1,:)) !vector from first to second atom
   uvector2=(dihedral_members(3,:)-dihedral_members(4,:)) !vector from fourth to third atom
   CALL normalize3D(connection_vector)
   CALL normalize3D(uvector1)
   CALL normalize3D(uvector2)
   !INITIALIZE ROTATION, RX=0
   uvector(1)=0.0d0 !X component is zero
   sina=SQRT(connection_vector(2)**2+connection_vector(3)**2)! 'sina' is here the length of the projection of the connection_vector along the x axis.
   IF (sina>=tolerance) THEN
     uvector(2)=connection_vector(3)/sina
     uvector(3)=-connection_vector(2)/sina!now, uvector is normalized.
     cosa=connection_vector(1)
     CALL make_rotation_matrix()
     CALL rotate(uvector1)
     CALL rotate(uvector2)
   ELSE!there are no components in y or z direction in the connection vector.
    IF (DEVELOPERS_VERSION) WRITE(*,*) "  ! dihedral_angle: already preoriented"
   ENDIF
   !Project to plane
   projection1=uvector1(2:3)
   projection2=uvector2(2:3)
   CALL normalize2D(projection1)
   CALL normalize2D(projection2)
   projection1=(projection1*projection2)
   length=(projection1(1)+projection1(2))
   IF (length<(-1.0d0)) THEN
    IF (length<(-1.0d0-tolerance)) CALL report_error(2)
    length=(-1.0d0)
   ELSEIF (length>(1.0d0)) THEN
    IF (length>(1.0d0+tolerance)) CALL report_error(2)
    length=(1.0d0)
   ENDIF
   clip=ACOS(length)
   ccw(1)=+COS(clip)*uvector1(2)-SIN(clip)*uvector1(3)
   ccw(2)=+SIN(clip)*uvector1(2)+COS(clip)*uvector1(3)
   CALL normalize2D(ccw)
   dihedral_angle=(clip*degrees)
   ccw=(ccw*projection2)
   clip=ccw(1)+ccw(2)
   IF (clip<(-1.0d0)) THEN
    IF (clip<(-1.0d0-tolerance)) CALL report_error(2)
    clip=(-1.0d0)
   ELSEIF (clip>(1.0d0)) THEN
    IF (clip>(1.0d0+tolerance)) CALL report_error(2)
    clip=(1.0d0)
   ENDIF
   IF ((ACOS(clip)*degrees)<=tolerance) THEN
    IF (PRESENT(minusplus)) THEN
     IF (minusplus) THEN 
      dihedral_angle=(-dihedral_angle)
     ELSE
      dihedral_angle=(360.0d0-dihedral_angle)
     ENDIF
    ELSE
     dihedral_angle=(360.0d0-dihedral_angle)
    ENDIF
   ENDIF
  END FUNCTION dihedral_angle

END MODULE ANGLES
!--------------------------------------------------------------------------------------------------------------------------------!
!This module is responsible for handling the trajectory and passing information to other modules.
MODULE MOLECULAR ! Copyright (C) 2020 Frederik Philippi
!Atomic masses are handled with single precision.
    USE SETTINGS
 IMPLICIT NONE
 PRIVATE
 !parameter
 REAL,PARAMETER :: default_mass_hydrogen=01.008
 REAL,PARAMETER :: default_mass_fluorine=18.998
 REAL,PARAMETER :: default_mass_boron=10.81
 REAL,PARAMETER :: default_mass_chlorine=35.45
 REAL,PARAMETER :: default_mass_bromine=79.904
 REAL,PARAMETER :: default_mass_iodine=126.904
 REAL,PARAMETER :: default_mass_nitrogen=14.007
 REAL,PARAMETER :: default_mass_oxygen=15.999
 REAL,PARAMETER :: default_mass_carbon=12.011
 REAL,PARAMETER :: default_mass_sulfur=32.066
 REAL,PARAMETER :: default_mass_phosphorus=30.974
 REAL,PARAMETER :: default_mass_lithium=6.94
 !variables
 REAL :: mass_hydrogen=default_mass_hydrogen
 REAL :: mass_fluorine=default_mass_fluorine
 REAL :: mass_boron=default_mass_boron
 REAL :: mass_chlorine=default_mass_chlorine
 REAL :: mass_bromine=default_mass_bromine
 REAL :: mass_iodine=default_mass_iodine
 REAL :: mass_nitrogen=default_mass_nitrogen
 REAL :: mass_oxygen=default_mass_oxygen
 REAL :: mass_carbon=default_mass_carbon
 REAL :: mass_sulfur=default_mass_sulfur
 REAL :: mass_phosphorus=default_mass_phosphorus
 REAL :: mass_lithium=default_mass_lithium
 LOGICAL :: fragments_initialised=.FALSE.!Status boolean, is true if the fragment_list has been initialised.
 !fragment lists: store the atom_indices of the fragments.
 INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_base(:) !List of centre-of-mass fragments (defined as atom_indices) for base atom
 INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_tip(:) !List of centre-of-mass fragments (defined as atom_indices) for tip atom
 INTEGER :: number_of_tip_atoms,number_of_base_atoms !extent of the two fragment_lists
 REAL(KIND=WORKING_PRECISION) :: mass_of_tip_fragment,mass_of_base_fragment !total masses of the two fragments
 INTEGER :: molecule_type_index_for_fragments !molecule type index of the molecule to which tip and base atoms belong to
 LOGICAL :: dihedrals_initialised=.FALSE.!Status boolean, is true if the dihedral_member_indices has been initialised.
 LOGICAL :: drudes_assigned=.FALSE.
 LOGICAL :: drudes_allocated=.FALSE.
 LOGICAL :: drude_details=.FALSE. ! will be TRUE when reduced mass, minimum_drude_distance and maximum_drude_distance are initialised.
 LOGICAL :: use_barycentre=.FALSE.
 INTEGER,DIMENSION(:,:),ALLOCATABLE :: dihedral_member_indices !list of atom indices used for reporting dihedral angles.
 INTEGER :: number_of_dihedrals,molecule_type_index_for_dihedrals!number of dihedrals to report, type of molecule they belong to.
 INTEGER :: headerlines_to_skip!header lines in trajectory file, e.g. 9 for a lammps file or 2 for xyz format.
 INTEGER :: lines_to_skip!the total lines that one timestep takes.
 TYPE,PRIVATE :: atom
        REAL(KIND=STORAGE_PRECISION) :: coordinates(3)=0.0d0
    END TYPE atom
 TYPE,PRIVATE :: drude_pair
  INTEGER :: drude_flag=-1 ! atom index of the drude particle attached to this core. Will be -1 if no drude particle is found, and 0 if this is a drude itself.
  REAL(KIND=GENERAL_PRECISION) :: reduced_mass=0.0d0 ! reduced mass of this drude/core pair
  REAL(KIND=STORAGE_PRECISION) :: minimum_drude_distance,maximum_drude_distance !minimum and maximum distances for this pair in the first step
 END TYPE drude_pair
 TYPE,PRIVATE :: molecule
  INTEGER :: constraints=0
  INTEGER :: charge=0
  INTEGER :: number_of_atoms=0 ! = extent of dimension one of trajectory, number of atoms PER SINGLE MOLECULE, not total!
  INTEGER :: total_molecule_count=0 ! = extent of dimension two of trajectory, number of molecules of this type in the box
  INTEGER :: number_of_drudes_in_molecule=0 ! = number of drude particles in this molecule (counter for successfully assigned drude pairs)
        REAL(KIND=GENERAL_PRECISION) :: mass=0.0d0
  CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE :: list_of_elements !--> Turned on support for  2-letter elements!
  TYPE(drude_pair),DIMENSION(:),ALLOCATABLE :: list_of_drude_pairs ! list of pairs of drude particles / drude cores / drudes
  REAL(KIND=GENERAL_PRECISION),DIMENSION(:),ALLOCATABLE :: list_of_atom_masses !corresponding masses for the atoms
  LOGICAL,DIMENSION(:),ALLOCATABLE :: manual_atom_mass_specified !any elements that are .TRUE. will not be initialised from the trajectory!
  TYPE(atom),DIMENSION(:,:),ALLOCATABLE :: snapshot !like trajectory, but for one timestep only. Required for READ_SEQUENTIAL.
  TYPE(atom),DIMENSION(:,:,:),ALLOCATABLE :: trajectory
  !trajectory is now in column major order!
  !first dimension: Atom index in molecule
  !second dimension: index of the molecule (not to be confused with molecule_type_index)
  !third dimension: timestep
  TYPE(atom),DIMENSION(:,:,:),ALLOCATABLE :: queue!the queue for circular parallel operation. conceptually the same as 'trajectory', but only big enough to act as a buffer.
 END TYPE molecule
 REAL(KIND=STORAGE_PRECISION) :: box_dimensions(2,3)!low and high for x,y,z
 REAL(KIND=STORAGE_PRECISION) :: box_size(3) !size of the box.
 REAL(KIND=STORAGE_PRECISION) :: maximum_distance !maximum possible distance in box.
 REAL(KIND=STORAGE_PRECISION) :: maximum_distance_squared !square of maximum possible distance in box.
    REAL(KIND=SP) :: drude_mass=0.0e0
 REAL(KIND=SP) :: COM_mass_list(IACHAR("a"):(IACHAR("a")+25)) !A list of user-specified masses for centre-of-mass trajectories.
 LOGICAL :: custom_default_masses !if 'T', then the user has specified his own default masses for lowercase letters or elements.
 LOGICAL :: custom_atom_masses !if 'T', then the user has manually specified atom masses.
 LOGICAL :: custom_constraints !if 'T', then the user has specified custom constraints on some molecules.
 !use_firstatom_as_com=.TRUE. ONLY with whole trajectory, NOT with read_sequential.
 LOGICAL :: use_firstatom_as_com=.FALSE. ! If 'T', then the first atom in any molecule is used instead of centre of mass.
 INTEGER :: file_position=-1!at which timestep the file to read is positioned. The first step is 1.
 INTEGER :: number_of_steps=0 !number of timesteps in whole trajectory
 INTEGER :: number_of_molecule_types=0 !number of different molecules, usually two (cation and anion)
 INTEGER :: total_number_of_atoms=0 !number of atoms per timestep
 INTEGER :: number_of_drude_particles=0 !number of drude particles per box
 INTEGER :: ndrudes_check=0 !number of drude particles specified in the molecular input file
 TYPE(molecule),DIMENSION(:),ALLOCATABLE :: molecule_list !list of molecules, usually two members: cation and anion. Has to be same order as in lammps trajectory. The different molecule types / members are usually referred to as 'molecule_type_index' in subroutines.
 !PRIVATE VARIABLES
 PRIVATE :: box_dimensions,drude_mass,number_of_molecule_types,total_number_of_atoms,number_of_steps,box_size
 PRIVATE :: dihedrals_initialised,dihedral_member_indices,number_of_dihedrals,wrap_snap
 PRIVATE :: file_position,goto_timestep,headerlines_to_skip,lines_to_skip,COM_mass_list,custom_default_masses,wrap_full
 PRIVATE :: number_of_drude_particles,allocate_drude_list,drudes_allocated,drudes_assigned,ndrudes_check,molecule_list
 PRIVATE :: fragments_initialised,fragment_list_base,fragment_list_tip,number_of_tip_atoms,number_of_base_atoms
 PRIVATE :: mass_of_tip_fragment,mass_of_base_fragment,molecule_type_index_for_fragments,custom_constraints
 PRIVATE :: custom_atom_masses,use_firstatom_as_com
 !PRIVATE/PUBLIC declarations
 PRIVATE :: report_trajectory_properties,reset_trajectory_file
 PUBLIC :: atomic_weight,load_trajectory,initialise_molecular,finalise_molecular,write_molecule,give_temperature
 PUBLIC :: give_dihedrals,initialise_dihedrals,give_number_of_molecule_types,give_number_of_atoms_per_molecule
 PUBLIC :: give_number_of_molecules_per_step,give_mass_of_molecule,show_molecular_settings,print_memory_requirement
 PUBLIC :: give_total_number_of_molecules_per_step,show_drude_settings,give_box_volume,give_box_limit
 PUBLIC :: assign_drudes,give_intramolecular_distances,give_total_temperature,give_sum_formula,constraints_available
 PUBLIC :: give_total_degrees_of_freedom,give_number_of_atoms_per_step,convert_parallel,compute_drude_temperature
 PUBLIC :: initialise_fragments,give_tip_fragment,give_base_fragment,give_number_of_drude_particles
 PUBLIC :: give_smallest_distance,give_number_of_neighbours,wrap_vector,give_smallest_distance_squared,compute_drude_properties
 PUBLIC :: compute_squared_radius_of_gyration,are_drudes_assigned,write_molecule_merged_drudes,set_cubic_box
 PUBLIC :: give_intermolecular_contact_distance,give_smallest_atom_distance,give_smallest_atom_distance_squared
 PUBLIC :: give_charge_of_molecule,give_number_of_timesteps,give_fragment_information,give_dihedral_member_indices
 PUBLIC :: report_element_lists,give_center_of_mass,write_header,give_maximum_distance_squared,set_default_masses
 PUBLIC :: print_atomic_masses,give_comboost,switch_to_barycenter
 CONTAINS

  SUBROUTINE set_default_masses()
  IMPLICIT NONE
   mass_hydrogen=default_mass_hydrogen
   mass_fluorine=default_mass_fluorine
   mass_boron=default_mass_boron
   mass_chlorine=default_mass_chlorine
   mass_bromine=default_mass_bromine
   mass_iodine=default_mass_iodine
   mass_nitrogen=default_mass_nitrogen
   mass_oxygen=default_mass_oxygen
   mass_carbon=default_mass_carbon
   mass_sulfur=default_mass_sulfur
   mass_phosphorus=default_mass_phosphorus
   mass_lithium=default_mass_lithium
  END SUBROUTINE set_default_masses

  SUBROUTINE subtract_drude_masses()
  IMPLICIT NONE
   IF (VERBOSE_OUTPUT) PRINT *," Subtracting drude masses from N,O,C,S,P,Li"
   mass_nitrogen=mass_nitrogen-drude_mass
   mass_oxygen=mass_oxygen-drude_mass
   mass_carbon=mass_carbon-drude_mass
   mass_sulfur=mass_sulfur-drude_mass
   mass_phosphorus=mass_phosphorus-drude_mass
   mass_lithium=mass_lithium-drude_mass
  END SUBROUTINE subtract_drude_masses

  SUBROUTINE print_atomic_masses()
  IMPLICIT NONE
  INTEGER :: molecule_type_index,atom_index
  REAL :: native_mass,atomic_mass
  CHARACTER(LEN=2) :: element_name
   WRITE(*,'(" The following lines can be used/altered to request custom atomic masses:")')
   WRITE(*,'("   atomic_masses ",I0)') SUM(molecule_list(:)%number_of_atoms)
   DO molecule_type_index=1,number_of_molecule_types,1
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
     atomic_mass=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
     native_mass=atomic_weight(molecule_list(molecule_type_index)%list_of_elements(atom_index))
     element_name=molecule_list(molecule_type_index)%list_of_elements(atom_index)
     WRITE(*,FMT='("   ",I0," ",I0," ",F0.3," (")',ADVANCE="NO") molecule_type_index,atom_index,atomic_mass
     IF (ABS(atomic_mass-native_mass)>0.001) THEN
      WRITE(*,'(A,", default mass was ",F0.3,")")') TRIM(element_name),native_mass
     ELSE
      WRITE(*,'("Element ",A,")")') TRIM(element_name)
     ENDIF
    ENDDO
   ENDDO
  END SUBROUTINE print_atomic_masses

  !The following subroutine sets the cubic box limits
  SUBROUTINE set_cubic_box(lower,upper)
  IMPLICIT NONE
  REAL(KIND=STORAGE_PRECISION),INTENT(IN) :: lower,upper
   box_dimensions(1,:)=lower
   box_dimensions(2,:)=upper
   !initialise box size
   box_size(:)=box_dimensions(2,:)-box_dimensions(1,:)
   maximum_distance_squared=box_size(2)**2+SQRT(box_size(1)**2+box_size(3)**2)
   maximum_distance=SQRT(maximum_distance_squared)
   BOX_VOLUME_GIVEN=.TRUE.
  END SUBROUTINE set_cubic_box

  !The following subroutine computes the squared radius of gyration (rgy_sq) of a given molecule,
  !as well as the furthest distance of any atom from the centre of mass of the molecule.
  SUBROUTINE compute_squared_radius_of_gyration(timestep,molecule_type_index,molecule_index,rgy_sq,maxdist)
  IMPLICIT NONE
  REAL(KIND=WORKING_PRECISION),INTENT(OUT) :: rgy_sq,maxdist
  REAL(KIND=WORKING_PRECISION) :: centre_of_mass(3),difference_vector(3),maxdist_squared,current_distance_squared
  INTEGER, INTENT(IN) :: timestep,molecule_type_index,molecule_index
  INTEGER :: atom_index
   IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
   maxdist_squared=0.0d0
   rgy_sq=0.0d0
   !centre of mass has to be computed first... because required for the distance
   centre_of_mass(:)=give_center_of_mass(timestep,molecule_type_index,molecule_index)
   DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
    IF (READ_SEQUENTIAL) THEN
     difference_vector(:)=DBLE(molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:))
    ELSE
     difference_vector(:)=DBLE(molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:))
    ENDIF
    difference_vector(:)=difference_vector(:)-centre_of_mass(:)
    current_distance_squared=SUM(difference_vector(:)**2)
    !check if that's a new record, amend if necessary
    IF (current_distance_squared>maxdist_squared) maxdist_squared=current_distance_squared
    !then, this *squared* position is weighted with the atom's mass
    rgy_sq=rgy_sq+(molecule_list(molecule_type_index)%list_of_atom_masses(atom_index))*current_distance_squared
   ENDDO
   rgy_sq=rgy_sq/molecule_list(molecule_type_index)%mass
   maxdist=SQRT(maxdist_squared)
  END SUBROUTINE compute_squared_radius_of_gyration

  !The following subroutine computes equation (13), (14), and (15) in 10.1021/acs.jpclett.9b02983.
  !These are assigned as TCM, TR, and TD, respectively
  SUBROUTINE compute_drude_temperature(timestep,TCM,TR,TD,Nf,Nmol,ND)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep
  INTEGER,INTENT(OUT) :: Nf,Nmol,ND
  REAL(KIND=WORKING_PRECISION),INTENT(OUT) :: TCM,TR,TD
   CALL compute_TCM() !compute_TCM also computes Nmol
   CALL compute_TR()  !compute_TR also computes Nf
   CALL compute_TD()  !compute_TD also computes ND
  
  CONTAINS

   SUBROUTINE compute_TCM()
   IMPLICIT NONE
   REAL(KIND=WORKING_PRECISION) :: mass_clipboard
   INTEGER :: molecule_type_index,molecule_index
    !initialise output variables
    Nmol=give_total_number_of_molecules_per_step()
    TCM=0.0d0
    !calculate the temperature, iterate over every molecule:
    DO molecule_type_index=1,number_of_molecule_types,1
     mass_clipboard=molecule_list(molecule_type_index)%mass
     DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory    
      !sum up all the mv²
      TCM=TCM+mass_clipboard*SUM((give_center_of_mass(timestep,molecule_type_index,molecule_index))**2)
     ENDDO
    ENDDO
    !divide temperatures by boltzmann constant as well as degrees of freedom.
    !For TCM, the degrees of freedom are 3*Nmol
    TCM=TCM/(3.0d0*Nmol*boltzmann)
    !boltzmann was given in J/K. m*v² is in (g*angstroms²)/(mol*fs²). Change that.
    TCM=1.0d7*TCM/avogadro
   END SUBROUTINE compute_TCM

   SUBROUTINE compute_TR()
   IMPLICIT NONE
   REAL(KIND=WORKING_PRECISION) :: core_velocity(3),drude_velocity(3),mass_clipboard,core_mass,drude_mass
   INTEGER :: drude_flag,molecule_type_index,molecule_index,atom_index,nmol_per_step
    !initialise output variables
    TR=0.0d0
    Nf=give_total_degrees_of_freedom()
    IF (READ_SEQUENTIAL) CALL goto_timestep(timestep)
    DO molecule_type_index=1,number_of_molecule_types,1
     nmol_per_step=molecule_list(molecule_type_index)%total_molecule_count
     DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
      !for a certain atom in this molecule type index, check whether it is NOT a drude particle.
      drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
      IF (drude_flag/=0) THEN
       !This atom contributes to TR.
       IF (drude_flag==-1) THEN
        !It is a non-polarisable atom!
        mass_clipboard=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
        DO molecule_index=1,nmol_per_step,1
         IF (READ_SEQUENTIAL) THEN
          core_velocity(:)=&
          &molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
         ELSE
          core_velocity(:)=&
          &molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
         ENDIF
         !calculate the temperature from the velocity
         TR=TR+mass_clipboard*SUM((core_velocity(:))**2)
        ENDDO
       ELSE
        !it is a drude core --> get centre of mass velocity
        core_mass=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
        drude_mass=molecule_list(molecule_type_index)%list_of_atom_masses(drude_flag)
        mass_clipboard=core_mass+drude_mass
        DO molecule_index=1,nmol_per_step,1
         IF (READ_SEQUENTIAL) THEN
          core_velocity(:)=&
          &molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
          drude_velocity(:)=&
          &molecule_list(molecule_type_index)%snapshot(drude_flag,molecule_index)%coordinates(:)
         ELSE
          core_velocity(:)=&
          &molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
          drude_velocity(:)=&
          &molecule_list(molecule_type_index)%trajectory(drude_flag,molecule_index,timestep)%coordinates(:)
         ENDIF
         core_velocity(:)=(core_mass*core_velocity(:)+drude_mass*drude_velocity(:))/mass_clipboard
         !calculate the temperature from the velocity
         TR=TR+mass_clipboard*SUM((core_velocity(:))**2)
        ENDDO
       ENDIF
      ENDIF
     ENDDO
    ENDDO
    !divide temperatures by boltzmann constant as well as degrees of freedom.
    !For TR, the degrees of freedom are Nf
    TR=TR/(Nf*boltzmann)
    !boltzmann was given in J/K. m*v² is in (g*angstroms²)/(mol*fs²). Change that.
    TR=1.0d7*TR/avogadro
   END SUBROUTINE compute_TR

   SUBROUTINE compute_TD()
   IMPLICIT NONE
   REAL(KIND=WORKING_PRECISION) :: relative_drude_velocity(3)
   INTEGER :: drude_flag,molecule_type_index,molecule_index,atom_index,nmol_per_step
    !initialise output variables
    TD=0.0d0
    ND=0
    IF (READ_SEQUENTIAL) CALL goto_timestep(timestep)
    DO molecule_type_index=1,number_of_molecule_types,1
     DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
      !for a certain atom in this molecule type index, check whether it is a core with drude attached to it.
      drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
      IF (drude_flag>0) THEN
       !If so, compute the relative drude velocity and sum them.
       nmol_per_step=molecule_list(molecule_type_index)%total_molecule_count
       ND=ND+nmol_per_step
       DO molecule_index=1,nmol_per_step,1
        IF (READ_SEQUENTIAL) THEN
         relative_drude_velocity(:)=&
         &molecule_list(molecule_type_index)%snapshot(drude_flag,molecule_index)%coordinates(:)-&
         &molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
        ELSE
         relative_drude_velocity(:)=&
         &molecule_list(molecule_type_index)%trajectory(drude_flag,molecule_index,timestep)%coordinates(:)-&
         &molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
        ENDIF
        TD=TD+molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%reduced_mass*&
        &SUM((relative_drude_velocity(:))**2)
       ENDDO
      ENDIF
     ENDDO
    ENDDO
    !divide temperatures by boltzmann constant as well as degrees of freedom.
    !For TD, the degrees of freedom are 3*ND
    TD=TD/(3.0d0*ND*boltzmann)
    !boltzmann was given in J/K. m*v² is in (g*angstroms²)/(mol*fs²). Change that.
    TD=1.0d7*TD/avogadro
   END SUBROUTINE compute_TD

  END SUBROUTINE compute_drude_temperature

  SUBROUTINE allocate_drude_list()
  IMPLICIT NONE
  INTEGER :: allocstatus,m,molecule_type_index
   IF (drudes_allocated) RETURN
   !allocate memory for detailed drude pair list
   DO molecule_type_index=1,number_of_molecule_types,1
    ALLOCATE(molecule_list(molecule_type_index)%&
    &list_of_drude_pairs(molecule_list(molecule_type_index)%number_of_atoms),STAT=allocstatus)
    IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
    molecule_list(molecule_type_index)%number_of_drudes_in_molecule=0
    DO m=1,molecule_list(molecule_type_index)%number_of_atoms,1
     molecule_list(molecule_type_index)%list_of_drude_pairs(m)%drude_flag=-1
     molecule_list(molecule_type_index)%list_of_drude_pairs(m)%reduced_mass=0.0d0
     molecule_list(molecule_type_index)%list_of_drude_pairs(m)%minimum_drude_distance=maximum_distance
     molecule_list(molecule_type_index)%list_of_drude_pairs(m)%maximum_drude_distance=0.0d0
    ENDDO
   ENDDO
   drudes_allocated=.TRUE.
  END SUBROUTINE allocate_drude_list

  SUBROUTINE assign_drudes()
  IMPLICIT NONE
  INTEGER :: molecule_type_index,current_atom_index,atom_index,atom_index_observed,drude_flag
  REAL(KIND=STORAGE_PRECISION) :: smallest_distance_found,current_distance
   drude_details=.FALSE.
   IF (drudes_assigned) THEN
    !drude particles are already assigned, but the reduced masses need to be calculated.
    IF (number_of_drude_particles/=0) THEN
     DO molecule_type_index=1,number_of_molecule_types,1
      DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
       drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
       IF (drude_flag>0) CALL compute_drude_properties(molecule_type_index,drude_flag,atom_index,skip_position=.TRUE.)
      ENDDO
     ENDDO
    ENDIF
    RETURN
   ENDIF
   IF (INFORMATION_IN_TRAJECTORY=="VEL") THEN
    PRINT *,"Cannot assign drude particles from velocity information."
    RETURN
   ENDIF
   PRINT *,"Assigning drudes from trajectory file."
   IF (VERBOSE_OUTPUT) THEN
    PRINT *,"Use keyword 'show_drude' to print detailed information."
    PRINT *,"(including the appropriate input section for the molecular input file)"
   ENDIF
   IF (READ_SEQUENTIAL) CALL goto_timestep(1)
   CALL allocate_drude_list()
   !Fill drude pair list:
   ! - drude particles will take the value 0
   ! - drude cores will have the index of their respective drude
   ! - non-polarisable atoms remain -1.
   DO molecule_type_index=1,number_of_molecule_types,1 !iterate over all molecule types.
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !for each molecule type, look at the atoms one by one and assign their values.
     !checking the atom represented by atom_index:
     IF ((TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))=="X").OR.&
     &(TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))=="D")) THEN
      !this atom is a drude particle. First, change flag to '0'
      molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag=0
      !THEN, check which atom this drude belongs to by searching for the closest distance.
      smallest_distance_found=maximum_distance
      !initialise current_atom_index to -1, if it stays like that then there was no core available!
      current_atom_index=-1
      DO atom_index_observed=1,molecule_list(molecule_type_index)%number_of_atoms,1
       IF (((TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index_observed))/="X").AND.&
       &(TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index_observed))/="D")).AND.&
       &(molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_observed)%drude_flag==-1)) THEN
        !The observed atom is not a drude particle, and has not yet been assigned a drude particle.
        ! No need to check for (atom_index_observed/=atom_index) because of /= "X"
        !check the first molecule in the first timestep.
        current_distance=give_smallest_atom_distance&
        &(1,1,molecule_type_index,molecule_type_index,1,1,atom_index,atom_index_observed)
        IF (current_distance<smallest_distance_found) THEN
         !new 'record', i.e. a new smallest distance encountered.
         smallest_distance_found=current_distance
         current_atom_index=atom_index_observed
        ENDIF
       ENDIF
      ENDDO
      !current_atom_index should now store the appropriate drude core for the drude with index 'atom_index'.
      !three checks will be done:
      ! check 1) Almost trivial: check if current_atom_index is -1! This means that this drude particle will be ignored.
      IF (current_atom_index==-1) THEN
       CALL report_error(86,exit_status=atom_index)
      ELSE
       !increase counter for successfully assigned drude pairs
       molecule_list(molecule_type_index)%number_of_drudes_in_molecule=&
       &molecule_list(molecule_type_index)%number_of_drudes_in_molecule+1
       ! check 2) It is technically not necessary to check for whether the core is already assigned. Thus, error 0.
       !          The same applies to drude-to-drude assignments - they just shouldn't happen, because of the conditions above.
       IF (molecule_list(molecule_type_index)%list_of_drude_pairs(current_atom_index)%drude_flag/=-1) THEN
        CALL report_error(0)
       ENDIF
       !Change drude flag *OF THE CORE ATOM* accordingly.
       molecule_list(molecule_type_index)%list_of_drude_pairs(current_atom_index)%drude_flag=atom_index
       ! check 3) get the smallest and highest distance between this drude and its alleged core (globally), as well as their reduced mass.
       !The following subroutine is called for every pair of drude (atom_index) and core (current_atom_index) found.
       CALL compute_drude_properties(molecule_type_index,atom_index,current_atom_index)
      ENDIF
     ENDIF
    ENDDO
   ENDDO
   drudes_assigned=.TRUE.
  END SUBROUTINE assign_drudes

  !The following subroutine is called for every found pair of drude and core particle.
  SUBROUTINE compute_drude_properties(molecule_type_index,atom_index_drude,atom_index_core,skip_position)
  IMPLICIT NONE
  REAL(KIND=GENERAL_PRECISION) :: mass_a,mass_b,current_distance
  INTEGER :: molecule_index
  INTEGER,INTENT(IN) :: molecule_type_index,atom_index_drude,atom_index_core
  LOGICAL,INTENT(IN),OPTIONAL :: skip_position
   !Compute reduced mass
   mass_a=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index_drude)
   mass_b=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index_core)
   molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%reduced_mass=&
   &(mass_a*mass_b)/(mass_a+mass_b)
   IF (PRESENT(skip_position)) THEN
    IF (skip_position) RETURN
   ENDIF
   IF (INFORMATION_IN_TRAJECTORY=="POS") THEN
    !Compute smallest distance and highest distance in the first step
    DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
     current_distance=give_smallest_atom_distance&
     &(1,1,molecule_type_index,molecule_type_index,molecule_index,molecule_index,atom_index_drude,atom_index_core)
     IF (current_distance<molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%minimum_drude_distance) THEN
      molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%minimum_drude_distance=current_distance
     ENDIF
     !use two separate IF's - there could be only one molecule!
     IF (current_distance>molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%maximum_drude_distance) THEN
      molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%maximum_drude_distance=current_distance
     ENDIF
    ENDDO
   ELSE
    molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%maximum_drude_distance=0.0d0
    molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index_core)%minimum_drude_distance=0.0d0
   ENDIF
   drude_details=.TRUE.
  END SUBROUTINE compute_drude_properties

  INTEGER FUNCTION give_number_of_drude_particles()
  IMPLICIT NONE
   give_number_of_drude_particles=number_of_drude_particles
  END FUNCTION give_number_of_drude_particles

  LOGICAL FUNCTION give_comboost()
  IMPLICIT NONE
   give_comboost=use_firstatom_as_com
  END FUNCTION give_comboost

  REAL(KIND=GENERAL_PRECISION) FUNCTION give_smallest_distance&
  &(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2
   !pass on to the function that computes the squared distance (less SQRTs to evaluate)
   give_smallest_distance=SQRT(give_smallest_distance_squared&
   &(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2))
  END FUNCTION give_smallest_distance

  !This subroutine computes the smallest and largest distance within all molecules of the given type and timestep.
  SUBROUTINE give_intramolecular_distances(timestep,molecule_type_index,smallest,largest)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep,molecule_type_index
  REAL,INTENT(OUT) :: smallest,largest
  REAL :: current_squared
  INTEGER :: molecule_index,atom_index_1,atom_index_2
   !two atoms can be no further apart than the diagonale of the box... that's what I initialise to
   smallest=maximum_distance_squared
   largest=0.0
   DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
    DO atom_index_1=1,molecule_list(molecule_type_index)%number_of_atoms-1,1
     DO atom_index_2=atom_index_1+1,molecule_list(molecule_type_index)%number_of_atoms,1
      current_squared=give_smallest_atom_distance_squared&
      &(timestep,timestep,molecule_type_index,molecule_type_index,molecule_index,molecule_index,atom_index_1,atom_index_2)
      IF (current_squared>largest) largest=current_squared
      IF (current_squared<smallest) smallest=current_squared
     ENDDO
    ENDDO
   ENDDO
   smallest=SQRT(smallest)
   largest=SQRT(largest)
  END SUBROUTINE give_intramolecular_distances

  !This subroutine computes the smallest intermolecular distance of the given type and timestep to any other molecule.
  SUBROUTINE give_intermolecular_contact_distance(timestep,molecule_type_index_1,smallest)
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_thread_num()
  !$  INTEGER :: OMP_get_thread_num
  !$  END FUNCTION OMP_get_thread_num
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$ END INTERFACE
  INTEGER,INTENT(IN) :: timestep,molecule_type_index_1
  REAL,INTENT(OUT) :: smallest
  INTEGER :: molecule_index_1
  REAL :: smallest_local
   !two atoms can be no further apart than the diagonale of the box... that's what I initialise to
   !$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(smallest_local)
   !$OMP SINGLE
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE (*,'(A,I0,A)') "     ### Parallel execution on ",OMP_get_num_threads()," threads (intermolecular contact distance)"
   !$  CALL timing_parallel_sections(.TRUE.)
   !$ ENDIF
   !$OMP END SINGLE
   smallest=maximum_distance_squared
   smallest_local=maximum_distance_squared
   !outer loop: goes over all the molecules of this molecule type...
   !$OMP DO
   DO molecule_index_1=1,molecule_list(molecule_type_index_1)%total_molecule_count,1
    CALL particular_molecule_contact_distance(molecule_index_1,smallest_local)
   ENDDO
   !$OMP END DO
   !$ IF (DEVELOPERS_VERSION) WRITE(*,'("  ! thread ",I0,", smallest value: ",F0.3,A)') OMP_get_thread_num(),SQRT(smallest_local)
   !CRITICAL directive to properly update the final value
   !$OMP CRITICAL
   IF (smallest_local<smallest) smallest=smallest_local
   !$OMP END CRITICAL
   !$OMP END PARALLEL
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE(*,ADVANCE="NO",FMT='("     ### End of parallelised section, took ")')
   !$  CALL timing_parallel_sections(.FALSE.)
   !$ ENDIF
   smallest=SQRT(smallest)
   CONTAINS

    SUBROUTINE particular_molecule_contact_distance(molecule_index,smallest_inout)
    IMPLICIT NONE
    INTEGER :: molecule_index,molecule_index_2,atom_index_1,atom_index_2,molecule_type_index_2
    REAL :: current_squared
    REAL,INTENT(INOUT) :: smallest_inout
     !and in the inner loop, over all the atoms of this particular molecule.
     DO atom_index_1=1,molecule_list(molecule_type_index_1)%number_of_atoms,1
      !thus here, molecule_index and atom_index_1 will take every possible value for any atom in molecule_type_index_1!
      !For each of these atoms, check all molecules *other* than the current one:
      DO molecule_type_index_2=1,number_of_molecule_types,1
       !This loop iterates over all molecules of the other type...
       DO molecule_index_2=1,molecule_list(molecule_type_index_2)%total_molecule_count,1
        IF ((molecule_type_index_2==molecule_type_index_1).AND.(molecule_index_2==molecule_index)) THEN
         !they are the same! abort!
         EXIT
        ELSE
         !now, finally, go through the atoms of that second molecule.
         DO atom_index_2=1,molecule_list(molecule_type_index_2)%number_of_atoms,1
          current_squared=give_smallest_atom_distance_squared&
          &(timestep,timestep,molecule_type_index_1,molecule_type_index_2,&
          &molecule_index,molecule_index_2,atom_index_1,atom_index_2)
          IF (current_squared<smallest_inout) smallest_inout=current_squared
         ENDDO  
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    END SUBROUTINE particular_molecule_contact_distance

  END SUBROUTINE give_intermolecular_contact_distance

  REAL(KIND=GENERAL_PRECISION) FUNCTION give_smallest_atom_distance&
  &(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2,atom_index_1,atom_index_2)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep1,timestep2,molecule_type_index_1,molecule_type_index_2
  INTEGER,INTENT(IN) :: molecule_index_1,molecule_index_2,atom_index_1,atom_index_2
   !pass on to the function that computes the squared distance (less SQRTs to evaluate)
   give_smallest_atom_distance=SQRT(give_smallest_atom_distance_squared&
   &(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2,atom_index_1,atom_index_2))
  END FUNCTION give_smallest_atom_distance

  !This FUNCTION returns the smallest squared distance of 2 atoms considering all PBCs - as well as the corresponding translation vector.
  REAL(KIND=GENERAL_PRECISION) FUNCTION give_smallest_atom_distance_squared&
  &(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2,atom_index_1,atom_index_2)
  IMPLICIT NONE
  INTEGER :: a,b,c
  INTEGER,INTENT(IN) :: timestep1,timestep2,molecule_type_index_1,molecule_type_index_2
  INTEGER,INTENT(IN) :: molecule_index_1,molecule_index_2,atom_index_1,atom_index_2
  REAL(KIND=WORKING_PRECISION) :: pos_1(3),pos_2(3),shift(3),distance_clip
   IF (READ_SEQUENTIAL) THEN
    CALL goto_timestep(timestep1)
    pos_1(:)=DBLE(molecule_list(molecule_type_index_1)%snapshot(atom_index_1,molecule_index_1)%coordinates(:))
    CALL goto_timestep(timestep2)
    pos_2(:)=DBLE(molecule_list(molecule_type_index_2)%snapshot(atom_index_2,molecule_index_2)%coordinates(:))
   ELSE
    pos_1(:)=DBLE(molecule_list(molecule_type_index_1)%trajectory(atom_index_1,molecule_index_1,timestep1)%coordinates(:))
    pos_2(:)=DBLE(molecule_list(molecule_type_index_2)%trajectory(atom_index_2,molecule_index_2,timestep2)%coordinates(:))
   ENDIF
   !two atoms can be no further apart than the diagonale of the box... that's what I initialise to
   give_smallest_atom_distance_squared=maximum_distance_squared
   IF (.NOT.(WRAP_TRAJECTORY)) THEN
    CALL wrap_vector(pos_1)
    CALL wrap_vector(pos_2)
   ENDIF
   !Now, check all mirror images
   DO a=-1,1,1! a takes the values (-1, 0, 1)
    DO b=-1,1,1! b takes the values (-1, 0, 1)
     DO c=-1,1,1! c takes the values (-1, 0, 1)
      shift(1)=FLOAT(a)
      shift(2)=FLOAT(b)
      shift(3)=FLOAT(c)
      shift=shift*box_size
      !shift is now the translation vector to the mirror image.
      distance_clip=SUM(((pos_2+shift)-pos_1)**2)
      IF (distance_clip<give_smallest_atom_distance_squared) THEN
       !a distance has been found that's closer than the current best - amend that.
       give_smallest_atom_distance_squared=distance_clip
      ENDIF
       ENDDO
    ENDDO
   ENDDO
  END FUNCTION give_smallest_atom_distance_squared

  !This FUNCTION returns the smallest squared distance of centres of mass considering all PBCs - as well as the corresponding translation vector.
  REAL(KIND=GENERAL_PRECISION) FUNCTION give_smallest_distance_squared&
  &(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2,translation)
  IMPLICIT NONE
  INTEGER :: a,b,c
  INTEGER,INTENT(IN) :: timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,molecule_index_1,molecule_index_2
  REAL(KIND=WORKING_PRECISION) :: pos_1(3),pos_2(3),shift(3),distance_clip,wrapshift(3)
  REAL(KIND=WORKING_PRECISION),INTENT(OUT),OPTIONAL :: translation(3)
   pos_1(:)=give_center_of_mass(timestep1,molecule_type_index_1,molecule_index_1)
   pos_2(:)=give_center_of_mass(timestep2,molecule_type_index_2,molecule_index_2)
   !two molecules can be no further apart than the diagonale of the box... that's what I initialise to
   give_smallest_distance_squared=maximum_distance_squared
   IF (.NOT.(WRAP_TRAJECTORY)) THEN
    CALL wrap_vector(pos_1,shift(:))
    CALL wrap_vector(pos_2,wrapshift(:))
    !"shift" is now the vector to translate the reference molecule into the box by wrapping.
    !"wrapshift" is the same for the second, observed molecule.
    !now, store in wrapshift the vector to bring the second *into the same box* as the first molecule:
    wrapshift(:)=wrapshift(:)-shift(:)
   ENDIF
   !Now, check all mirror images
   DO a=-1,1,1! a takes the values (-1, 0, 1)
    DO b=-1,1,1! b takes the values (-1, 0, 1)
     DO c=-1,1,1! c takes the values (-1, 0, 1)
      shift(1)=FLOAT(a)
      shift(2)=FLOAT(b)
      shift(3)=FLOAT(c)
      shift=shift*box_size
      !shift is now the translation vector to the mirror image.
      distance_clip=SUM(((pos_2+shift)-pos_1)**2)
      IF (distance_clip<give_smallest_distance_squared) THEN
       !a distance has been found that's closer than the current best - amend that.
       give_smallest_distance_squared=distance_clip
       IF (PRESENT(translation)) translation(:)=shift(:)
      ENDIF
       ENDDO
    ENDDO
   ENDDO
   IF (PRESENT(translation).AND.(.NOT.(WRAP_TRAJECTORY))) translation(:)=translation(:)+wrapshift(:)
  END FUNCTION give_smallest_distance_squared

  SUBROUTINE give_number_of_neighbours&
  &(timestep,molecule_type_index,molecule_index,neighbour_molecules,neighbour_atoms,cutoff,output_unit)
  IMPLICIT NONE
  INTEGER,INTENT(OUT),OPTIONAL :: neighbour_atoms,neighbour_molecules
  REAL(KIND=WORKING_PRECISION) :: cutoff_squared,distance_squared,shift(3)
  REAL(KIND=WORKING_PRECISION),INTENT(IN) :: cutoff
  INTEGER,INTENT(IN),OPTIONAL :: output_unit !this unit (if specified) will be filled with the neighbours
  INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
  INTEGER :: molecule_type_index_counter,molecule_index_counter,natoms
   cutoff_squared=cutoff**2
   IF (PRESENT(neighbour_atoms)) neighbour_atoms=0
   IF (PRESENT(neighbour_molecules)) neighbour_molecules=0
   DO molecule_type_index_counter=1,number_of_molecule_types,1
    !How many atoms are in this molecule?
    natoms=molecule_list(molecule_type_index_counter)%number_of_atoms
    DO molecule_index_counter=1,molecule_list(molecule_type_index_counter)%total_molecule_count,1 !gives dimension 2 of trajectory    
     !get the smallest distance (as its square, because who needs SQRT?)
     distance_squared=give_smallest_distance_squared&
     &(timestep,timestep,molecule_type_index,molecule_type_index_counter,molecule_index,molecule_index_counter,shift)
     IF (distance_squared<=cutoff_squared) THEN
      IF ((molecule_type_index/=molecule_type_index_counter).OR.(molecule_index/=molecule_index_counter)) THEN
       !update output variables
       IF (PRESENT(neighbour_atoms)) neighbour_atoms=neighbour_atoms+natoms
       IF (PRESENT(neighbour_molecules)) neighbour_molecules=neighbour_molecules+1
       !save the encountered atoms in the appropriate unit
       IF (PRESENT(output_unit)) THEN
        CALL write_molecule&
        &(output_unit,timestep,molecule_type_index_counter,molecule_index_counter,include_header=.FALSE.,translate_by=shift)
       ENDIF
      ENDIF
     ENDIF
    ENDDO
   ENDDO
  END SUBROUTINE give_number_of_neighbours

  !The task of convert_parallel is to reduce a trajectory to centre-of-mass. Not part of the official version. Doesn't really work, because File I/O is the slow part.
  SUBROUTINE convert_parallel()
  IMPLICIT NONE
  INTEGER :: queuesize,first_in,last_in,stepcounter,nprocs,number_of_molecules
  INTEGER,DIMENSION(:),ALLOCATABLE :: element_status !0=to convert, 1=converting currently, 2=converted, 3=written
  CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE :: list_of_elements_output
  REAL(KIND=WORKING_PRECISION),DIMENSION(:,:,:),ALLOCATABLE :: coordinates_output
  LOGICAL :: completed,skip
  !$ INTERFACE
  !$  FUNCTION OMP_get_thread_num()
  !$  INTEGER :: OMP_get_thread_num
  !$  END FUNCTION OMP_get_thread_num
  !$  FUNCTION OMP_get_max_threads()
  !$  INTEGER :: OMP_get_max_threads
  !$  END FUNCTION OMP_get_max_threads
  !$  SUBROUTINE OMP_set_nested(enable)
  !$  LOGICAL,INTENT(IN) :: enable
  !$  END SUBROUTINE OMP_set_nested
  !$ END INTERFACE
   !$ nprocs=OMP_get_max_threads()
   !$ IF (nprocs<3) THEN
   !$  WRITE(*,*) " ! less than 3 threads available. Returning control to main module."
   !$  RETURN
   !$ ENDIF
   !$ CALL OMP_set_num_threads(nprocs)
   !$ WRITE(*,'("  ! number of threads set to ",I0)') nprocs
   !omp: if false then return.
   !$ IF (.FALSE.) THEN
    WRITE(*,*) " ! -fopenmp flag not set. Returning control to main module."
    RETURN
   !$ ENDIF
   first_in=0
   last_in=0
   completed=.FALSE.
   queuesize=100
   CALL initialise_queue()
   REWIND 9
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(stepcounter,skip)
   !$OMP SECTIONS
   !$OMP SECTION
   !$ WRITE(*,'("  ! Thread number ",I0," is reading the trajectory file.")') OMP_get_thread_num()
   !This thread is reading the trajectory step by step, writing into the queue, and marking as 'to convert'.
   DO stepcounter=1,number_of_steps,1
    CALL enqueue()
   ENDDO
   !$ WRITE(*,'("  ! Thread number ",I0," is done reading the trajectory.")') OMP_get_thread_num()
   !$OMP SECTION
   !$ WRITE(*,'("  ! Thread number ",I0," writes to output trajectory.")') OMP_get_thread_num()
   DO stepcounter=1,number_of_steps,1
   1 IF (element_status((MOD(stepcounter-1,queuesize)+1))==2) THEN
     !write this step.
     CALL write_from_queue(stepcounter*TIME_SCALING_FACTOR,(MOD(stepcounter-1,queuesize)+1))
     !set flag to 'written'.
     element_status(MOD(stepcounter-1,queuesize)+1)=3
    ELSE
     GOTO 1
    ENDIF
   ENDDO
   !$ WRITE(*,'("  ! Thread number ",I0," is done writing to output trajectory.")') OMP_get_thread_num()
   completed=.TRUE.
   !$OMP END SECTIONS NOWAIT
   IF (.NOT.completed) THEN
   !$ WRITE(*,'("  ! Thread ",I0," is ready to convert.")') OMP_get_thread_num()
   !These threads are constantly converting everything that needs converting.
    DO stepcounter=1,number_of_steps,1
    2 skip=.FALSE.
     IF (completed) THEN
     !$  WRITE(*,'("  ! Thread number ",I0," caught termination signal.")') OMP_get_thread_num()
      EXIT
     ENDIF
     !$OMP CRITICAL
     !Check if this step is to be converted, and if yes, change flag to '1' and start converting.
     IF (element_status(MOD(stepcounter-1,queuesize)+1)==0) THEN
      element_status(MOD(stepcounter-1,queuesize)+1)=1
     ELSE
      skip=.TRUE.
     ENDIF
     !$OMP END CRITICAL
     IF (skip) GOTO 2
     CALL convert_queue(MOD(stepcounter-1,queuesize)+1)
     !set flag to 'converted'
     element_status(MOD(stepcounter-1,queuesize)+1)=2
    ENDDO
   ENDIF
   !$OMP END PARALLEL
   CALL finalise_queue()
   CALL reset_trajectory_file()
   CONTAINS

    SUBROUTINE initialise_queue()
    IMPLICIT NONE
    INTEGER :: allocstatus,n,m,counter
    LOGICAL :: connected
    CHARACTER(LEN=1024) :: fstring
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     WRITE(fstring,'(2A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"COM_parallel.lmp"
     OPEN(UNIT=3,FILE=TRIM(fstring))
     number_of_molecules=0
     ALLOCATE(element_status(queuesize),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     !converted elements must be initialised to 'written', to keep every thread but the read thread idle at the beginning.
     element_status(:)=3
     DO n=1,number_of_molecule_types,1
      number_of_molecules=number_of_molecules+molecule_list(n)%total_molecule_count
      ALLOCATE(molecule_list(n)%queue(queuesize,molecule_list(n)%total_molecule_count,molecule_list(n)%number_of_atoms)&
      &,STAT=allocstatus)
      IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     ENDDO
     ALLOCATE(list_of_elements_output(number_of_molecules),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     ALLOCATE(coordinates_output(queuesize,number_of_molecules,3),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     counter=0
     DO n=1,number_of_molecule_types,1
      DO m=1,molecule_list(n)%total_molecule_count,1
       counter=counter+1
       list_of_elements_output(counter)=CHAR(ALPHABET_small(MOD((n-1),26)+1))
      ENDDO
     ENDDO
    END SUBROUTINE initialise_queue

    SUBROUTINE enqueue() !reads a step and puts it in the queue
    IMPLICIT NONE
    INTEGER :: dummy_iterations,molecule_type_index,atom_index,molecule_index,counter
    CHARACTER(LEN=2) :: dummystring
     dummy_iterations=0
     last_in=MOD(last_in,queuesize)+1
     DO
      IF ((first_in/=last_in).AND.(element_status(last_in)==3)) THEN
       !read/skip the header
       DO counter=1,headerlines_to_skip,1
        READ(9,*)
       ENDDO
       !read the body into queue(last_in)
       DO molecule_type_index=1,number_of_molecule_types,1
        !For each molecule type, read the corresponding number of molecules:
        DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of queue
         !Finally, iterate over the atoms in that particular molecule:
         DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives third dimension of queue
          !LOOP VARIABLES:
          !molecule_type_index: current molecule type, e.g. 1 (only 1 and 2 for a binary IL)
          !molecule_index: current explicit molecule, e.g. molecule number 231
          !atom_index: current atom in that molecule, e.g. atom 2 (in molecule number 231 of type 1 in timestep 1234...)
          READ(9,*) dummystring,&
          &molecule_list(molecule_type_index)%queue(last_in,molecule_index,atom_index)%coordinates
         ENDDO
        ENDDO 
       ENDDO
       !set flag to convert. Only one thread can change back to zero - this one.
       element_status(last_in)=0
       EXIT
      ELSE
       dummy_iterations=dummy_iterations+1
      ENDIF
     ENDDO
     IF (dummy_iterations>0) PRINT *,"DUMMY ITERATIONS ",dummy_iterations
    END SUBROUTINE enqueue

    SUBROUTINE convert_queue(position_input) !convert queue to coordinates_output, which is centre of mass
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: weighted_pos(3)
    INTEGER :: counter,molecule_type_index,molecule_index,atom_index,position_input
     coordinates_output(position_input,:,:)=0.0d0
     counter=1
     DO molecule_type_index=1,number_of_molecule_types,1
      DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of queue
       DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives third dimension of queue
        !molecule_type_index: current molecule type, e.g. 1 (only 1 and 2 for a binary IL)
        !molecule_index: current explicit molecule, e.g. molecule number 231
        !atom_index: current atom in that molecule, e.g. atom 2 (in molecule number 231 of type 1 in timestep 1234...)
        !store current atom position in weighted_pos
        weighted_pos=DBLE(molecule_list(molecule_type_index)%queue(position_input,molecule_index,atom_index)%coordinates(:))
        !weigh position with mass...
        weighted_pos(:)=weighted_pos(:)*DBLE(molecule_list(molecule_type_index)%list_of_atom_masses(atom_index))
        !... and add to centre of mass.
        coordinates_output(position_input,counter,:)=coordinates_output(position_input,counter,:)+weighted_pos(:)
       ENDDO
       !normalise, increase counter for the next molecule.
       coordinates_output(position_input,counter,:)=&
       &coordinates_output(position_input,counter,:)/DBLE(molecule_list(molecule_type_index)%mass)
       counter=counter+1
      ENDDO 
     ENDDO
    END SUBROUTINE convert_queue

    SUBROUTINE write_from_queue(step_number,position_input)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: step_number,position_input
    INTEGER :: n
     !WRITE header
     WRITE(3,'("ITEM: TIMESTEP")')
     WRITE(3,'(I0)') step_number
     WRITE(3,'("ITEM: NUMBER OF ATOMS")')
     WRITE(3,'(I0)') number_of_molecules
     WRITE(3,'("ITEM: BOX BOUNDS pp pp pp")')
     WRITE(3,*) box_dimensions(:,1)
     WRITE(3,*) box_dimensions(:,2)
     WRITE(3,*) box_dimensions(:,3)
     !Append the line that tells the user about the content!
     SELECT CASE (INFORMATION_IN_TRAJECTORY)
     CASE("UNK")
      WRITE(3,'("ITEM: ATOMS element x? y? z?")')
     CASE("VEL")
      WRITE(3,'("ITEM: ATOMS element vx vy vz")')
     CASE("POS")
      WRITE(3,'("ITEM: ATOMS element xu yu zu")')
     CASE DEFAULT
      CALL report_error(0)
     END SELECT
     !WRITE body (centre of mass)
     DO n=1,number_of_molecules,1
      WRITE(3,'(A2,3E19.10)') list_of_elements_output(n),coordinates_output(position_input,n,:)
     ENDDO
    END SUBROUTINE write_from_queue

    SUBROUTINE finalise_queue()
    IMPLICIT NONE
    INTEGER :: deallocstatus,n
     CLOSE(UNIT=3)
     DEALLOCATE(element_status,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(list_of_elements_output,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(coordinates_output,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DO n=1,number_of_molecule_types,1
      DEALLOCATE(molecule_list(n)%queue,STAT=deallocstatus)
      IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     ENDDO
    END SUBROUTINE finalise_queue

  END SUBROUTINE convert_parallel

  !This subroutine rewinds the trajectory file and adjusts the 'pointer' accordingly. (Hell, I'm not starting with real pointers here)
  SUBROUTINE reset_trajectory_file()
  IMPLICIT NONE
   IF (READ_SEQUENTIAL) THEN
    REWIND 9
    file_position=0
    CALL goto_timestep(1)
   ENDIF
  END SUBROUTINE reset_trajectory_file

  !Wrap a single position so its inside the box - could be a centre of mass, for example.
  SUBROUTINE wrap_vector(input_vector,wrapshift)
  IMPLICIT NONE
  REAL(KIND=WORKING_PRECISION),INTENT(INOUT) :: input_vector(3)
  REAL(KIND=WORKING_PRECISION),INTENT(OUT),OPTIONAL :: wrapshift(3)
  INTEGER :: xyzcounter
   IF (PRESENT(wrapshift)) wrapshift(:)=0.0d0
   DO xyzcounter=1,3,1
    DO
     IF (input_vector(xyzcounter)<box_dimensions(1,xyzcounter)) THEN
      !input vector is outside of box (smaller)
      input_vector(xyzcounter)=input_vector(xyzcounter)+box_size(xyzcounter)
      IF (PRESENT(wrapshift)) wrapshift(xyzcounter)=wrapshift(xyzcounter)+box_size(xyzcounter)
      CYCLE
     ELSE
      IF (input_vector(xyzcounter)>box_dimensions(2,xyzcounter)) THEN
       !input vector is outside of box (bigger)
       input_vector(xyzcounter)=input_vector(xyzcounter)-box_size(xyzcounter)
       IF (PRESENT(wrapshift)) wrapshift(xyzcounter)=wrapshift(xyzcounter)-box_size(xyzcounter)
       CYCLE
      ELSE
       !input vector is inside box!
       EXIT
      ENDIF
     ENDIF
    ENDDO
   ENDDO
  END SUBROUTINE wrap_vector

  !wrap the molecules into the box - full trajectory.
  !This routine benefits from parallelisation.
  SUBROUTINE wrap_full()
  IMPLICIT NONE
  INTEGER :: molecule_type_index,molecule_index,stepcounter
  REAL :: centre_of_mass(3)
  INTEGER :: xyzcounter
  !$ INTERFACE
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$ END INTERFACE
   !$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(molecule_type_index,molecule_index,centre_of_mass,xyzcounter)
   !$OMP SINGLE
   !start the timer for the parallel section.
   !$ CALL timing_parallel_sections(.TRUE.)
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION))&
   !$ &WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (wrapping)"
   !$OMP END SINGLE
   !$OMP DO SCHEDULE(STATIC,1)
   DO stepcounter=1,number_of_steps,1
    DO molecule_type_index=1,number_of_molecule_types,1
     DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
      centre_of_mass=give_center_of_mass(stepcounter,molecule_type_index,molecule_index)
      DO xyzcounter=1,3,1
      !See? I got rid of it! (The jump)
       DO
        IF (centre_of_mass(xyzcounter)<=box_dimensions(1,xyzcounter)) THEN
         !centre of mass is outside of box (smaller)
         centre_of_mass(xyzcounter)=centre_of_mass(xyzcounter)+box_size(xyzcounter)
         molecule_list(molecule_type_index)%trajectory(:,molecule_index,stepcounter)%coordinates(xyzcounter)=&
         &molecule_list(molecule_type_index)%trajectory(:,molecule_index,stepcounter)%coordinates(xyzcounter)&
         &+box_size(xyzcounter)
         CYCLE
        ELSE
         IF (centre_of_mass(xyzcounter)>=box_dimensions(2,xyzcounter)) THEN
          !centre of mass is outside of box (bigger)
          centre_of_mass(xyzcounter)=centre_of_mass(xyzcounter)-box_size(xyzcounter)
          molecule_list(molecule_type_index)%trajectory(:,molecule_index,stepcounter)%coordinates(xyzcounter)=&
          &molecule_list(molecule_type_index)%trajectory(:,molecule_index,stepcounter)%coordinates(xyzcounter)&
          &-box_size(xyzcounter)
          CYCLE
         ELSE
          !centre of mass is inside box!
          EXIT
         ENDIF
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
   !$  CALL timing_parallel_sections(.FALSE.)
   !$ ENDIF
  END SUBROUTINE wrap_full

  !wrap the molecules into the box - just the current snapshot.
  SUBROUTINE wrap_snap()
  IMPLICIT NONE
  INTEGER :: molecule_type_index,molecule_index
  REAL :: centre_of_mass(3)
  INTEGER :: xyzcounter
   !Parallelisation is not available anyway...
   DO molecule_type_index=1,number_of_molecule_types,1
    DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
     !timestep has no meaning because snapshot.
     centre_of_mass=give_center_of_mass(file_position,molecule_type_index,molecule_index)
     DO xyzcounter=1,3,1
     !Apologies for the jump. We are all weak sometimes.
     10  IF (centre_of_mass(xyzcounter)<=box_dimensions(1,xyzcounter)) THEN
       centre_of_mass(xyzcounter)=centre_of_mass(xyzcounter)+box_size(xyzcounter)
       molecule_list(molecule_type_index)%snapshot(:,molecule_index)%coordinates(xyzcounter)=&
       &molecule_list(molecule_type_index)%snapshot(:,molecule_index)%coordinates(xyzcounter)+box_size(xyzcounter)
       GOTO 10
      ELSE
       IF (centre_of_mass(xyzcounter)>=box_dimensions(2,xyzcounter)) THEN
        centre_of_mass(xyzcounter)=centre_of_mass(xyzcounter)-box_size(xyzcounter)
        molecule_list(molecule_type_index)%snapshot(:,molecule_index)%coordinates(xyzcounter)=&
        &molecule_list(molecule_type_index)%snapshot(:,molecule_index)%coordinates(xyzcounter)-box_size(xyzcounter)
        GOTO 10
       ENDIF
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  END SUBROUTINE wrap_snap

  !The following subroutine calculates the instantaneous temperature for a given timestep and molecule type.
  !Temperatures are given resolved in x,y,z - 'corrected_temperature' is drift-corrected!
  SUBROUTINE give_temperature(timestep,drift,molecule_type_index,temperature,corrected_temperature,&
  &kinetic_energy,constraints_correction)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep,molecule_type_index
  INTEGER :: atom_index,degrees_of_freedom,molecule_index,internal_constraints
  REAL(KIND=WORKING_PRECISION),INTENT(OUT) :: drift(3),temperature(3),corrected_temperature(3)
  REAL(KIND=WORKING_PRECISION) :: TL(3),T0(3),T1(3),T2(3) !LAMMPS temperature, desired corrected temperature, first approximation, second approximation.
  REAL(KIND=WORKING_PRECISION) :: atom_clipboard(3),mass_clipboard
  REAL(KIND=WORKING_PRECISION),INTENT(OUT),OPTIONAL :: kinetic_energy,constraints_correction
   IF (READ_SEQUENTIAL) CALL goto_timestep(timestep)
   !initialise output variables
   drift(:)=0.0d0
   temperature(:)=0.0d0
   corrected_temperature(:)=0.0d0
   !The drift is obtained as the average of the center of mass, assuming that all molecules of one type have the same mass.
   DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory
    drift(:)=drift(:)+give_center_of_mass(timestep,molecule_type_index,molecule_index)
   ENDDO
   !Normalise the drift
   drift(:)=drift(:)/FLOAT(molecule_list(molecule_type_index)%total_molecule_count)
   !At this point, the drift is correct. Now, calculate the temperatures:
   DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory  
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives first dimension of trajectory
     !Get the current 'atom_clipboard', which usually would be the velocity.
     IF (READ_SEQUENTIAL) THEN
      atom_clipboard(:)=molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
     ELSE
      atom_clipboard(:)=molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
     ENDIF
     mass_clipboard=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
     !first, the normal, uncorrected temperature.
     temperature(:)=temperature(:)+mass_clipboard*(atom_clipboard(:)**2) !in every direction, T=m*v²
     !now, correct for the drift, and then compute the drift-corrected temperature from that.
     atom_clipboard(:)=atom_clipboard(:)-drift(:)
     corrected_temperature(:)=corrected_temperature(:)+mass_clipboard*(atom_clipboard(:)**2)
    ENDDO
   ENDDO
   IF (DEVELOPERS_VERSION) THEN
    TL(:)=0.0d0
    T0(:)=0.0d0
    T1(:)=0.0d0
    T2(:)=0.0d0
    DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory  
     DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives first dimension of trajectory
      !Get the current 'atom_clipboard', which usually would be the velocity.
      IF (READ_SEQUENTIAL) THEN
       atom_clipboard(:)=molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
      ELSE
       atom_clipboard(:)=molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
      ENDIF
      mass_clipboard=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
      !first, the normal, uncorrected temperature.
      TL(:)=TL(:)+mass_clipboard*(atom_clipboard(:))**2
      T0(:)=T0(:)+mass_clipboard*(atom_clipboard(:)-drift(:))**2
      T1(:)=T1(:)+mass_clipboard*(drift(:))**2
      T2(:)=T2(:)+mass_clipboard*2.0d0*(atom_clipboard(:)-drift(:))*(drift(:))
     ENDDO
    ENDDO
   ENDIF
   IF (PRESENT(kinetic_energy)) kinetic_energy=SUM(temperature(:))/2.0d0
   !divide temperatures by boltzmann constant as well as degrees of freedom...
   degrees_of_freedom=molecule_list(molecule_type_index)%number_of_atoms*molecule_list(molecule_type_index)%total_molecule_count
   !So far, degrees_of_freedom just contains the number of particles of that type. Which is sufficient, as we work in one dimension here.
   temperature(:)=temperature(:)/(DFLOAT(degrees_of_freedom)*boltzmann)
   corrected_temperature(:)=corrected_temperature(:)/(DFLOAT(degrees_of_freedom)*boltzmann)
   !account for constraints imposed on the molecule.
   IF (PRESENT(constraints_correction)) THEN
    !the internal_constraints are in three dimensions!
    internal_constraints=&
    &molecule_list(molecule_type_index)%total_molecule_count*molecule_list(molecule_type_index)%constraints
    !one additional degree of freedom per dimension is lost when the centre of mass is subtracted.
    constraints_correction=(DFLOAT(3*degrees_of_freedom))/(DFLOAT(3*degrees_of_freedom-3-internal_constraints))
   ENDIF
   !boltzmann was given in J/K. m*v² is in (g*angstroms²)/(mol*fs²). Change that.
   temperature(:)=1.0d7*temperature(:)/(avogadro)
   corrected_temperature(:)=1.0d7*corrected_temperature(:)/(avogadro)
   IF (DEVELOPERS_VERSION) THEN
    TL(:)=(1.0d7*TL(:))/(DFLOAT(degrees_of_freedom)*boltzmann*avogadro)
    T0(:)=(1.0d7*T0(:))/(DFLOAT(degrees_of_freedom)*boltzmann*avogadro)
    T1(:)=(1.0d7*T1(:))/(DFLOAT(degrees_of_freedom)*boltzmann*avogadro)
    T2(:)=(1.0d7*T2(:))/(DFLOAT(degrees_of_freedom)*boltzmann*avogadro)
    WRITE(*,'("  ! TL: ",3EN11.1)') SNGL(TL(:))
    WRITE(*,'("  ! T0: ",3EN11.1)') SNGL(T0(:))
    WRITE(*,'("  ! T1: ",3EN11.1)') SNGL(T1(:))
    WRITE(*,'("  ! T2: ",3EN11.1)') SNGL(T2(:))
    WRITE(*,'("  ! vd: ",3EN11.1)') SNGL(drift(:))
   ENDIF
  END SUBROUTINE give_temperature

  !The following subroutine calculates instantaneous Temperature for the whole box.
  !The result satisfies eq (14) in DOI 10.1021/acs.jpclett.9b02983
  REAL(KIND=WORKING_PRECISION) FUNCTION give_total_temperature(timestep)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep
  INTEGER :: molecule_type_index,molecule_index,atom_index
  REAL(KIND=WORKING_PRECISION) :: atom_clipboard(3),mass_clipboard
   IF (READ_SEQUENTIAL) CALL goto_timestep(timestep)
   !initialise output variables
   give_total_temperature=0.0d0
   !calculate the temperature, iterate over every molecule:
   DO molecule_type_index=1,number_of_molecule_types,1
    DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory    
     DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives first dimension of trajectory
      !Get the current 'atom_clipboard', which usually would be the velocity.
      IF (READ_SEQUENTIAL) THEN
       atom_clipboard(:)=molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
      ELSE
       atom_clipboard(:)=molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
      ENDIF
      mass_clipboard=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
      !sum up all the mv²
      give_total_temperature=give_total_temperature+mass_clipboard*SUM((atom_clipboard(:))**2)
     ENDDO
    ENDDO
   ENDDO
   !divide temperatures by boltzmann constant as well as degrees of freedom...
   give_total_temperature=give_total_temperature/(give_total_degrees_of_freedom()*boltzmann)
   !boltzmann was given in J/K. m*v² is in (g*angstroms²)/(mol*fs²). Change that.
   give_total_temperature=1.0d7*give_total_temperature/avogadro
  END FUNCTION give_total_temperature

  SUBROUTINE show_molecular_settings()
  IMPLICIT NONE
   WRITE(*,*) "Printing current molecular settings."
   WRITE(*,'("    ",A," ",I0)') "headerlines_to_skip       ",headerlines_to_skip
   WRITE(*,'("    ",A," ",I0)') "total_lines_to_skip       ",lines_to_skip
   WRITE(*,'("    ",A," ",I0)') "file_position             ",file_position
   WRITE(*,'("    ",A," ",L1)') "dihedrals_initialised     ",dihedrals_initialised
   WRITE(*,'("    ",A," ",L1)') "custom_constraints        ",custom_constraints
   WRITE(*,'("    ",A," ",L1)') "custom_default_masses     ",custom_default_masses
   WRITE(*,'("    ",A," ",L1)') "custom_atom_masses        ",custom_atom_masses
   WRITE(*,'("    ",A," ",L1)') "drudes_assigned           ",drudes_assigned
   WRITE(*,'("    ",A," ",L1)') "drude_details             ",drude_details
   WRITE(*,'("    ",A," ",L1)') "drudes_allocated          ",drudes_allocated
   WRITE(*,'("    ",A," ",L1)') "use_firstatom_as_com      ",use_firstatom_as_com
   WRITE(*,'("    ",A," ",L1)') "use_barycentre            ",use_barycentre
   WRITE(*,'("    ",A," ",I0)') "total_number_of_atoms     ",total_number_of_atoms
   WRITE(*,'("    ",A," ",I0)') "number_of_molecule_types  ",number_of_molecule_types
   WRITE(*,'("    ",A," ",I0)') "number_of_drude_particles ",number_of_drude_particles
   WRITE(*,'("    ",A," ",I0)') "number_of_steps           ",number_of_steps
   WRITE(*,'("    ",A," ",F5.3)')"drude_mass (Dalton)       ",drude_mass
   IF (drudes_assigned) WRITE(*,*) "invoke 'show_drude' to print detailed drude particle assignments."
   WRITE(*,*) "Memory requirement for storage of entire trajectory in RAM:"
   !getting a memory requirement estimate
   CALL print_memory_requirement()
  END SUBROUTINE show_molecular_settings

  SUBROUTINE show_drude_settings()
  IMPLICIT NONE
  INTEGER :: drude_flag,molecule_type_index,atom_index,number_of_assigned_DC_pairs
   IF (drudes_assigned) THEN
    WRITE(*,*) "Printing detailed drude information."
   ELSE
    CALL report_error(91)
    RETURN
   ENDIF
   number_of_assigned_DC_pairs=0 !counter for the number of assigned drude pairs
   WRITE(*,*) "Printing detailed drude information."
   DO molecule_type_index=1,number_of_molecule_types,1 !iterate over all molecule types.
    WRITE(*,'("   Molecule type ",I0," has ",I0," drude particles.")') &
    &molecule_type_index,molecule_list(molecule_type_index)%number_of_drudes_in_molecule
    IF (VERBOSE_OUTPUT) THEN
     WRITE(*,FMT='("     Indices of non-polarisable atoms:")')
     DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
      drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
      IF (drude_flag==-1) THEN
       WRITE(*,FMT='("       ",I0," (",A,")")')&
       &atom_index,TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))
      ENDIF
     ENDDO
    ENDIF
    WRITE(*,ADVANCE="NO",FMT='("     Indices of polarisable atoms (drude cores)")')
    IF (drude_details) THEN
     WRITE(*,'(" and detailed information from first step:")')
     WRITE(*,'("     (minimum and maximum core-drude distance, reduced mass, atom index of drude)")')
    ELSE
     WRITE(*,'(", atom index of drude:")')
    ENDIF
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
     drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
     IF (drude_flag>0) THEN
      number_of_assigned_DC_pairs=number_of_assigned_DC_pairs+1
      WRITE(*,FMT='("       ",I0," (",A,")")',ADVANCE="NO")&
      &atom_index,TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))
      IF (drude_details) THEN
       WRITE(*,'(2E9.2," 0",F0.4," ",I0)')&
       &molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%minimum_drude_distance,&
       &molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%maximum_drude_distance,&
       &molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%reduced_mass,&
       &drude_flag
      ELSE
       WRITE(*,'(" ",I0)') drude_flag
      ENDIF
     ENDIF
    ENDDO
    IF (VERBOSE_OUTPUT) THEN
     WRITE(*,FMT='("     Indices of drude particles:")')
     DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
      drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
      IF (drude_flag==0) THEN
        WRITE(*,FMT='("       ",I0," (",A,")")')&
       &atom_index,TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))
      ENDIF
     ENDDO
    ENDIF
   ENDDO
   WRITE(*,*) "To manually request the above assignment, add the "
   WRITE(*,'(" following ",I0," lines to your molecular input file:")') number_of_assigned_DC_pairs+1
   WRITE(*,'("   drudes ",I0)') number_of_assigned_DC_pairs 
   DO molecule_type_index=1,number_of_molecule_types,1
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
     drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
     IF (drude_flag>0) THEN
      WRITE(*,'("   ",I0," ",I0," ",I0)') molecule_type_index,atom_index,drude_flag
     ENDIF
    ENDDO
   ENDDO
  END SUBROUTINE show_drude_settings

  SUBROUTINE print_memory_requirement(required_storage_manual)
  IMPLICIT NONE
  REAL :: required_storage
  REAL(KIND=(WORKING_PRECISION)),OPTIONAL,INTENT(IN) :: required_storage_manual
  CHARACTER(LEN=2) :: memory_suffix
   IF (PRESENT(required_storage_manual)) THEN
    required_storage=required_storage_manual
   ELSE
    required_storage=DFLOAT(total_number_of_atoms)*DFLOAT(number_of_steps)*(12.0/1024.0d0)!That's Kibibytes KiB (just KB because whatever)
   ENDIF
   IF (required_storage<1000.0) THEN
    memory_suffix="KB"
   ELSE
    required_storage=required_storage/1024.0d0!Now, we're at Mebibytes MiB (who uses these symbols anyway?)
    IF (required_storage<1000.0) THEN
     memory_suffix="MB"
    ELSE
     required_storage=required_storage/1024.0d0!You get the idea.
     IF (required_storage<1000.0) THEN
      memory_suffix="GB"
     ELSE
      required_storage=required_storage/1024.0d0!and another one. I added that one later because it actually happens.
      memory_suffix="TB"
      IF (required_storage>100.0) THEN !100 TB is the capacity of "ephemereal". That's disk space, not RAM though...
       WRITE(*,*) "just... wow."
       RETURN
      ENDIF
     ENDIF
    ENDIF
   ENDIF
   IF (PRESENT(required_storage_manual)) THEN
    WRITE(*,ADVANCE="NO",FMT='(F6.1,A2)') required_storage,memory_suffix
   ELSE
    WRITE(*,'(" 3*",I0,"*",I0,"*4Byte =",F6.1,A2," (single precision)")')& !printing memory requirement
    &total_number_of_atoms,number_of_steps,required_storage,memory_suffix
   ENDIF
  END SUBROUTINE print_memory_requirement

  !This procedure is a serious bottleneck.
  SUBROUTINE goto_timestep(timestep)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep
  INTEGER :: counter,molecule_type_index,atom_index,molecule_index
  CHARACTER(LEN=2) :: dummystring
   !Important: Note that the file is positioned just AFTER the step specified in file_position! That saves backspacing every time.
   IF (timestep==file_position) THEN!check first if changing the timestep is necessary.
    RETURN
   ENDIF
   !The following section should only be executed if a change is necessary.
   IF (timestep>file_position) THEN!Have to advance
    !Skip over the unnecessary (including body) lines, i.e. (header+natoms)*(steps to skip)
    DO counter=1,lines_to_skip*(timestep-file_position-1),1!DO loop is not executed if the 'next step' is to be read!
     READ(9,*)
    ENDDO
    CALL read_snapshot_body()
   ELSE!timestep MUST be smaller than the file position - have to go back!
    !How far back is that new, requested timestep?
    IF (timestep>(file_position/2)) THEN
     !Use backspacing.
     !Dear future self, the following line might look weird, but it works. Trust me.
     DO counter=-1,lines_to_skip*(timestep-file_position-1),-1!
      BACKSPACE 9
     ENDDO
     CALL read_snapshot_body()
    ELSE
     !complete REWIND is advisable.
     REWIND 9
     !in principle, one could now set the file_position flag to zero and let the subroutine call itself recursively.
     !However, I prefer iteration:
     !Skip over the unnecessary (including body) lines, i.e. (header+natoms)*(steps to skip)
     DO counter=1,lines_to_skip*(timestep-1),1!'file_position' is zero here due to REWIND
      READ(9,*)
     ENDDO
     CALL read_snapshot_body()
    ENDIF
   ENDIF
   !set flag to new position.
   file_position=timestep
   !If necessary, wrap into box
   IF (WRAP_TRAJECTORY) CALL wrap_snap()
   !If necessary, change to barycentric reference frame.
   IF (use_barycentre) CALL remove_barycenter_sequential()
   CONTAINS

    SUBROUTINE read_snapshot_body()
    IMPLICIT NONE
    INTEGER :: ios
     !Read the required part, skipping only over the headerlines_to_skip
     DO counter=1,headerlines_to_skip,1
      READ(9,IOSTAT=ios,FMT=*)
      IF (ios/=0) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,'("stopped at step ",I0,", Headerline ",I0,".")') timestep,counter
       !This is a severe error - stop execution.
       CALL report_error(85,exit_status=ios)
      ENDIF
     ENDDO
     !THEN, read one molecule type after the other:
     DO molecule_type_index=1,number_of_molecule_types,1
      !For each molecule type, read the corresponding number of molecules:
      DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of snapshot (would be "2" for trajectory)
       !Finally, iterate over the atoms in that particular molecule:
       DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives first dimension of snapshot (would be "1" for trajectory)
        !LOOP VARIABLES:
        !molecule_type_index: current molecule type, e.g. 1 (only 1 and 2 for a binary IL)
        !molecule_index: current explicit molecule, e.g. molecule number 231
        !atom_index: current atom in that molecule, e.g. atom 2 (in molecule number 231 of type 1 in timestep 1234...)
        READ(9,*) dummystring,&
        &molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates
       ENDDO
      ENDDO 
     ENDDO
     !And now, the file is positioned just after the timestep that has been read.
    END SUBROUTINE read_snapshot_body

  END SUBROUTINE goto_timestep

  !The following set of functions provides the values of important variables to other routines. This serves the purpose of keeping variables local.
  CHARACTER(LEN=1024) FUNCTION give_sum_formula(molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index
  INTEGER :: outer,inner,n
  LOGICAL :: element_unused(molecule_list(molecule_type_index)%number_of_atoms)!has this atom been used up yet?
  CHARACTER(LEN=2) :: current_element
   element_unused(:)=.TRUE.
   give_sum_formula=""
   DO outer=1,molecule_list(molecule_type_index)%number_of_atoms,1
    current_element=TRIM(molecule_list(molecule_type_index)%list_of_elements(outer))
    !Print the element in outer, if not yet done:
    IF (element_unused(outer)) THEN
     !append the new element
     give_sum_formula=TRIM(give_sum_formula)//TRIM(current_element)
     !count how many are there, and label them as used
     n=1
     DO inner=(outer+1),molecule_list(molecule_type_index)%number_of_atoms,1
      IF (TRIM(current_element)==TRIM(molecule_list(molecule_type_index)%list_of_elements(inner))) THEN
       element_unused(inner)=.FALSE.
       n=n+1
      ENDIF
     ENDDO
     !append the number
     IF (n>1) THEN
      WRITE(give_sum_formula,'(A,I0)') TRIM(give_sum_formula),n
     ENDIF
    ENDIF
   ENDDO
  END FUNCTION give_sum_formula

  INTEGER FUNCTION give_number_of_molecule_types()
  IMPLICIT NONE
   give_number_of_molecule_types=number_of_molecule_types
  END FUNCTION give_number_of_molecule_types

  LOGICAL FUNCTION are_drudes_assigned()
  IMPLICIT NONE
   are_drudes_assigned=drudes_assigned
  END FUNCTION are_drudes_assigned

  LOGICAL FUNCTION constraints_available()
  IMPLICIT NONE
   constraints_available=custom_constraints
  END FUNCTION constraints_available

  REAL(KIND=WORKING_PRECISION)  FUNCTION give_box_volume()
  IMPLICIT NONE
   IF (BOX_VOLUME_GIVEN) THEN
    give_box_volume=box_size(1)*box_size(2)*box_size(3)
   ELSE
    give_box_volume=-1.0d0
    CALL report_error(41)
   ENDIF
  END FUNCTION give_box_volume

  REAL(KIND=WORKING_PRECISION)  FUNCTION give_box_limit()
  IMPLICIT NONE
   IF (BOX_VOLUME_GIVEN) THEN
    give_box_limit=MINVAL(box_size(:))
   ELSE
    give_box_limit=-1.0d0
    CALL report_error(41)
   ENDIF
  END FUNCTION give_box_limit

  INTEGER FUNCTION give_number_of_atoms_per_molecule(molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index
   give_number_of_atoms_per_molecule=molecule_list(molecule_type_index)%number_of_atoms
  END FUNCTION give_number_of_atoms_per_molecule

  INTEGER FUNCTION give_maximum_distance_squared()
  IMPLICIT NONE
   give_maximum_distance_squared=maximum_distance_squared
  END FUNCTION give_maximum_distance_squared

  INTEGER FUNCTION give_number_of_molecules_per_step(molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index
   give_number_of_molecules_per_step=molecule_list(molecule_type_index)%total_molecule_count
  END FUNCTION give_number_of_molecules_per_step

  INTEGER FUNCTION give_total_number_of_molecules_per_step()
  IMPLICIT NONE
  INTEGER :: molecule_type_index
   give_total_number_of_molecules_per_step=0
   DO molecule_type_index=1,number_of_molecule_types,1
    give_total_number_of_molecules_per_step=give_total_number_of_molecules_per_step+&
    &molecule_list(molecule_type_index)%total_molecule_count
   ENDDO
  END FUNCTION give_total_number_of_molecules_per_step

  REAL(KIND=GENERAL_PRECISION) FUNCTION give_mass_of_molecule(molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index
   give_mass_of_molecule=molecule_list(molecule_type_index)%mass
  END FUNCTION give_mass_of_molecule

  INTEGER FUNCTION give_charge_of_molecule(molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index
   give_charge_of_molecule=molecule_list(molecule_type_index)%charge
  END FUNCTION give_charge_of_molecule

  INTEGER FUNCTION give_number_of_atoms_per_step()
  IMPLICIT NONE
   give_number_of_atoms_per_step=total_number_of_atoms
  END FUNCTION give_number_of_atoms_per_step

  !The following function gives Nf in 10.1021/acs.jpclett.9b02983
  REAL(KIND=GENERAL_PRECISION) FUNCTION give_total_degrees_of_freedom()
  IMPLICIT NONE
  INTEGER :: degrees_of_freedom,m
   degrees_of_freedom=0
   DO m=1,number_of_molecule_types,1
    degrees_of_freedom=degrees_of_freedom-molecule_list(m)%constraints*molecule_list(m)%total_molecule_count
   ENDDO
   give_total_degrees_of_freedom=DFLOAT(degrees_of_freedom)+(total_number_of_atoms-number_of_drude_particles)*3.0d0
  END FUNCTION give_total_degrees_of_freedom

  INTEGER FUNCTION give_number_of_timesteps()
  IMPLICIT NONE
   give_number_of_timesteps=number_of_steps
  END FUNCTION give_number_of_timesteps

  !Writes the dihedrals defined in initialise_dihedrals for the specified timestep and molecule_index into dihedral_list.
  SUBROUTINE give_dihedrals(dihedral_list,timestep,molecule_index,dump_xyz)
  USE ANGLES
  IMPLICIT NONE
  REAL(KIND=GENERAL_PRECISION),INTENT(OUT) :: dihedral_list(:)
  REAL(KIND=GENERAL_PRECISION) :: dihedral_members(4,3)
  INTEGER :: n,m
  INTEGER,INTENT(IN) :: timestep,molecule_index
  LOGICAL,INTENT(IN),OPTIONAL :: dump_xyz
  LOGICAL :: writexyz,connected
  CHARACTER(LEN=1024) :: fstring
   IF (READ_SEQUENTIAL) CALL goto_timestep(timestep)
   writexyz=.FALSE.
   IF (PRESENT(dump_xyz)) THEN
    IF (dump_xyz) THEN
     writexyz=.TRUE. !if requested, then the skeleton of the dihedral will be written into an xyz file.
    ENDIF
   ENDIF
   IF (SIZE(dihedral_list)/=number_of_dihedrals) CALL report_error(13)
   !The following part is responsible for writing the dihedral angles (in degrees) in the output list.
   INQUIRE(UNIT=3,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=3)
   DO m=1,number_of_dihedrals,1 ! m iterates over all the dihedrals to dump
    IF (writexyz) THEN
     WRITE(fstring,'(2A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"dihedral_",m,".xyz"
     OPEN(UNIT=3,FILE=TRIM(fstring))
     WRITE(3,*) 4
     WRITE(3,*)
    ENDIF
    DO n=1,4,1 ! n iterates over the four atom indices required for a dihedral
     IF (READ_SEQUENTIAL) THEN
      dihedral_members(n,:)=molecule_list(molecule_type_index_for_dihedrals)%&
      &snapshot(dihedral_member_indices(m,n),molecule_index)%coordinates(:)
     ELSE
      dihedral_members(n,:)=molecule_list(molecule_type_index_for_dihedrals)%&
      &trajectory(dihedral_member_indices(m,n),molecule_index,timestep)%coordinates(:)
     ENDIF
     IF (writexyz) WRITE(3,*) molecule_list(molecule_type_index_for_dihedrals)%&
      &list_of_elements(dihedral_member_indices(m,n)),dihedral_members(n,:)
    ENDDO
    IF (writexyz) CLOSE(UNIT=3)
    dihedral_list(m)=dihedral_angle(dihedral_members)
   ENDDO
  END SUBROUTINE give_dihedrals

  !report information about the fragments
  SUBROUTINE give_fragment_information(tip_fragment)
  IMPLICIT NONE
  LOGICAL,INTENT(IN) :: tip_fragment
  INTEGER :: natoms,first_atom_index
  REAL(KIND=WORKING_PRECISION) ::  mass
  CHARACTER(LEN=1024) :: fragment_formula
   !for fragment: report mass, number of atoms, sum formula
   IF (tip_fragment) THEN
    natoms=number_of_tip_atoms
    WRITE(*,FMT='(" tip  fragment ")',ADVANCE="NO")
    first_atom_index=fragment_list_tip(1)
    mass=mass_of_tip_fragment
   ELSE
    natoms=number_of_base_atoms
    WRITE(*,FMT='(" base fragment ")',ADVANCE="NO")
    first_atom_index=fragment_list_base(1)
    mass=mass_of_base_fragment
   ENDIF
   CALL write_fragment_formula()
   IF (natoms==1) THEN
    WRITE(*,'("consists of one ",A," atom with index ",I0,".")') TRIM(fragment_formula),first_atom_index
   ELSE
    WRITE(*,'(" (",A,") consists of ",I0," atoms and has a molecular weight of ",F0.4,".")') TRIM(fragment_formula),natoms,mass
   ENDIF

  CONTAINS

   SUBROUTINE write_fragment_formula()
   IMPLICIT NONE
   INTEGER :: outer,inner,n,maxcount
   LOGICAL :: element_unused(MAXVAL((/ number_of_base_atoms,number_of_tip_atoms /)))!has this atom been used up yet?
   CHARACTER(LEN=2) :: current_element,list_element
    element_unused(:)=.TRUE.
    fragment_formula=""
    IF (tip_fragment) THEN
     maxcount=number_of_tip_atoms
    ELSE
     maxcount=number_of_base_atoms
    ENDIF
    IF (maxcount==1) THEN
     !There is just one atom in the molecule.
     IF (tip_fragment) THEN
      fragment_formula=TRIM(molecule_list(molecule_type_index_for_fragments)%list_of_elements(fragment_list_tip(1)))
     ELSE
      fragment_formula=TRIM(molecule_list(molecule_type_index_for_fragments)%list_of_elements(fragment_list_base(1)))
     ENDIF
     RETURN
    ENDIF
    DO outer=1,maxcount,1
     IF (tip_fragment) THEN
      current_element=TRIM(molecule_list(molecule_type_index_for_fragments)%list_of_elements(fragment_list_tip(outer)))
     ELSE
      current_element=TRIM(molecule_list(molecule_type_index_for_fragments)%list_of_elements(fragment_list_base(outer)))
     ENDIF
     !Print the element in outer, if not yet done:
     IF (element_unused(outer)) THEN
      !append the new element
      fragment_formula=TRIM(fragment_formula)//TRIM(current_element)
      !count how many are there, and label them as used
      n=1
      DO inner=(outer+1),maxcount,1
       IF (tip_fragment) THEN
        list_element=TRIM(molecule_list(molecule_type_index_for_fragments)%list_of_elements(fragment_list_tip(inner)))
       ELSE
        list_element=TRIM(molecule_list(molecule_type_index_for_fragments)%list_of_elements(fragment_list_base(inner)))
       ENDIF
       IF (TRIM(current_element)==TRIM(list_element)) THEN
        element_unused(inner)=.FALSE.
        n=n+1
       ENDIF
      ENDDO
      !append the number
      IF (n>1) THEN
       WRITE(fragment_formula,'(A,I0)') TRIM(fragment_formula),n
      ENDIF
     ENDIF
    ENDDO
   END SUBROUTINE write_fragment_formula
  
  END SUBROUTINE give_fragment_information

  !Calculate center of mass of tip fragment
  FUNCTION give_tip_fragment(timestep,molecule_index)
  IMPLICIT NONE
  REAL(KIND=WORKING_PRECISION) :: give_tip_fragment(3),element_mass,unweighted_pos(3)
  INTEGER,INTENT(IN) :: timestep,molecule_index
  INTEGER :: counter
   IF (number_of_tip_atoms==1) THEN
    !centre of mass doesn't make a difference for just one atom.
    IF (READ_SEQUENTIAL) THEN
     IF (timestep/=file_position) CALL goto_timestep(timestep)
     give_tip_fragment(:)=molecule_list(molecule_type_index_for_fragments)%&
     &snapshot(fragment_list_tip(1),molecule_index)%coordinates(:)
    ELSE
     give_tip_fragment(:)=molecule_list(molecule_type_index_for_fragments)%&
     &trajectory(fragment_list_tip(1),molecule_index,timestep)%coordinates(:)
    ENDIF
    RETURN
   ENDIF
   give_tip_fragment(:)=0.0d0
   DO counter=1,number_of_tip_atoms,1
    element_mass=molecule_list(molecule_type_index_for_fragments)%list_of_atom_masses(fragment_list_tip(counter))
    IF (READ_SEQUENTIAL) THEN
     IF (timestep/=file_position) CALL goto_timestep(timestep)
     unweighted_pos(:)=DBLE(molecule_list(molecule_type_index_for_fragments)%&
     &snapshot(fragment_list_tip(counter),molecule_index)%coordinates(:))
    ELSE
     unweighted_pos(:)=DBLE(molecule_list(molecule_type_index_for_fragments)%&
     &trajectory(fragment_list_tip(counter),molecule_index,timestep)%coordinates(:))
    ENDIF
    give_tip_fragment(:)=give_tip_fragment(:)+element_mass*unweighted_pos(:)
   ENDDO
   give_tip_fragment(:)=give_tip_fragment(:)/mass_of_tip_fragment
  END FUNCTION give_tip_fragment

  !Calculate center of mass of base fragment
  FUNCTION give_base_fragment(timestep,molecule_index)
  IMPLICIT NONE
  REAL(KIND=WORKING_PRECISION) :: give_base_fragment(3),element_mass,unweighted_pos(3)
  INTEGER,INTENT(IN) :: timestep,molecule_index
  INTEGER :: counter
   IF (number_of_base_atoms==1) THEN
    !centre of mass doesn't make a difference for just one atom.
    IF (READ_SEQUENTIAL) THEN
     IF (timestep/=file_position) CALL goto_timestep(timestep)
     give_base_fragment(:)=molecule_list(molecule_type_index_for_fragments)%&
     &snapshot(fragment_list_base(1),molecule_index)%coordinates(:)
    ELSE
     give_base_fragment(:)=molecule_list(molecule_type_index_for_fragments)%&
     &trajectory(fragment_list_base(1),molecule_index,timestep)%coordinates(:)
    ENDIF
    RETURN
   ENDIF
   give_base_fragment(:)=0.0d0
   DO counter=1,number_of_base_atoms,1
    element_mass=molecule_list(molecule_type_index_for_fragments)%list_of_atom_masses(fragment_list_base(counter))
    IF (READ_SEQUENTIAL) THEN
     IF (timestep/=file_position) CALL goto_timestep(timestep)
     unweighted_pos(:)=DBLE(molecule_list(molecule_type_index_for_fragments)%&
     &snapshot(fragment_list_base(counter),molecule_index)%coordinates(:))
    ELSE
     unweighted_pos(:)=DBLE(molecule_list(molecule_type_index_for_fragments)%&
     &trajectory(fragment_list_base(counter),molecule_index,timestep)%coordinates(:))
    ENDIF
    give_base_fragment(:)=give_base_fragment(:)+element_mass*unweighted_pos(:)
   ENDDO
   give_base_fragment(:)=give_base_fragment(:)/mass_of_base_fragment
  END FUNCTION give_base_fragment

  !This subroutine gives the four dihedral member indices for a specified dihedral index
  FUNCTION give_dihedral_member_indices(dihedral_index)
  IMPLICIT NONE
  INTEGER :: give_dihedral_member_indices(4)
  INTEGER,INTENT(IN) :: dihedral_index
   give_dihedral_member_indices(:)=dihedral_member_indices(dihedral_index,:)
  END FUNCTION give_dihedral_member_indices

  !this subroutine provides an interface to choose the dihedral angles to be reported by give_dihedrals in a controlled way.
  !set_number_of_dihedrals is usually two, for example for NTf2 
  SUBROUTINE initialise_dihedrals(set_dihedral_member_indices,molecule_type_index,set_number_of_dihedrals)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: set_dihedral_member_indices(:,:)!list of dihedrals. First dimension = index/number of dihedral, second dimension = the atom indices of the dihedral members.
  INTEGER :: allocstatus,deallocstatus,m,n
  INTEGER,INTENT(IN) :: molecule_type_index,set_number_of_dihedrals !which molecule type the dihedrals belong to and how many dihedrals there are.
   !re-initialise the dihedral list, if necessary.
   IF (dihedrals_initialised) THEN
    !error 12 is treated as warning if no exit status is passed, and as severe error otherwise.
    CALL report_error(12)
    DEALLOCATE(dihedral_member_indices,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(12,exit_status=deallocstatus)
   ENDIF
   number_of_dihedrals=set_number_of_dihedrals
   ALLOCATE(dihedral_member_indices(number_of_dihedrals,4),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(11,exit_status=allocstatus)
   DO m=1,number_of_dihedrals,1! m iterates over all the dihedrals
    DO n=1,4,1 ! n iterates over the four atom indices required for a dihedral
     dihedral_member_indices(m,n)=set_dihedral_member_indices(m,n)
    ENDDO
   ENDDO
   molecule_type_index_for_dihedrals=molecule_type_index
   dihedrals_initialised=.TRUE.
  END SUBROUTINE initialise_dihedrals

  SUBROUTINE initialise_fragments&
  &(set_fragments_tip,set_fragments_base,set_number_of_tip_atoms,set_number_of_base_atoms,molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: set_fragments_tip(:),set_fragments_base(:),set_number_of_tip_atoms,set_number_of_base_atoms
  INTEGER,INTENT(IN) :: molecule_type_index
  INTEGER :: allocstatus,deallocstatus,m
   !check molecule type index before passing it to this routine!
   IF ((molecule_type_index<1).OR.(molecule_type_index>number_of_molecule_types)) CALL report_error(0)
   !re-initialise the fragment list, if necessary
   IF (fragments_initialised) THEN
    !error 78 is treated as warning if no exit status is passed, and as severe error otherwise.
    CALL report_error(78)
    DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(78,exit_status=deallocstatus)
    DEALLOCATE(fragment_list_base,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(78,exit_status=deallocstatus)
   ENDIF
   molecule_type_index_for_fragments=molecule_type_index
   number_of_base_atoms=set_number_of_base_atoms
   number_of_tip_atoms=set_number_of_tip_atoms
   ALLOCATE(fragment_list_tip(set_number_of_tip_atoms),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(79,exit_status=allocstatus)
   ALLOCATE(fragment_list_base(set_number_of_base_atoms),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(79,exit_status=allocstatus)
   mass_of_base_fragment=0.0d0
   mass_of_tip_fragment=0.0d0
   DO m=1,number_of_base_atoms,1
    fragment_list_base(m)=set_fragments_base(m)
    mass_of_base_fragment=mass_of_base_fragment+&
    &molecule_list(molecule_type_index_for_fragments)%list_of_atom_masses(fragment_list_base(m))
    IF ((set_fragments_base(m)<1).OR.(set_fragments_base(m)>molecule_list(molecule_type_index)%number_of_atoms)) THEN
     !out of bounds for number of atoms
     CALL report_error(80,exit_status=set_fragments_base(m))
    ENDIF
   ENDDO
   DO m=1,number_of_tip_atoms,1
    fragment_list_tip(m)=set_fragments_tip(m)
    mass_of_tip_fragment=mass_of_tip_fragment+&
    &molecule_list(molecule_type_index_for_fragments)%list_of_atom_masses(fragment_list_tip(m))
    IF ((set_fragments_tip(m)<1).OR.(set_fragments_tip(m)>molecule_list(molecule_type_index)%number_of_atoms)) THEN
     !out of bounds for number of atoms
     CALL report_error(80,exit_status=set_fragments_tip(m))
    ENDIF
   ENDDO
   molecule_type_index_for_fragments=molecule_type_index
   fragments_initialised=.TRUE.
  END SUBROUTINE initialise_fragments

  !Writes the element lists to the standard output
  SUBROUTINE report_element_lists
  IMPLICIT NONE
  INTEGER :: molecule_type_index,atom_index,natoms
   WRITE(*,*) " ! ELEMENT LISTS:"
   DO molecule_type_index=1,number_of_molecule_types,1
    natoms=molecule_list(molecule_type_index)%number_of_atoms
    WRITE(*,'("  ! molecule #",I0,":")') molecule_type_index
    WRITE(*,*) " ! ",(TRIM(molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1))&
    &,atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1)
    CALL print_atom_properties(molecule_type_index)
   ENDDO
  END SUBROUTINE report_element_lists

  SUBROUTINE print_atom_properties(molecule_type_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index
  INTEGER :: atom_index
   !replace
   !atomic_weight(molecule_list(molecule_type_index)%list_of_elements(atom_index))
   !with
   !molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
   WRITE(*,'("  ! number element atomic_mass native_mass")')
   DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
    WRITE(*,'("  !   ",I0," ",A2," ",F0.3," ",F0.3)') atom_index,&
    &molecule_list(molecule_type_index)%list_of_elements(atom_index),&
    &molecule_list(molecule_type_index)%list_of_atom_masses(atom_index),&
    &atomic_weight(molecule_list(molecule_type_index)%list_of_elements(atom_index))
   ENDDO
  END SUBROUTINE print_atom_properties

  FUNCTION give_center_of_mass(timestep,molecule_type_index,molecule_index)
  IMPLICIT NONE
  REAL(KIND=WORKING_PRECISION) :: give_center_of_mass(3),weighted_pos(3)!higher precision, because intermediate result.
  INTEGER :: atom_index
  INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
   IF (use_firstatom_as_com) THEN
    give_center_of_mass(:)=DBLE(molecule_list(molecule_type_index)%trajectory(1,molecule_index,timestep)%coordinates(:))
    RETURN
   ENDIF
   IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
   give_center_of_mass(:)=0.0d0
   DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
    !first, the current atom's position is stored in weighted_pos.
    !added support for sequential read.
    IF (READ_SEQUENTIAL) THEN
     weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:))
    ELSE
     weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:))
    ENDIF
    !then, this position is weighted with the atom's mass
    weighted_pos(:)=weighted_pos(:)*molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
    !this weighted position is now added to the center of mass.
    give_center_of_mass(:)=give_center_of_mass(:)+weighted_pos(:)
   ENDDO
   !finally, the center of mass has to be normalised by the total mass.
   give_center_of_mass(:)=give_center_of_mass(:)/DBLE(molecule_list(molecule_type_index)%mass)
  END FUNCTION give_center_of_mass

  FUNCTION give_box_center_of_mass(timestep)
  IMPLICIT NONE
  REAL :: give_box_center_of_mass(3)
  REAL(KIND=WORKING_PRECISION) :: pos_atom(3),weighted_pos(3),total_mass
  INTEGER :: atom_index
  INTEGER,INTENT(IN) :: timestep
  INTEGER :: molecule_type_index,molecule_index
   IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
   weighted_pos(:)=0.0d0
   total_mass=0.0d0
   DO molecule_type_index=1,number_of_molecule_types,1
    total_mass=total_mass+DBLE(molecule_list(molecule_type_index)%total_molecule_count)*&
    &molecule_list(molecule_type_index)%mass
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
     pos_atom(:)=0.0d0
     DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
      !first, the current atom's position is added to pos_atom.
      !added support for sequential read.
      IF (READ_SEQUENTIAL) THEN
       pos_atom(:)=pos_atom(:)+&
       &DBLE(molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:))
      ELSE
       pos_atom(:)=pos_atom(:)+&
       &DBLE(molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:))
      ENDIF
      !then, this position is weighted with the atom's mass
     ENDDO
     !since every atom 'atom_index' in a molecule type is the same, we can take advantage of:
     !m1*r1+m1*r2+m2*r3+m2*r4=m1*(r1+r2)+m2*(r3+r4)
     weighted_pos(:)=weighted_pos(:)+&
     &pos_atom(:)*molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
    ENDDO
   ENDDO
   !finally, the center of mass has to be normalised by the total mass.
   give_box_center_of_mass(:)=(weighted_pos(:)/total_mass)
  END FUNCTION give_box_center_of_mass

  SUBROUTINE remove_barycenter_sequential()
  IMPLICIT NONE
  INTEGER :: molecule_type_index,molecule_index,atom_index
  REAL :: box_COM(3)
  !needs file_position to be up to date!
  box_COM(:)=give_box_center_of_mass(file_position)
  DO molecule_type_index=1,number_of_molecule_types,1
   DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
    DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
     molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)=&
     &molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)-box_COM(:)
    ENDDO
   ENDDO
  ENDDO
  END SUBROUTINE remove_barycenter_sequential

  SUBROUTINE remove_barycenter_fulltraj()
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$ END INTERFACE
  INTEGER :: stepcounter,molecule_type_index,molecule_index,atom_index
  REAL :: box_COM(3)
   IF (READ_SEQUENTIAL) CALL report_error(0)
   IF (VERBOSE_OUTPUT) WRITE(*,*) "Removing box centre-of-mass from every step."
   !$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
   !$OMP PRIVATE(box_COM,molecule_type_index,molecule_index,atom_index)
   !$OMP SINGLE
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (remove barycentre)"
   !$  CALL timing_parallel_sections(.TRUE.)
   !$ ENDIF
   !$OMP END SINGLE
   !$OMP DO SCHEDULE(STATIC,1)
   DO stepcounter=1,number_of_steps,1
    box_COM(:)=give_box_center_of_mass(stepcounter)
    DO molecule_type_index=1,number_of_molecule_types,1
     DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
      DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
       molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,stepcounter)%coordinates(:)=&
       &molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,stepcounter)%coordinates(:)-box_COM(:)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
   !$  CALL timing_parallel_sections(.FALSE.)
   !$ ENDIF
  END SUBROUTINE remove_barycenter_fulltraj

  SUBROUTINE switch_to_barycenter()
  IMPLICIT NONE
   IF (use_barycentre) THEN
    PRINT *,"Already using barycentric reference frame."
    RETURN
   ENDIF
   use_barycentre=.TRUE.
   PRINT *,"Switching to barycentric reference frame."
   IF (READ_SEQUENTIAL) THEN
    CALL reset_trajectory_file()
    CALL remove_barycenter_sequential()
   ELSE
    CALL remove_barycenter_fulltraj()
   ENDIF
  END SUBROUTINE switch_to_barycenter

  !writes the specified molecule in xyz format into the specified unit. unit is not overwritten if opened with APPEND!
  !the blanks are not added!
  !If include_header=.TRUE. (which is the default), then the first two lines for an xyz file are added.
  SUBROUTINE write_molecule(unit_number,timestep_in,molecule_type_index,molecule_index,include_header,&
  &custom_header,translate_by)
  IMPLICIT NONE
  LOGICAL,INTENT(IN),OPTIONAL :: include_header
  INTEGER,INTENT(IN) :: unit_number,molecule_type_index,molecule_index,timestep_in
  INTEGER :: atom_index,natoms,timestep
  CHARACTER(LEN=*),OPTIONAL :: custom_header
  REAL(KIND=WORKING_PRECISION),INTENT(IN),OPTIONAL :: translate_by(3)
  REAL(KIND=WORKING_PRECISION) :: shift(3)
   IF (PRESENT(translate_by)) THEN
    shift(:)=translate_by(:)
   ELSE
    shift(:)=0.0d0
   ENDIF
   timestep=timestep_in
   IF (READ_SEQUENTIAL) THEN
    IF (timestep==-1) THEN
     !not really required, but included for safety here in case I change something in the procedures' body that requires the timestep.
     timestep=file_position
    ELSE
     CALL goto_timestep(timestep)
    ENDIF
   ELSE
    IF (timestep==-1) timestep=1
   ENDIF
   natoms=molecule_list(molecule_type_index)%number_of_atoms
   IF (PRESENT(include_header)) THEN
    IF (include_header) THEN
     WRITE(unit_number,'(I0)') molecule_list(molecule_type_index)%number_of_atoms
     IF (PRESENT(custom_header)) THEN
      WRITE(unit_number,*) TRIM(ADJUSTL(custom_header))
     ELSE
      WRITE(unit_number,'(" Molecular weight = ",F0.4)') molecule_list(molecule_type_index)%mass
     ENDIF
    ENDIF
   ELSE
    WRITE(unit_number,'(I0)') molecule_list(molecule_type_index)%number_of_atoms
    IF (PRESENT(custom_header)) THEN
     WRITE(unit_number,*) TRIM(ADJUSTL(custom_header))
    ELSE
     WRITE(unit_number,'(" Molecular weight = ",F0.4)') molecule_list(molecule_type_index)%mass
    ENDIF
   ENDIF
   DO atom_index=1,natoms,1
    IF (READ_SEQUENTIAL) THEN
     WRITE(unit_number,*) molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),&
     &SNGL(molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)+shift(:))
    ELSE
     WRITE(unit_number,*) molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),&
     &SNGL(molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)+shift(:))
    ENDIF
   ENDDO
  END SUBROUTINE write_molecule

  !writes the specified molecule in xyz format into the specified unit, merging drudes into cores.
  !no blanks are added. header is included only if include_header is set to .TRUE. - the default is that no header is added here!
  SUBROUTINE write_molecule_merged_drudes(unit_number,timestep_in,molecule_type_index,molecule_index,include_header)
  IMPLICIT NONE
  LOGICAL,INTENT(IN),OPTIONAL :: include_header
  INTEGER,INTENT(IN) :: unit_number,molecule_type_index,molecule_index,timestep_in
  INTEGER :: atom_index,natoms,timestep,drude_flag
  REAL :: drudepos(3),corepos(3)
   timestep=timestep_in
   IF (READ_SEQUENTIAL) THEN
    IF (timestep==-1) THEN
     !not really required, but included for safety here in case I change something in the procedures' body that requires the timestep.
     timestep=file_position
    ELSE
     CALL goto_timestep(timestep)
    ENDIF
   ELSE
    IF (timestep==-1) timestep=1
   ENDIF
   natoms=molecule_list(molecule_type_index)%number_of_atoms
   IF (PRESENT(include_header)) THEN
    IF (include_header) THEN
     WRITE(unit_number,'(I0)') molecule_list(molecule_type_index)%number_of_atoms
     WRITE(unit_number,'(" Molecular weight = ",F0.4)') molecule_list(molecule_type_index)%mass
    ENDIF
   ENDIF
   DO atom_index=1,natoms,1
    !drude_flag is the atom index of the drude particle attached to this core.
    !Will be -1 if no drude particle is found, and 0 if this is a drude itself.
    drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
    SELECT CASE (drude_flag)
    CASE (-1) !non-polarisable atom - print just as it is.
     IF (READ_SEQUENTIAL) THEN
      WRITE(unit_number,*) molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),&
      &molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
     ELSE
      WRITE(unit_number,*) molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),&
      &molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
     ENDIF
    CASE (0) !this is a drude particle - skip this one by cycling.
     CYCLE
    CASE DEFAULT
     !everything else should be drude cores - merge with the drude particle.
     IF (READ_SEQUENTIAL) THEN
      corepos(:)=molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
      drudepos(:)=molecule_list(molecule_type_index)%snapshot(drude_flag,molecule_index)%coordinates(:)
     ELSE
      corepos(:)=molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
      drudepos(:)=molecule_list(molecule_type_index)%trajectory(drude_flag,molecule_index,timestep)%coordinates(:)
     ENDIF
     !merge drudepos into corepos:
     corepos(:)=(corepos(:)*molecule_list(molecule_type_index)%list_of_atom_masses(MODULO(atom_index-1,natoms)+1)+& !core position * its mass
     &drudepos(:)*molecule_list(molecule_type_index)%list_of_atom_masses(MODULO(drude_flag-1,natoms)+1))/& !drude position * its mass
     !The following the lines give the total mass, by which will be divided:
     &(molecule_list(molecule_type_index)%list_of_atom_masses(MODULO(atom_index-1,natoms)+1)& !for the core...
     &+molecule_list(molecule_type_index)%list_of_atom_masses(MODULO(drude_flag-1,natoms)+1)) !and the drude.
     !Then, write to output.
     WRITE(unit_number,*)&
     &molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),corepos(:)
    END SELECT
   ENDDO
  END SUBROUTINE write_molecule_merged_drudes

  !initialises the molecular module by reading the input file 'molecular.inp'
  SUBROUTINE initialise_molecular()
  IMPLICIT NONE
  LOGICAL :: file_exists,connected
  INTEGER :: ios,n,allocstatus,a,b,c,totalcharge,headerlines_molecular,m
  CHARACTER(LEN=16) :: inputstring
  REAL(KIND=SP) :: mass_input
   !Turn COMboost off. Use ONLY with whole trajectory, NOT with read_sequential.
   use_firstatom_as_com=.FALSE.
   use_barycentre=.FALSE.
   file_position=-1
   dihedrals_initialised=.FALSE.
   ! first, check if file exists. If not, switch to user input for this part.
   INQUIRE(FILE=TRIM(FILENAME_MOLECULAR_INPUT),EXIST=file_exists)!no input path is added for the molecular file!
   IF (file_exists) THEN
    IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(FILENAME_MOLECULAR_INPUT),"'"!no input path is added for the molecular file!
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    OPEN(UNIT=3,FILE=TRIM(FILENAME_MOLECULAR_INPUT),&!no input path is added for the molecular file!
    &ACTION='READ',IOSTAT=ios)
    IF (ios/=0) CALL report_error(7,exit_status=ios)
    !Read header!
    CALL read_molecular_input_file_header()
    IF (totalcharge/=0) CALL report_error(52,exit_status=totalcharge)
    !Read default masses!
    CALL read_molecular_input_file_default_masses()
    IF (drude_mass>0.001d0) CALL subtract_drude_masses()
    !Read atomic masses
    CALL read_molecular_input_file_atomic_masses()
    !Read constraints!
    CALL read_molecular_input_file_constraints()
    !Read drudes!
    CALL read_molecular_input_file_drudes()
    CLOSE(UNIT=3)
   ELSE
    CALL report_error(17)
   ENDIF
   CONTAINS

    SUBROUTINE read_molecular_input_file_header()
    IMPLICIT NONE
     headerlines_molecular=2
     READ(3,IOSTAT=ios,FMT=*) number_of_steps
     IF (ios/=0) CALL report_error(7,exit_status=ios)
     READ(3,IOSTAT=ios,FMT=*) number_of_molecule_types
     IF (ios/=0) CALL report_error(7,exit_status=ios)
     total_number_of_atoms=0
     totalcharge=0
     !allocate memory for list of molecules
     ALLOCATE(molecule_list(number_of_molecule_types),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
     headerlines_molecular=headerlines_molecular+number_of_molecule_types
     !Iterate over all the molecule types - n is the molecule index here.
     DO n=1,number_of_molecule_types,1
      READ(3,IOSTAT=ios,FMT=*) a,b,c
      IF (ios/=0) CALL report_error(7,exit_status=ios)
      molecule_list(n)%charge=a
      molecule_list(n)%number_of_atoms=b
      IF (b<1) CALL report_error(80,exit_status=b)
      molecule_list(n)%total_molecule_count=c
      IF (c<1) CALL report_error(119,exit_status=c)
      total_number_of_atoms=total_number_of_atoms+b*c
      totalcharge=totalcharge+a*c
      ALLOCATE(molecule_list(n)%list_of_elements(b),STAT=allocstatus)
      IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
      ALLOCATE(molecule_list(n)%list_of_atom_masses(b),STAT=allocstatus)
      ALLOCATE(molecule_list(n)%manual_atom_mass_specified(b),STAT=allocstatus)
      IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
      IF (READ_SEQUENTIAL) THEN
       ALLOCATE(molecule_list(n)%snapshot(b,c),STAT=allocstatus)
       IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
      ELSE
       ALLOCATE(molecule_list(n)%trajectory(b,c,number_of_steps),STAT=allocstatus)
       IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
      ENDIF
      !at this point, the list_of_elements and mass are not specified. This information will be read from the lammps trajectory.
     ENDDO
     !Check if COM trajectory
     use_firstatom_as_com=.TRUE.
     DO n=1,number_of_molecule_types,1
      IF (molecule_list(n)%number_of_atoms>1) THEN
       use_firstatom_as_com=.FALSE.
       EXIT
      ENDIF
     ENDDO
     IF ((use_firstatom_as_com).AND.(VERBOSE_OUTPUT)) THEN
      PRINT *,"All molecules have only 1 atom - turn on com boost."
      PRINT *,"(centre-of-mass position = position of this atom)"
     ENDIF
    END SUBROUTINE read_molecular_input_file_header

    SUBROUTINE read_molecular_input_file_default_masses()
    IMPLICIT NONE
    CHARACTER(LEN=2) :: shortstring
     COM_mass_list(:)=0.0d0
     custom_default_masses=.FALSE.
     WRITE(*,ADVANCE="NO",FMT='(" Searching for ",A," statement...")') "'default_masses'"
     !Skip over headerlines in molecular input file
     REWIND 3
     DO n=1,headerlines_molecular,1
      READ(3,*)
     ENDDO
     !search for masses section.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF (ios<0) THEN
       !end of file encountered
       WRITE(*,'("done (none found, end of file encountered).")')
       EXIT
      ENDIF
      !support for synonyms
      IF (TRIM(inputstring)=="default_masses") inputstring="masses"
      IF (ios==0) THEN
       IF (TRIM(inputstring)=="masses") THEN
        WRITE(*,'("found in line ",I0,".")') n+headerlines_molecular
        BACKSPACE 3
        READ(3,IOSTAT=ios,FMT=*) inputstring,a
        IF (ios/=0) THEN
         !something went wrong
         custom_default_masses=.FALSE.
        ELSE
         !keyword ok - read the section.
         IF (a>0) THEN
          custom_default_masses=.TRUE.
          WRITE(*,FMT='(" Trying to read ",I0," custom default masses...")',ADVANCE="NO") a
          DO m=1,a,1
           READ(3,IOSTAT=ios,FMT=*) shortstring,mass_input
           IF (ios/=0) THEN
            !wrong format... abort.
            custom_default_masses=.FALSE.
            EXIT
           ELSE
            IF (ANY(ALPHABET_small==IACHAR(shortstring(1:1)))) THEN
             !lowercase letter found! Check if already specified.
             IF (COM_mass_list(IACHAR(shortstring(1:1)))>0.001d0) THEN
              WRITE(*,*)
              CALL report_error(58,exit_status=IACHAR(shortstring(1:1)))
             ENDIF
             IF (mass_input<0.0d0) THEN
              WRITE(*,*)
              CALL report_error(59,exit_status=IACHAR(shortstring(1:1)))
             ENDIF
             COM_mass_list(IACHAR(shortstring(1:1)))=mass_input
            ELSE
             IF ((shortstring(1:1)=="X").OR.(shortstring(1:1)=="D")) THEN
              IF (drude_mass>0.001d0) THEN
               WRITE(*,*)
               CALL report_error(58,exit_status=IACHAR(shortstring(1:1)))
              ENDIF
             ENDIF
             CALL change_default_mass(shortstring,mass_input)
            ENDIF
           ENDIF
          ENDDO
          !test if custom_default_masses still true. If not, print error message.
          IF (custom_default_masses) THEN
           WRITE(*,'(" done.")') 
          ELSE
           WRITE(*,'(" failed.")') 
           CALL report_error(61)
           COM_mass_list(:)=0.0d0
          ENDIF
         ENDIF
        ENDIF
        EXIT
       ELSEIF (TRIM(inputstring)=="quit") THEN
        WRITE(*,'("done (none found before ",A,").")') "'quit'"
        EXIT
       ENDIF
      ENDIF
     ENDDO
    END SUBROUTINE read_molecular_input_file_default_masses

    SUBROUTINE change_default_mass(element_name_input,new_mass)
    IMPLICIT NONE
    REAL,INTENT(IN) :: new_mass
    REAL :: old_mass
    LOGICAL,SAVE :: print_blank=.TRUE.
    !If you change this part, then also change Module_SETTINGS and atomic_weight
    CHARACTER(LEN=*),INTENT(IN) :: element_name_input
    CHARACTER(LEN=32) :: element_name_full
     IF (print_blank) THEN
      IF (VERBOSE_OUTPUT) WRITE(*,*)
      print_blank=.FALSE.
     ENDIF
     SELECT CASE (TRIM(element_name_input))
     CASE ("H")
      old_mass=mass_hydrogen
      mass_hydrogen=new_mass
      element_name_full="Hydrogen"
     CASE ("F")
      old_mass=mass_fluorine
      mass_fluorine=new_mass
      element_name_full="Fluorine"
     CASE ("B")
      old_mass=mass_boron
      mass_boron=new_mass
      element_name_full="Boron"
     CASE ("Cl")
      old_mass=mass_chlorine
      mass_chlorine=new_mass
      element_name_full="Chlorine"
     CASE ("Br")
      old_mass=mass_bromine
      mass_bromine=new_mass
      element_name_full="Bromine"
     CASE ("I")
      old_mass=mass_iodine
      mass_iodine=new_mass
      element_name_full="Iodine"
     CASE ("N")
      old_mass=mass_nitrogen
      mass_nitrogen=new_mass
      element_name_full="Nitrogen"
     CASE ("O")
      old_mass=mass_oxygen
      mass_oxygen=new_mass
      element_name_full="Oxygen"
     CASE ("C")
      old_mass=mass_carbon
      mass_carbon=new_mass
      element_name_full="Carbon"
     CASE ("S")
      old_mass=mass_sulfur
      mass_sulfur=new_mass
      element_name_full="Sulfur"
     CASE ("P")
      old_mass=mass_phosphorus
      mass_phosphorus=new_mass
      element_name_full="Phosphorus"
     CASE ("Li")
      old_mass=mass_lithium
      mass_lithium=new_mass
      element_name_full="Lithium"
     CASE ("D","X")
      old_mass=drude_mass
      drude_mass=new_mass
      element_name_full="drude particle"
     CASE DEFAULT
      CALL report_error(60,exit_status=IACHAR(element_name_input(1:1)))
      RETURN
     END SELECT
     IF (VERBOSE_OUTPUT) WRITE(*,'(" Changing mass of ",A," from ",F0.3," to ",F0.4,"")')&
     &TRIM(element_name_full),old_mass,new_mass
    END SUBROUTINE change_default_mass

    SUBROUTINE read_molecular_input_file_atomic_masses()
    IMPLICIT NONE
    INTEGER :: molecule_type_index,atom_index
     !initialise logicals
     DO molecule_type_index=1,number_of_molecule_types,1
      molecule_list(molecule_type_index)%manual_atom_mass_specified(:)=.FALSE.
     ENDDO
     custom_atom_masses=.FALSE.
     WRITE(*,ADVANCE="NO",FMT='(" Searching for ",A," statement...")') "'atomic_masses'"
     !Skip over headerlines in molecular input file
     REWIND 3
     DO n=1,headerlines_molecular,1
      READ(3,*)
     ENDDO
     !search for masses section.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF (ios<0) THEN
       !end of file encountered
       WRITE(*,'("done (none found, end of file encountered).")')
       EXIT
      ENDIF
      IF (ios==0) THEN
       !support for synonyms
       IF (TRIM(inputstring)=="atom_masses") inputstring="atomic_masses"
       IF (TRIM(inputstring)=="atomic_masses") THEN
        WRITE(*,'("found in line ",I0,".")') n+headerlines_molecular
        BACKSPACE 3
        READ(3,IOSTAT=ios,FMT=*) inputstring,a
        IF (ios/=0) THEN
         !something went wrong
         CALL report_error(109)
         custom_atom_masses=.FALSE.
        ELSE
         !keyword ok - read the section.
         IF (a>0) THEN
          custom_atom_masses=.TRUE.
          WRITE(*,FMT='(" Trying to read ",I0," custom atomic masses...")',ADVANCE="NO") a
          DO m=1,a,1
           READ(3,IOSTAT=ios,FMT=*) molecule_type_index,atom_index,mass_input
           IF (ios/=0) THEN
            !wrong format... abort.
            custom_atom_masses=.FALSE.
            EXIT
           ELSE
            !check molecule_type_index and atom_index
            IF ((molecule_type_index>0).AND.(molecule_type_index<=number_of_molecule_types)) THEN
             !valid molecule type - check atom_index
             IF ((atom_index>0).AND.(atom_index<=molecule_list(molecule_type_index)%number_of_atoms)) THEN
              IF (molecule_list(molecule_type_index)%manual_atom_mass_specified(atom_index)) CALL report_error(112)
              molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)=mass_input
              molecule_list(molecule_type_index)%manual_atom_mass_specified(atom_index)=.TRUE.
             ELSE
              WRITE(*,*)
              CALL report_error(111,exit_status=atom_index)
             ENDIF
            ELSE
             WRITE(*,*)
             CALL report_error(110,exit_status=molecule_type_index)
            ENDIF
           ENDIF
          ENDDO
          !test if custom_atom_masses still true. If not, print error message.
          IF (custom_atom_masses) THEN
           WRITE(*,'("done.")') 
          ELSE
           WRITE(*,'("failed.")') 
           CALL report_error(109)
           DO molecule_type_index=1,number_of_molecule_types,1
            molecule_list(molecule_type_index)%manual_atom_mass_specified(:)=.FALSE.
           ENDDO
          ENDIF
         ENDIF
        ENDIF
        EXIT
       ELSEIF (TRIM(inputstring)=="quit") THEN
        WRITE(*,'("done (none found before ",A,").")') "'quit'"
        EXIT
       ENDIF
      ENDIF
     ENDDO
    END SUBROUTINE read_molecular_input_file_atomic_masses

    SUBROUTINE read_molecular_input_file_constraints()
    IMPLICIT NONE
    INTEGER :: number_of_constraints
     !initialise constraints to zero.
     DO n=1,number_of_molecule_types,1
      molecule_list(n)%constraints=0
     ENDDO
     number_of_constraints=0
     custom_constraints=.FALSE.
     WRITE(*,ADVANCE="NO",FMT='(" Searching for ",A," statement...")') "'constraints'"
     !Skip over headerlines in molecular input file
     REWIND 3
     DO n=1,headerlines_molecular,1
      READ(3,*)
     ENDDO
     !search for constraints section.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF (ios<0) THEN
       !end of file encountered
       WRITE(*,'("done (none found, end of file encountered).")')
       EXIT
      ENDIF
      IF (ios==0) THEN
       IF (TRIM(inputstring)=="constraints") THEN
        WRITE(*,'("found in line ",I0,".")') n+headerlines_molecular
        BACKSPACE 3
        READ(3,IOSTAT=ios,FMT=*) inputstring,a
        IF (ios/=0) THEN
         !something went wrong
         custom_constraints=.FALSE.
        ELSE
         !keyword ok - read the section.
         IF (a>0) THEN
          custom_constraints=.TRUE.
          WRITE(*,'(" Trying to read ",I0," custom constraints...")') a
          DO m=1,a,1
           READ(3,IOSTAT=ios,FMT=*) b,c
           !do some fools proof checks
           IF (ios/=0) THEN
            !wrong format... abort.
            custom_constraints=.FALSE.
            EXIT
           ELSE!format is formally correct. Check for sensible values.
            IF ((b>0).AND.(b<=number_of_molecule_types)) THEN
             IF (molecule_list(b)%constraints/=0) THEN
              CALL report_error(65,exit_status=b)
             ELSE
              IF (c/=0) number_of_constraints=number_of_constraints+1
             ENDIF
             molecule_list(b)%constraints=c
            ELSE
             CALL report_error(64,exit_status=(n+headerlines_molecular+m))
            ENDIF
           ENDIF
          ENDDO
          !test if custom_constraints still true. If not, print error message.
          IF (custom_constraints) THEN
           WRITE(*,'(" ...done.")') 
          ELSE
           WRITE(*,'(" ...failed.")') 
           CALL report_error(63)
           DO m=1,number_of_molecule_types,1
            molecule_list(m)%constraints=0
           ENDDO
          ENDIF
         ENDIF
        ENDIF
        EXIT
       ELSEIF (TRIM(inputstring)=="quit") THEN
        WRITE(*,'("done (none found before ",A,").")') "'quit'"
        EXIT
       ENDIF
      ENDIF
     ENDDO
     IF ((VERBOSE_OUTPUT).AND.(number_of_constraints>0)) THEN
      WRITE(*,'(" List of ",I0," (nonzero) constraints:")') number_of_constraints
      DO n=1,number_of_molecule_types,1
       IF (molecule_list(n)%constraints/=0) WRITE(*,'("   Molecule ",I0," has ",I0," constraints.")')&
       &n,molecule_list(n)%constraints
      ENDDO
     ENDIF
    END SUBROUTINE read_molecular_input_file_constraints

    SUBROUTINE read_molecular_input_file_drudes()
    IMPLICIT NONE
    INTEGER :: inputinteger
     drudes_assigned=.FALSE.
     ndrudes_check=0
     WRITE(*,ADVANCE="NO",FMT='(" Searching for ",A," statement...")') "'drudes'"
     !Skip over headerlines in molecular input file
     REWIND 3
     DO n=1,headerlines_molecular,1
      READ(3,*)
     ENDDO
     !search for drudes section.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF (ios<0) THEN
       !end of file encountered
       WRITE(*,'("done (none found, end of file encountered).")')
       EXIT
      ENDIF
      IF (ios==0) THEN
       IF (TRIM(inputstring)=="drudes") THEN
        WRITE(*,'("found in line ",I0,".")') n+headerlines_molecular
        BACKSPACE 3
        READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
        IF (ios/=0) THEN
         !something went wrong
         drudes_assigned=.FALSE.
         !couldn't read integer
         WRITE(*,'("failed.")') 
        ELSE
         !keyword ok - read the section.
         CALL allocate_drude_list()
         IF (inputinteger>0) THEN
          drudes_assigned=.TRUE.
          IF (VERBOSE_OUTPUT) WRITE(*,FMT='(" Trying to read ",I0," custom drude assignments...")',ADVANCE="NO") inputinteger
          DO m=1,inputinteger,1
           !read in:
           ! a - molecule type index
           ! b - atom index (core)
           ! c - atom index (drude)
           READ(3,IOSTAT=ios,FMT=*) a,b,c
           IF (ios/=0) THEN
            !wrong format... abort.
            drudes_assigned=.FALSE.
            EXIT
           ELSE
            !Add core/drude pair to molecule type index 'a'
            molecule_list(a)%list_of_drude_pairs(b)%drude_flag=c
            !The drude_flag of the drude particle itself is set to zero.
            molecule_list(a)%list_of_drude_pairs(c)%drude_flag=0
            molecule_list(a)%number_of_drudes_in_molecule=molecule_list(a)%number_of_drudes_in_molecule+1
            ndrudes_check=ndrudes_check+molecule_list(a)%total_molecule_count
           ENDIF
          ENDDO
          !test if drudes_assigned still true. If not, print error message.
          IF (drudes_assigned) THEN
           WRITE(*,'("done.")')
          ELSE
           WRITE(*,'("failed.")') 
           CALL report_error(87)
          ENDIF
         ELSE
          WRITE(*,'("failed.")') 
         ENDIF
        ENDIF
        EXIT
       ELSEIF (TRIM(inputstring)=="quit") THEN
        WRITE(*,'("done (none found before ",A,").")') "'quit'"
        EXIT
       ENDIF
      ENDIF
     ENDDO
    END SUBROUTINE read_molecular_input_file_drudes

  END SUBROUTINE initialise_molecular

  !finalises the molecular module.
  SUBROUTINE finalise_molecular()
  IMPLICIT NONE
  INTEGER :: deallocstatus,n
   DO n=1,number_of_molecule_types,1
    IF (drudes_assigned) THEN
     !deallocate memory for detailed drude pair list
     DEALLOCATE(molecule_list(n)%list_of_drude_pairs,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    ENDIF
    DEALLOCATE(molecule_list(n)%list_of_elements,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    DEALLOCATE(molecule_list(n)%list_of_atom_masses,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    DEALLOCATE(molecule_list(n)%manual_atom_mass_specified,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    IF (READ_SEQUENTIAL) THEN
     DEALLOCATE(molecule_list(n)%snapshot,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    ELSE
     DEALLOCATE(molecule_list(n)%trajectory,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    ENDIF
   ENDDO
   drudes_assigned=.FALSE.
   drudes_allocated=.FALSE.
   COM_mass_list(:)=0.0d0
   drude_mass=0.0d0
   DEALLOCATE(molecule_list,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
   IF (dihedrals_initialised) THEN
    DEALLOCATE(dihedral_member_indices,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
   ENDIF
   IF (fragments_initialised) THEN
    DEALLOCATE(fragment_list_base,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
    DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
   ENDIF
   IF (READ_SEQUENTIAL) THEN
    IF (VERBOSE_OUTPUT) WRITE(*,*) "closing file '",TRIM(FILENAME_MOLECULAR_INPUT),"'"!no input path is added for the molecular file!
    CLOSE(UNIT=9)
   ENDIF
  END SUBROUTINE finalise_molecular

  !reports properties: Box size, density, molecule types and masses, formulae, charge...
  SUBROUTINE report_trajectory_properties()
  IMPLICIT NONE
  INTEGER :: molecule_type_index
  REAL(KIND=GENERAL_PRECISION) :: volume,box_weight
  CHARACTER(LEN=8) :: chargestring
  LOGICAL :: comtraj !test if the trajectory format satisfies the centre-of-mass output...
   comtraj=.TRUE.
   box_weight=0.0d0
   WRITE(*,*) "General information read from the trajectory:"
   DO molecule_type_index=1,number_of_molecule_types,1
    IF (comtraj) THEN !format not yet violated - test further!
     IF (molecule_list(molecule_type_index)%number_of_atoms==1) THEN
      IF (.NOT.(ANY(ALPHABET_small==IACHAR(molecule_list(molecule_type_index)%list_of_elements(1)(1:1))))) THEN
       comtraj=.FALSE. !cannot be COM output, because uppercase letter!
      ENDIF
     ELSE
      comtraj=.FALSE. !cannot be COM output, because more (or less) then one 'atom'!
     ENDIF
    ENDIF
    SELECT CASE (molecule_list(molecule_type_index)%charge)
    CASE (0)
     chargestring="neutral "
    CASE (-1)
     chargestring="an anion"
    CASE (1)
     chargestring="a cation"
    CASE DEFAULT
     chargestring="an ion"
    END SELECT
    IF (molecule_list(molecule_type_index)%mass<999.0d0) THEN
     WRITE(*,'(3A,I0,3A,I0,A,F7.3,A)') "   Molecule ",TRIM(give_sum_formula(molecule_type_index)),&
     &" (#",molecule_type_index, ") is ",TRIM(chargestring),&
     &" with ",molecule_list(molecule_type_index)%number_of_atoms," atoms and mass = ",&
     &molecule_list(molecule_type_index)%mass," Da."
    ELSE
     WRITE(*,'(3A,I0,3A,I0,A,E12.6,A)') "   Molecule ",TRIM(give_sum_formula(molecule_type_index)),&
     &" (#",molecule_type_index, ") is ",TRIM(chargestring),&
     &" with ",molecule_list(molecule_type_index)%number_of_atoms," atoms and mass = ",&
     &molecule_list(molecule_type_index)%mass," Da."
    ENDIF
    box_weight=box_weight+molecule_list(molecule_type_index)%mass*molecule_list(molecule_type_index)%total_molecule_count
   ENDDO
   IF (DEVELOPERS_VERSION) CALL report_element_lists()
   !Test for comtraj complete!
   IF ((comtraj).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  Trajectory is in centre-of-mass layout. Check if above masses are correct!"
   IF (TRAJECTORY_TYPE=="lmp") THEN
    WRITE(*,*) "  Box dimensions:"
    IF ((MINVAL(box_dimensions(:,:))<-99.9).OR.(MAXVAL(box_dimensions(:,:))>+99.9)) THEN
     WRITE(*,'(A,2E14.6)') "     x: ",box_dimensions(:,1)
     WRITE(*,'(A,2E14.6)') "     y: ",box_dimensions(:,2)
     WRITE(*,'(A,2E14.6)') "     z: ",box_dimensions(:,3)
    ELSE
     WRITE(*,'(A,2F8.3)') "     x: ",box_dimensions(:,1)
     WRITE(*,'(A,2F8.3)') "     y: ",box_dimensions(:,2)
     WRITE(*,'(A,2F8.3)') "     z: ",box_dimensions(:,3)
    ENDIF
    volume=give_box_volume()
    WRITE(*,'(A,E9.3,A)') "   Box volume is ",volume," cubic Angströms"
    WRITE(*,'(A,F5.3,A)') "   Density is ",(box_weight/volume)*(1d24/avogadro)," g/mL"
    IF (number_of_drude_particles/=0) WRITE(*,'(A,I0,A)') "   Detected ",number_of_drude_particles," drude particles"
   ENDIF
  END SUBROUTINE report_trajectory_properties

  !This subroutine is responsible for loading the whole (lammps) trajectory.
  SUBROUTINE load_trajectory()
  IMPLICIT NONE
  LOGICAL :: file_exists,connected
  INTEGER :: ios,stepcounter,dummy,molecule_type_index,atom_index,molecule_index
  CHARACTER(LEN=2) :: element_name
   INQUIRE(FILE=TRIM(PATH_TRAJECTORY)//TRIM(FILENAME_TRAJECTORY),EXIST=file_exists)
   IF (file_exists) THEN
    IF (VERBOSE_OUTPUT) THEN
     IF (READ_SEQUENTIAL) THEN
      WRITE(*,*) "trajectory file will be read sequentially (needs less RAM, but slow)."
      WRITE(*,*)
      WRITE(*,*) "opening file '",TRIM(PATH_TRAJECTORY)//TRIM(FILENAME_TRAJECTORY),"'"
     ELSE
      WRITE(*,*) "load complete trajectory into RAM. Very fast for some analyses, like diffusion."
      WRITE(*,*)
      WRITE(*,*) "reading file '",TRIM(PATH_TRAJECTORY)//TRIM(FILENAME_TRAJECTORY),"'"
     ENDIF
    ENDIF
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    OPEN(UNIT=3,FILE=TRIM(PATH_TRAJECTORY)//TRIM(FILENAME_TRAJECTORY),ACTION='READ',IOSTAT=ios)
    IF (ios/=0) CALL report_error(38,exit_status=ios)
    !Here starts the part where the trajectory is actually read in. The rest is error handling etc.
    !first, read the header of the trajectory file to get box sizes.
    CALL load_trajectory_header_information()
    IF (READ_SEQUENTIAL) THEN
     CLOSE(UNIT=3)
     IF ((WRAP_TRAJECTORY).AND.(VERBOSE_OUTPUT)) THEN
      WRITE(*,*) "centres of mass of each specified molecule will be wrapped into box."
     ENDIF
     !The trajectory will be read step-wise. Thus, we have to initialise the sequential access.
     CALL initialise_sequential_access()
    ELSE
     !Read the whole trajectory. This is quite some effort.
     CALL load_trajectory_body()
     !IF necessary, wrap trajectory
     IF (WRAP_TRAJECTORY) THEN
      IF (VERBOSE_OUTPUT) WRITE(*,*) "Wrapping centres of mass of each specified molecule into box *now*."
      CALL wrap_full()
     ENDIF
     !At this point, the trajectory should be available.
     CLOSE(UNIT=3)
    ENDIF
    !Important: keep assign_drudes here, because it needs the first step to be read.
    IF (number_of_drude_particles/=0) CALL assign_drudes()
   ELSE
    CALL report_error(9)
   ENDIF
   CONTAINS

    SUBROUTINE load_trajectory_header_information()
    IMPLICIT NONE
    CHARACTER(LEN=32) :: dummystring,edump,xdump,ydump,zdump
    REAL(KIND=SP) :: element_mass
    INTEGER :: current_atom_index,n,skipped_masses
     SELECT CASE (TRAJECTORY_TYPE)
     CASE ("lmp")
      BOX_VOLUME_GIVEN=.TRUE.
      headerlines_to_skip=9
      READ(3,IOSTAT=ios,FMT=*)
      IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
      READ(3,IOSTAT=ios,FMT=*)
      IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
      READ(3,IOSTAT=ios,FMT=*)
      IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
      READ(3,IOSTAT=ios,FMT=*) dummy
      IF (ios/=0) THEN
       IF (ERROR_CODE/=71) CALL report_error(71)
      ELSE
       IF (total_number_of_atoms/=dummy) CALL report_error(10)
      ENDIF
      READ(3,IOSTAT=ios,FMT=*)
      IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
      DO n=1,3,1
       READ(3,IOSTAT=ios,FMT=*) box_dimensions(:,n)
       IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
      ENDDO
      !initialise box size
      box_size(:)=box_dimensions(2,:)-box_dimensions(1,:)
      maximum_distance_squared=box_size(2)**2+SQRT(box_size(1)**2+box_size(3)**2)
      maximum_distance=SQRT(maximum_distance_squared)
      READ(3,IOSTAT=ios,FMT=*) dummystring,dummystring,edump,xdump,ydump,zdump
      IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
      IF (TRIM(edump)/="element") CALL report_error(53)
      IF ((TRIM(xdump)=="vx").AND.(TRIM(ydump)=="vy").AND.(TRIM(zdump)=="vz")) THEN
       INFORMATION_IN_TRAJECTORY="VEL"
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory seems to contain velocities."
      ELSEIF ((TRIM(xdump)=="xu").AND.(TRIM(ydump)=="yu").AND.(TRIM(zdump)=="zu")) THEN
       INFORMATION_IN_TRAJECTORY="POS"
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Trajectory seems to contain Cartesian coordinates."
      ELSE
       INFORMATION_IN_TRAJECTORY="UNK"
       CALL report_error(54)
      ENDIF
     CASE ("xyz")
      BOX_VOLUME_GIVEN=.FALSE.
      headerlines_to_skip=2
      READ(3,IOSTAT=ios,FMT=*) dummy
      IF (ios/=0) THEN
       IF (ERROR_CODE/=71) CALL report_error(71)
      ELSE
       IF (total_number_of_atoms/=dummy) CALL report_error(10)
      ENDIF
      READ(3,IOSTAT=ios,FMT=*)
      IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
     CASE DEFAULT
      CALL report_error(0)!unknown trajectory format, which should never be passed to this subroutine.
     END SELECT
     !define the number of lines to skip when advancing through the trajectory file
     lines_to_skip=headerlines_to_skip+total_number_of_atoms
     number_of_drude_particles=0
     !skipped_masses counts how many masses were already specified manually.
     skipped_masses=0
     !Get the elements - assumes consistent ordering
     DO molecule_type_index=1,number_of_molecule_types,1
      IF ((ERROR_CODE)==70)  ERROR_CODE=ERROR_CODE_DEFAULT
      DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
       READ(3,IOSTAT=ios,FMT=*) element_name
       IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
       IF ((TRIM(element_name)=="X").OR.(TRIM(element_name)=="D")) number_of_drude_particles=number_of_drude_particles+1
       molecule_list(molecule_type_index)%list_of_elements(atom_index)=element_name
       element_mass=atomic_weight(element_name)
       IF (molecule_list(molecule_type_index)%manual_atom_mass_specified(atom_index)) THEN
        skipped_masses=skipped_masses+1
        element_mass=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
       ELSE
        molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)=element_mass
       ENDIF
       IF ((element_mass<0.001d0).AND.(ERROR_CODE/=62)) CALL report_error(62,exit_status=atom_index)
       molecule_list(molecule_type_index)%mass=molecule_list(molecule_type_index)%mass+element_mass
      ENDDO
      IF (molecule_list(molecule_type_index)%total_molecule_count/=0) THEN
       dummy=molecule_list(molecule_type_index)%number_of_atoms*(molecule_list(molecule_type_index)%total_molecule_count-1)
       DO atom_index=1,dummy,1!skip over the remaining ones
        READ(3,IOSTAT=ios,FMT=*) element_name
        IF ((ios/=0).AND.(ERROR_CODE/=71)) CALL report_error(71)
        IF ((TRIM(element_name)=="X").OR.(TRIM(element_name)=="D")) number_of_drude_particles=number_of_drude_particles+1
        !While 'skipping', also check if there are some violations so far.
        current_atom_index=(MOD(atom_index-1,molecule_list(molecule_type_index)%number_of_atoms)+1)
        IF (TRIM(molecule_list(molecule_type_index)%list_of_elements(current_atom_index))/=TRIM(element_name)) THEN
         !element mismatch - probably the wrong trajectory or molecular input file!
         IF ((ERROR_CODE)/=70) CALL report_error(70,exit_status=molecule_type_index)
        ENDIF
       ENDDO
      ENDIF
     ENDDO
     IF ((ndrudes_check/=0).AND.(number_of_drude_particles>0)) THEN
      IF (ndrudes_check/=number_of_drude_particles) CALL report_error(88)
     ENDIF
     REWIND 3
     IF (DEVELOPERS_VERSION) WRITE(*,'("  ! skipped ",I0," previously specified atomic masses.")') skipped_masses
     IF ((skipped_masses>0).AND.(.NOT.(custom_atom_masses))) CALL report_error(0)
     CALL report_trajectory_properties() !report all the useful general properties
    END SUBROUTINE load_trajectory_header_information

    SUBROUTINE load_trajectory_body()
    IMPLICIT NONE
    INTEGER :: item_timestep,current,previous,step
    LOGICAL :: use_scaling
     use_scaling=.FALSE.
     previous=-1
     current=-1
     step=-1
     IF (VERBOSE_OUTPUT) THEN
      WRITE(*,FMT='(A26)',ADVANCE="NO") " reading the trajectory..."
      IF (number_of_steps>100) WRITE(*,*)
      CALL print_progress(number_of_steps)
     ENDIF
     !Flush I/O to ease identification of bottlenecks
     CALL refresh_IO()
     REWIND 3
     !The outer loop iterates over the timesteps.
     DO stepcounter=1,number_of_steps,1 !gives third dimension of trajectory
      !first, skip the head of the trajectory
      DO dummy=1,headerlines_to_skip,1
       IF ((dummy==2).AND.(TRAJECTORY_TYPE=="lmp")) THEN
        READ(3,IOSTAT=ios,FMT=*) item_timestep
        IF (ios/=0) THEN
         IF (ERROR_CODE/=107) THEN
          IF (VERBOSE_OUTPUT) WRITE(*,'("issue at step ",I0,", Headerline ",I0,".")') stepcounter,dummy
          CALL report_error(107,exit_status=stepcounter)
         ENDIF
        ELSE
         !check timestep item for consistency
         previous=current
         current=item_timestep
         IF (((previous+step)/=current).AND.(previous>0)) THEN
          !at least the first two steps have been read, and there has been a change in 'step'.
          IF (step>0) THEN
           !'step' has already been initialised - there must have been an issue!
           use_scaling=.FALSE.
           IF (ERROR_CODE/=108) THEN
            IF (VERBOSE_OUTPUT) WRITE(*,'("issue at step ",I0,", Headerline ",I0,".")')&
            &stepcounter,dummy
            CALL report_error(108)
           ENDIF
          ELSE
           !step not yet initialised. ideally, this situation is encountered exactly once.
           use_scaling=.TRUE.
           step=current-previous
          ENDIF
         ENDIF
        ENDIF
        CYCLE
       ENDIF
       READ(3,IOSTAT=ios,FMT=*)
       IF (ios/=0) THEN
        IF (VERBOSE_OUTPUT) WRITE(*,'("stopped at step ",I0,", Headerline ",I0,".")') stepcounter,dummy
        CALL report_error(84,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,*) "Note that the memory for the trajectory has already been allocated."
        WRITE(*,'(" Only ",I0," steps will be used (specified: ",I0," steps)")') stepcounter-1,number_of_steps
        number_of_steps=stepcounter-1
        WRITE(*,'(" number_of_steps reset to ",I0,"!")') number_of_steps
        RETURN
       ENDIF
      ENDDO
      !THEN, read one molecule type after the other:
      DO molecule_type_index=1,number_of_molecule_types,1
       !For each molecule type, read the corresponding number of molecules:
       DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory
        !Finally, iterate over the atoms in that particular molecule:
        DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1 !gives first dimension of trajectory
         !LOOP VARIABLES:
         !stepcounter: current timestep, e.g. 1234
         !molecule_type_index: current molecule type, e.g. 1 (only 1 and 2 for a binary IL)
         !molecule_index: current explicit molecule, e.g. molecule number 231
         !atom_index: current atom in that molecule, e.g. atom 2 (in molecule number 231 of type 1 in timestep 1234...)
         READ(3,*) element_name,&
         &molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,stepcounter)%coordinates
        ENDDO
       ENDDO 
      ENDDO
      !Show progress
      CALL print_progress()
     ENDDO
     IF (VERBOSE_OUTPUT) THEN
      IF (number_of_steps>100) THEN
       WRITE(*,*)
       WRITE(*,FMT='(" ")',ADVANCE="NO")
      ENDIF
      WRITE(*,'("done.")')
     ENDIF
     !set TIME_SCALING_FACTOR
     IF ((use_scaling).AND.(step>1)) THEN
      WRITE(*,'(" Setting time scaling factor to ",I0," based on ITEM: TIMESTEP")') step
      TIME_SCALING_FACTOR=step
     ENDIF
    END SUBROUTINE load_trajectory_body

    SUBROUTINE initialise_sequential_access()
    IMPLICIT NONE
     !Since sequential access has been requested, the trajectory file will be kept open for a long time (maybe days).
     !There is one and only one unit reserved for that, which is number 9.
     INQUIRE(UNIT=9,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=9)
     OPEN(UNIT=9,FILE=TRIM(PATH_TRAJECTORY)//TRIM(FILENAME_TRAJECTORY),ACTION="READ",STATUS="OLD",IOSTAT=ios,POSITION="REWIND")!Just to be on the safe side. You never know.
     IF (ios/=0) CALL report_error(38,exit_status=ios)
     file_position=0!set flag for subroutine. Has to be zero so the first step is read.
     CALL goto_timestep(1)
    END SUBROUTINE initialise_sequential_access

  END SUBROUTINE load_trajectory

  SUBROUTINE write_header(unit_number,step_number,natoms,output_format)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: unit_number,step_number,natoms
  CHARACTER(LEN=3),INTENT(IN) :: output_format
   !Write head, depending on which type the trajectory has...
   SELECT CASE (output_format)
   CASE ("lmp")
    WRITE(unit_number,'("ITEM: TIMESTEP")')
    WRITE(unit_number,'(I0)') step_number
    WRITE(unit_number,'("ITEM: NUMBER OF ATOMS")')
    WRITE(unit_number,'(I0)') natoms
    WRITE(unit_number,'("ITEM: BOX BOUNDS pp pp pp")')
    WRITE(unit_number,*) box_dimensions(:,1)
    WRITE(unit_number,*) box_dimensions(:,2)
    WRITE(unit_number,*) box_dimensions(:,3)
    !Append the line that tells the user about the content!
    SELECT CASE (INFORMATION_IN_TRAJECTORY)
    CASE("UNK")
     WRITE(unit_number,'("ITEM: ATOMS element x? y? z?")')
    CASE("VEL")
     WRITE(unit_number,'("ITEM: ATOMS element vx vy vz")')
    CASE("POS")
     WRITE(unit_number,'("ITEM: ATOMS element xu yu zu")')
    CASE DEFAULT
     CALL report_error(0)
    END SELECT
   CASE ("xyz")
    WRITE(unit_number,'(I0)') natoms
    !Append the line that tells the user about the content!
    SELECT CASE (INFORMATION_IN_TRAJECTORY)
    CASE("UNK")
     WRITE(unit_number,'("!Unknown content!")')
    CASE("VEL")
     WRITE(unit_number,'("Center-of-Mass velocities:")')
    CASE("POS")
     WRITE(unit_number,'("Center-of-Mass positions:")')
    CASE DEFAULT
     CALL report_error(0)
    END SELECT
   CASE DEFAULT
    CALL report_error(0)!unknown trajectory output format, which should never be passed to this subroutine.
   END SELECT
  END SUBROUTINE write_header

  REAL(KIND=SP) FUNCTION atomic_weight(element_name) !this function returns the atomic weight for a given element.
  IMPLICIT NONE
  !If you change this part, then also change Module_SETTINGS and change_atomic_mass!
  CHARACTER(LEN=*),INTENT(IN) :: element_name
   SELECT CASE (TRIM(element_name))
   CASE ("H")
    atomic_weight=mass_hydrogen
   CASE ("F")
    atomic_weight=mass_fluorine
   CASE ("B")
    atomic_weight=mass_boron
   CASE ("Cl")
    atomic_weight=mass_chlorine
   CASE ("Br")
    atomic_weight=mass_bromine
   CASE ("I")
    atomic_weight=mass_iodine
   CASE ("N")
    atomic_weight=(mass_nitrogen) !IF you change this part, THEN change Module_Main, too!
   CASE ("O")
    atomic_weight=(mass_oxygen) !IF you change this part, THEN change Module_Main, too!
   CASE ("C")
    atomic_weight=(mass_carbon) !IF you change this part, THEN change Module_Main, too!
   CASE ("S")
    atomic_weight=(mass_sulfur) !IF you change this part, THEN change Module_Main, too!
   CASE ("P")
    atomic_weight=(mass_phosphorus) !IF you change this part, THEN change Module_Main, too!
   CASE ("Li")
    atomic_weight=(mass_lithium) !IF you change this part, THEN change Module_Main, too!
   CASE ("X")
    atomic_weight=(drude_mass) !IF you change this part, THEN change Module_Main, too!
   CASE ("D")
    atomic_weight=(drude_mass) !IF you change this part, THEN change Module_Main, too!
   CASE DEFAULT
    !the 'convert' keyword produces a trajectory with a,b,c,...,z as element names.
    IF (ANY(ALPHABET_small==IACHAR(element_name(1:1)))) THEN
     atomic_weight=COM_mass_list(IACHAR(element_name(1:1)))
    ELSE
     CALL report_error(4)
    ENDIF
   END SELECT
  END FUNCTION atomic_weight

END MODULE MOLECULAR
!--------------------------------------------------------------------------------------------------------------------------------!

!This module contains procedures for debugging and technical purposes.
MODULE DEBUG ! Copyright (C) 2020 Frederik Philippi
    USE SETTINGS
 USE MOLECULAR
 IMPLICIT NONE
 PRIVATE
 !variables
 !$ REAL(8) :: timeline_begin_real=0.0d0
 INTEGER :: timeline_begin
 !PRIVATE/PUBLIC declarations
 PUBLIC center_xyz,timing,dump_example,dump_snapshot,dump_split,convert,report_temperature,dump_single,report_drude_temperature
 PUBLIC remove_drudes,dump_dimers,contact_distance,dump_cut,report_gyradius,check_timesteps,jump_analysis,check_timestep
 PUBLIC dump_neighbour_traj
 PRIVATE test_dihedrals
 CONTAINS

  !checks whether timesteps are reasonable.
  SUBROUTINE check_timesteps(startstep,endstep)
  IMPLICIT NONE
  INTEGER,INTENT(INOUT) :: startstep,endstep
  INTEGER :: step_clip
   IF (startstep>endstep) THEN
    CALL report_error(105)
    step_clip=startstep
    startstep=endstep
    endstep=step_clip
   ENDIF
   IF (startstep<1) THEN
    CALL report_error(57,exit_status=startstep)
    startstep=1
   ELSEIF (startstep>give_number_of_timesteps()) THEN
    CALL report_error(57,exit_status=startstep)
    startstep=give_number_of_timesteps()
   ENDIF
   IF (endstep<1) THEN
    CALL report_error(57,exit_status=endstep)
    endstep=1
   ELSEIF (endstep>give_number_of_timesteps()) THEN
    CALL report_error(57,exit_status=endstep)
    endstep=give_number_of_timesteps()
   ENDIF
  END SUBROUTINE check_timesteps

  !checks whether timesteps are reasonable.
  SUBROUTINE check_timestep(timestep)
  IMPLICIT NONE
  INTEGER,INTENT(INOUT) :: timestep
   IF (timestep<1) THEN
    CALL report_error(57,exit_status=timestep)
    timestep=1
   ELSEIF (timestep>give_number_of_timesteps()) THEN
    CALL report_error(57,exit_status=timestep)
    timestep=give_number_of_timesteps()
   ENDIF
  END SUBROUTINE check_timestep

  !Constructs a histogram of jump velocity with increasing jump time.
  !For use_velocity=FALSE, see 10.1021/acs.jpcb.5b01093
  SUBROUTINE jump_analysis(delta_steps,bin_count_in,molecule_type_index_in,startstep_in,endstep_in,use_velocity,maximum_in)
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$ END INTERFACE
  INTEGER,INTENT(IN) :: delta_steps,bin_count_in,molecule_type_index_in,startstep_in,endstep_in
  REAL,INTENT(IN),OPTIONAL :: maximum_in
  LOGICAL,INTENT(IN) :: use_velocity
  REAL :: maximum,stepsize
  INTEGER :: allocstatus,deallocstatus,ios,bin_count,type_counter
  INTEGER(KIND=DP) :: within_bin,out_of_bounds
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: jump_histogram
   !First, do the fools-proof checks
   IF (molecule_type_index_in>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_in)
    RETURN
   ENDIF
   IF (give_number_of_timesteps()<(delta_steps+1)) THEN
    CALL report_error(103,exit_status=delta_steps+1)
    RETURN
   ENDIF
   IF (bin_count_in<10) THEN
    bin_count=10
    CALL report_error(104,exit_status=10)
   ELSEIF (bin_count_in>1000) THEN
    bin_count=1000
    CALL report_error(104,exit_status=1000)
   ELSE
    bin_count=bin_count_in
   ENDIF
   ALLOCATE(jump_histogram(bin_count,delta_steps),STAT=allocstatus)
   IF (allocstatus/=0) THEN
    CALL report_error(22,exit_status=ios)
    RETURN
   ENDIF
   IF (molecule_type_index_in>0) THEN
    WRITE(*,'(" Molecule type ",I0," (",A,"):")')&
    &molecule_type_index_in,TRIM(give_sum_formula(molecule_type_index_in))
    CALL optimize_stepsize(molecule_type_index_in)
    jump_histogram(:,:)=0
    CALL fill_histogram(molecule_type_index_in)
    CALL write_histogram(molecule_type_index_in)
   ELSE
    DO type_counter=1,give_number_of_molecule_types()
     within_bin=0
     out_of_bounds=0
     WRITE(*,'(" Molecule type ",I0," out of ",I0," (",A,"):")')&
     &type_counter,give_number_of_molecule_types(),TRIM(give_sum_formula(type_counter))
     CALL optimize_stepsize(type_counter)
     jump_histogram(:,:)=0
     CALL fill_histogram(type_counter)
     CALL write_histogram(type_counter)
    ENDDO
   ENDIF
   DEALLOCATE(jump_histogram,STAT=deallocstatus)
   IF (deallocstatus/=0) THEN
    CALL report_error(23,exit_status=deallocstatus)
    RETURN
   ENDIF

   CONTAINS

    SUBROUTINE optimize_stepsize(molecule_type_index)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index
    INTEGER :: molecule_index
    REAL :: jump_velocity,jump_vector(3),jump_gyradius
    LOGICAL :: optimize
     optimize=.TRUE.
     IF (PRESENT(maximum_in)) THEN
      IF (maximum_in>0.0) THEN
       optimize=.FALSE.
       maximum=maximum_in
      ENDIF
     ENDIF
     IF (optimize) THEN
      IF (use_velocity) THEN
       !get maximum velocity from first step
       maximum=0.0
       DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
        jump_vector(:)=give_center_of_mass(2,molecule_type_index,molecule_index)-&
        &give_center_of_mass(1,molecule_type_index,molecule_index)
        jump_velocity=SQRT(SUM(jump_vector(:)**2))/FLOAT(TIME_SCALING_FACTOR)
        IF (jump_velocity>maximum) maximum=jump_velocity
       ENDDO
      ELSE
       !get maximum gyradius from last step
       maximum=0.0
       DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
        jump_gyradius=chain_gyradius(1,delta_steps,molecule_type_index,molecule_index)
        IF (jump_gyradius>maximum) maximum=jump_gyradius
       ENDDO
      ENDIF
     ENDIF
     IF (use_velocity) THEN
      IF (VERBOSE_OUTPUT) WRITE(*,'("   Highest binned velocity "ES9.3", largest time difference ",I0," steps.")')&
      &maximum,delta_steps
     ELSE
      IF (VERBOSE_OUTPUT) WRITE(*,'("   Highest binned gyradius "ES9.3", largest time difference ",I0," steps.")')&
      &maximum,delta_steps
     ENDIF
     stepsize=(maximum)/FLOAT(bin_count)
    END SUBROUTINE optimize_stepsize

    SUBROUTINE fill_histogram(molecule_type_index)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index
    INTEGER :: molecule_index,timestep,time_shift,number_of_timesteps,binpos,jump_histogram_local(bin_count,delta_steps)
    REAL :: jump_velocity,jump_vector(3),jump_time,jump_gyradius
     number_of_timesteps=give_number_of_timesteps()
     !$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
     !$OMP PRIVATE(time_shift,jump_time,molecule_index,jump_vector,jump_velocity,binpos,jump_histogram_local,jump_gyradius)
     !$OMP SINGLE
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,'(A,I0,A)') "   ### Parallel execution on ",OMP_get_num_threads()," threads (jump histogram)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     jump_histogram_local(:,:)=0
     !$OMP DO SCHEDULE(STATIC,1)
     DO timestep=startstep_in,endstep_in-1,1
      !the jump distance will be evaluated between timestep and timestep+time_shift.
      !timestep+time_shift must not exceed the number of timesteps.
      DO time_shift=1,MIN((number_of_timesteps-timestep),delta_steps),1
       jump_time=FLOAT(time_shift*TIME_SCALING_FACTOR)
       DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
        IF (use_velocity) THEN
         !compute jump velocity
         jump_vector(:)=give_center_of_mass(timestep+time_shift,molecule_type_index,molecule_index)-&
         &give_center_of_mass(timestep,molecule_type_index,molecule_index)
         jump_velocity=SQRT(SUM(jump_vector(:)**2))/jump_time
         binpos=(INT(jump_velocity/stepsize)+1)
        ELSE
         !compute radius of gyration of chain
         jump_gyradius=chain_gyradius(timestep,time_shift,molecule_type_index,molecule_index)
         binpos=(INT(jump_gyradius/stepsize)+1)
        ENDIF
        IF (.NOT.((binpos>bin_count).OR.(binpos<0))) THEN
         !$OMP ATOMIC
         within_bin=within_bin+1
         jump_histogram_local(binpos,time_shift)=jump_histogram_local(binpos,time_shift)+1
        ELSE
         !$OMP ATOMIC
         out_of_bounds=out_of_bounds+1
        ENDIF
       ENDDO
      ENDDO
     ENDDO
     !$OMP END DO
     !$OMP CRITICAL
     jump_histogram(:,:)=jump_histogram(:,:)+jump_histogram_local(:,:)
     !$OMP END CRITICAL
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='("   ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
    END SUBROUTINE fill_histogram

    REAL FUNCTION chain_gyradius(timestep,time_shift,molecule_type_index,molecule_index)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: timestep,time_shift,molecule_type_index,molecule_index
    REAL :: line_segment(0:time_shift,3),rgc(3)
    INTEGER :: n
     !first, get the member of the trajectory segment
     DO n=0,time_shift,1
      line_segment(n,:)=give_center_of_mass(timestep+n,molecule_type_index,molecule_index)
     ENDDO
     !Compute radius of gyration of trajectory segment
     DO n=1,3,1
      !first, compute rgc=geometric centre of trajectory segment (equal to centre of mass in this case)
      rgc(n)=SUM(line_segment(:,n))/FLOAT(time_shift+1)
      !Then, subtract rgc to give the vector (rj-rgc)
      line_segment(:,n)=line_segment(:,n)-rgc(n)
     ENDDO
     chain_gyradius=SUM(line_segment(:,1)**2+line_segment(:,2)**2+line_segment(:,3)**2)
     chain_gyradius=SQRT(chain_gyradius/FLOAT(time_shift+1))
    END FUNCTION chain_gyradius

    SUBROUTINE write_histogram(molecule_type_index)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index
    INTEGER :: time_shift,binpos
    REAL :: jump_quantity,jump_probability_histogram(bin_count,delta_steps)
    CHARACTER(LEN=1024) :: fstring
    LOGICAL :: connected
     !normalise histogram
     DO time_shift=1,delta_steps,1
      jump_probability_histogram(:,time_shift)=&
      &FLOAT(jump_histogram(:,time_shift))/FLOAT(SUM(jump_histogram(:,time_shift)))
     ENDDO
     !Print results
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (use_velocity) THEN
      WRITE(fstring,'(2A,I0,A,I0,A,I0)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),&
      &"jump_velocity_probability_type_",molecule_type_index,"_step_",startstep_in,"-",endstep_in
      IF (VERBOSE_OUTPUT) THEN
       WRITE(*,'("   within bin ",I0,", out of bounds ",I0,".")') within_bin,out_of_bounds
       WRITE(*,'(3A)') "   Writing jump velocity histogram (probability) to file '",TRIM(fstring),"'"
      ENDIF
      OPEN(UNIT=3,FILE=TRIM(fstring))
      WRITE(3,'(A,I0,A)') "This file contains the probabilities of any molecule of type ",molecule_type_index&
      &," to jump with the indicated velocity."
      WRITE(3,*) "jump_velocity jump_time probability"
      DO binpos=1,bin_count,1
       jump_quantity=FLOAT(binpos)*stepsize-0.5*stepsize
       DO time_shift=1,delta_steps,1
        WRITE(3,*) jump_quantity,time_shift*TIME_SCALING_FACTOR,jump_probability_histogram(binpos,time_shift)
       ENDDO
      ENDDO
     ELSE
      WRITE(fstring,'(2A,I0,A,I0,A,I0)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),&
      &"jump_gyradius_probability_type_",molecule_type_index,"_step_",startstep_in,"-",endstep_in
      IF (VERBOSE_OUTPUT) THEN
       WRITE(*,'("   within bin ",I0,", out of bounds ",I0,".")') within_bin,out_of_bounds
       WRITE(*,'(3A)') "   Writing jump gyradius histogram (probability) to file '",TRIM(fstring),"'"
      ENDIF
      OPEN(UNIT=3,FILE=TRIM(fstring))
      WRITE(3,'(A,I0,A)') "This file contains the probabilities of any molecule of type ",molecule_type_index&
      &," to perform a jump with the given gyradius."
      WRITE(3,*) "jump_gyradius jump_time probability"
      DO binpos=1,bin_count,1
       jump_quantity=FLOAT(binpos)*stepsize-0.5*stepsize
       DO time_shift=1,delta_steps,1
        WRITE(3,*) jump_quantity,time_shift*TIME_SCALING_FACTOR,jump_probability_histogram(binpos,time_shift)
       ENDDO
      ENDDO
     ENDIF
     CLOSE(UNIT=3)
    END SUBROUTINE write_histogram

  END SUBROUTINE jump_analysis

  !This SUBROUTINE writes a trajectory with a certain number of closest neighbours.
  !all closest molecules of type molecule_type_index_2 around all molecules of type molecule_type_index_1.
  SUBROUTINE dump_neighbour_traj(update_com,startstep_in,endstep_in,&
  &molecule_type_index_1,molecule_index_1,molecule_type_index_2,neighbour_num)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: startstep_in,endstep_in,molecule_type_index_1,molecule_type_index_2,molecule_index_1,neighbour_num
  LOGICAL,INTENT(IN) :: update_com
  INTEGER :: natoms,stepcounter
  CHARACTER(LEN=1024) :: fstring,header
  REAL(KIND=WORKING_PRECISION) :: current_centre(3)
  TYPE :: molecule
   INTEGER :: molecule_index
   REAL :: distance
   REAL(KIND=WORKING_PRECISION) :: shift(3)
  END TYPE molecule
  TYPE(molecule),DIMENSION(:),ALLOCATABLE :: neighbour_list
  LOGICAL :: connected
   !First, do the fools-proof checks
   IF (molecule_type_index_1>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_1)
    RETURN
   ELSEIF (molecule_type_index_1<1) THEN
    CALL report_error(33,exit_status=molecule_type_index_1)
    RETURN
   ENDIF
   IF (molecule_type_index_2>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_2)
    RETURN
   ENDIF
   IF ((give_number_of_molecules_per_step(molecule_type_index_2)<neighbour_num)&
   &.OR.((molecule_type_index_2==molecule_type_index_1).AND.&
   &(give_number_of_molecules_per_step(molecule_type_index_2)<(neighbour_num-1)))) THEN
    CALL report_error(97)
    RETURN
   ENDIF
   CALL initialise_neighbourtraj()
   DO stepcounter=startstep_in,endstep_in,1
    IF (update_com) current_centre(:)=give_center_of_mass(stepcounter,molecule_type_index_1,molecule_index_1)
    CALL make_neighbour_list()
    CALL write_neighbours()
   ENDDO
   CALL finalise_neighbourtraj()

  CONTAINS

   SUBROUTINE initialise_neighbourtraj()
   IMPLICIT NONE
   INTEGER :: allocstatus,ios
    ALLOCATE(neighbour_list(give_number_of_molecules_per_step(molecule_type_index_2)),STAT=allocstatus)
    IF (allocstatus/=0) THEN
     CALL report_error(22,exit_status=ios)
     RETURN
    ENDIF
    natoms=give_number_of_atoms_per_molecule(molecule_type_index_1)+&
    &give_number_of_atoms_per_molecule(molecule_type_index_2)*neighbour_num
    INQUIRE(UNIT=4,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=4)
    WRITE(fstring,'(A,I0,A,I0,A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"neighbourtraj_type_",&
    &molecule_type_index_2,"_around_type_",molecule_type_index_1,"_step_",startstep_in,"-",endstep_in,".xyz"
    OPEN(UNIT=4,FILE=TRIM(fstring))
    current_centre(:)=give_center_of_mass(startstep_in,molecule_type_index_1,molecule_index_1)
   END SUBROUTINE initialise_neighbourtraj

   SUBROUTINE finalise_neighbourtraj()
   IMPLICIT NONE
   INTEGER deallocstatus
    DEALLOCATE(neighbour_list,STAT=deallocstatus)
    IF (deallocstatus/=0) THEN
     CALL report_error(23,exit_status=deallocstatus)
     RETURN
    ENDIF
    CLOSE(UNIT=4)
   END SUBROUTINE finalise_neighbourtraj

   SUBROUTINE make_neighbour_list()
   IMPLICIT NONE
   REAL(KIND=WORKING_PRECISION) :: current_shift(3),pivot
   INTEGER :: counter!counter iterates over all possible neighbours.
    !first, get all the neighbours...
    DO counter=1,give_number_of_molecules_per_step(molecule_type_index_2),1
     neighbour_list(counter)%molecule_index=counter
     neighbour_list(counter)%distance=give_smallest_distance_squared&
     &(stepcounter,stepcounter,molecule_type_index_1,molecule_type_index_2,&
     &molecule_index_1,counter,neighbour_list(counter)%shift)
    ENDDO
    !Then, sort it!
    CALL sort_neighbour_list(1,give_number_of_molecules_per_step(molecule_type_index_2))
   END SUBROUTINE make_neighbour_list

   !The following subroutine is a quicksort algorithm to sort the neighbour list.
   RECURSIVE SUBROUTINE sort_neighbour_list(left,right)
   IMPLICIT NONE
   INTEGER left,right,a,b
   REAL pivot
   TYPE(molecule) :: clipboard_element
    IF (left<right) THEN
     pivot=neighbour_list(left)%distance
     a=left
     b=right
     DO
      DO WHILE (neighbour_list(a)%distance<pivot)
       a=a+1
      ENDDO
      DO WHILE (neighbour_list(b)%distance>pivot)
       b=b-1
      ENDDO
      IF (a>=b) EXIT
      !swap elements, unless they're the same
      IF (neighbour_list(a)%distance==neighbour_list(b)%distance) THEN
       b=b-1
       IF (a==b) EXIT
      ELSE
       clipboard_element=neighbour_list(a)
       neighbour_list(a)=neighbour_list(b)
       neighbour_list(b)=clipboard_element
      ENDIF
     ENDDO
     CALL sort_neighbour_list(left,b)
     CALL sort_neighbour_list(b+1,right)
    ENDIF
   END SUBROUTINE sort_neighbour_list

   !When this subroutine is invoked, the neighbour list is ready and sorted.
   SUBROUTINE write_neighbours()
   IMPLICIT NONE
   INTEGER :: counter,ignore
    !write the header to the output trajectory
    WRITE(4,'(I0)') natoms
    IF (stepcounter==startstep_in) THEN
     WRITE(4,'("Neighbours: closest ",I0," of type ",I0," around type ",I0," / index ",I0,", step ",I0,".")')&
     &neighbour_num,molecule_type_index_2,molecule_type_index_1,molecule_index_1,stepcounter
    ELSE
     WRITE(4,'("Step: ",I0)') stepcounter
    ENDIF
    !Check if the molecules are *exactly* the same! could happen because we might also want to dump dimers in, e.g., pure methanol.
    IF (molecule_type_index_1==molecule_type_index_2) THEN
     ignore=1
    ELSE
     ignore=0
    ENDIF
    CALL write_molecule(4,stepcounter,molecule_type_index_1,molecule_index_1,include_header=.FALSE.,translate_by=-current_centre(:))
    DO counter=1+ignore,neighbour_num+ignore,1
     CALL write_molecule(4,stepcounter,molecule_type_index_2,neighbour_list(counter)%molecule_index&
     &,include_header=.FALSE.,translate_by=neighbour_list(counter)%shift(:)-current_centre(:))
    ENDDO
   END SUBROUTINE write_neighbours

  END SUBROUTINE dump_neighbour_traj

  !This SUBROUTINE writes the dimers for a given timestep_in, i.e.
  !all closest molecules of type molecule_type_index_2 around all molecules of type molecule_type_index_1.
  SUBROUTINE dump_dimers(timestep_in,combine_to_single_traj,molecule_type_index_1,molecule_type_index_2)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep_in,molecule_type_index_1,molecule_type_index_2
  LOGICAL,INTENT(IN) :: combine_to_single_traj
  INTEGER :: molecule_index_1,natoms,molecule_index_2
  CHARACTER(LEN=1024) :: fstring,header
  REAL(KIND=WORKING_PRECISION) :: shift(3)
  LOGICAL :: connected
   !First, do the fools-proof checks
   IF (molecule_type_index_1>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_1)
    RETURN
   ELSEIF (molecule_type_index_1<1) THEN
    CALL report_error(33,exit_status=molecule_type_index_1)
    RETURN
   ENDIF
   IF (molecule_type_index_2>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_2)
    RETURN
   ELSEIF (molecule_type_index_2<1) THEN
    CALL report_error(33,exit_status=molecule_type_index_2)
    RETURN
   ENDIF
   IF ((give_number_of_molecules_per_step(molecule_type_index_2)==1)&
   &.AND.(molecule_type_index_2==molecule_type_index_1)) THEN
    CALL report_error(97)
    RETURN
   ENDIF
   natoms=give_number_of_atoms_per_molecule(molecule_type_index_1)+give_number_of_atoms_per_molecule(molecule_type_index_2)
   !First, do the fools-proof checks
   INQUIRE(UNIT=4,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=4)
   !If merged into one file: open the unit here.
   IF (combine_to_single_traj) THEN
    WRITE(fstring,'(A,I0,A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"dimers_type_",&
    &molecule_type_index_2,"_around_type_",molecule_type_index_1,"_step_",timestep_in,".xyz"
    OPEN(UNIT=4,FILE=TRIM(fstring))
    INQUIRE(UNIT=10,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=10)
    OPEN(UNIT=10,STATUS="SCRATCH")
   ENDIF
   DO molecule_index_1=1,give_number_of_molecules_per_step(molecule_type_index_1),1
    molecule_index_2=0
    shift(:)=0.0
    CALL find_closest_partner()
    CALL write_closest()
   ENDDO
   IF (combine_to_single_traj) THEN
    CLOSE(UNIT=4)
    CLOSE(UNIT=10)
   ENDIF

  CONTAINS

   SUBROUTINE find_closest_partner()
   IMPLICIT NONE
   REAL :: smallest_distance_squared,current_distance_squared
   REAL(KIND=WORKING_PRECISION) :: current_shift(3)
   INTEGER :: counter
   smallest_distance_squared=give_maximum_distance_squared()
    DO counter=1,give_number_of_molecules_per_step(molecule_type_index_2),1
     !Check if the molecules are *exactly* the same! could happen because we might also want to dump dimers in, e.g., pure methanol.
     IF ((molecule_index_1==counter).AND.(molecule_type_index_1==molecule_type_index_2)) CYCLE
     current_distance_squared=give_smallest_distance_squared&
     &(timestep_in,timestep_in,molecule_type_index_1,molecule_type_index_2,molecule_index_1,counter,current_shift(:))
     IF (current_distance_squared<smallest_distance_squared) THEN
      !found new best
      molecule_index_2=counter
      smallest_distance_squared=current_distance_squared
      shift(:)=current_shift(:)
     ENDIF
    ENDDO
   END SUBROUTINE find_closest_partner

   !When this subroutine is invoked, the two molecule type indices and molecule indices *are* the dimer pair!
   SUBROUTINE write_closest()
   IMPLICIT NONE
    WRITE(header,'("Dimer: index ",I0," of type ",I0," around index ",I0," of type ",I0,".")')&
    &molecule_index_2,molecule_type_index_2,molecule_index_1,molecule_type_index_1
    IF (combine_to_single_traj) THEN
     REWIND 10
     WRITE(10,'(I0)') natoms
     WRITE(10,*)
     CALL write_molecule(10,timestep_in,molecule_type_index_1,molecule_index_1,include_header=.FALSE.)
     CALL write_molecule(10,timestep_in,molecule_type_index_2,molecule_index_2,include_header=.FALSE.,translate_by=shift(:))
     CALL center_xyz(10,addhead=.TRUE.,outputunit=4,custom_header=header)
    ELSE
     !Write separate files... might be quite a few.
     WRITE(fstring,'(A,I0,A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"dimer_around_type-index_",&
     &molecule_type_index_1,"-",molecule_index_1,"_step_",timestep_in,".xyz"
     OPEN(UNIT=4,FILE=TRIM(fstring))
     WRITE(4,'(I0)') natoms
     WRITE(4,*)
     CALL write_molecule(4,timestep_in,molecule_type_index_1,molecule_index_1,include_header=.FALSE.)
     CALL write_molecule(4,timestep_in,molecule_type_index_2,molecule_index_2,include_header=.FALSE.,translate_by=shift(:))
     CALL center_xyz(4,addhead=.TRUE.,custom_header=header)
     CLOSE(UNIT=4)
    ENDIF
   END SUBROUTINE write_closest

  END SUBROUTINE dump_dimers

  !This SUBROUTINE reports the smallest and largest intramolecular distance and the smallest intermolecular distance for all molecule types in the given timestep
  SUBROUTINE contact_distance(startstep_in,molecule_type_index_in)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: startstep_in,molecule_type_index_in
  INTEGER :: molecule_type_index
  REAL :: smallest,largest
   !availability of box volume check in calling routine.
   IF (molecule_type_index_in>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_in)
    RETURN
   ENDIF
   IF (molecule_type_index_in<1) THEN
    WRITE(*,'(" Calculating intra- and intermolecular contact distances at timestep ",I0," for all molecule types.")') startstep_in
    DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
     WRITE(*,'("   Molecule type index ",I0," out of ",I0,".")') molecule_type_index,give_number_of_molecule_types()
     CALL print_distances()
    ENDDO
   ELSE
    molecule_type_index=molecule_type_index_in
    WRITE(*,FMT='(" Calculating intra- and intermolecular contact distances for molecule type ",I0," at timestep ",I0,".")')&
    &,molecule_type_index,startstep_in
    CALL print_distances()
   ENDIF

   CONTAINS

    SUBROUTINE print_distances()
    IMPLICIT NONE
     CALL give_intramolecular_distances(startstep_in,molecule_type_index,smallest,largest)
     IF (molecule_type_index_in<1) WRITE(*,FMT='("  ")',ADVANCE="NO")
     WRITE(*,'("   Largest intramolecular distance:  ",F0.3)') largest
     IF (molecule_type_index_in<1) WRITE(*,FMT='("  ")',ADVANCE="NO")
     WRITE(*,'("   Smallest intramolecular distance: ",F0.3)') smallest
     CALL give_intermolecular_contact_distance(startstep_in,molecule_type_index,smallest)
     IF (molecule_type_index_in<1) WRITE(*,FMT='("  ")',ADVANCE="NO")
     WRITE(*,'("   Smallest intermolecular distance: ",F0.3)') smallest
    END SUBROUTINE print_distances

  END SUBROUTINE contact_distance

  !SUBROUTINE to center the molecule provided in the specified unit in xyz format.
  !If addhead is .TRUE. THEN a line with the number of atoms and a blank line are added. note that 'addhead' defaults to .FALSE.!
  !If the outputunit is present, THEN it is opened with "append"!
  SUBROUTINE center_xyz(unit_number,addhead,outputunit,custom_header,geometric_center)
  IMPLICIT NONE
  INTEGER,INTENT(IN),OPTIONAL :: outputunit
  INTEGER,INTENT(IN) :: unit_number
  INTEGER :: number_of_atoms,ios
  LOGICAL,OPTIONAL :: addhead,geometric_center
  CHARACTER(LEN=*),OPTIONAL :: custom_header
  TYPE :: atom
   CHARACTER (LEN=1) :: atom_type='X'
   REAL(KIND=WORKING_PRECISION) :: atom_position(3)
   REAL(KIND=WORKING_PRECISION) :: mass
  END TYPE atom
  TYPE(atom) :: center_of_mass
  TYPE(atom),DIMENSION(:),ALLOCATABLE :: molecule
   !first, initialise everything and allocate memory.
   CALL initialize_xyz()
   IF ((ERROR_CODE/=29).AND.(ERROR_CODE/=22)) THEN
    !If no problems were encountered, continue with reading the unit and centering the molecule
    CALL read_center_write()
   ENDIF
   IF (ERROR_CODE/=29) THEN
    CALL finalize()
   ENDIF
   REWIND unit_number
   CONTAINS

    !allocates memory, initializes center of mass, reads number_of_atoms from unit.
    SUBROUTINE initialize_xyz()
    IMPLICIT NONE
    INTEGER allocstatus
     REWIND unit_number
     !First line should contain the number of atoms (at least for proper .xyz format)
     READ(unit_number,IOSTAT=ios,FMT=*) number_of_atoms
     IF (ios/=0) THEN
      CALL report_error(29,exit_status=ios)
      RETURN
     ENDIF
     !THEN a blank line (error handling, because it could as well be EOF)
     READ(unit_number,IOSTAT=ios,FMT=*)
     IF (ios/=0) THEN
      CALL report_error(29,exit_status=ios)
      RETURN
     ENDIF
     !allocate the memory for the molecule in the file.
     ALLOCATE(molecule(number_of_atoms),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(22,exit_status=ios)
      RETURN
     ENDIF
     !initialise the variable for center of mass
     center_of_mass%mass=0.0d0
     center_of_mass%atom_position(1)=0.0d0
     center_of_mass%atom_position(2)=0.0d0
     center_of_mass%atom_position(3)=0.0d0
    END SUBROUTINE initialize_xyz

       !reads molecule from unit, subtracts center of mass, writes to unit
    SUBROUTINE read_center_write()
    IMPLICIT NONE
    INTEGER :: n
    REAL :: weigh
     !Iterate over the body of the XYZ file, read data into 
     DO n=1,number_of_atoms,1
      READ(unit_number,IOSTAT=ios,FMT=*) molecule(n)%atom_type,molecule(n)%atom_position(:)
      IF (ios/=0) THEN
       CALL report_error(29,exit_status=ios)
       RETURN
      ENDIF
      molecule(n)%mass=atomic_weight(molecule(n)%atom_type)
      center_of_mass%mass=center_of_mass%mass+molecule(n)%mass
      weigh=molecule(n)%mass
      IF (PRESENT(geometric_center)) THEN
       !for the centroid, the atomic position are weighed equally.
       IF (geometric_center) weigh=1.0
      ENDIF
      center_of_mass%atom_position(:)=center_of_mass%atom_position(:)+weigh*molecule(n)%atom_position(:)
     END DO
     IF (PRESENT(geometric_center)) THEN
      !for the centroid, normalisation is by the number of atoms
      IF (geometric_center) center_of_mass%mass=FLOAT(number_of_atoms)
     ENDIF
     !normalize sum of positions by total mass, so that the center of mass is obtained
     center_of_mass%atom_position(:)=center_of_mass%atom_position(:)/center_of_mass%mass
     !subtract the centre of mass from atomic coordinates:
     DO n=1,number_of_atoms,1
      molecule(n)%atom_position(:)=molecule(n)%atom_position(:)-center_of_mass%atom_position(:)
     END DO
     !write molecule to the specified unit
     !TO DO not elegant, change order of if statements
     IF (PRESENT(outputunit)) THEN
      IF (PRESENT(addhead)) THEN
       IF (addhead) THEN
        WRITE(outputunit,'(I0)') number_of_atoms
        IF (PRESENT(custom_header)) THEN
         WRITE(outputunit,'(A)') TRIM(ADJUSTL(custom_header))
        ELSE
         WRITE(outputunit,*)
        ENDIF
       ENDIF
      ENDIF
     ELSE
      REWIND unit_number
      IF (PRESENT(addhead)) THEN
       IF (addhead) THEN
        WRITE(unit_number,'(I0)') number_of_atoms
        IF (PRESENT(custom_header)) THEN
         WRITE(unit_number,'(A)') TRIM(ADJUSTL(custom_header))
        ELSE
         WRITE(unit_number,*)
        ENDIF
       ENDIF
      ENDIF
     ENDIF
     DO n=1,number_of_atoms,1
      IF (PRESENT(outputunit)) THEN
       WRITE(outputunit,*) molecule(n)%atom_type,SNGL(molecule(n)%atom_position(:))
      ELSE
       WRITE(unit_number,*) molecule(n)%atom_type,SNGL(molecule(n)%atom_position(:))
      ENDIF
     END DO
     IF (.NOT.(PRESENT(outputunit))) THEN
      WRITE(unit_number,*)
      WRITE(unit_number,*)
      ENDFILE unit_number
      REWIND unit_number
     ENDIF
    END SUBROUTINE read_center_write

    !deallocates memory.
    SUBROUTINE finalize()
    IMPLICIT NONE
    INTEGER allocstatus
     DEALLOCATE(molecule,STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(23,exit_status=allocstatus)
      RETURN
     ENDIF
    END SUBROUTINE finalize

  END SUBROUTINE center_xyz

  SUBROUTINE timing(total)
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_wtime()
  !$  REAL(8) :: OMP_get_wtime
  !$  END FUNCTION OMP_get_wtime
  !$ END INTERFACE
  LOGICAL,INTENT(IN),OPTIONAL :: total
  INTEGER :: clipboard
  !$ REAL(8) :: clipboard_real
  INTEGER,SAVE :: timeline=0
  !$ REAL(8),SAVE :: timeline_real=0.0d0
   !$ clipboard_real=OMP_get_wtime()
   !$ IF (PRESENT(total).AND.(TIME_OUTPUT)) THEN
   !$  IF (total) THEN
   !$   WRITE(*,ADVANCE="NO",FMT='(" TOTAL elapsed time: ")')
   !$   CALL user_friendly_time_output(clipboard_real-timeline_begin_real)
   !$  ENDIF
   !$  WRITE(*,*)
   !$  RETURN
   !$ ENDIF
   !$ IF ((timeline_real>0.0d0).AND.(TIME_OUTPUT)) THEN
   !$  WRITE(*,ADVANCE="NO",FMT='(A)') " elapsed time: "
   !$  CALL user_friendly_time_output(clipboard_real-timeline_real)
   !$ ELSE
   !$  timeline_begin_real=clipboard_real
   !$ ENDIF
   !$ timeline_real=clipboard_real
   !$ WRITE(*,*)
   !If the -fopenmp flag has not been set, THEN OMP_get_wtime is not available, only SYSTEM_CLOCK.
   !$ IF (.FALSE.) THEN
    CALL SYSTEM_CLOCK(clipboard)
    IF ((timeline/=0).AND.(TIME_OUTPUT)) THEN
     WRITE(*,'(A,EN9.1)') " elapsed time: ",DFLOAT(clipboard-timeline)
    ELSE
     timeline_begin=clipboard
    ENDIF
    WRITE(*,*)
    timeline=clipboard
    IF (PRESENT(total)) THEN
     IF (total) WRITE(*,'(A,EN9.1)') " total execution time: ",DFLOAT(clipboard-timeline_begin)
     WRITE(*,*)
    ENDIF
   !$ ENDIF
   !Flush I/O to ease identification of bottlenecks
   CALL refresh_IO()
  END SUBROUTINE timing

  !dumps a snapshot of the given timestep in .xyz format in either separate files or in one file per molecule type.
  SUBROUTINE dump_snapshot(timestep,separate_files)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: timestep
  INTEGER :: n,molecule_type_index
  CHARACTER(LEN=1024) :: fstring
  LOGICAL,INTENT(IN) :: separate_files
  LOGICAL :: connected
   !First, do the fools-proof check
   INQUIRE(UNIT=4,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=4)
   !If merged into one file: open the unit here
   IF (.NOT.(separate_files)) THEN
    WRITE(fstring,'(2A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",timestep,".xyz"
    OPEN(UNIT=4,FILE=TRIM(fstring))
    WRITE(4,'(I0)') give_number_of_atoms_per_step()
    WRITE(4,*)
   ENDIF
   DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
    IF (separate_files) THEN
     DO n=1,give_number_of_molecules_per_step(molecule_type_index),1
      WRITE(fstring,'(2A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"type_",molecule_type_index,"_nr_",n,".xyz"
      OPEN(UNIT=4,FILE=TRIM(fstring),STATUS="REPLACE")
      CALL write_molecule(4,timestep,molecule_type_index,n,include_header=.TRUE.)
      CALL center_xyz(4,addhead=.TRUE.)
      CLOSE(UNIT=4)
     ENDDO
    ELSE
     DO n=1,give_number_of_molecules_per_step(molecule_type_index),1
      CALL write_molecule(4,timestep,molecule_type_index,n,include_header=.FALSE.)
     ENDDO
    ENDIF
   ENDDO
   IF (.NOT.(separate_files)) THEN
    CALL center_xyz(4,addhead=.TRUE.)
    CLOSE(UNIT=4)
   ENDIF
  END SUBROUTINE dump_snapshot

  !dumps a trajectory of a single molecule plus its surrounding neighbours.
  SUBROUTINE dump_cut(use_com,startstep_in,endstep_in,molecule_type_index,molecule_index,cutoff)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: startstep_in,endstep_in,molecule_type_index,molecule_index
  REAL(KIND=WORKING_PRECISION),INTENT(IN) :: cutoff
  LOGICAL,INTENT(IN) :: use_com
  CHARACTER(LEN=1024) :: fstring
  REAL(KIND=WORKING_PRECISION),SAVE :: origin(3)=0.0d0
  LOGICAL :: connected
  INTEGER :: stepcounter
  INTEGER :: number_neighbours,number_of_neighbouring_atoms
   !First, do the fools-proof checks
   IF (molecule_type_index>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index)
    RETURN
   ENDIF
   IF (molecule_type_index<1) THEN
    CALL report_error(33,exit_status=molecule_type_index)
    RETURN
   ENDIF
   IF (molecule_index>give_number_of_molecules_per_step(molecule_type_index)) THEN
    CALL report_error(69,exit_status=molecule_index)
    RETURN
   ENDIF
   IF (molecule_index<1) THEN
    CALL report_error(69,exit_status=molecule_index)
    RETURN
   ENDIF
   !While searching the neighbours, the found atoms will be written into the scratch file in unit 10.
   !THEN, final output will be unit 4 - which is filled directly from unit 10.
   INQUIRE(UNIT=10,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=10)
   OPEN(UNIT=10,STATUS="SCRATCH")
   !open the output file
   INQUIRE(UNIT=4,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=4)
   WRITE(fstring,'(2A,I0,A,I0,A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),&
   &"molecule_",molecule_index,"_type_",molecule_type_index,"_step_",startstep_in,"-",endstep_in,"_neighbours.xyz"
   OPEN(UNIT=4,FILE=TRIM(fstring))
   REWIND 4
   !find suitable origin
   origin(:)=give_center_of_mass(startstep_in,molecule_type_index,molecule_index)
   !iterate over the specified timesteps
   DO stepcounter=startstep_in,endstep_in,1
    !First, add the reference molecule to the xyz file.
    REWIND 10
    IF (use_com) origin(:)=give_center_of_mass(stepcounter,molecule_type_index,molecule_index)
    CALL write_molecule(10,stepcounter,molecule_type_index,molecule_index,include_header=.FALSE.)
    !Search for neighbours and write them into unit 10
    CALL give_number_of_neighbours&
    &(stepcounter,molecule_type_index,molecule_index,number_neighbours,number_of_neighbouring_atoms,cutoff,10)
    !Write the string to pass as custom header later
    WRITE(fstring,'("Timestep nr. ",I0," with cutoff ",F0.2)') stepcounter,cutoff
    CALL transfer_to_output()
   ENDDO
   WRITE(4,*)
   WRITE(4,*)
   ENDFILE 4
   CLOSE(UNIT=4)
   CLOSE(UNIT=10)
   CONTAINS

    !This SUBROUTINE writes the trajectory including neighbours into unit 4. It also wraps and centers, if necessary.
    SUBROUTINE transfer_to_output()
    IMPLICIT NONE
    INTEGER :: atomcounter,natoms
    REAL(KIND=WORKING_PRECISION) :: position_clipboard(3)
    CHARACTER(LEN=2) :: element
     natoms=number_of_neighbouring_atoms+give_number_of_atoms_per_molecule(molecule_type_index)
     !Put reference molecule into origin. Directly transfer from unit 10 to unit 4.
     REWIND 10
     WRITE(4,'(I0)') natoms
     WRITE(4,'(A)') TRIM(fstring)
     DO atomcounter=1,natoms,1
      READ(10,*) element,position_clipboard(:)
      WRITE(4,*) element,SNGL(position_clipboard(:)-origin(:))
     ENDDO
    END SUBROUTINE transfer_to_output

  END SUBROUTINE dump_cut

  !dumps a trajectory of a single molecule.
  SUBROUTINE dump_single(use_com,startstep_in,endstep_in,molecule_type_index,molecule_index)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: startstep_in,endstep_in,molecule_type_index,molecule_index
  LOGICAL,INTENT(IN) :: use_com
  CHARACTER(LEN=1024) :: fstring
  LOGICAL :: connected
  INTEGER :: stepcounter
   !First, do the fools-proof checks
   IF (molecule_type_index>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index)
    RETURN
   ENDIF
   IF (molecule_type_index<1) THEN
    CALL report_error(33,exit_status=molecule_type_index)
    RETURN
   ENDIF
   IF (molecule_index>give_number_of_molecules_per_step(molecule_type_index)) THEN
    CALL report_error(69,exit_status=molecule_index)
    RETURN
   ENDIF
   IF (molecule_index<1) THEN
    CALL report_error(69,exit_status=molecule_index)
    RETURN
   ENDIF
   IF (use_com) THEN
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    OPEN(UNIT=3,STATUS="SCRATCH")
   ENDIF
   !open the output file
   INQUIRE(UNIT=4,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=4)
   WRITE(fstring,'(2A,I0,A,I0,A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),&
   &"molecule_",molecule_index,"_type_",molecule_type_index,"_step_",startstep_in,"-",endstep_in,".xyz"
   OPEN(UNIT=4,FILE=TRIM(fstring))
   !iterate over the specified timesteps
   DO stepcounter=startstep_in,endstep_in,1
    !Write the string to pass as custom header later
    WRITE(fstring,'("Timestep nr. ",I0)') stepcounter
    IF (use_com) THEN
     !if centre-of-mass is desired, THEN the scratch file will be used.
     CALL write_molecule(3,stepcounter,molecule_type_index,molecule_index,include_header=.TRUE.,custom_header=TRIM(fstring))
     CALL center_xyz(3,addhead=.TRUE.,outputunit=4)
    ELSE
     CALL write_molecule(4,stepcounter,molecule_type_index,molecule_index,include_header=.TRUE.,custom_header=TRIM(fstring))
    ENDIF
   ENDDO
   WRITE(4,*)
   WRITE(4,*)
   ENDFILE 4
   CLOSE(UNIT=4)
   IF (use_com) THEN
    CLOSE(UNIT=3)
   ENDIF
  END SUBROUTINE dump_single

  SUBROUTINE dump_example()
  IMPLICIT NONE
  INTEGER :: molecule_type_index
  LOGICAL :: connected
  CHARACTER(LEN=1024) :: fstring
   DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
    WRITE(fstring,'(2A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"type_",molecule_type_index,".xyz"
    INQUIRE(UNIT=4,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=4)
    OPEN(UNIT=4,FILE=TRIM(fstring),STATUS="REPLACE")
    CALL write_molecule(4,-1,molecule_type_index,1,include_header=.TRUE.)
    CALL center_xyz(4,addhead=.TRUE.)
    CLOSE(UNIT=4)
   ENDDO
  END SUBROUTINE dump_example

  SUBROUTINE dump_split(startstep_in,endstep_in)
  IMPLICIT NONE
  INTEGER :: startstep,endstep
  INTEGER :: molecule_type_index,stepcounter,moleculecounter,atomcount
  INTEGER,INTENT(IN) :: startstep_in,endstep_in
  LOGICAL :: connected
  CHARACTER(LEN=1024) :: fstring
   !First, do the fools-proof checks
   startstep=startstep_in
   endstep=endstep_in
   IF (startstep<1) THEN
    CALL report_error(57,exit_status=startstep)
    startstep=1
   ENDIF
   IF (endstep>give_number_of_timesteps()) THEN
    CALL report_error(57,exit_status=endstep)
    endstep=give_number_of_timesteps()
   ENDIF
   IF (endstep<startstep) THEN
    CALL report_error(57,exit_status=endstep)
    endstep=startstep
   ENDIF
   OPEN(UNIT=3,STATUS="SCRATCH")
   DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
    WRITE(*,ADVANCE="NO",FMT='(" Writing separate trajectory for molecule number ",I0," ...")')&
    & molecule_type_index
    WRITE(fstring,'(2A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_",molecule_type_index,".xyz"
    INQUIRE(UNIT=4,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=4)
    OPEN(UNIT=4,FILE=TRIM(fstring),STATUS="REPLACE")
    atomcount=give_number_of_atoms_per_molecule(molecule_type_index)*give_number_of_molecules_per_step(molecule_type_index)
    DO stepcounter=startstep,endstep,1
     WRITE(4,'(I0)') atomcount
     WRITE(4,'("Timestep number ",I0)') stepcounter
     !Reset unit 3
     REWIND 3
     !Add the header
     WRITE(3,*) atomcount
     WRITE(3,*)
     DO moleculecounter=1,give_number_of_molecules_per_step(molecule_type_index),1
      !the following line appends one molecule to the scratch file.
      CALL write_molecule(3,stepcounter,molecule_type_index,moleculecounter,include_header=.FALSE.)
     ENDDO
     !Here, all the molecules for the current timestep have been appended. Thus, transfer to output:
     CALL center_xyz(3,addhead=.FALSE.,outputunit=4)
    ENDDO
    ENDFILE 4
    CLOSE(UNIT=4)
    WRITE(*,*)"done"
   ENDDO
   CLOSE(UNIT=3)
  END SUBROUTINE dump_split

  SUBROUTINE convert(writemolecularinputfile,output_format)
  IMPLICIT NONE
  CHARACTER(LEN=3),INTENT(IN) :: output_format
  INTEGER :: molecule_type_index,stepcounter,moleculecounter,output(3),ncentres
  LOGICAL,INTENT(IN) :: writemolecularinputfile
  LOGICAL :: connected
  CHARACTER(LEN=1024) :: fstring
  CHARACTER(LEN=1) :: element
   IF (DEVELOPERS_VERSION) THEN
    PRINT *," ! CAREFUL: 'PARALLELISED' CONVERTER USED"
    CALL convert_parallel()
    RETURN
   ELSE
    !first, get the number of lines that will be added.
    ncentres=0
    DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
     ncentres=ncentres+give_number_of_molecules_per_step(molecule_type_index)
    ENDDO
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    WRITE(fstring,'(3A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_COM.",output_format
    OPEN(UNIT=3,FILE=TRIM(fstring))
    WRITE(*,ADVANCE="NO",FMT='(" Writing new trajectory to file ",A,"...")') "'"//TRIM(fstring)//"'"
    DO stepcounter=1,give_number_of_timesteps(),1
     !Write head, depending on which type the trajectory has...
     CALL write_header(3,stepcounter,ncentres,output_format)
     DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
      element=CHAR(ALPHABET_small(MOD((molecule_type_index-1),26)+1)) !assign the element names a,b,c,... to the centred molecules.
      DO moleculecounter=1,give_number_of_molecules_per_step(molecule_type_index),1
       !Sort of high accuracy should be kept here because of the way I use this routine.
       !reduced to 16.8 from 18.10
       WRITE(3,'(A1,3E16.8)') element,give_center_of_mass(stepcounter,molecule_type_index,moleculecounter)
      ENDDO
     ENDDO
    ENDDO
    CLOSE(UNIT=3)
    WRITE(*,*) "done"
   ENDIF
   IF (writemolecularinputfile) THEN
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    WRITE(fstring,'(3A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"COM_molecular.inp"
    OPEN(UNIT=3,FILE=TRIM(fstring))
    WRITE(*,ADVANCE="NO",FMT='(" Writing new molecular input file in ",A,"...")') "'"//TRIM(fstring)//"'"
    WRITE(3,'(I0," ### number of timesteps")') give_number_of_timesteps()
    WRITE(3,'(I0," ### number of different types of molecules. Followed by list of molecules.")') give_number_of_molecule_types()
    output(2)=1
    DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
     output(1)=give_charge_of_molecule(molecule_type_index)
     output(3)=give_number_of_molecules_per_step(molecule_type_index)
     WRITE(3,ADVANCE="NO",FMT='(SP,I2,SS," ",I0," ",I0," ### ")') output(:) !write the crucial part
     WRITE(3,'("There are ",I0," molecules per step with charge ",SP,I2,SS,", given as centre of mass.")')& !write the comments
     & output(3),output(1)
    ENDDO
    !write the custom masses section.
    WRITE(3,'("masses ",I0," ### The following lines contain the masses of every molecule type.")')&
    &give_number_of_molecule_types()
    DO molecule_type_index=1,give_number_of_molecule_types(),1
     WRITE(3,'(A1," ",F9.3)') CHAR(ALPHABET_small(MOD((molecule_type_index-1),26)+1)),give_mass_of_molecule(molecule_type_index)
    ENDDO
    CLOSE(UNIT=3)
    WRITE(*,*) "done"
   ENDIF
  END SUBROUTINE convert

  SUBROUTINE report_gyradius(molecule_type_index_in,startstep,endstep)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: molecule_type_index_in,startstep,endstep
  INTEGER :: timestep,molecule_type_index
  REAL(KIND=GENERAL_PRECISION),DIMENSION(:),ALLOCATABLE :: averages_maxdist,averages_rgysquared,averages_rgy
   !first, do all the annoying fools-proof tests...
   IF (molecule_type_index_in>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_in)
    RETURN
   ENDIF
   IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
   CALL initialise_gyradius()
   IF (molecule_type_index_in<1) THEN
    WRITE(*,FMT='(" Calculating radius of gyration for all molecule types.")')
    WRITE(*,'(" Taking ensemble average from step ",I0," to ",I0,".")') startstep,endstep
    DO molecule_type_index=1,give_number_of_molecule_types(),1
     WRITE(*,'("   Molecule type index ",I0," out of ",I0,".")') molecule_type_index,give_number_of_molecule_types()
     CALL print_gyradius()
    ENDDO
   ELSE
    molecule_type_index=molecule_type_index_in
    WRITE(*,FMT='(" Calculating radius of gyration for molecule type ",I0,".")'),molecule_type_index
    WRITE(*,'(" Taking ensemble average from step ",I0," to ",I0,".")') startstep,endstep
    CALL print_gyradius()
   ENDIF
   CALL finalise_gyradius()
   
   CONTAINS

    SUBROUTINE initialise_gyradius()
    IMPLICIT NONE
    INTEGER :: allocstatus,ios
     ALLOCATE(averages_maxdist(startstep:endstep),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(22,exit_status=ios)
      RETURN
     ENDIF
     ALLOCATE(averages_rgysquared(startstep:endstep),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(22,exit_status=ios)
      RETURN
     ENDIF
     ALLOCATE(averages_rgy(startstep:endstep),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(22,exit_status=ios)
      RETURN
     ENDIF
    END SUBROUTINE initialise_gyradius

    SUBROUTINE finalise_gyradius()
    IMPLICIT NONE
    INTEGER deallocstatus
     DEALLOCATE(averages_maxdist,STAT=deallocstatus)
     IF (deallocstatus/=0) THEN
      CALL report_error(23,exit_status=deallocstatus)
      RETURN
     ENDIF
     DEALLOCATE(averages_rgysquared,STAT=deallocstatus)
     IF (deallocstatus/=0) THEN
      CALL report_error(23,exit_status=deallocstatus)
      RETURN
     ENDIF
     DEALLOCATE(averages_rgy,STAT=deallocstatus)
     IF (deallocstatus/=0) THEN
      CALL report_error(23,exit_status=deallocstatus)
      RETURN
     ENDIF
    END SUBROUTINE finalise_gyradius

    !This SUBROUTINE writes <Rgy²>, <Rgy>, and <maxdist> to the screen, including standard deviations.
    SUBROUTINE print_gyradius()
    IMPLICIT NONE
    REAL(KIND=GENERAL_PRECISION) :: rgy_sq,maxdist,globalav_rgy_sq,globalav_rgy,globalav_maxdist
    REAL(KIND=GENERAL_PRECISION) :: stdev_rgy_sq,stdev_rgy,stdev_maxdist
    INTEGER :: molecule_index,timestep
    INTEGER(KIND=WORKING_PRECISION) :: N
    LOGICAL :: connected
    CHARACTER(LEN=1024) :: fstring
    !the local variables in this routine must be declared as private if parallelised
     averages_maxdist(:)=0.0d0
     averages_rgysquared(:)=0.0d0
     averages_rgy(:)=0.0d0
     globalav_maxdist=0.0d0
     globalav_rgy=0.0d0
     globalav_rgy_sq=0.0d0
     stdev_rgy_sq=0.0d0
     stdev_rgy=0.0d0
     stdev_maxdist=0.0d0
     N=give_number_of_molecules_per_step(molecule_type_index)*(endstep-startstep+1)
     IF (DEVELOPERS_VERSION) THEN
      WRITE(fstring,'(2A,I0)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"gyradius_type_",molecule_type_index
      PRINT *," ! dumping values to file ",TRIM(fstring)
      INQUIRE(UNIT=3,OPENED=connected)
      IF (connected) CALL report_error(27,exit_status=3)
      OPEN(UNIT=3,FILE=TRIM(fstring))
      WRITE(3,'("rgy_sq maxdist")')
     ENDIF
     !TO DO Parallelise this? Actually not necessary - is quite fast, and converges quickly
     !First, compute global average of maxdist, rgy, and rgy_squared.
     DO timestep=startstep,endstep,1
      !iterate over molecule indices
      DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
       CALL compute_squared_radius_of_gyration(timestep,molecule_type_index,molecule_index,rgy_sq,maxdist)
       IF (DEVELOPERS_VERSION) WRITE(3,*) rgy_sq,maxdist
       averages_maxdist(timestep)=averages_maxdist(timestep)+maxdist
       averages_rgysquared(timestep)=averages_rgysquared(timestep)+rgy_sq
       averages_rgy(timestep)=averages_rgy(timestep)+SQRT(rgy_sq)
      ENDDO
     ENDDO
     !update memory / OMP flush if parallelised
     IF (DEVELOPERS_VERSION) CLOSE(UNIT=3)
     !correct every element by number of molecules
     averages_maxdist(:)=averages_maxdist(:)/FLOAT(give_number_of_molecules_per_step(molecule_type_index))
     averages_rgysquared(:)=averages_rgysquared(:)/FLOAT(give_number_of_molecules_per_step(molecule_type_index))
     averages_rgy(:)=averages_rgy(:)/FLOAT(give_number_of_molecules_per_step(molecule_type_index))
     !take global averages
     globalav_maxdist=SUM(averages_maxdist(:))/FLOAT(endstep-startstep+1)
     globalav_rgy=SUM(averages_rgy(:))/FLOAT(endstep-startstep+1)
     globalav_rgy_sq=SUM(averages_rgysquared(:))/FLOAT(endstep-startstep+1)
     averages_maxdist(:)=0.0d0
     averages_rgysquared(:)=0.0d0
     averages_rgy(:)=0.0d0
     !Second, compute standard deviation from this average.
     DO timestep=startstep,endstep,1
      !iterate over molecule indices
      DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
       CALL compute_squared_radius_of_gyration(timestep,molecule_type_index,molecule_index,rgy_sq,maxdist)
       !use 'averages' arrays to store the standard deviations.
       averages_maxdist(timestep)=averages_maxdist(timestep)+(maxdist-globalav_maxdist)**2
       averages_rgysquared(timestep)=averages_rgysquared(timestep)+(rgy_sq-globalav_rgy_sq)**2
       averages_rgy(timestep)=averages_rgy(timestep)+(SQRT(rgy_sq)-globalav_rgy)**2
      ENDDO
     ENDDO
     !Sum everything up
     stdev_rgy_sq=SUM(averages_rgysquared(:))
     stdev_rgy=SUM(averages_rgy(:))
     stdev_maxdist=SUM(averages_maxdist(:))
     !Divide by total number (N-1)
     stdev_maxdist=stdev_maxdist/FLOAT(N-1)
     stdev_rgy=stdev_rgy/FLOAT(N-1)
     stdev_rgy_sq=stdev_rgy_sq/FLOAT(N-1)
     !take the square root to arrive at standard deviations.
     stdev_maxdist=SQRT(stdev_maxdist)
     stdev_rgy=SQRT(stdev_rgy)
     stdev_rgy_sq=SQRT(stdev_rgy_sq)
     !Print ensemble average (=global averages) and the standard deviations
     IF (molecule_type_index_in<1) WRITE(*,FMT='("  ")',ADVANCE="NO")
     WRITE(*,'(" ",I0," averages taken. Results:")') N
     CALL formatted_print("   <gyradius> ",globalav_rgy,stdev_rgy)
     CALL formatted_print("   <gyrad**2> ",globalav_rgy_sq,stdev_rgy_sq)
     CALL formatted_print("   <maxdist>  ",globalav_maxdist,stdev_maxdist)
    END SUBROUTINE print_gyradius

    SUBROUTINE formatted_print(inputstring,real1,real2)
    IMPLICIT NONE
    CHARACTER(LEN=14),INTENT(IN) :: inputstring
    REAL(KIND=GENERAL_PRECISION),INTENT(IN) :: real1,real2
     IF (molecule_type_index_in<1) WRITE(*,FMT='("  ")',ADVANCE="NO")
     WRITE(*,'(A14,E11.4,", stdev =",E10.3)') inputstring,real1,real2
    END SUBROUTINE formatted_print

  END SUBROUTINE report_gyradius

  !This SUBROUTINE writes a SUBROUTINE with removed drude particles. This requires assigned drude particles!
  SUBROUTINE remove_drudes(startstep_in,endstep_in)
  IMPLICIT NONE
  INTEGER :: stepcounter,molecule_type_index,molecule_index
  INTEGER,INTENT(IN) :: startstep_in,endstep_in
  LOGICAL :: connected
  CHARACTER(LEN=1024) :: fstring
   !open the output file
   INQUIRE(UNIT=4,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=4)
   IF (startstep_in==endstep_in) THEN
    WRITE(fstring,'(2A,I0,A)') &
    &TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_nodrudes_step_",startstep_in,".xyz"
   ELSE
    WRITE(fstring,'(2A,I0,A,I0,A)') &
    &TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_nodrudes_step_",startstep_in,"-",endstep_in,".xyz"
   ENDIF
   OPEN(UNIT=4,FILE=TRIM(fstring))
   !iterate over the specified timesteps
   DO stepcounter=startstep_in,endstep_in,1
    WRITE(4,'(I0)') (give_number_of_atoms_per_step()-give_number_of_drude_particles())
    WRITE(4,'("Timestep number ",I0)') stepcounter
    DO molecule_type_index=1,give_number_of_molecule_types(),1
     DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
      !the following line appends one molecule to the output trajectory.
      CALL write_molecule_merged_drudes(4,stepcounter,molecule_type_index,molecule_index,include_header=.FALSE.)
     ENDDO
    ENDDO
   ENDDO
   ENDFILE 4
   CLOSE(UNIT=4)
  END SUBROUTINE

  !This SUBROUTINE reports the temperatures as given in Eq. (13), (14) and (15) in 10.1021/acs.jpclett.9b02983
  SUBROUTINE report_drude_temperature(startstep,endstep)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: startstep,endstep
  INTEGER :: Nf,Nmol,ND
  REAL(KIND=WORKING_PRECISION) :: TCM,TR,TD
   !fools-proof tests
   IF (INFORMATION_IN_TRAJECTORY=="POS") CALL report_error(56)
   IF ((startstep==1).AND.(endstep==1)) THEN
    !do the dirty quick print.
    IF (VERBOSE_OUTPUT) WRITE(*,*) "Quick print - no output file will be produced."
    PRINT *,"temperatures in first step, based on equation (13), (14), and (15) in 10.1021/acs.jpclett.9b02983"
    PRINT *,"degrees of freedoms are given in brackets."
    CALL compute_drude_temperature(1,TCM,TR,TD,Nf,Nmol,ND)
    WRITE(*,ADVANCE="NO",FMT='(" TCM: ")')
    CALL print_simple_temperature_output(TCM)
    WRITE(*,'(" (",I0,")")') Nmol
    WRITE(*,ADVANCE="NO",FMT='(" TR:  ")')
    CALL print_simple_temperature_output(TR)
    WRITE(*,'(" (",I0,")")') Nf
    WRITE(*,ADVANCE="NO",FMT='(" TD:  ")')
    CALL print_simple_temperature_output(TD)
    WRITE(*,'(" (",I0,")")') ND
   ELSE
    CALL write_temp_with_drudes()
   ENDIF
   CONTAINS

    SUBROUTINE print_simple_temperature_output(temperature)
    IMPLICIT NONE
    54 FORMAT (EN9.1," K")
    55 FORMAT (F5.1," K")
    REAL(KIND=WORKING_PRECISION),INTENT(IN) :: temperature
     !This routine is just for nice output.
     IF ((temperature<999.0d0).AND.(temperature>1.0d0)) THEN
      WRITE(*,FMT=55,ADVANCE="NO") temperature
     ELSE
      WRITE(*,FMT=54,ADVANCE="NO") temperature
     ENDIF
    END SUBROUTINE print_simple_temperature_output

    SUBROUTINE write_temp_with_drudes()
    IMPLICIT NONE
    INTEGER :: timestep,nsteps
    REAL(KIND=WORKING_PRECISION) :: TCM_average,TR_average,TD_average
    CHARACTER(LEN=1024) :: fstring
    LOGICAL :: connected
     TCM_average=0.0d0
     TR_average=0.0d0
     TD_average=0.0d0
     WRITE(fstring,'(2A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"drude_temp"
     IF (VERBOSE_OUTPUT) WRITE(*,FMT='(3A)',ADVANCE="NO") "   Writing file '",TRIM(fstring),"'..."
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     OPEN(UNIT=3,FILE=TRIM(fstring))
     WRITE(3,*) "This file contains effective temperatures (centre of mass, total, drude temperature)."
     WRITE(3,*) "Reference: equations (13), (14), and (15) in 10.1021/acs.jpclett.9b02983"
     WRITE(3,'(A15,2A11,A12)') "timeline","TCM","TR","TD"
     DO timestep=startstep,endstep,1
      CALL compute_drude_temperature(timestep,TCM,TR,TD,Nf,Nmol,ND)
      WRITE(3,'(I15,2EN11.1,EN12.2)') timestep*TIME_SCALING_FACTOR,TCM,TR,TD
      TCM_average=TCM_average+TCM
      TR_average=TR_average+TR
      TD_average=TD_average+TD
     ENDDO
     CLOSE(UNIT=3)
     IF (VERBOSE_OUTPUT) WRITE(*,'("done.")')
     nsteps=(endstep-startstep+1)
     TCM_average=TCM_average/FLOAT(nsteps)
     TR_average=TR_average/FLOAT(nsteps)
     TD_average=TD_average/FLOAT(nsteps)
     WRITE(*,'(" Average Temperatures over ",I0," steps:")') nsteps
     PRINT *,"degrees of freedoms are given in brackets."
     WRITE(*,ADVANCE="NO",FMT='("   TCM: ")')
     CALL print_simple_temperature_output(TCM)
     WRITE(*,'(" (",I0,")")') Nmol
     WRITE(*,ADVANCE="NO",FMT='("   TR:  ")')
     CALL print_simple_temperature_output(TR)
     WRITE(*,'(" (",I0,")")') Nf
     WRITE(*,ADVANCE="NO",FMT='("   TD:  ")')
     CALL print_simple_temperature_output(TD)
     WRITE(*,'(" (",I0,")")') ND
    END SUBROUTINE write_temp_with_drudes

  END SUBROUTINE report_drude_temperature

  SUBROUTINE report_temperature(molecule_type_index_in,startstep,endstep)
  IMPLICIT NONE
 50 FORMAT ("   T",A1,": ",EN9.1," K (",EN9.1," K)")
 51 FORMAT ("   T",A1,": ",F5.1," K (",F5.1," K)")
  INTEGER,INTENT(IN) :: molecule_type_index_in,startstep,endstep
  INTEGER :: counter
  REAL(KIND=WORKING_PRECISION) :: drift(3),temperature(3),corrected_temperature(3),&
  &kinetic_energy,kinetic_energy_total,total_temp,constraints_correction
   kinetic_energy_total=0.0d0
   !first, do all the annoying fools-proof tests...
   IF (molecule_type_index_in>give_number_of_molecule_types()) THEN
    CALL report_error(33,exit_status=molecule_type_index_in)
    RETURN
   ENDIF
   IF (INFORMATION_IN_TRAJECTORY=="POS") CALL report_error(56)
   IF (molecule_type_index_in<1) THEN
    IF ((startstep==1).AND.(endstep==1)) THEN
     !do the dirty quick print.
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Quick print - no output file will be produced."
     IF (VERBOSE_OUTPUT) WRITE(*,*) "drift-corrected temperatures are given in brackets."
     DO counter=1,give_number_of_molecule_types(),1
      CALL give_temperature(1,drift,counter,temperature,corrected_temperature,&
      &kinetic_energy=kinetic_energy,constraints_correction=constraints_correction)
      kinetic_energy_total=kinetic_energy_total+kinetic_energy
      WRITE(*,'(" Molecule Type ",I0,":")') counter
      IF ((MAXVAL(temperature(:))<999.0d0).AND.(MINVAL(temperature(:))>1.0d0).AND.&
      &(MAXVAL(corrected_temperature(:))<999.0d0).AND.(MINVAL(corrected_temperature(:))>1.0d0)) THEN
       WRITE(*,51) "",SUM(temperature(:))/3.0d0,SUM(corrected_temperature(:))/3.0d0
       WRITE(*,51) "x",temperature(1),corrected_temperature(1)
       WRITE(*,51) "y",temperature(2),corrected_temperature(2)
       WRITE(*,51) "z",temperature(3),corrected_temperature(3)
      ELSE
       WRITE(*,50) "x",temperature(1),corrected_temperature(1)
       WRITE(*,50) "y",temperature(2),corrected_temperature(2)
       WRITE(*,50) "z",temperature(3),corrected_temperature(3)
      ENDIF
      IF (VERBOSE_OUTPUT) WRITE(*,'("   Total drift is ",EN9.1)')SQRT(SUM(drift(:)**2))
      IF (constraints_available()) THEN
       IF (constraints_correction<=0.0d0) THEN
        CALL report_error(74)
       ELSE
        WRITE(*,'(" To correct for constraints, multiply values by ",F6.4)') constraints_correction
        WRITE(*,*) "(Temperatures given above contain no constraints correction)"
       ENDIF
      ENDIF
     ENDDO
     total_temp=2.0d7*kinetic_energy_total/(boltzmann*avogadro*give_total_degrees_of_freedom())
    ELSE
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Iterating over all molecule types."
     DO counter=1,give_number_of_molecule_types(),1
      CALL write_temp_for_one_molecule_type(counter)
     ENDDO
     total_temp=2.0d7*kinetic_energy_total/(boltzmann*avogadro*give_total_degrees_of_freedom()*FLOAT(endstep-startstep+1))
    ENDIF
    IF ((total_temp<999.0d0).AND.(total_temp>1.0d0)) THEN
     WRITE(*,ADVANCE="NO",FMT='(" Total average Temperature is ",F5.1," K")') total_temp
    ELSE
     WRITE(*,ADVANCE="NO",FMT='(" Total average Temperature is ",EN9.1," K")') total_temp
    ENDIF
    IF ((give_number_of_drude_particles()/=0).OR.(constraints_available())) THEN
     WRITE(*,'(", including drudes and constraints.")')
    ELSE
     WRITE(*,'(".")')
    ENDIF
    IF (DEVELOPERS_VERSION) THEN
     WRITE(*,*) " ! KINETIC ENERGY IS ",kinetic_energy_total
     total_temp=give_total_temperature(1)
     WRITE(*,*) " ! TEMPERATURE WITHOUT CONTRAINTS / DRUDES IS ",total_temp
    ENDIF
   ELSE
    !Just one molecule type.
    CALL write_temp_for_one_molecule_type(molecule_type_index_in)
   ENDIF
   CONTAINS

    SUBROUTINE write_temp_for_one_molecule_type(molecule_type_index)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index
    INTEGER :: timestep
    REAL :: drift_av(3),temperature_av(3),corrected_temperature_av(3)
    CHARACTER(LEN=1024) :: fstring
    LOGICAL :: connected
     IF (VERBOSE_OUTPUT) WRITE(*,'(" Molecule Type ",I0,":")') molecule_type_index
     !initialise variables
     drift_av(:)=0.0d0
     temperature_av(:)=0.0d0
     corrected_temperature_av(:)=0.0d0
     WRITE(fstring,'(2A,I0)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"temp_type_",molecule_type_index
     IF (VERBOSE_OUTPUT) WRITE(*,FMT='(3A)',ADVANCE="NO") "   Writing file '",TRIM(fstring),"'..."
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     OPEN(UNIT=3,FILE=TRIM(fstring))
     WRITE(3,*) "This file contains instantaneous direction-resolved (drift-corrected) kinetic temperatures."
     WRITE(3,'(8A14)') "timeline","Tx","Ty","Tz","Tx_corr","Ty_corr","Tz_corr","T_isotropic"
     DO timestep=startstep,endstep,1
      CALL give_temperature(timestep,drift,molecule_type_index,temperature,corrected_temperature,&
      &kinetic_energy,constraints_correction)
      kinetic_energy_total=kinetic_energy_total+kinetic_energy
      IF ((MAXVAL(temperature(:))<999.0d0).AND.(MINVAL(temperature(:))>1.0d0)) THEN
       WRITE(3,'(I14,7F14.6)') timestep*TIME_SCALING_FACTOR,temperature(:),corrected_temperature(:),SUM(temperature(:))/3.0d0
      ELSE
       WRITE(3,'(I14,7E14.6)') timestep*TIME_SCALING_FACTOR,temperature(:),corrected_temperature(:),SUM(temperature(:))/3.0d0
      ENDIF
      drift_av(:)=drift_av(:)+drift(:)
      temperature_av(:)=temperature_av(:)+temperature(:)
      corrected_temperature_av(:)=corrected_temperature_av(:)+corrected_temperature(:)
     ENDDO
     CLOSE(UNIT=3)
     IF (VERBOSE_OUTPUT) WRITE(*,'("done. Statistics:")')
     drift_av(:)=drift_av(:)/FLOAT(endstep-startstep+1)
     temperature_av(:)=temperature_av(:)/FLOAT(endstep-startstep+1)
     corrected_temperature_av(:)=corrected_temperature_av(:)/FLOAT(endstep-startstep+1)
     WRITE(*,'("   Average Temperatures:")')
     CALL readable_temperature_output(temperature_av(:))
     WRITE(*,'("   Average drift-corrected Temperatures:")')
     CALL readable_temperature_output(corrected_temperature_av(:))
     IF (VERBOSE_OUTPUT) THEN
      WRITE(*,'("   Average drift velocity (centre of mass):")')
      WRITE(*,'("     v:  ",EN10.1)') SQRT(SUM(drift_av(:)**2))
      WRITE(*,'("     vx: ",EN10.1)') drift_av(1)
      WRITE(*,'("     vy: ",EN10.1)') drift_av(2)
      WRITE(*,'("     vz: ",EN10.1)') drift_av(3)
     ENDIF
     IF (constraints_available()) THEN
      IF (constraints_correction<=0.0d0) THEN
       CALL report_error(74)
      ELSE
       WRITE(*,'(" To correct for constraints, multiply values by ",F6.4)') constraints_correction
       WRITE(*,*) "(reported Temperatures contain no constraints correction)"
      ENDIF
     ENDIF
    END SUBROUTINE write_temp_for_one_molecule_type

    SUBROUTINE readable_temperature_output(input_average)
    IMPLICIT NONE
    52 FORMAT ("     T",A,": ",EN9.1," K")
    53 FORMAT ("     T",A,": ",F5.1," K")
    REAL,INTENT(IN) :: input_average(3)
     !This routine is just for nice output.
     IF ((MAXVAL(input_average(:))<999.0d0).AND.(MINVAL(input_average(:))>1.0d0)) THEN
      WRITE(*,53) "",SUM(input_average(:))/3.0d0
      WRITE(*,53) "x",input_average(1)
      WRITE(*,53) "y",input_average(2)
      WRITE(*,53) "z",input_average(3)
     ELSE
      WRITE(*,52) "",SUM(input_average(:))/3.0d0
      WRITE(*,52) "x",input_average(1)
      WRITE(*,52) "y",input_average(2)
      WRITE(*,52) "z",input_average(3)
     ENDIF
    END SUBROUTINE readable_temperature_output

  END SUBROUTINE report_temperature

  SUBROUTINE test_dihedrals
  IMPLICIT NONE
  INTEGER :: dihedral_member_indices(2,4)
  REAL(KIND=GENERAL_PRECISION) :: dihedral_list(2)
  !INITIALISE dihedral_member_indices
   dihedral_member_indices(1,:)=(/1,14,9,15/)
   dihedral_member_indices(2,:)=(/14,9,15,2/)
   CALL initialise_dihedrals(dihedral_member_indices,1,2)
   CALL give_dihedrals(dihedral_list,1,1,dump_xyz=.FALSE.)
   WRITE(*,*) dihedral_list(:)
  END SUBROUTINE test_dihedrals

END MODULE DEBUG
!--------------------------------------------------------------------------------------------------------------------------------!
!This Module is capable of calculating autocorrelation functions. Currently implemented:
!the intermittent autocorrelation function for a binary operator, based on an arbitrary set of dihedral constraints.
!these constraints can be only one (e.g. for simple chain conformer analyses), two (e.g. two dihedrals plus folding for cisoid/transoid transitions), or more.
!reorientational autocorrelation functions.
!relative mean molecular velocity correlation functions.
MODULE AUTOCORRELATION ! Copyright (C) 2020 Frederik Philippi
    USE SETTINGS
 USE MOLECULAR
 IMPLICIT NONE
 PRIVATE
 !default values
 LOGICAL,PARAMETER :: fold_default=.FALSE.
 LOGICAL,PARAMETER :: dump_verbose_default=.FALSE.
 LOGICAL,PARAMETER :: skip_autocorr_default=.FALSE.
 LOGICAL,PARAMETER :: export_dihedral_default=.FALSE.
 INTEGER,PARAMETER :: sampling_interval_default=1
 INTEGER,PARAMETER :: legendre_order_default=2
 INTEGER(KIND=GENERAL_PRECISION),PARAMETER :: tmax_default=1000
 INTEGER(KIND=GENERAL_PRECISION),PARAMETER :: bin_count_default=100
 !Variables.
 TYPE,PRIVATE :: vc_component
        INTEGER :: reference_type=-1 ! reference molecule_type_index for velocity correlation
  INTEGER :: observed_type=-1 ! observed molecule_type_index
  LOGICAL :: self=.FALSE. ! if true, then this component is a self-correlation - the types must be identical
    END TYPE vc_component
 TYPE(vc_component),ALLOCATABLE :: vc_components(:) ! list of velocity correlation components to be calculated.
 INTEGER :: legendre_order !the order of the legendre polynomial to use, usually 2, maybe 1.
 LOGICAL,ALLOCATABLE :: autocorr_array(:,:)!first dimension: number of timesteps. second dimension: number of molecules per step.
 LOGICAL :: fold=fold_default !when true, then on top of the values a,b,... specified in the dihedral_list, (360-b),(360-a)... will be considered, too.
 LOGICAL :: dump_verbose=dump_verbose_default!controls if additional information is dumped into separate files for not.
 LOGICAL :: export_dihedral=export_dihedral_default!Controls whether dihedrals are to be exported.
 LOGICAL :: skip_autocorr=skip_autocorr_default!skips the actual autocorrelation. The preparation is still done, which is useful when only PES subsets or dihedral shares are required.
 INTEGER :: molecule_type_index_b!molecule_type_indices for the second molecule, 'b'
 REAL(KIND=GENERAL_PRECISION),ALLOCATABLE :: boundaries(:,:)!first dimension: index of dihedral. second dimension: upper and lower boundary.
 INTEGER,ALLOCATABLE :: export_list(:) !the list of molecules to export.
 INTEGER :: export_total_number !how many molecules are to be exported?
 INTEGER(KIND=GENERAL_PRECISION),ALLOCATABLE :: PES_subset_independent(:,:)!first dimension: number of condition. second dimension: binned PES subset
 INTEGER(KIND=GENERAL_PRECISION),ALLOCATABLE :: PES_subset_dependent(:,:)!first dimension: condition 1, second dimension: condition 2. Only for two dimensional subsets, i.e. for number_of_dihedral_conditions=2.
 REAL(KIND=WORKING_PRECISION) :: average_h!average of the population operator <h>
 INTEGER(KIND=GENERAL_PRECISION) :: tmax=tmax_default!max number of timesteps into the future for the autocorrelation function. Default is 1000 (maximaler shift)
 INTEGER :: sampling_interval=sampling_interval_default!every so many steps will be sampled
 INTEGER(KIND=WORKING_PRECISION) :: global_incidence,number_of_entries_in_array
 CHARACTER (LEN=32) :: operation_mode="dihedral"!operation mode of the autocorrelation module.
 CHARACTER (LEN=32),ALLOCATABLE :: formatted_dihedral_names(:)
 INTEGER :: molecule_type_index,number_of_dihedral_conditions,bin_count=bin_count_default
 !PRIVATE/PUBLIC declarations
 PRIVATE :: autocorr_array,molecule_type_index,operation_mode,number_of_dihedral_conditions,boundaries,formatted_dihedral_names
 PRIVATE :: initialise_autocorrelation,dihedral_autocorrelation,finalise_autocorrelation,bin_count,skip_autocorr
 PRIVATE :: global_incidence,dump_verbose,number_of_entries_in_array,PES_subset_independent,PES_subset_dependent,average_h
 PRIVATE :: tmax,calculate_autocorrelation_function_from_binary_array,molecule_type_index_b,export_dihedral,export_list
 PRIVATE :: fold_default,dump_verbose_default,skip_autocorr_default,tmax_default,bin_count_default,legendre_order
 PUBLIC :: perform_autocorrelation,user_dihedral_input,user_vacf_input,user_reorientation_input,write_simple_conductivity
 PUBLIC :: user_vcf_components_input,user_conductivity_input

 CONTAINS

  !WRITING input file to unit 8, which shouldn't be open.
  !'header' is either vcf or cacf.
  !has to be compliant with 'read_vc_components' and 'read_velocity_correlation_body' in 'AUTOCORRELATION' module
  SUBROUTINE user_vcf_components_input&
  &(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_vc_components,header)
  IMPLICIT NONE
  CHARACTER (LEN=*) :: filename_vc_components,header
  LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
  INTEGER,INTENT(IN) :: number_of_molecules,nsteps
  INTEGER :: maxmol,allocstatus,deallocstatus,ios,ncomponents,n
  LOGICAL :: connected
   SELECT CASE (TRIM(header))
   CASE ("cacf")
    PRINT *,"Generating input for electric current autocorrelation function."
    PRINT *,"You now have to define the CACF components."
   CASE ("vcf")
    PRINT *,"Generating input for velocity correlation functions."
    PRINT *,"You now have to define the VCF components."
   END SELECT
   maxmol=number_of_molecules
   IF (number_of_molecules==-1) maxmol=10000!unknown molecule number... expect the worst.
   PRINT *,"How many of these custom components would you like to define?"
   ncomponents=user_input_integer(1,10000)
   ALLOCATE(vc_components(ncomponents),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(11,exit_status=allocstatus)
   PRINT *,"For each of these components, you need to specify:"
   PRINT *," - the molecule type indices of the species to be correlated"
   PRINT *," - whether you want the self or distinct contributions."
   DO n = 1,ncomponents,1
    WRITE(*,'(" Reading the information for custom component number ",I0,":")') n
    PRINT *,"Please enter the first molecule type index:"
    vc_components(n)%reference_type=user_input_integer(1,maxmol)
    PRINT *,"Please enter the second molecule type index:"
    vc_components(n)%observed_type=user_input_integer(1,maxmol)
    IF (vc_components(n)%reference_type==vc_components(n)%observed_type) THEN
     PRINT *,"Would you like to compute self (y) or distinct (n) contributions?"
     vc_components(n)%self=user_input_logical()
    ELSE
     vc_components(n)%self=.FALSE.
    ENDIF
   ENDDO
   PRINT *,"How many steps do you want the shift of the (auto)correlation functions to be?"
   PRINT *,"A good starting value is usually 2500 (with 2fs steps / time_scaling 2)."
   WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
   tmax=user_input_integer(1,(nsteps-1))
   !tmax is initialised now.
   PRINT *,"Every how many steps would you like to use?"
   WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
   sampling_interval=user_input_integer(1,nsteps)
   parallelisation_possible=.TRUE.
   IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
    PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
    IF (user_input_logical()) parallelisation_requested=.TRUE.
   ENDIF
   WRITE(*,FMT='(A)',ADVANCE="NO") " writing "//TRIM(header)//" input file..."
   INQUIRE(UNIT=8,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=8)
   OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_vc_components),IOSTAT=ios)!input path is added for the VCF file!
   IF (ios/=0) CALL report_error(46,exit_status=ios)
   WRITE(8,'(A,I0,A)') " "//TRIM(header)//" ",ncomponents," ### type of analysis + number of components"
   DO n = 1,ncomponents,1
    WRITE(8,'(" ",I0," ",I0," ",L1)') vc_components(n)%reference_type,vc_components(n)%observed_type,vc_components(n)%self
   ENDDO
   WRITE(8,'(" tmax ",I0," ### maximum time shift of the correlation function")') tmax
   WRITE(8,'(" sampling_interval ",I0," ### every so many timesteps will be used")') sampling_interval
   WRITE(8,*) "quit"
   WRITE(8,*)
   WRITE(8,*) "This is an input file for the calculation of "//TRIM(header)//" components."
   WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
   ENDFILE 8
   CLOSE(UNIT=8)
   DEALLOCATE(vc_components,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   WRITE(*,*) "done"
  END SUBROUTINE user_vcf_components_input

  !WRITING input file to unit 8, which shouldn't be open.
  !has to be compliant with 'read_velocity_correlation_body' in 'AUTOCORRELATION' module
  SUBROUTINE user_conductivity_input&
  &(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_conductivity)
  IMPLICIT NONE
  CHARACTER (LEN=*) :: filename_conductivity
  LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
  INTEGER,INTENT(IN) :: number_of_molecules,nsteps
  INTEGER :: maxmol,allocstatus,deallocstatus,ios,ncomponents,n
  LOGICAL :: connected
   PRINT *,"Generating input for electrical conductivity."
   PRINT *,"How many steps do you want the shift of the (auto)correlation functions to be?"
   PRINT *,"A good starting value is usually 10000 (with 20fs steps / time_scaling 20)."
   WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
   tmax=user_input_integer(1,(nsteps-1))
   !tmax is initialised now.
   PRINT *,"Every how many steps would you like to use?"
   WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
   sampling_interval=user_input_integer(1,nsteps)
   parallelisation_possible=.TRUE.
   IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
    PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
    IF (user_input_logical()) parallelisation_requested=.TRUE.
   ENDIF
   WRITE(*,FMT='(A)',ADVANCE="NO") " writing conductivity input file..."
   INQUIRE(UNIT=8,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=8)
   OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_conductivity),IOSTAT=ios)!input path is added for the VCF file!
   IF (ios/=0) CALL report_error(46,exit_status=ios)
   WRITE(8,'(A)') " conductivity ### type of analysis"
   WRITE(8,'(" tmax ",I0," ### maximum time shift of the correlation function")') tmax
   WRITE(8,'(" sampling_interval ",I0," ### every so many timesteps will be used")') sampling_interval
   WRITE(8,*) "quit"
   WRITE(8,*)
   WRITE(8,*) "This is an input file for the calculation of electrical conductivity."
   WRITE(8,*) "To actually perform the implied calculation, it has to be referenced in 'general.inp'."
   ENDFILE 8
   CLOSE(UNIT=8)
   WRITE(*,*) "done"
  END SUBROUTINE user_conductivity_input

  !WRITING input file to unit 8, which shouldn't be open.
  !has to be compliant with 'read_input_for_rmmvcf' in 'AUTOCORRELATION' module
  SUBROUTINE user_vacf_input&
  &(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_rmmvcf)
  IMPLICIT NONE
  CHARACTER (LEN=*) :: filename_rmmvcf
  LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
  INTEGER,INTENT(IN) :: number_of_molecules,nsteps
  INTEGER :: maxmol,ios
  LOGICAL :: connected
   PRINT *,"Generating VACF input. (in particular, RMM-VCF)"
   !the case 'number_of_molecules==1' has already been caught earlier, in the calling routine.
   IF (number_of_molecules==2) THEN
    PRINT *,"Only two molecule types are present, which will be used as input."
    molecule_type_index=1
    molecule_type_index_b=2
   ELSE
    PRINT *,"You need to specify two molecule types now."
    maxmol=number_of_molecules
    IF (number_of_molecules==-1) maxmol=10000!unknown molecule number... expect the worst.
    PRINT *,"Please enter the index of the first molecule as integer:"
    molecule_type_index=user_input_integer(1,maxmol)
    PRINT *,"Please enter the index of the second molecule as integer:"
    DO
     molecule_type_index_b=user_input_integer(1,maxmol)
     IF (molecule_type_index/=molecule_type_index_b) THEN
      EXIT!valid input.
     ELSE
      PRINT *,"Please give two *different* molecule types."
     ENDIF
    ENDDO
   ENDIF
   !molecule_type_index and molecule_type_index_b are initialised now.
   PRINT *,"How many steps do you want the shift of the (auto)correlation functions to be?"
   PRINT *,"A good starting value is usually 2500 (with 2fs steps / time_scaling 2)."
   WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
   tmax=user_input_integer(1,(nsteps-1))
   !tmax is initialised now (molecule_type_index and molecule_type_index_b too).
   PRINT *,"Should the self-contributions also be calculated? (y/n)"
   IF (user_input_logical()) THEN
    PRINT *,"Every how many steps would you like to use for the self-contributions?"
    WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
    sampling_interval=user_input_integer(1,nsteps)
    !include self-contributions. This means that parallelisation is possible... (and advisable!)
    skip_autocorr=.FALSE.
    parallelisation_possible=.TRUE.
    IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
     PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
     IF (user_input_logical()) parallelisation_requested=.TRUE.
    ENDIF
   ELSE
    skip_autocorr=.TRUE.
   ENDIF
   WRITE(*,FMT='(A30)',ADVANCE="NO") " writing RMM-VCF input file..."
   INQUIRE(UNIT=8,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=8)
   OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_rmmvcf),IOSTAT=ios)!input path is added for the RMMVCF file!
   IF (ios/=0) CALL report_error(46,exit_status=ios)
   WRITE(8,'(" ",I0," ",I0," ### Molecule type indices (i.e. the molecules to observe)")')&
   &molecule_type_index,molecule_type_index_b
   WRITE(8,*) "rmm-vcf ### the type of analysis"
   WRITE(8,'(" tmax ",I0," ### maximum time shift of the correlation function")') tmax
   WRITE(8,'(" sampling_interval ",I0," ### every so many timesteps will be used for the self-contributions")') sampling_interval
   WRITE(8,FMT='(" skip_autocorrelation ",L1)',ADVANCE="NO") skip_autocorr
   IF (skip_autocorr) THEN
    WRITE(8,*) " ### don't calculate self-contributions"
   ELSE
    WRITE(8,*) " ### calculate self-contributions, print related information"
   ENDIF
   WRITE(8,*) "quit"
   WRITE(8,*)
   WRITE(8,*) "This is an input file for the calculation of velocity correlation coefficients."
   WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
   ENDFILE 8
   CLOSE(UNIT=8)
   WRITE(*,*) "done"
  END SUBROUTINE user_vacf_input

  !WRITING input file to unit 8, which shouldn't be open.
  !has to be compliant with 'read_input_for_dihedral_mode'
  SUBROUTINE user_dihedral_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_dihedral)
  IMPLICIT NONE
  CHARACTER (LEN=*) :: filename_dihedral
  LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
  INTEGER,INTENT(IN) :: number_of_molecules,nsteps
  INTEGER :: maxmol,n,allocstatus,deallocstatus,ios,inputinteger
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: dihedral_member_indices !list of atom indices used to generate input
  LOGICAL :: connected
   PRINT *,"Generating dihedral condition input."
   PRINT *,"Please enter the number of the molecule type you would like to observe."
   maxmol=number_of_molecules
   IF (number_of_molecules<1) maxmol=10000!unknown molecule number... expect the worst.
   molecule_type_index=user_input_integer(1,maxmol)
   PRINT *,"You now have to define the set of dihedral conditions to be fulfilled simultaneously."
   PRINT *,"How many conditions would you like to define?"
   number_of_dihedral_conditions=user_input_integer(1,10000)
   !allocate the temporary memory.
   ALLOCATE(dihedral_member_indices(number_of_dihedral_conditions,4),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(11,exit_status=allocstatus)
   ALLOCATE(boundaries(number_of_dihedral_conditions,2),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   PRINT *,"For each of these conditions, a dihedral is required,"
   PRINT *,"as well as the lower and upper bounds this dihedral is allowed"
   PRINT *,"to take for the condition to be fulfilled."
   DO n = 1,number_of_dihedral_conditions,1 !reading in the specified dihedral conditions that have to be fulfilled simultaneously
    WRITE(*,'(" Reading the information for the dihedral condition number ",I0,":")') n
    PRINT *,"Please enter the index of the first atom:"
    dihedral_member_indices(n,1)=user_input_integer(1,10000)
    PRINT *,"Please enter the index of the second atom:"
    dihedral_member_indices(n,2)=user_input_integer(1,10000)
    PRINT *,"Please enter the index of the third atom:"
    dihedral_member_indices(n,3)=user_input_integer(1,10000)
    PRINT *,"Please enter the index of the fourth atom:"
    dihedral_member_indices(n,4)=user_input_integer(1,10000)
    WRITE(*,'(" You specified the dihedral ",I0,"-",I0,"-",I0,"-",I0,".")') dihedral_member_indices(n,:)
    PRINT *,"What would you like the lower boundary to be?"
    PRINT *,"(Bear in mind that the dihedrals are defined from 0.0° to 360.0°)"
    boundaries(n,1)=user_input_real(0.0,360.0)
    PRINT *,"Please enter the value for the upper boundary."
    boundaries(n,2)=user_input_real(SNGL(boundaries(n,1)),360.0)
   ENDDO
   !Boundaries and dihedral members are initialised now. (Also: molecule_type_index, number_of_dihedral_conditions)
   PRINT *,"It is possible to consider the 'folded' values."
   PRINT *,"'(360-upper) to (360-lower)' then also fulfills the condition, not only 'lower to upper'."
   PRINT *,"Would you like to use these extended boundaries? (y/n)"
   fold=user_input_logical()
   PRINT *,"Please enter the bin count, e.g. '36' equals to binning in steps of 10°."
   WRITE(*,'(" The default is currently set to ",I0,".")') bin_count_default
   bin_count=user_input_integer(10,360)
   PRINT *,"Do you want to compute static properties? (y/n)"
   IF (number_of_dihedral_conditions==2) THEN
    PRINT *,"(Independent and dependent incidences, share of fulfilled conditions)"
   ELSE
    PRINT *,"(Independent incidences ('counts') and share of fulfilled conditions)"
   ENDIF
   dump_verbose=user_input_logical()
   PRINT *,"Would you like to print the values for the dihedrals of a particular molecule? (y/n)"
   PRINT *,"(Written into a file containing a column for the timestep and one for each dihedral)"
   export_dihedral=user_input_logical()
   IF (export_dihedral) THEN
    IF (number_of_molecules==1) THEN
     PRINT *,"Just one molecule, which will be exported."
     export_total_number=1
    ELSE
     PRINT *,"Of how many molecules would you like to export the dihedral values?"
     export_total_number=user_input_integer(1,maxmol)
    ENDIF
    !I mean, I have that variable. Might as well use it.
    ALLOCATE(export_list(export_total_number),STAT=allocstatus)
    IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
    WRITE(*,'(A,I0,A)') "You have to enter the indices of the molecules (of type ",molecule_type_index,") now."
    DO n=1,export_total_number,1
     WRITE(*,'(A,I0,A,I0,A)') "Please enter molecule index ",n," out of ",export_total_number," you would like to export."
     inputinteger=user_input_integer(1,maxmol)
     !check if sensible, i.e. no double specifications. Otherwise, complain.
     IF (n>1) THEN
      IF (ANY(export_list(1:n-1)==inputinteger)) THEN
       WRITE(*,'(A,I0,A)') "The molecule with index ",n," has already been specified."
       CYCLE
      ENDIF
     ENDIF
     export_list(n)=inputinteger
    ENDDO
   ENDIF
   skip_autocorr=.FALSE.
   IF (dump_verbose) THEN
    PRINT *,"Do you want to skip the autocorrelation and just compute the static properties? (y/n)"
    skip_autocorr=user_input_logical()
   ENDIF
   IF (.NOT.(skip_autocorr)) THEN
    parallelisation_possible=.TRUE.
    PRINT *,"How many steps do you want the shift of the autocorrelation functions to be?"
    PRINT *,"A good starting value is usually 2000 (with 2fs steps / time_scaling 2)."
    WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
    tmax=user_input_integer(1,(nsteps-1))
    IF (.NOT.(parallelisation_requested)) THEN!parallelisation no requested?
     PRINT *,"Parallelisation is available for the autocorrelation function. Would you like to turn it on? (y/n)"
     IF (user_input_logical()) parallelisation_requested=.TRUE.
    ENDIF
   ENDIF
   !sufficient information collected.
   WRITE(*,FMT='(A30)',ADVANCE="NO") " writing dihedral input file..."
   INQUIRE(UNIT=8,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=8)
   OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_dihedral),IOSTAT=ios)!input path is added for the dihedral file!
   IF (ios/=0) CALL report_error(46,exit_status=ios)
   WRITE(8,'(" ",I0," ### Molecule type index (i.e. the molecule to observe)")') molecule_type_index
   WRITE(8,'(" dihedral ",I0," ### dihedral autocorrelation + number of dihedral conditions to be fulfilled ")')&
   & number_of_dihedral_conditions
   DO n = 1,number_of_dihedral_conditions,1
    WRITE(8,'(" ",I0," ",I0," ",I0," ",I0," ",2F6.1)') dihedral_member_indices(n,:),boundaries(n,:)
   ENDDO
   WRITE(8,'(" tmax ",I0," ### maximum time shift of the correlation function")') tmax
   WRITE(8,FMT='(" skip_autocorrelation ",L1)',ADVANCE="NO") skip_autocorr
   IF (skip_autocorr) THEN
    WRITE(8,*) " ### no autocorrelation, just the PES & indicences."
   ELSE
    WRITE(8,*) " ### compute the autocorrelation function."
   ENDIF
   IF (fold) THEN
    WRITE(8,'(" fold ",L1," ### also check for the range (360-b) to (360-a), not just a to b")') fold
   ELSE
    WRITE(8,'(" fold ",L1," ### just use the range from a to b, no folding")') fold
   ENDIF
   IF (export_dihedral) THEN
    DO n=1,export_total_number,1
     WRITE(8,'(" export ",I0," ### write the dihedrals for molecule index ",I0," into separate file.")')&
     & export_list(n),export_list(n)
    ENDDO
    !DEALLOCATE memory again
    DEALLOCATE(export_list,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
    !reset variables to avoid interference
    export_dihedral=.FALSE.
    export_total_number=0
   ENDIF
   IF (dump_verbose) THEN
    WRITE(8,'(" dump_verbose ",L1," ### dump verbose information such as the PES subset population.")') dump_verbose
   ELSE
    WRITE(8,'(" dump_verbose ",L1," ### do not dump verbose information.")') dump_verbose
   ENDIF
   WRITE(8,'(" bin_count ",I0," ### Setting the bin count to ",I0," (default is ",I0,")")')&
   & bin_count,bin_count,bin_count_default
   WRITE(8,*) "quit"
   WRITE(8,*)
   WRITE(8,*) "This is an input file for the dihedral conditions analysis."
   WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
   ENDFILE 8
   CLOSE(UNIT=8)
   DEALLOCATE(boundaries,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   DEALLOCATE(dihedral_member_indices,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   WRITE(*,*) "done"
  END SUBROUTINE user_dihedral_input

  !WRITING input file to unit 8, which shouldn't be open.
  !has to be compliant with 'read_input_for_rmmvcf' in 'AUTOCORRELATION' module
  SUBROUTINE user_reorientation_input&
  &(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_reorient)
  IMPLICIT NONE
  CHARACTER (LEN=*) :: filename_reorient
  LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
  INTEGER,INTENT(IN) :: number_of_molecules,nsteps
  INTEGER :: maxmol,ios,allocstatus,deallocstatus,n,number_of_tip_atoms,number_of_base_atoms
  LOGICAL :: connected
  INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_base(:) !temporary list of centre-of-mass fragments (defined as atom_indices) for base atom
  INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_tip(:) !temporary list of centre-of-mass fragments (defined as atom_indices) for tip atom
   PRINT *,"Generating input for reorientational time correlation function."
   maxmol=number_of_molecules
   IF (number_of_molecules<1) maxmol=10000!unknown molecule number... expect the worst.
   PRINT *,"Please specify the molecule type index of the molecule to observe as an integer:"
   molecule_type_index=user_input_integer(1,maxmol)
   PRINT *,"How many steps do you want the shift of the (auto)correlation functions to be?"
   WRITE(*,'(" The default is currently set to ",I0,".")') tmax_default
   tmax=user_input_integer(1,(nsteps-1))
   PRINT *,"Please enter the order of the Legendre polynomial you wish to use as an integer:"
   WRITE(*,'(A,I0,A)') " If in doubt, enter '2'. The default is currently set to '",legendre_order_default,"'."
   legendre_order=user_input_integer(0,4)
   !tmax, molecule_type_index, legendre_order are initialised now.
   PRINT *,"You have to specify the vector whose reorientation dynamics are to be computed."
   PRINT *,"For this vector, you have to define a base fragment and a tip fragment."
   PRINT *,"The centres of mass of these fragments define the base and tip points of the vector."
   PRINT *,"It is possible to just specify two atoms (e.g. for a N-H vector)."
   PRINT *,"How many atoms are in the base fragment?"
   number_of_base_atoms=user_input_integer(1,100)
   ALLOCATE(fragment_list_base(number_of_base_atoms),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(79,exit_status=allocstatus)
   PRINT *,"You have to enter the atom indices of the (base) fragment atoms now."
   DO n=1,number_of_base_atoms,1
    WRITE(*,'(A,I0,A,I0,A)') " Please enter the atom index for fragment atom ",n," out of ",number_of_base_atoms,"."
    fragment_list_base(n)=user_input_integer(1,10000)
    !check for double specifications, and print warning
    IF (n>1) THEN
     IF (ANY(fragment_list_base(1:n-1)==fragment_list_base(n))) &
     &PRINT *,"This atom has already been specified. (Which is allowed, but not necessarily sensible.)"
    ENDIF
   ENDDO
   PRINT *,"How many atoms are in the tip fragment?"
   number_of_tip_atoms=user_input_integer(1,100)
   ALLOCATE(fragment_list_tip(number_of_tip_atoms),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(79,exit_status=allocstatus)
   PRINT *,"You have to enter the atom indices of the (tip) fragment atoms now."
   DO n=1,number_of_tip_atoms,1
    WRITE(*,'(A,I0,A,I0,A)') " Please enter the atom index for fragment atom ",n," out of ",number_of_tip_atoms,"."
    fragment_list_tip(n)=user_input_integer(1,10000)
    !check for double specifications, and print warning
    IF (n>1) THEN
     IF (ANY(fragment_list_tip(1:n-1)==fragment_list_tip(n))) &
     &PRINT *,"This atom has already been specified. (Which is allowed, but not necessarily sensible.)"
    ENDIF
   ENDDO
   PRINT *,"Every how many steps would you like to use for the time correlation function?"
   WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
   sampling_interval=user_input_integer(1,nsteps)
   !parallelisation is possible, because autocorrelation.
   parallelisation_possible=.TRUE.
   IF (.NOT.(parallelisation_requested)) THEN! ask for parallelisation, if not yet requested.
    PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
    IF (user_input_logical()) parallelisation_requested=.TRUE.
   ENDIF
   WRITE(*,FMT='(A30)',ADVANCE="NO") " writing input file for reorientational time correlation function..."
   INQUIRE(UNIT=8,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=8)
   OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_reorient),IOSTAT=ios)!input path is added for the reorientational tcf file!
   IF (ios/=0) CALL report_error(46,exit_status=ios)
   WRITE(8,'(" ",I0," ### Molecule type index (i.e. the molecule to observe)")')&
   &molecule_type_index
   WRITE(8,*) "reorientation ### the type of analysis"
   WRITE(8,'(" tmax ",I0," ### maximum time shift of the time correlation function")') tmax
   WRITE(8,'(" sampling_interval ",I0," ### every so many timesteps will be used for the tcf")') sampling_interval
   WRITE(8,'(" legendre ",I0," ### use legendre polynomial of order ",I0)') legendre_order,legendre_order
   WRITE(8,FMT='(" base ",I0)',ADVANCE="NO") number_of_base_atoms
   IF (number_of_base_atoms==1) THEN
    WRITE(8,'(" ### the atom with index ",I0," is used as base point")') fragment_list_base(1)
   ELSE
    WRITE(8,'(" ### ",I0," atoms are used to define the base point (as centre of mass)")') number_of_base_atoms
   ENDIF
   DO n=1,number_of_base_atoms,1
    WRITE(8,FMT='(" ",I0)',ADVANCE="NO") fragment_list_base(n)
   ENDDO
   WRITE(8,*)
   WRITE(8,FMT='(" tip ",I0)',ADVANCE="NO") number_of_tip_atoms
   IF (number_of_tip_atoms==1) THEN
    WRITE(8,'(" ### the atom with index ",I0," is used as tip point")') fragment_list_tip(1)
   ELSE
    WRITE(8,'(" ### ",I0," atoms are used to define the tip point (as centre of mass)")') number_of_tip_atoms
   ENDIF
   DO n=1,number_of_tip_atoms,1
    WRITE(8,FMT='(" ",I0)',ADVANCE="NO") fragment_list_tip(n)
   ENDDO
   WRITE(8,*)
   WRITE(8,*) "quit"
   WRITE(8,*)
   WRITE(8,*) "This is an input file for the calculation of a vector reorientational time correlation function."
   WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
   ENDFILE 8
   CLOSE(UNIT=8)
   WRITE(*,*) "done"
   DEALLOCATE(fragment_list_base,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
   DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
  END SUBROUTINE user_reorientation_input

  !initialises the autocorrelation module by reading the specified input file.
  SUBROUTINE initialise_autocorrelation()
  IMPLICIT NONE
  LOGICAL :: file_exists,connected
  INTEGER :: ios,allocstatus
   ! first, check if file exists.
   INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_AUTOCORRELATION_INPUT),EXIST=file_exists)
   IF (file_exists) THEN
    CALL set_defaults()
    IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(PATH_INPUT)//TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    OPEN(UNIT=3,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_AUTOCORRELATION_INPUT),&
    &ACTION='READ',IOSTAT=ios)
    IF (ios/=0) CALL report_error(14,exit_status=ios)
    READ(3,IOSTAT=ios,FMT=*) operation_mode
    IF (ios/=0) operation_mode="UNK"
    IF (TRIM(operation_mode)=="eccf") operation_mode="ecaf"!support for synonyms
    IF (TRIM(operation_mode)=="cacf") operation_mode="ecaf"!support for synonyms
    IF (TRIM(operation_mode)=="conductivity_components") operation_mode="ecaf"!support for synonyms
    IF ((TRIM(operation_mode)=="rmm-vcf")&
    &.OR.(TRIM(operation_mode)=="conductivity")&
    &.OR.(TRIM(operation_mode)=="vcf")&
    &.OR.(TRIM(operation_mode)=="ecaf")) THEN
     molecule_type_index=1
     molecule_type_index_b=2
    ELSE
     !The part in this level is kept for backwards compatibility.
     REWIND 3
     READ(3,IOSTAT=ios,FMT=*) molecule_type_index
     IF (ios/=0) THEN
      READ(3,IOSTAT=ios,FMT=*) operation_mode
      IF (ios/=0) CALL report_error(14,exit_status=ios)
      IF (TRIM(operation_mode)=="rmm-vcf") THEN
       BACKSPACE 3 !rmm-vcf can handle 'no' input for the indices, it then just takes the first two.
       molecule_type_index=1
       molecule_type_index_b=2
      ELSE
       CALL report_error(14,exit_status=ios)!ERROR 14: incorrect format in autocorrelation.inp
      ENDIF
     ENDIF
     IF ((molecule_type_index>give_number_of_molecule_types()).OR.(molecule_type_index<1)) THEN
      !the specified molecule type doesn't exist.
      CALL report_error(33,exit_status=molecule_type_index)
      CLOSE(UNIT=3)
      RETURN
     ENDIF
     READ(3,IOSTAT=ios,FMT=*) operation_mode!read the operation mode.
     IF (ios/=0) CALL report_error(14,exit_status=ios)
    ENDIF
    !Now read the body of the autocorrelation input file in line with the requested operation mode:
    SELECT CASE (TRIM(operation_mode))
    CASE ("dihedral")
     WRITE(*,*) "Performing autocorrelation analysis of dihedral subspace."
     CALL read_input_for_dihedral_mode()!uses unit 3!!
     IF (ERROR_CODE==121) RETURN
     number_of_entries_in_array=give_number_of_timesteps()*give_number_of_molecules_per_step(molecule_type_index)
     IF (dump_verbose) THEN!verbose output - array required to bin into. Has to be allocated now:
      ALLOCATE(PES_subset_independent(number_of_dihedral_conditions,0:bin_count),STAT=allocstatus)
      IF (SIZE(autocorr_array)/=number_of_entries_in_array) CALL report_error(0)
      IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
      PES_subset_independent(:,:)=0
      IF (number_of_dihedral_conditions==2) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,*) "2D PES subset - reporting dependent subset as well."
       ALLOCATE(PES_subset_dependent(0:bin_count,0:bin_count),STAT=allocstatus)
       IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
       PES_subset_dependent(:,:)=0
      ENDIF
     ENDIF
    CASE ("reorientation")
     WRITE(*,*) "Performing autocorrelation analysis of vector reorientation."
     WRITE(*,*) "See for example dx.doi.org/10.1016/B978-0-12-387032-2.00011-8"
     WRITE(*,*) "or, for ionic liquids, dx.doi.org/10.1016/j.cplett.2007.03.084"
     WRITE(*,'(" Molecule type index is ",I0)'),molecule_type_index
     CALL read_input_for_reorientation()
    CASE ("rmm-vcf")!correlation module!
     IF (give_number_of_molecule_types()==1) THEN
      CALL report_error(39)
      RETURN
     ELSE
      WRITE(*,*) "Performing rmm-vcf analysis. See dx.doi.org/10.1103/PhysRevE.50.1162"
      CALL report_error(40)!Warn the user about the unusual format of the input file (velocities).
      CALL read_input_for_rmmvcf()!uses unit 3!!
     ENDIF
    CASE ("conductivity")
     WRITE(*,*) "Calculating electrical conductivity from electric current autocorrelation function."
     WRITE(*,*) "See, for example, equation (10) in J. Chem. Phys., 2002, 116, 3018–3026."
     CALL report_error(40)!Warn the user about the unusual format of the input file (velocities).
     CALL read_velocity_correlation_body()
    CASE ("vcf")
     WRITE(*,*) "Calculating custom velocity correlation functions."
     WRITE(*,*) "See, for example, equations (A6) to (A10) in J. Phys. Chem. B, 2011, 115, 13212–13221."
     IF (VERBOSE_OUTPUT) THEN
      WRITE(*,*) "Normalisation follows the convention given in the above paper."
      WRITE(*,*) "Thus, what will be reported are the quantities enclosed in <...>."
      WRITE(*,*) "(i.e. including Nan*(Nan-1), but not N or 1/3 in (A9))"
     ENDIF
     CALL report_error(40)
     CALL read_vc_components()
     IF (ERROR_CODE==118) RETURN
     CALL read_velocity_correlation_body()
    CASE ("ecaf")
     WRITE(*,*) "Calculating custom components of the microscopic charge current autocorrelation function."
     WRITE(*,*) "See also: THEORY OF SIMPLE LIQUIDS (Hansen / McDonald), fourth edition, chapter 7.7 and 10.5."
     IF (VERBOSE_OUTPUT) THEN
      WRITE(*,*) "Reported quantities are not normalised (apart from the number of ensemble averages)."
      WRITE(*,*) "Thus, the SUM(zi*zj*<vi*vj>) will be reported (not including elementary charge)."
     ENDIF
     CALL report_error(40)
     CALL read_vc_components()
     IF (ERROR_CODE==118) RETURN
     CALL read_velocity_correlation_body()
    CASE DEFAULT
     CALL report_error(14)
    END SELECT
    CLOSE(UNIT=3)
    !Finally, set global_incidence to zero.
    global_incidence=0
   ELSE
    CALL report_error(21)!No input - no output. easy as that.
   ENDIF
   CONTAINS

    SUBROUTINE read_vc_components()
    IMPLICIT NONE
    INTEGER :: ncomponents,n,allocstatus,a,b,c,d
    CHARACTER(LEN=32) :: inputstring
     REWIND 3
     READ(3,IOSTAT=ios,FMT=*) inputstring,ncomponents
     IF (ios/=0) CALL report_error(14) !no compromises!
     IF (ncomponents<1) THEN
      CALL report_error(118,ncomponents)
      RETURN
     ENDIF
     IF (VERBOSE_OUTPUT) WRITE(*,'(" Reading ",I0," components:")') ncomponents
     !allocate memory
     IF (ALLOCATED(vc_components)) CALL report_error(0)
     ALLOCATE(vc_components(ncomponents),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(115,exit_status=allocstatus)
     DO n=1,ncomponents,1
      READ(3,IOSTAT=ios,FMT=*) a,b
      IF (ios/=0) CALL report_error(14)
      IF (a==b) THEN
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) c,d,vc_components(n)%self
       IF (ios/=0) CALL report_error(14)
      ELSE
       vc_components(n)%self=.FALSE.
       !not the same molecule type, must be distinct.
      ENDIF
      vc_components(n)%reference_type=a
      vc_components(n)%observed_type=b
      IF (VERBOSE_OUTPUT) THEN
       IF (vc_components(n)%self) THEN
        WRITE(*,FMT='("   ",I0," is <v",I0,"*v",I0,"> (self)")',ADVANCE="NO")n,a,b
       ELSE
        WRITE(*,FMT='("   ",I0," is <v",I0,"*v",I0,"> (distinct)")',ADVANCE="NO")n,a,b
       ENDIF
       IF ((a<1).OR.(a>give_number_of_molecule_types())&
       &.OR.(b<1).OR.(b>give_number_of_molecule_types())&
       &.OR.((a/=b).AND.(vc_components(n)%self))) THEN
        WRITE(*,'(" ! INVALID !")')
       ELSE
        WRITE(*,*)
       ENDIF
      ENDIF
     ENDDO
    END SUBROUTINE read_vc_components

    SUBROUTINE read_velocity_correlation_body()
    IMPLICIT NONE
    INTEGER :: n
    CHARACTER(LEN=32) :: inputstring
     !File should be positioned correctly already.
     !Here, the indices are ready and the body of the input file can be read.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
      IF (ios/=0) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      ENDIF
      SELECT CASE (TRIM(inputstring))
      CASE ("tmax")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tmax' to default (=",tmax_default,")"
        tmax=tmax_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
       ENDIF
      CASE ("sampling_interval")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') &
        &"setting 'sampling_interval' to default (=",sampling_interval_default,")"
        sampling_interval=sampling_interval_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'sampling_interval' to ",sampling_interval
       ENDIF
      CASE ("quit")
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      CASE DEFAULT
       IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
      END SELECT
     ENDDO
    END SUBROUTINE read_velocity_correlation_body

    SUBROUTINE read_input_for_velocity_correlations()
    IMPLICIT NONE
    INTEGER :: n,a,b
    CHARACTER(LEN=32) :: inputstring
     REWIND 3
     READ(3,IOSTAT=ios,FMT=*) a,b
     IF (ios/=0) THEN
      molecule_type_index=1!equal to molecule 'a'
      molecule_type_index_b=2!equal to molecule 'b'
     ELSE
      !check if the molecule indices are sensible
      IF ((((a>0).AND.(b>0)).AND.(a/=b)).AND.&
       &((a<=give_number_of_molecule_types()).AND.(b<=give_number_of_molecule_types()))) THEN
       !the indices are valid.
       molecule_type_index=a
       molecule_type_index_b=b
      ELSE
       IF (VERBOSE_OUTPUT) WRITE(*,'(" invalid input for the indices (read ",I0," and ",I0,")")') a,b
       molecule_type_index=1
       molecule_type_index_b=2
      ENDIF
     ENDIF
     WRITE(*,'(" Using molecule_type ",I0," (for a) and ",I0," (for b).")') &
     &molecule_type_index,molecule_type_index_b
     !skip over line with operation_mode
     READ(3,*)
     !Here, the indices are ready and the body of the input file can be read.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
      IF (ios/=0) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      ENDIF
      SELECT CASE (TRIM(inputstring))
      CASE ("tmax")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tmax' to default (=",tmax_default,")"
        tmax=tmax_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
       ENDIF
      CASE ("sampling_interval")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') &
        &"setting 'sampling_interval' to default (=",sampling_interval_default,")"
        sampling_interval=sampling_interval_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'sampling_interval' to ",sampling_interval
       ENDIF
      CASE ("skip_autocorrelation")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,skip_autocorr
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') "setting 'skip_autocorr' to default (=",skip_autocorr_default,")"
        skip_autocorr=skip_autocorr_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'skip_autocorr' to ",skip_autocorr
       ENDIF
      CASE ("quit")
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      CASE DEFAULT
       IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
      END SELECT
     ENDDO
    END SUBROUTINE read_input_for_velocity_correlations

    SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
    IMPLICIT NONE
     bin_count=bin_count_default
     tmax=tmax_default
     fold=fold_default
     dump_verbose=dump_verbose_default
     skip_autocorr=skip_autocorr_default
     sampling_interval=sampling_interval_default
     legendre_order=legendre_order_default
    END SUBROUTINE set_defaults

    SUBROUTINE read_input_for_reorientation()
    IMPLICIT NONE
    LOGICAL :: tip_read,base_read
    INTEGER :: number_of_tip_atoms,number_of_base_atoms,counter,allocstatus,deallocstatus,n
    INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_base(:) !temporary list of centre-of-mass fragments (defined as atom_indices) for base atom
    INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_tip(:) !temporary list of centre-of-mass fragments (defined as atom_indices) for tip atom
    CHARACTER(LEN=32) :: inputstring
     tip_read=.FALSE.
     base_read=.FALSE.
     legendre_order=legendre_order_default
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
      IF (ios/=0) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      ENDIF
      SELECT CASE (TRIM(inputstring))
      CASE ("tmax")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') " setting 'tmax' to default (=",tmax_default,")"
        tmax=tmax_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
       ENDIF
      CASE ("quit")
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      CASE ("legendre")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,legendre_order
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') " setting 'legendre_order' to default (=",legendre_order_default,")"
        legendre_order=legendre_order_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(" ",A,I0)') "requesting legendre polynomial of order ",legendre_order
       ENDIF
      CASE ("sampling_interval")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') &
        &"setting 'sampling_interval' to default (=",sampling_interval_default,")"
        sampling_interval=sampling_interval_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'sampling_interval' to ",sampling_interval
       ENDIF
      CASE ("tip")
       IF (tip_read) THEN
        DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
        IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
       ENDIF
       !start with the assumption that something will go wrong:
       tip_read=.FALSE.
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,number_of_tip_atoms
       IF (ios==0) THEN
        !Check if positive number of fragments has been specified
        IF (number_of_tip_atoms<1) THEN
         CALL report_error(82,exit_status=number_of_tip_atoms)
         CYCLE
        ENDIF
        !allocate memory for fragments
        ALLOCATE(fragment_list_tip(number_of_tip_atoms),STAT=allocstatus)
        IF (allocstatus/=0) CALL report_error(79,exit_status=allocstatus)
        !Read fragment record
        READ(3,IOSTAT=ios,FMT=*) (fragment_list_tip(counter),counter=1,number_of_tip_atoms,1)
        IF (ios==0) THEN
         tip_read=.TRUE.
         IF (VERBOSE_OUTPUT) WRITE(*,'(" ",I0,A)') number_of_tip_atoms," tip atoms read from fragment record."
        ELSE
         CALL report_error(81,exit_status=n+2)
         DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
         IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
        ENDIF
       ELSE
        CALL report_error(81,exit_status=n+2)
       ENDIF
      CASE ("base")
       IF (base_read) THEN
        DEALLOCATE(fragment_list_base,STAT=deallocstatus)
        IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
       ENDIF
       !start with the assumption that something will go wrong:
       base_read=.FALSE.
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,number_of_base_atoms
       IF (ios==0) THEN
        !Check if positive number of fragments has been specified
        IF (number_of_base_atoms<1) THEN
         CALL report_error(82,exit_status=number_of_base_atoms)
         CYCLE
        ENDIF
        !allocate memory for fragments
        ALLOCATE(fragment_list_base(number_of_base_atoms),STAT=allocstatus)
        IF (allocstatus/=0) CALL report_error(79,exit_status=allocstatus)
        !Read fragment record
        READ(3,IOSTAT=ios,FMT=*) (fragment_list_base(counter),counter=1,number_of_base_atoms,1)
        IF (ios==0) THEN
         base_read=.TRUE.
         IF (VERBOSE_OUTPUT) WRITE(*,'(" ",I0,A)') number_of_base_atoms," base atoms read from fragment record."
        ELSE
         CALL report_error(81,exit_status=n+2)
         DEALLOCATE(fragment_list_base,STAT=deallocstatus)
         IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
        ENDIF
       ELSE
        CALL report_error(81,exit_status=n+2)
       ENDIF
      CASE DEFAULT
       IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
      END SELECT
     ENDDO
     IF ((base_read).AND.(tip_read)) THEN
      !both fragment lists should be filled now, call the initialisation in MODULE MOLECULAR
      CALL initialise_fragments(fragment_list_tip,fragment_list_base,number_of_tip_atoms,number_of_base_atoms,molecule_type_index)
     ELSE
      CALL report_error(83)
     ENDIF
     IF (base_read) THEN
      DEALLOCATE(fragment_list_base,STAT=deallocstatus)
      IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     ENDIF
     IF (tip_read) THEN
      DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
      IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     ENDIF
    END SUBROUTINE read_input_for_reorientation

    SUBROUTINE read_input_for_rmmvcf()
    IMPLICIT NONE
    INTEGER :: n,a,b
    CHARACTER(LEN=32) :: inputstring
     REWIND 3
     READ(3,IOSTAT=ios,FMT=*) a,b
     IF (ios/=0) THEN
      molecule_type_index=1!equal to molecule 'a'
      molecule_type_index_b=2!equal to molecule 'b'
     ELSE
      !check if the molecule indices are sensible
      IF ((((a>0).AND.(b>0)).AND.(a/=b)).AND.&
       &((a<=give_number_of_molecule_types()).AND.(b<=give_number_of_molecule_types()))) THEN
       !the indices are valid.
       molecule_type_index=a
       molecule_type_index_b=b
      ELSE
       IF (VERBOSE_OUTPUT) WRITE(*,'(" invalid input for the indices (read ",I0," and ",I0,")")') a,b
       molecule_type_index=1
       molecule_type_index_b=2
      ENDIF
     ENDIF
     WRITE(*,'(" Using molecule_type ",I0," (for a) and ",I0," (for b).")') &
     &molecule_type_index,molecule_type_index_b
     !skip over line with operation_mode
     READ(3,*)
     !Here, the indices are ready and the body of the input file can be read.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
      IF (ios/=0) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      ENDIF
      SELECT CASE (TRIM(inputstring))
      CASE ("tmax")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "setting 'tmax' to default (=",tmax_default,")"
        tmax=tmax_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
       ENDIF
      CASE ("sampling_interval")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') &
        &"setting 'sampling_interval' to default (=",sampling_interval_default,")"
        sampling_interval=sampling_interval_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'sampling_interval' to ",sampling_interval
       ENDIF
      CASE ("skip_autocorrelation")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,skip_autocorr
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') "setting 'skip_autocorr' to default (=",skip_autocorr_default,")"
        skip_autocorr=skip_autocorr_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'skip_autocorr' to ",skip_autocorr
       ENDIF
      CASE ("quit")
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      CASE DEFAULT
       IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
      END SELECT
     ENDDO
    END SUBROUTINE read_input_for_rmmvcf

    SUBROUTINE read_input_for_dihedral_mode()!This subroutine is responsible for reading the body of the autocorrelation input file.
    IMPLICIT NONE
    INTEGER :: n,deallocstatus,inputinteger
    CHARACTER(LEN=32) :: inputstring
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: dihedral_member_indices !list of atom indices used for reporting dihedral angles to be passed on to the module MOLECULAR
     BACKSPACE 3
     READ(3,IOSTAT=ios,FMT=*) operation_mode,number_of_dihedral_conditions
     IF (ios/=0) CALL report_error(14,exit_status=ios)
     IF (number_of_dihedral_conditions<1) THEN
      CALL report_error(121,number_of_dihedral_conditions)
      RETURN
     ENDIF
     !formatted_dihedral_names have to be allocated here, because number_of_dihedral_conditions is not available before.
     ALLOCATE(formatted_dihedral_names(number_of_dihedral_conditions),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(11,exit_status=allocstatus)
     !allocate memory for autocorr_array
     ALLOCATE(autocorr_array(give_number_of_timesteps(),give_number_of_molecules_per_step(molecule_type_index)),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(16,exit_status=allocstatus)
     !allocate memory for the dihedral members list
     ALLOCATE(dihedral_member_indices(number_of_dihedral_conditions,4),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(11,exit_status=allocstatus)
     ALLOCATE(boundaries(number_of_dihedral_conditions,2),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     WRITE(*,*) "Dihedrals / lower boundary / upper boundary:"
     DO n = 1,number_of_dihedral_conditions,1 !reading in the specified dihedral conditions that have to be fulfilled simultaneously
      READ(3,IOSTAT=ios,FMT=*) dihedral_member_indices(n,:),boundaries(n,:)
      WRITE(inputstring,'(I0,"-",I0,"-",I0,"-",I0)') dihedral_member_indices(n,:)
      formatted_dihedral_names(n)=TRIM(ADJUSTL(inputstring))
      WRITE(*,'(" ",A,2F6.1)') TRIM(inputstring),boundaries(n,:)
      IF (ios/=0) CALL report_error(14,exit_status=ios)
     ENDDO
     CALL initialise_dihedrals(dihedral_member_indices,molecule_type_index,number_of_dihedral_conditions)
     DEALLOCATE(dihedral_member_indices,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
     !initialise variables for dihedral export.
     export_total_number=0
     export_dihedral=.FALSE.
     DO n=1,MAXITERATIONS,1
      READ(3,IOSTAT=ios,FMT=*) inputstring
      IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
      IF (ios/=0) THEN
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      ENDIF
      SELECT CASE (TRIM(inputstring))
      CASE ("fold")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,fold
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') " setting 'fold' to default (=",fold_default,")"
        fold=.TRUE.
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'fold' to ",fold
       ENDIF
      CASE ("export")
       !prepare the list of molecule indices whose dihedrals are to be exported in a separate file.
       !molecule_type_index as well as the dihedral list are initialised elsewhere.
       IF (export_dihedral) THEN !export_dihedral is already active.
        !another one has been requested. Check if total number is exceeded...
        IF (export_total_number>=give_number_of_molecules_per_step(molecule_type_index)) THEN
         !too many molecules to export
         CALL report_error(75)
        ELSE
         BACKSPACE 3
         READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
         IF (ios/=0) THEN
          CALL report_error(24,exit_status=ios)
          inputinteger=1
         ENDIF
         !Check if input is sensible
         IF ((inputinteger<1).OR.(inputinteger>give_number_of_molecules_per_step(molecule_type_index))) THEN
          CALL report_error(69) !unavailable molecule_index - abort.
         ELSE
          IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A,I0,A)') &
          &" Dihedrals for molecule number '",inputinteger,"' of type '",molecule_type_index,"' will be exported."
          export_total_number=export_total_number+1
          !add molecule_index to list
          export_list(export_total_number)=inputinteger
         ENDIF
        ENDIF
       ELSE !first molecule_index in list!
        BACKSPACE 3
        READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
        IF (ios/=0) THEN
         CALL report_error(24,exit_status=ios)
         inputinteger=1
        ENDIF
        !Check if input is sensible
        IF ((inputinteger<1).OR.(inputinteger>give_number_of_molecules_per_step(molecule_type_index))) THEN
         CALL report_error(69) !unavailable molecule_index - abort.
        ELSE
         IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A,I0,A)') &
         &" Dihedrals for molecule number '",inputinteger,"' of type '",molecule_type_index,"' will be exported."
         export_total_number=1
         export_dihedral=.TRUE.
         !first occurrence of "export". Allocate memory for list.
         ALLOCATE(export_list(give_number_of_molecules_per_step(molecule_type_index)),STAT=allocstatus)
         IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
         !add molecule_index to list
         export_list(export_total_number)=inputinteger
        ENDIF
       ENDIF
      CASE ("dump_verbose")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,dump_verbose
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') " setting 'dump_verbose' to default (=",dump_verbose_default,")"
        dump_verbose=.FALSE.
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'dump_verbose' to ",dump_verbose
       ENDIF
      CASE ("tmax")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') " setting 'tmax' to default (=",tmax_default,")"
        tmax=tmax_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') " setting 'tmax' to ",tmax
       ENDIF
      CASE ("quit")
       IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
       EXIT
      CASE ("skip_autocorrelation")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,skip_autocorr
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') " setting 'skip_autocorr' to default (=",skip_autocorr_default,")"
        skip_autocorr=skip_autocorr_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'skip_autocorr' to ",skip_autocorr
       ENDIF
      CASE ("bin_count")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,bin_count
       IF (ios/=0) THEN
        CALL report_error(24,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') " setting 'bin_count' to default (=",bin_count_default,")"
        bin_count=bin_count_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(" ",A,I0)') "setting 'bin_count' to ",bin_count
       ENDIF
       IF (bin_count<10) THEN
        bin_count=10
        CALL report_error(104,exit_status=10)
       ELSEIF (bin_count>1000) THEN
        bin_count=1000
        CALL report_error(104,exit_status=1000)
       ENDIF
      CASE DEFAULT
       IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
      END SELECT
     ENDDO
     !If necessary, check for issues with export_list
     IF (export_dihedral) THEN
      DO n=1,export_total_number-1,1
       IF (ANY(export_list((n+1):(export_total_number))==export_list(n))) THEN
        CALL report_error(76)
        EXIT
       ENDIF
      ENDDO
     ENDIF
    END SUBROUTINE read_input_for_dihedral_mode

  END SUBROUTINE initialise_autocorrelation

  !finalises the autocorrelation module.
  SUBROUTINE finalise_autocorrelation()
  IMPLICIT NONE
  INTEGER :: deallocstatus
   IF (TRIM(operation_mode)=="dihedral") THEN
    DEALLOCATE(autocorr_array,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
    DEALLOCATE(formatted_dihedral_names,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
    DEALLOCATE(boundaries,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
    IF (dump_verbose) THEN
     DEALLOCATE(PES_subset_independent,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
     IF (number_of_dihedral_conditions==2) THEN
      DEALLOCATE(PES_subset_dependent,STAT=deallocstatus)
      IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
     ENDIF
    ENDIF
    IF (export_dihedral) THEN
     DEALLOCATE(export_list,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
    ENDIF
   ENDIF
   IF (ALLOCATED(vc_components)) THEN
    DEALLOCATE(vc_components,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
   ENDIF
  END SUBROUTINE finalise_autocorrelation

  !The following SUBROUTINE is tailored toward the calculation of RMM-VCFs.
  SUBROUTINE cross_correlation()
  IMPLICIT NONE
  REAL(KIND=WORKING_PRECISION),ALLOCATABLE :: average_velocities(:,:,:)!ua and ub for rmm-vcf's. first dimension: timestep index, second dimension: 1 (=a) or 2 (=b), third dimension: vector with velocity.
  REAL(WORKING_PRECISION),ALLOCATABLE :: correlation_function(:),autocorrelation_function(:,:)!second dimension of the autocorrelation_function is the two particles.
  INTEGER :: na,nb !the number of molecules for the two types
  REAL(KIND=WORKING_PRECISION) :: nareal,nbreal
  INTEGER :: allocstatus,deallocstatus,nsteps
  REAL(KIND=WORKING_PRECISION) :: xa,xb,integral_cross,integral_a,integral_b,ma,mb,temperature,temperature_b
  REAL(KIND=WORKING_PRECISION) :: firstvalue_a,area_a,area_b,firstvalue_b!these are initialised/computed by 'report_autocorrelation_function'
  REAL(KIND=WORKING_PRECISION) :: firstvalue,area!these are initialised/computed by 'report_correlation_function'
  REAL(KIND=WORKING_PRECISION) :: delta,D0,D_distinct!these are initialised/computed by 'report_summary'
  INTEGER,ALLOCATABLE :: x_num(:)!number of averages taken for self-contributions, not including the averages over the particles in one snapshot.
   nsteps=give_number_of_timesteps()
   IF ((tmax>(nsteps-1)).OR.(tmax<1)) THEN
    tmax=(nsteps-1)
    CALL report_error(28,exit_status=INT(tmax))
   ENDIF
   na=give_number_of_molecules_per_step(molecule_type_index)
   nb=give_number_of_molecules_per_step(molecule_type_index_b)
   nareal=DFLOAT(na)
   nbreal=DFLOAT(nb)
   !xa and xb are the molar fractions
   xa=nareal/(DFLOAT(na+nb))
   xb=nbreal/(DFLOAT(na+nb))
   !also request the masses
   ma=give_mass_of_molecule(molecule_type_index)
   mb=give_mass_of_molecule(molecule_type_index_b)
   !allocate average velocity memory.
   ALLOCATE(average_velocities(nsteps,2,3),STAT=allocstatus)!needs ~ 50 MB for 1ns@1fs with DP
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   !first, fill average_velocities array with the required information.
   CALL compute_average_velocities()
   !allocate memory for the correlation_function (from t=0 to t=tmax)
   ALLOCATE(correlation_function(tmax+1),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   correlation_function(:)=0.0d0
   !then, calculate the lambda quantity
   CALL correlate_average_velocities()
   DEALLOCATE(average_velocities,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
   !Print the intermediate information.
   CALL report_correlation_function()
   IF (.NOT.(skip_autocorr)) THEN
    !allocate memory for the autocorrelation_function (self-diffusion)
    ALLOCATE(autocorrelation_function(tmax+1,2),STAT=allocstatus)
    IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
    ALLOCATE(x_num(tmax+1),STAT=allocstatus)
    IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
    !compute the autocorrelation function / self-diffusion / VACF
    IF (PARALLEL_OPERATION) THEN
     CALL compute_self_contribution_parallel(sampling_interval)
    ELSE
     CALL compute_self_contribution(sampling_interval)
    ENDIF
   ENDIF
   IF (.NOT.(skip_autocorr)) THEN
    CALL report_autocorrelation_function()
    CALL report_summary()
    IF (BOX_VOLUME_GIVEN) CALL report_conductivity()
    DEALLOCATE(autocorrelation_function,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
    DEALLOCATE(x_num,STAT=deallocstatus)
    IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   ENDIF
   DEALLOCATE(correlation_function,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   CONTAINS

    !This subroutine computes eq (5) from https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162
    SUBROUTINE compute_average_velocities()
    IMPLICIT NONE
    INTEGER :: timestep_counter,molecule_counter
     !initialise average velocities.
     average_velocities(:,:,:)=0.0d0
     !THEN, add all the molecular velocities together, for all the timesteps.
     DO timestep_counter=1,nsteps,1
      !first molecule (particle 'a')
      DO molecule_counter=1,na,1
       average_velocities(timestep_counter,1,:)=average_velocities(timestep_counter,1,:)&
       &+give_center_of_mass(timestep_counter,molecule_type_index,molecule_counter)
      ENDDO
      !second molecule (particle 'b')
      DO molecule_counter=1,nb,1
       average_velocities(timestep_counter,2,:)=average_velocities(timestep_counter,2,:)&
       &+give_center_of_mass(timestep_counter,molecule_type_index_b,molecule_counter)
      ENDDO
     ENDDO
     !normalise by number of members to arrive at equation (5):
     average_velocities(:,1,:)=average_velocities(:,1,:)/nareal
     average_velocities(:,2,:)=average_velocities(:,2,:)/nbreal
    END SUBROUTINE compute_average_velocities

    !This subroutine computes eq (8) from https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162
    SUBROUTINE compute_self_contribution(sampling)
    IMPLICIT NONE
    INTEGER :: molecule_counter,local_tmax,startstep,timeline
    REAL(KIND=WORKING_PRECISION) :: temp_value
    REAL(WORKING_PRECISION),ALLOCATABLE :: initial_velocities_a(:,:),initial_velocities_b(:,:)!variables to store initial velocities, i.e. 'ua,n(t0)'
    INTEGER,INTENT(IN) :: sampling
     !Compute velocity autocorrelation functions
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Computing self contributions"
     !allocate memory for the initial velocities of every starting timestep
     ALLOCATE(initial_velocities_a(na,3),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     ALLOCATE(initial_velocities_b(nb,3),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     !first, initialise
     autocorrelation_function(:,:)=0.0d0
     x_num(:)=0
     !Here, the outer loop is chosen to be the starting step, i.e. 't0'
     DO startstep=1,nsteps,sampling
      !Thus, first get the initial velocities of every particle for the two molecule types
      DO molecule_counter=1,na,1
       initial_velocities_a(molecule_counter,:)=give_center_of_mass(startstep,molecule_type_index,molecule_counter)
      ENDDO
      DO molecule_counter=1,nb,1
       initial_velocities_b(molecule_counter,:)=give_center_of_mass(startstep,molecule_type_index_b,molecule_counter)
      ENDDO
      !These velocities now have to be correlated with 'themselves' at a later time step.
      !Careful: startstep+timeline can of course not exceed the number of available steps.
      IF ((startstep+tmax)>nsteps) THEN!This check is the necessary price for switching the two loops, which in turn was introduced to speed up sequential read.
       local_tmax=(nsteps-startstep)
      ELSE
       local_tmax=tmax
      ENDIF
      DO timeline=0,local_tmax,1
       !increment number of averages taken, is the same for both functions even if they are based on different numbers of molecules.
       x_num(timeline+1)=x_num(timeline+1)+1
       temp_value=0.0d0
       DO molecule_counter=1,na,1
        !the center of mass is acutally a good quantity here (gives the velocity of the center of mass)
        temp_value=temp_value+DOT_PRODUCT(&
        &give_center_of_mass(startstep+timeline,molecule_type_index,molecule_counter),&
        &initial_velocities_a(molecule_counter,:))
       ENDDO
       autocorrelation_function(timeline+1,1)=autocorrelation_function(timeline+1,1)+(temp_value/nareal)
       temp_value=0.0d0
       DO molecule_counter=1,nb,1
        !you don't believe me? Take the derivative of the centre of mass with respect to time...
        !The code might be wrong, but the intention was right.
        temp_value=temp_value+DOT_PRODUCT(&
        &give_center_of_mass(startstep+timeline,molecule_type_index_b,molecule_counter),&
        &initial_velocities_b(molecule_counter,:))
       ENDDO
       autocorrelation_function(timeline+1,2)=autocorrelation_function(timeline+1,2)+(temp_value/nbreal)
      ENDDO
     ENDDO
     !deallocate initial velocity memory
     DEALLOCATE(initial_velocities_a,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(initial_velocities_b,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !normalise by the number of startsteps / averages taken.
     DO timeline=0,tmax,1
      !both molecule a and b have to be normalised, hence ':'. They are based on the same number of starting steps.
      autocorrelation_function(timeline+1,:)=autocorrelation_function(timeline+1,:)/DFLOAT(x_num(timeline+1))
     ENDDO
     !multiply with 1/3 to arrive at equation (8):
     autocorrelation_function(:,:)=autocorrelation_function(:,:)/3.0d0
    END SUBROUTINE compute_self_contribution

    !This subroutine computes eq (8) from https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162
    !like compute_self_contribution, but parallelised! The other one is kept since I don't trust my own parallelisation skills.
    !results seem to be the same though (after a lot of pain admittedly)
    SUBROUTINE compute_self_contribution_parallel(sampling)
    IMPLICIT NONE
    !$ INTERFACE
    !$  FUNCTION OMP_get_num_threads()
    !$  INTEGER :: OMP_get_num_threads
    !$  END FUNCTION OMP_get_num_threads
    !$ END INTERFACE
    INTEGER,INTENT(IN) :: sampling
    INTEGER :: molecule_counter,startstep,timeline
    REAL(WORKING_PRECISION),ALLOCATABLE :: temp_function(:,:)
    REAL(WORKING_PRECISION),ALLOCATABLE :: initial_velocities_a(:,:),initial_velocities_b(:,:)!variables to store initial velocities, i.e. 'ua,n(t0)'
    INTEGER,ALLOCATABLE :: x_num_temp(:)
     !Compute velocity autocorrelation functions
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Computing self contributions"
     autocorrelation_function(:,:)=0.0d0
     x_num(:)=0
     !$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
     !$OMP PRIVATE(temp_function,x_num_temp,initial_velocities_a,initial_velocities_b,molecule_counter)
     !$OMP SINGLE
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (autocorrelation)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     !allocate memory for temporary functions (used for parallelisation)
     ALLOCATE(temp_function(tmax+1,2),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     ALLOCATE(x_num_temp(tmax+1),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     !allocate memory for the initial velocities of every starting timestep
     ALLOCATE(initial_velocities_a(na,3),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     ALLOCATE(initial_velocities_b(nb,3),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     x_num_temp(:)=0
     !first, initialise
     temp_function(:,:)=0.0d0
     !Here, the outer loop is chosen to be the starting step, i.e. 't0'
     !$OMP DO SCHEDULE(STATIC,1)
     DO startstep=1,nsteps,sampling
      !Thus, first get the initial velocities of every particle for the two molecule types
      DO molecule_counter=1,na,1
       initial_velocities_a(molecule_counter,:)=give_center_of_mass(startstep,molecule_type_index,molecule_counter)
      ENDDO
      DO molecule_counter=1,nb,1
       initial_velocities_b(molecule_counter,:)=give_center_of_mass(startstep,molecule_type_index_b,molecule_counter)
      ENDDO
      CALL iterate_timesteps_self_contributions(&
      &startstep,temp_function,initial_velocities_a,initial_velocities_b,x_num_temp)
     ENDDO
     !$OMP END DO
     !CRITICAL directive to properly update the autocorrelation_function
     !$OMP CRITICAL
     autocorrelation_function(:,:)=autocorrelation_function(:,:)+temp_function(:,:)
     !$OMP END CRITICAL
     !CRITICAL directive to properly update x_num
     !$OMP CRITICAL
     x_num(:)=x_num(:)+x_num_temp(:)
     !$OMP END CRITICAL
     !deallocate private memory used for parallelisation
     DEALLOCATE(temp_function,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(x_num_temp,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !deallocate initial velocity memory
     DEALLOCATE(initial_velocities_a,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(initial_velocities_b,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
     !normalise by the number of startsteps / averages taken.
     DO timeline=0,tmax,1
      !both molecule a and b have to be normalised, hence ':'. They are based on the same number of starting steps.
      autocorrelation_function(timeline+1,:)=autocorrelation_function(timeline+1,:)/DFLOAT(x_num(timeline+1))
     ENDDO
     !multiply with 1/3 to arrive at equation (8):
     autocorrelation_function(:,:)=autocorrelation_function(:,:)/3.0d0
    END SUBROUTINE compute_self_contribution_parallel

    SUBROUTINE iterate_timesteps_self_contributions(startstep,temp_function,vel_a,vel_b,x_num_temp)
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: temp_value
    INTEGER :: molecule_counter,local_tmax
    INTEGER :: timeline
    INTEGER,INTENT(IN) :: startstep
    REAL(WORKING_PRECISION),ALLOCATABLE,INTENT(INOUT) :: temp_function(:,:)
    REAL(WORKING_PRECISION),ALLOCATABLE,INTENT(IN) :: vel_a(:,:)
    REAL(WORKING_PRECISION),ALLOCATABLE,INTENT(IN) :: vel_b(:,:)
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: x_num_temp(:)
     !The velocities now have to be correlated with 'themselves' at a later time step.
     !Careful: startstep+timeline can of course not exceed the number of available steps.
     IF ((startstep+tmax)>nsteps) THEN!still ok this way round, because the average velocities should be allocated in the calling routine anyway.
      local_tmax=(nsteps-startstep)
     ELSE
      local_tmax=tmax
     ENDIF
     DO timeline=0,local_tmax,1
      !increment number of averages taken, is the same for both functions even if they are based on different numbers of molecules.
      x_num_temp(timeline+1)=x_num_temp(timeline+1)+1
      temp_value=0.0d0
      DO molecule_counter=1,na,1
       !the center of mass is acutally a good quantity here (gives the velocity of the center of mass)
       temp_value=temp_value+DOT_PRODUCT(&
       &give_center_of_mass(startstep+timeline,molecule_type_index,molecule_counter),&
       &vel_a(molecule_counter,:))
      ENDDO
      temp_function(timeline+1,1)=temp_function(timeline+1,1)+(temp_value/nareal)
      temp_value=0.0d0
      DO molecule_counter=1,nb,1
       !you don't believe me? Take the derivative of the centre of mass with respect to time...
       !The code might be wrong, but the intention was right.
       temp_value=temp_value+DOT_PRODUCT(&
       &give_center_of_mass(startstep+timeline,molecule_type_index_b,molecule_counter),&
       &vel_b(molecule_counter,:))
      ENDDO
      temp_function(timeline+1,2)=temp_function(timeline+1,2)+(temp_value/nbreal)
     ENDDO
    END SUBROUTINE iterate_timesteps_self_contributions

    !This subroutine computes eq (4) from https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162 without the prefactors (just the part in {...})
    SUBROUTINE correlate_average_velocities()
    IMPLICIT NONE
    INTEGER :: timeline,timesteps
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Computing cross contributions"
     !'timeline' is the argument of the correlation_function, i.e. the time shift of the original function.
     !for example, if the current shift is timeline=1000, and there are 10000 timesteps in total,
     !then the argument has to be evaluated from (h(0+0)-<h>)(h(0+1000)-<h>) up to (h(9000+0)-<h>)(h(9000+1000)-<h>)
     DO timeline=0,tmax,1
      !inner loop iterates over the subset of the autocorr_array
      DO timesteps=1,give_number_of_timesteps()-timeline,1
       !this is the central part of the whole autocorrelation process.
       correlation_function(timeline+1)=correlation_function(timeline+1)&
       &+DOT_PRODUCT((average_velocities((timesteps+timeline),1,:)-average_velocities((timesteps+timeline),2,:)),&
       &(average_velocities(timesteps,1,:)-average_velocities(timesteps,2,:)))
       !equation (4), that's the part [ua(t)-ub(t)][ua(t0)-ub(t0)]
      ENDDO
      !Normalise result. For the above example, this would be division by 9000
      correlation_function(timeline+1)=correlation_function(timeline+1)/DFLOAT(give_number_of_timesteps()-timeline)
      !equation (4), until here the part <[ua(t)-ub(t)][ua(t0)-ub(t0)]>
     ENDDO! End of outer loop over time shifts
     !multiply with N:
     correlation_function(:)=correlation_function(:)*DFLOAT(na+nb)
     !divide by 3, multiply with xa xb:
     correlation_function(:)=correlation_function(:)*(1.0d0/3.0d0)*xa*xb
    END SUBROUTINE correlate_average_velocities

    !reporting the correlation_function, integrating on the fly.
    SUBROUTINE report_correlation_function()
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing mean molecular relative VCFs into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"RMM_VCF","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"RMM_VCF",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains velocity cross-correlation coefficients based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "reference: https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162"
     WRITE(3,'(A15,3A25)') "timeline","lambda_ab(t)","integral","C_ab(t)"
     firstvalue=correlation_function(1)
     xa=nareal/(DFLOAT(na+nb))
     xb=nbreal/(DFLOAT(na+nb))
     IF (VERBOSE_OUTPUT) THEN
      IF (na==nb) THEN
       WRITE(*,*) "same number of molecules for both types, xa=xb=1/2"
      ELSE
       WRITE(*,'(" xa = ",F4.2," xb = ",F4.2)') xa,xb
      ENDIF
      WRITE(*,'(" Normalise cross-correlation function: dividing by",E10.3,",")') firstvalue
      !Based on equation (13), calculate the Temperature!
      temperature=firstvalue*ma*mb/(boltzmann*(xa*ma+xb*mb))
      !correct to Kelvin
      temperature=(temperature*1.0d7)/avogadro
      !print the temperature as a valuable check!
      WRITE(*,'(" which corresponds to a temperature of ",F7.2," K.")') temperature
      WRITE(*,'(" Check your results if this is far off from what you would expect!")') 
     ENDIF
     integral_cross=0.0d0
     DO timeline=0,tmax-1,1
      WRITE(3,'(I15,E25.16,E25.16,F25.16)') timeline*TIME_SCALING_FACTOR,SNGL(correlation_function(timeline+1)),&
      &SNGL(integral_cross),SNGL(correlation_function(timeline+1)/firstvalue)
      area=correlation_function(timeline+2)+correlation_function(timeline+1)
      area=area*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
      integral_cross=integral_cross+area
     ENDDO
     IF (VERBOSE_OUTPUT) WRITE(*,'(" Integral of lambda_ab over ",I0," timesteps is",E13.6)') tmax,integral_cross
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE report_correlation_function

    !reporting the autocorrelation_function, integrating on the fly.
    SUBROUTINE report_autocorrelation_function()
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline
     !Opening output file for VACFs
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing VACFs into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"VACFs","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"VACFs",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains velocity autocorrelation coefficients based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "reference: https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162"
     WRITE(3,'(A15,7A25)') "timeline","lambdas_a(t)","integral","Cs_a(t)","lambdas_b(t)","integral","Cs_b(t)"
     firstvalue_a=autocorrelation_function(1,1)
     firstvalue_b=autocorrelation_function(1,2)
     IF (VERBOSE_OUTPUT) THEN
      WRITE(*,'(" Normalise autocorrelation functions: dividing by ",E9.3," and ",E9.3," (molecules a/",I0," and b/",I0,"),")')&
      & firstvalue_a,firstvalue_b,molecule_type_index,molecule_type_index_b
      !Based on equation (15), calculate the Temperature!
      temperature=firstvalue_a*ma/boltzmann
      temperature_b=firstvalue_b*mb/boltzmann
      !correct to Kelvin
      temperature=(temperature*1.0d7)/avogadro
      temperature_b=(temperature_b*1.0d7)/avogadro
      !print the temperature as a valuable check!
      WRITE(*,'(" which corresponds to temperatures of ",F7.2," K (a/",I0,") and ",F7.2," K (b/",I0,").")')&
      &temperature,molecule_type_index,temperature_b,molecule_type_index_b
      WRITE(*,'(" Check your results if this is far off from what you would expect!")') 
     ENDIF
     integral_a=0.0d0
     integral_b=0.0d0
     DO timeline=0,tmax-1,1
      WRITE(3,'(I15,2E25.16,F25.16,2E25.16,F25.20)') timeline*TIME_SCALING_FACTOR,& !I25: time variable "timeline"
      &SNGL(autocorrelation_function((timeline+1),1)),SNGL(integral_a),& !autocorrelation function for molecule a "lambdas_a(t)" and its "integral"
      &SNGL(autocorrelation_function((timeline+1),1)/firstvalue_a),& !The "normalised" function for molecule a
      &SNGL(autocorrelation_function((timeline+1),2)),SNGL(integral_b),& !autocorrelation function for molecule b "lambdas_b(t)" and its "integral"
      &SNGL(autocorrelation_function((timeline+1),2)/firstvalue_b) !The "normalised" function for molecule b
      !integrating the trapezoids:
      area_a=autocorrelation_function((timeline+2),1)+autocorrelation_function((timeline+1),1)
      area_a=area_a*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
      integral_a=integral_a+area_a
      area_b=autocorrelation_function((timeline+2),2)+autocorrelation_function((timeline+1),2)
      area_b=area_b*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
      integral_b=integral_b+area_b
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE report_autocorrelation_function

    !reporting a summary with the 'C' functions.
    SUBROUTINE report_summary()
    IMPLICIT NONE
    !Output formats for diffusion coefficients
   13 FORMAT ("   ",A5," ",E15.6)
    LOGICAL :: connected
    INTEGER :: ios,timeline
    REAL(KIND=WORKING_PRECISION) :: C0,firstvalue_C0
     !Opening output file for the summary
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing summary into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"VCF_summary","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"VCF_summary",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the core quantities based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "reference: https://journals.aps.org/pre/pdf/10.1103/PhysRevE.50.1162"
     WRITE(3,'(6A15)') "timeline","Cs_a(t)","Cs_b(t)","C0(t)","C_ab(t)","delta(t)"
     D0=xb*integral_a+xa*integral_b
     delta=(integral_cross-D0)/D0
     D_distinct=(integral_cross-D0)/(xa*xb)
     WRITE(*,*) "Reference frame independent diffusion coefficients:"
     WRITE(*,13) "Ds_a ",integral_a!equation (8)
     WRITE(*,13) "Ds_b ",integral_b!equation (8)
     WRITE(*,13) "D0_ab",D0!equation (7)
     WRITE(*,13) "D_ab ",integral_cross!equation (4)
     WRITE(*,13) "Dd_ab",D_distinct!equation (19)
     WRITE(*,13) "delta",delta!equation (20)
     WRITE(*,'("   Distinct contributions: ",F0.2,"%")') 1.0d2*(D0-integral_cross)/D0
     WRITE(*,*) "If input units are femtoseconds and Angströms, then divide values by 10 to arrive at cm²/s."
     firstvalue_C0=xb*autocorrelation_function(1,1)+xa*autocorrelation_function(1,2)
     DO timeline=0,tmax,1
      C0=(xb*autocorrelation_function((timeline+1),1)+xa*autocorrelation_function((timeline+1),2))/firstvalue_C0
      WRITE(3,'(I15,5F15.10)') timeline*TIME_SCALING_FACTOR,& !I15: time variable "timeline"
      &SNGL(autocorrelation_function((timeline+1),1)/firstvalue_a),& !C1s in Figure 1
      &SNGL(autocorrelation_function((timeline+1),2)/firstvalue_b),& !C2s in Figure 1
      &SNGL(C0),& !C0
      &SNGL(correlation_function(timeline+1)/firstvalue),& !C12 in Figure 1
      &SNGL((correlation_function(timeline+1)/firstvalue)-C0)!The delta function... But defined as follows: delta=C12(t)-C0(t)
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE report_summary

    !calculates conductivity
    SUBROUTINE report_conductivity()
    IMPLICIT NONE
    !Output format for conductivities
   14 FORMAT ("   ",A8," ",E15.6," (",EN13.4," Scm²/mol)")!molar
   15 FORMAT ("   ",A8," ",E15.6," (",EN13.4," S/cm)")!specific
   16 FORMAT ("   ",A5," ",E15.6,A14)
    REAL(KIND=WORKING_PRECISION) :: prefactor!conversion factor for conductivities
    REAL(KIND=WORKING_PRECISION) :: molar_conductivity_self,molar_conductivity_distinct!molar conductivities, divided into self- and distinct contributions
    REAL(KIND=WORKING_PRECISION) :: ca,cb!the charges of particles a and b
    INTEGER :: ios
    LOGICAL :: connected
     !Writing into file, because why not.
     !Opening output file for the conductivities
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing conductivity and diffusion data into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"conductivity_rmmvcf","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"conductivity_rmmvcf",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains transport quantities based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*)
     WRITE(3,*) "Diffusion properties from integrated velocity correlation functions:"
     WRITE(3,*) "reference: dx.doi.org/10.1103/PhysRevE.50.1162"
     WRITE(3,16) "Ds_a ",integral_a," equation  (8)"
     WRITE(3,16) "Ds_b ",integral_b," equation  (8)"
     WRITE(3,16) "D0_ab",D0," equation  (7)"
     WRITE(3,16) "D_ab ",integral_cross," equation  (4)"
     WRITE(3,16) "Dd_ab",D_distinct," equation (19)"
     WRITE(3,16) "delta",delta," equation (20)"
     WRITE(3,*) "If input units are femtoseconds and Angströms, then divide values by 10 to arrive at cm²/s."
     WRITE(3,*)
     WRITE(3,*) "VACF electrical conductivity:"
     WRITE(3,*) "reference: https://aip.scitation.org/doi/10.1063/1.466191"
     WRITE(3,*) "molar conductivity, assuming T=298.15K:"
     ca=DFLOAT(give_charge_of_molecule(molecule_type_index))
     cb=DFLOAT(give_charge_of_molecule(molecule_type_index_b))
     prefactor=((elementary_charge**2)/(boltzmann*298.15))
     !the following conductivities are calculated WITHOUT prefactor.
     molar_conductivity_self=((ca**2)*xa*integral_a+(cb**2)*xb*integral_b)!equation (29) in JCP 99, page 3983 (1993)
     molar_conductivity_distinct=-ca*cb*xa*xb*D_distinct!equation (30) in JCP 99, page 3983 (1993)
     WRITE(*,*)"calculating electrical conductivity. Reference: https://aip.scitation.org/doi/10.1063/1.466191"
     IF (VERBOSE_OUTPUT) WRITE(*,*)&
     &"be aware of the definition of your concentration. If in doubt, use the specific conductivities."
     WRITE(*,*)"molar conductivity, assuming T=298.15K:"
     WRITE(*,14) "self    ",prefactor*molar_conductivity_self,prefactor*molar_conductivity_self*(avogadro/10.0)
     WRITE(*,14) "distinct",prefactor*molar_conductivity_distinct,prefactor*molar_conductivity_distinct*(avogadro/10.0)
     WRITE(*,14) "total   ",(molar_conductivity_self+molar_conductivity_distinct)*prefactor,&
     &(molar_conductivity_self+molar_conductivity_distinct)*prefactor*(avogadro/10.0)
     !same into file:
     WRITE(3,14) "self    ",prefactor*molar_conductivity_self,prefactor*molar_conductivity_self*(avogadro/10.0)
     WRITE(3,14) "distinct",prefactor*molar_conductivity_distinct,prefactor*molar_conductivity_distinct*(avogadro/10.0)
     WRITE(3,14) "total   ",(molar_conductivity_self+molar_conductivity_distinct)*prefactor,&
     &(molar_conductivity_self+molar_conductivity_distinct)*prefactor*(avogadro/10.0)
     !This part is only called when BOX_VOLUME_GIVEN, so no further action is required.
     prefactor=prefactor*((na+nb)/(give_box_volume()))!multiply with concentration to arrive at specific conductivity
     WRITE(3,*) "specific conductivity, assuming T=298.15K:"
     WRITE(3,15) "self    ",prefactor*molar_conductivity_self,prefactor*molar_conductivity_self*1.0d23
     WRITE(3,15) "distinct",prefactor*molar_conductivity_distinct,prefactor*molar_conductivity_distinct*1.0d23
     WRITE(3,15) "total   ",(molar_conductivity_self+molar_conductivity_distinct)*prefactor,&
     &(molar_conductivity_self+molar_conductivity_distinct)*prefactor*1.0d23
     WRITE(*,*)"specific conductivity, assuming T=298.15K:"
     WRITE(*,15) "self    ",prefactor*molar_conductivity_self,prefactor*molar_conductivity_self*1.0d23
     WRITE(*,15) "distinct",prefactor*molar_conductivity_distinct,prefactor*molar_conductivity_distinct*1.0d23
     WRITE(*,15) "total   ",(molar_conductivity_self+molar_conductivity_distinct)*prefactor,&
     &(molar_conductivity_self+molar_conductivity_distinct)*prefactor*1.0d23
     WRITE(*,'("   Predicted Haven Ratio: ",F0.2,"%")') 1.0d2*&
     &((molar_conductivity_self+molar_conductivity_distinct)/molar_conductivity_self)!Reporting the Haven Ratio.
     WRITE(3,'("   Predicted Haven Ratio: ",F0.2,"%")') 1.0d2*&
     &((molar_conductivity_self+molar_conductivity_distinct)/molar_conductivity_self)!Reporting the Haven Ratio.
     WRITE(3,*) "be aware of the definition of your concentration. If in doubt, use the specific conductivities."
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE report_conductivity

  END SUBROUTINE cross_correlation

  !The following SUBROUTINE computes the rotational reorientation function, as described in:
  !'Theory of Simple Liquids, Hansen/McDonald, Chapter 11' dx.doi.org/10.1016/B978-0-12-387032-2.00011-8
  SUBROUTINE reorientational_autocorrelation()
  IMPLICIT NONE
  INTEGER :: na !the number of molecules for the observed molecule type
  REAL(KIND=WORKING_PRECISION) :: nareal
  INTEGER :: nsteps,allocstatus,deallocstatus
  REAL(WORKING_PRECISION),ALLOCATABLE :: time_correlation_function(:)
  INTEGER,ALLOCATABLE :: x_num(:)!number of averages taken for autocorrelation, not including the averages over the particles in one snapshot.
   nsteps=give_number_of_timesteps()
   IF ((tmax>(nsteps-1)).OR.(tmax<1)) THEN
    tmax=(nsteps-1)
    CALL report_error(28,exit_status=INT(tmax))
   ENDIF
   na=give_number_of_molecules_per_step(molecule_type_index)
   nareal=DFLOAT(na)
   !allocate memory for the correlation_function (from t=0 to t=tmax)
   ALLOCATE(time_correlation_function(tmax+1),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   !allocate memory for the number of averages taken
   ALLOCATE(x_num(tmax+1),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   !Report fragment information
   CALL give_fragment_information(tip_fragment=.FALSE.)
   CALL give_fragment_information(tip_fragment=.TRUE.)
   !Main part: compute the tcf
   CALL compute_reorientational_tcf_parallel(sampling_interval)
   !Report the tcf
   CALL report_time_correlation_function()
   !Deallocate memory for tcf and average counter
   DEALLOCATE(time_correlation_function,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   DEALLOCATE(x_num,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   CONTAINS

    !This subroutine computes eq (11.11.1) from "THEORY OF SIMPLE LIQUIDS" (Hansen / McDonald)
    !'sampling' is the sampling interval.
    SUBROUTINE compute_reorientational_tcf_parallel(sampling)
    IMPLICIT NONE
    !$ INTERFACE
    !$  FUNCTION OMP_get_num_threads()
    !$  INTEGER :: OMP_get_num_threads
    !$  END FUNCTION OMP_get_num_threads
    !$ END INTERFACE
    INTEGER,INTENT(IN) :: sampling
    INTEGER :: molecule_counter,startstep,timeline
    REAL(WORKING_PRECISION),ALLOCATABLE :: temp_function(:)
    REAL(WORKING_PRECISION),ALLOCATABLE :: initial_orientation(:,:)!variable to store initial position
    INTEGER,ALLOCATABLE :: x_num_temp(:)
     !Compute reorientational autocorrelation functions
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Computing reorientational correlation function"
     time_correlation_function(:)=0.0d0
     x_num(:)=0
     !$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
     !$OMP PRIVATE(temp_function,x_num_temp,initial_orientation,molecule_counter)
   !  !$OMP PRIVATE(iterate_timesteps_self_contributions)
     !$OMP SINGLE
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (time correlation function)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     !allocate memory for temporary functions (used for parallelisation)
     ALLOCATE(temp_function(tmax+1),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     ALLOCATE(x_num_temp(tmax+1),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     !allocate memory for the initial positions of every starting timestep
     ALLOCATE(initial_orientation(na,3),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     x_num_temp(:)=0
     !first, initialise
     temp_function(:)=0.0d0
     !Here, the outer loop is chosen to be the starting step, i.e. 't0'
     !$OMP DO SCHEDULE(STATIC,1)
     DO startstep=1,nsteps,sampling
      DO molecule_counter=1,na,1
       !initial orientation = normalised vector from 'base atom' to 'tip atom' at the initial timestep
       initial_orientation(molecule_counter,:)=&
       &give_tip_fragment(startstep,molecule_counter)-give_base_fragment(startstep,molecule_counter)
       CALL normalize3D(initial_orientation(molecule_counter,:))
      ENDDO
      CALL iterate_timesteps_tcf(startstep,temp_function,initial_orientation,x_num_temp)
     ENDDO
     !$OMP END DO
     !CRITICAL directive to properly update the time_correlation_function
     !$OMP CRITICAL
     time_correlation_function(:)=time_correlation_function(:)+temp_function(:)
     !$OMP END CRITICAL
     !CRITICAL directive to properly update x_num
     !$OMP CRITICAL
     x_num(:)=x_num(:)+x_num_temp(:)
     !$OMP END CRITICAL
     !deallocate private memory used for parallelisation
     DEALLOCATE(temp_function,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(x_num_temp,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !deallocate initial position memory
     DEALLOCATE(initial_orientation,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
     !normalise by the number of startsteps / averages taken.
     DO timeline=0,tmax,1
      !both molecule a and b have to be normalised, hence ':'. They are based on the same number of starting steps.
      time_correlation_function(timeline+1)=time_correlation_function(timeline+1)/DFLOAT(x_num(timeline+1))
     ENDDO
     !is already normalised, because Pl(cos0)=Pl(1)=1 (Pl = legendre polynomial of order 'l')
    END SUBROUTINE compute_reorientational_tcf_parallel

    SUBROUTINE iterate_timesteps_tcf(startstep,temp_function,initial_orientation,x_num_temp)
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: temp_value,second_vector(3)
    INTEGER :: molecule_counter,local_tmax
    INTEGER :: timeline
    INTEGER,INTENT(IN) :: startstep
    REAL(WORKING_PRECISION),ALLOCATABLE,INTENT(INOUT) :: temp_function(:)
    REAL(WORKING_PRECISION),ALLOCATABLE,INTENT(IN) :: initial_orientation(:,:)
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: x_num_temp(:)
     !Careful: startstep+timeline can of course not exceed the number of available steps.
     IF ((startstep+tmax)>nsteps) THEN!still ok this way round, because the average positions should be allocated in the calling routine anyway.
      local_tmax=(nsteps-startstep)
     ELSE
      local_tmax=tmax
     ENDIF
     DO timeline=0,local_tmax,1
      !increment number of averages taken:
      x_num_temp(timeline+1)=x_num_temp(timeline+1)+1
      temp_value=0.0d0
      DO molecule_counter=1,na,1
       !The first vector is already stored in initial_orientation as unit vector.
       !get the second vector:
       second_vector(:)=&
       &give_tip_fragment(startstep+timeline,molecule_counter)-give_base_fragment(startstep+timeline,molecule_counter)
       !then, normalise:
       CALL normalize3D(second_vector(:))
       !take the two normalised values and compute the legendre polynomial of the dot product.
       !the dot product is equal to the cosine of the angle between the two vectors, which is the desired quantity.
       temp_value=temp_value+&
       &legendre_polynomial(DOT_PRODUCT(second_vector(:),initial_orientation(molecule_counter,:)),legendre_order)
      ENDDO
      temp_function(timeline+1)=temp_function(timeline+1)+(temp_value/nareal)
     ENDDO
    END SUBROUTINE iterate_timesteps_tcf

    !reporting the time correlation function, integrating on the fly.
    SUBROUTINE report_time_correlation_function()
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline
    CHARACTER(LEN=25) :: tcf_name
    REAL(KIND=WORKING_PRECISION) :: integral,area,firstvalue
     !Opening output file for rotcorr
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing TCF into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"rotcorr","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"rotcorr",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the reorientational time correlation function based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "reference: THEORY OF SIMPLE LIQUIDS (Hansen / McDonald), fourth edition, chapter 11.11."
     WRITE(tcf_name,'("C",I0,"(t)")') legendre_order
     WRITE(3,'(A15,2A25)') "timeline",TRIM(tcf_name),"integral"
     firstvalue=time_correlation_function(1)
     integral=0.0d0
     DO timeline=0,tmax-1,1
      WRITE(3,'(I15,F25.16,E25.16)') timeline*TIME_SCALING_FACTOR,& !I25: time variable "timeline"
      &SNGL(time_correlation_function(timeline+1)),integral !time correlation function and its integral
      !integrating the trapezoids:
      area=time_correlation_function(timeline+2)+time_correlation_function(timeline+1)
      area=area*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
      integral=integral+area
     ENDDO
     IF (VERBOSE_OUTPUT) WRITE(*,'(" last area value is",E8.1," percent of total integral.")') 100.0*ABS(area)/integral
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE report_time_correlation_function

  END SUBROUTINE reorientational_autocorrelation

  SUBROUTINE dihedral_autocorrelation()
  IMPLICIT NONE
  INTEGER :: timestep_counter,local_incidence,molecule_counter,n,m,ios,unit_number
  REAL(KIND=WORKING_PRECISION) :: standard_deviation,local_average_h
  CHARACTER(LEN=32) :: filename_export
  LOGICAL :: connected
   !open files for dihedral export, if necessary
   IF (export_dihedral) THEN
    !Files will be opened in unit 10,11,12,...
    DO n=1,export_total_number,1
     unit_number=n+9
     INQUIRE(UNIT=unit_number,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=unit_number)
     WRITE(filename_export,'("dihedral_export",I0,"_",I0)') n,export_list(n)
     OPEN(UNIT=unit_number,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//filename_export,IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(unit_number,'(A,I0,A,I0)') "This file contains the dihedrals for the molecule number "&
     &,export_list(n)," of type ",molecule_type_index
     WRITE(unit_number,FMT='(A15)',ADVANCE="NO") "step"
     DO m=1,number_of_dihedral_conditions,1 
      WRITE(unit_number,FMT='(A15)',ADVANCE="NO") TRIM(formatted_dihedral_names(m))
     ENDDO
     WRITE(unit_number,*)
    ENDDO
   ENDIF
   CALL fill_dihedral_array
   !close units
   IF (export_dihedral) THEN
    DO n=1,export_total_number,1
     unit_number=n+9
     CLOSE(UNIT=unit_number)
    ENDDO
   ENDIF
   !compute average of population operator
   average_h=(DFLOAT(global_incidence)/DFLOAT(number_of_entries_in_array))
   !compute standard deviation, print information
   WRITE(*,'(" ",I0," out of ",I0," molecules were within the specified boundaries. <h>=",F0.2,"%")')&
   &global_incidence,number_of_entries_in_array,average_h*100.0
   IF (give_number_of_timesteps()>5) THEN
    DO timestep_counter=1,give_number_of_timesteps(),1
     local_incidence=0
     DO molecule_counter=1,give_number_of_molecules_per_step(molecule_type_index),1
      IF (autocorr_array(timestep_counter,molecule_counter)) local_incidence=local_incidence+1
     ENDDO
     local_average_h=(DFLOAT(local_incidence)/DFLOAT(give_number_of_molecules_per_step(molecule_type_index)))
     standard_deviation=standard_deviation+(local_average_h-average_h)**2!sum over the squares of the deviations from the global average
    ENDDO
    standard_deviation=standard_deviation/(DFLOAT(give_number_of_timesteps()-1))!Correct for number of observations
    standard_deviation=SQRT(standard_deviation)
    WRITE(*,'(" absolute value of <h>: ",E11.5,", standard sample deviation: ",E11.5)')average_h,standard_deviation
    IF (VERBOSE_OUTPUT) WRITE(*,*) "(standard deviation calculated from the snapshot-wise ensemble averages)"
   ENDIF
   IF (dump_verbose) THEN
    CALL dump_PES
    CALL dump_incidence
   ENDIF
   CONTAINS

    SUBROUTINE fill_dihedral_array()
    IMPLICIT NONE
    INTEGER :: molecule_index,export_counter,unit_number
    INTEGER :: bin_position,bin_position2
    INTEGER :: deallocstatus,allocstatus,condition_counter
    LOGICAL :: within_boundary
    REAL(KIND=GENERAL_PRECISION),ALLOCATABLE :: dihedral_list(:)
     ALLOCATE(dihedral_list(number_of_dihedral_conditions),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     !There is some potential for parallelisation here. Maybe. Probably not much gained though.
     DO timestep_counter=1,give_number_of_timesteps(),1
      !Is there a molecule whose dihedrals are to be exported?
      IF (export_dihedral) THEN
       DO export_counter=1,export_total_number,1
        unit_number=export_counter+9
        WRITE(unit_number,ADVANCE="NO",FMT='(I15)') timestep_counter
        CALL give_dihedrals(dihedral_list,timestep_counter,export_list(export_counter))
        DO condition_counter = 1,number_of_dihedral_conditions,1
         WRITE(unit_number,ADVANCE="NO",FMT='(F15.3)') dihedral_list(condition_counter)
        ENDDO
        WRITE(unit_number,*)
       ENDDO
      ENDIF
      DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
       CALL give_dihedrals(dihedral_list,timestep_counter,molecule_index)
       !First of all, bin into PES_subsets if requested.
       IF (dump_verbose) THEN
        IF (number_of_dihedral_conditions==2) THEN !two-dimensional PES subset - also sort into the dependent one.
         bin_position=INT(dihedral_list(1)/(360.0/bin_count))
         IF ((bin_position<0).OR.(bin_position>bin_count)) CALL report_error(25,exit_status=bin_position)
         bin_position2=INT(dihedral_list(2)/(360.0/bin_count))
         IF ((bin_position2<0).OR.(bin_position2>bin_count)) CALL report_error(25,exit_status=bin_position2)
         PES_subset_dependent(bin_position,bin_position2)=&
         &PES_subset_dependent(bin_position,bin_position2)+1
        ENDIF
        DO condition_counter = 1,number_of_dihedral_conditions,1
         bin_position=INT(dihedral_list(condition_counter)/(360.0/bin_count))
         IF ((bin_position<0).OR.(bin_position>bin_count)) CALL report_error(25,exit_status=bin_position)
         PES_subset_independent(condition_counter,bin_position)=&
         &PES_subset_independent(condition_counter,bin_position)+1
        ENDDO
       ENDIF
       within_boundary=.TRUE.!start with the assumption that all conditions are fulfilled.
       DO condition_counter = 1,number_of_dihedral_conditions,1 !checking the *unfolded* conditions that have to be fulfilled simultaneously
        !stay in loop as long as there are conditions left AND none of them being violated so far!
        IF ((dihedral_list(condition_counter)<boundaries(condition_counter,1))&
        &.OR.((dihedral_list(condition_counter)>boundaries(condition_counter,2)))) THEN!seems to be out of range - check for folded values!
         within_boundary=.FALSE.
         EXIT !The exit statement is correct here, because the folded part is a separate loop.
        ENDIF
       ENDDO
       IF (fold.AND.(.NOT.within_boundary)) THEN!dihedral_list is not within specified boundary, but maybe with folding?
        ! if folding was not requested, then everything stays as it is. Otherwise, it has to be checked that the folded condition is violated as well.
        within_boundary=.TRUE.!same assumption as before.
        !Since we already know at this point that the unfolded condition is NOT fulfilled (otherwise this part would not have been entered),
        !it is correct and necessary to reset within_boundary to TRUE. A distinct working variable is not required.
        dihedral_list(:)=360.0-dihedral_list(:)
        !reasoning behind this line: start with a<x<b, with a and b being the lower and upper boundaries.
        !Thus, (360-b)<x<(360-a) is the folded condition.
        !(360-b)<x<(360-a) <=> (-b)<(x-360)<(-a) <=> b>(360-x)>a, and x is essentially the dihedral_list.
        ! The loop itself stays the same.
        DO condition_counter = 1,number_of_dihedral_conditions,1 !checking the *folded* conditions that have to be fulfilled simultaneously
         !stay in loop as long as there are conditions left AND none of them being violated so far!
         IF ((dihedral_list(condition_counter)<boundaries(condition_counter,1))&
         &.OR.((dihedral_list(condition_counter)>boundaries(condition_counter,2)))) THEN!seems to be out of range - check for folded values!
          within_boundary=.FALSE.
          EXIT
         ENDIF
        ENDDO
       ENDIF
       !Finally, switch the corresponding entry in the boolean array.
       autocorr_array(timestep_counter,molecule_index)=within_boundary
       IF (within_boundary) global_incidence=global_incidence+1
      ENDDO
     ENDDO
     DEALLOCATE(dihedral_list,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
    END SUBROUTINE fill_dihedral_array

    SUBROUTINE dump_PES()
    IMPLICIT NONE
    !Output formats - AUTOCORRELATION module
   3 FORMAT (2F15.3,I15) !Module AUTOCORRELATION - PES_subset_dependent
    INTEGER :: condition_counter,bin_position,ios,bin_position2
    LOGICAL :: connected
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     !First, write the File with the independent PES subset.
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing independent PES subset population into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"PES_subset_independent","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"PES_subset_independent",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the independent incidence (=counts) for the dihedrals specified in '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "bin_position angle_start angle_end ",&
     &(TRIM(formatted_dihedral_names(condition_counter))//"  ",condition_counter=1,number_of_dihedral_conditions,1)
     DO bin_position=0,bin_count,1
      WRITE(3,*) bin_position,(360.0/bin_count)*bin_position,(360.0/bin_count)*(bin_position+1),&
      &(PES_subset_independent(condition_counter,bin_position),condition_counter=1,number_of_dihedral_conditions,1)
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
     IF (number_of_dihedral_conditions==2) THEN !two-dimensional PES - also write dependent subset.
      IF (VERBOSE_OUTPUT) WRITE(*,*) "2D PES - writing dependent PES subset population into file '",&
      &TRIM(ADJUSTL(OUTPUT_PREFIX))//"PES_subset_dependent","'"
      OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"PES_subset_dependent",IOSTAT=ios)
      IF (ios/=0) CALL report_error(26,exit_status=ios)
      WRITE(3,*) "This file contains the 2D dependent incidence (=counts) for the dihedrals specified in '"&
      &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"' (angle_start given)"
      WRITE(3,'(3A15)') TRIM(formatted_dihedral_names(1)),TRIM(formatted_dihedral_names(2)),"counts"
      !Iterating over the two dihedral angles / bin positions
      DO bin_position=0,bin_count,1
       DO bin_position2=0,bin_count,1
        WRITE(3,3) (360.0/bin_count)*bin_position,(360.0/bin_count)*bin_position2,&
        &PES_subset_dependent(bin_position,bin_position2)
       ENDDO
      ENDDO
      ENDFILE 3
      CLOSE(UNIT=3)
     ENDIF
    END SUBROUTINE dump_PES
    
    SUBROUTINE dump_incidence()
    IMPLICIT NONE
    !Output formats - AUTOCORRELATION module
   5 FORMAT (2I15,F15.3) !Module AUTOCORRELATION - local_incidence
    LOGICAL :: connected
    INTEGER :: molecule_counter,timestep_counter,local_incidence,ios
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing local incidences ('shares') timestep-wise into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"local_incidences","'"
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"local_incidences",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the local incidences for the dihedrals specified in '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "(i.e. counts or 'shares' within boundaries per timestep)"
     WRITE(3,'(3A15)') "timestep","count","share"
     DO timestep_counter=1,give_number_of_timesteps(),1
      local_incidence=0
      DO molecule_counter=1,give_number_of_molecules_per_step(molecule_type_index),1
       IF (autocorr_array(timestep_counter,molecule_counter)) local_incidence=local_incidence+1
      ENDDO
      WRITE(3,5) timestep_counter,local_incidence,&
      &FLOAT(local_incidence*100)/FLOAT(give_number_of_molecules_per_step(molecule_type_index))
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE dump_incidence

  END SUBROUTINE dihedral_autocorrelation

  SUBROUTINE calculate_autocorrelation_function_from_binary_array()!This subroutine converts the binary autocorr_array to the correlation function.
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$ END INTERFACE
  REAL(WORKING_PRECISION),ALLOCATABLE :: autocorrelation_function(:),temp_function(:)!the quantity C(t), which is at first formed uncorrected.
  INTEGER :: n,allocstatus,deallocstatus
   !first, check for sensible input.
   IF ((tmax>(give_number_of_timesteps()-1)).OR.(tmax<100)) THEN
    tmax=(give_number_of_timesteps()-1)
    CALL report_error(28,exit_status=INT(tmax))
   ENDIF
   !allocate memory for the autocorrelation_function (from t=0 to t=tmax)
   ALLOCATE(autocorrelation_function(tmax+1),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   !initialise autocorrelation_function.
   autocorrelation_function(:)=0.0_DP
   !prepare parallelisation, if required.
   !temp_function is a local variable to avoid racing conditions
   !$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(temp_function)
   !$OMP SINGLE
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (intermittent autocorrelation function)"
   !$  CALL timing_parallel_sections(.TRUE.)
   !$ ENDIF
   !$OMP END SINGLE
   !allocate memory and initialise temp_function for every member of the team.
   ALLOCATE(temp_function(tmax+1),STAT=allocstatus)
   IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
   temp_function(:)=0.0_DP
   !iterate over molecules and pass chunks to the subroutine that iterates over timesteps
   !$OMP DO SCHEDULE(STATIC,1)
   DO n=1,give_number_of_molecules_per_step(molecule_type_index),1
    temp_function(:)=temp_function(:)+iterate_timesteps(autocorr_array(:,n))
   ENDDO
   !$OMP END DO
   !CRITICAL directive to properly update the autocorrelation_function
   !$OMP CRITICAL
   autocorrelation_function(:)=autocorrelation_function(:)+temp_function(:)
   !$OMP END CRITICAL
   DEALLOCATE(temp_function,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   !$OMP END PARALLEL
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
   !$  CALL timing_parallel_sections(.FALSE.)
   !$ ENDIF
   !normalise autocorrelation_function by number of molecules
   autocorrelation_function(:)=autocorrelation_function(:)/DFLOAT(give_number_of_molecules_per_step(molecule_type_index))
   !print / report autocorrelation function
   CALL report_autocorrelation_function()
   DEALLOCATE(autocorrelation_function,STAT=deallocstatus)
   IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
   CONTAINS

    FUNCTION iterate_timesteps(chunk)!This function calculates the autocorrelation function for ONE molecule.
    IMPLICIT NONE
    REAL(WORKING_PRECISION) :: iterate_timesteps(tmax+1)
    LOGICAL,INTENT(IN) :: chunk(:)
    INTEGER :: timeline,timesteps
     !initialise function
     iterate_timesteps(:)=0.0d0
     !'timeline' is the argument of the autocorrelation_function, i.e. the time shift of the original function.
     !for example, if the current shift is timeline=1000, and there are 10000 timesteps in total,
     !then the argument has to be evaluated from (h(0+0)-<h>)(h(0+1000)-<h>) up to (h(9000+0)-<h>)(h(9000+1000)-<h>)
     DO timeline=0,tmax,1
      !inner loop iterates over the whole chunk, i.e. the subset of the autocorr_array for the molecule in question.
      DO timesteps=1,give_number_of_timesteps()-timeline,1
       !this is the central part of the whole autocorrelation process.
       IF (chunk(timesteps)) THEN!h(t0) is one
        IF (chunk(timesteps+timeline)) THEN!h(t0+t) is one
         iterate_timesteps(timeline+1)=iterate_timesteps(timeline+1)&
         &+(1.0d0-average_h)*(1.0d0-average_h)
        ELSE!h(t0+t) is zero
         iterate_timesteps(timeline+1)=iterate_timesteps(timeline+1)&
         &+(1.0d0-average_h)*(0.0d0-average_h)
        ENDIF
       ELSE!h(t0) is zero
        IF (chunk(timesteps+timeline)) THEN!h(t0+t) is one
         iterate_timesteps(timeline+1)=iterate_timesteps(timeline+1)&
         &+(0.0d0-average_h)*(1.0d0-average_h)
        ELSE!h(t0+t) is zero
         iterate_timesteps(timeline+1)=iterate_timesteps(timeline+1)&
         &+(0.0d0-average_h)*(0.0d0-average_h)
        ENDIF
       ENDIF
      ENDDO
      !Normalise result. For the above example, this would be division by 9000
      iterate_timesteps(timeline+1)=iterate_timesteps(timeline+1)/DFLOAT(give_number_of_timesteps()-timeline)
     ENDDO! End of outer loop over time shifts
    END FUNCTION iterate_timesteps

    SUBROUTINE report_autocorrelation_function
    !Output formats - AUTOCORRELATION module
   1 FORMAT (I15,3F15.10) !Module AUTOCORRELATION - autocorrelation_function file
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline
    REAL(KIND=WORKING_PRECISION) :: firstvalue
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing autocorrelation function into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"autocorrelation_function","'"
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"autocorrelation_function",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the intermittent autocorrelation function based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     firstvalue=autocorrelation_function(1)!is already squared.
     IF (VERBOSE_OUTPUT) WRITE(*,'(" normalise autocorrelation function: dividing by ",F5.3)') firstvalue
     WRITE(3,*) "C(t)=<(h(t0+t)-<h>)(h(t0)-<h>)>/<(h(t0)-<h>)**2>"!up until here, autocorrelation_function is not normalised
     WRITE(3,'(4A15)') "timeline","C(t)","log(C(t))","C(t)uncorr"
     DO timeline=0,tmax,1
      WRITE(3,1) timeline*TIME_SCALING_FACTOR,SNGL(autocorrelation_function(timeline+1)/firstvalue),&
      &SNGL(LOG(autocorrelation_function(timeline+1)/firstvalue)),SNGL(autocorrelation_function(timeline+1))
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
    END SUBROUTINE report_autocorrelation_function

  END SUBROUTINE calculate_autocorrelation_function_from_binary_array

  !The following subroutine calculates *only* the microscopic charge current time autocorrelation function.
  !Handles the operation mode "conductivity".
  SUBROUTINE overall_conductivity()
  IMPLICIT NONE
  REAL(WORKING_PRECISION),ALLOCATABLE :: correlation_function(:)
  INTEGER :: nsteps
  CALL initialise_overall_conductivity()
  IF (.NOT.(ERROR_CODE==120)) THEN
   CALL compute_overall_conductivity(sampling_interval)
   CALL report_overall_conductivity()
  ENDIF
  CALL finalise_overall_conductivity()

   CONTAINS

    SUBROUTINE report_overall_conductivity()
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline
    REAL(KIND=WORKING_PRECISION) :: integral,area,conversion_factor
     !Opening output file for conductivity
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing CACF into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"CACF","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"conductivity_CACF",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the microscopic CACF and its running integral based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "reference: e.g. equation (10) in J. Chem. Phys., 2002, 116, 3018–3026."
     WRITE(3,'("    timeline      <j(t)j(0)>        integral")')
     integral=0.0d0
     DO timeline=0,tmax-1,1
      WRITE(3,*) timeline*TIME_SCALING_FACTOR,& !time variable "timeline"
      &SNGL(correlation_function(timeline+1)),SNGL(integral) !time correlation function and its integral
      !integrating the trapezoids:
      area=correlation_function(timeline+2)+correlation_function(timeline+1)
      area=area*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
      integral=integral+area
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
     IF (BOX_VOLUME_GIVEN) THEN
      conversion_factor=1.0d23/(3.0d0*boltzmann*298.15*give_box_volume())
      WRITE(*,'(A,E16.8,A,E16.8,A)')" specific conductivity is",integral*conversion_factor," S/cm (integral =",integral,")"
      WRITE(*,*)"(Calculated assuming T=298.15K, Angströms, and femtoseconds)"
     ENDIF
     IF (VERBOSE_OUTPUT) WRITE(*,'(" last area value is",E8.1," percent of total integral.")') 100.0*ABS(area)/integral
    END SUBROUTINE report_overall_conductivity

    SUBROUTINE initialise_overall_conductivity()
    IMPLICIT NONE
    INTEGER :: allocstatus,charge,molecule_type_counter
    LOGICAL :: charged_particles_present
     nsteps=give_number_of_timesteps()
     IF ((tmax>(nsteps-1)).OR.(tmax<1)) THEN
      tmax=(nsteps-1)
      CALL report_error(28,exit_status=INT(tmax))
     ENDIF
     !allocate memory for the correlation_function (from t=0 to t=tmax)
     ALLOCATE(correlation_function(tmax+1),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     correlation_function(:)=0.0d0
     !check that there are charged particles!
     charged_particles_present=.FALSE.
     DO molecule_type_counter=1,give_number_of_molecule_types()
      charge=give_charge_of_molecule(molecule_type_counter)
      IF (charge/=0) THEN
       charged_particles_present=.TRUE.
       EXIT
      ENDIF
     ENDDO
     IF (.NOT.(charged_particles_present)) CALL report_error(120)
    END SUBROUTINE initialise_overall_conductivity

    SUBROUTINE finalise_overall_conductivity()
    IMPLICIT NONE
    INTEGER :: deallocstatus
     DEALLOCATE(correlation_function,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
    END SUBROUTINE finalise_overall_conductivity

    !This subroutine computes <vi,vj>, and relevant sums thereof, as given in vc_components.
    !'sampling' is the sampling interval.
    SUBROUTINE compute_overall_conductivity(sampling)
    IMPLICIT NONE
    !$ INTERFACE
    !$  FUNCTION OMP_get_num_threads()
    !$  INTEGER :: OMP_get_num_threads
    !$  END FUNCTION OMP_get_num_threads
    !$ END INTERFACE
    INTEGER,INTENT(IN) :: sampling
    INTEGER :: startstep,timeline
    INTEGER :: allocstatus,deallocstatus
    REAL(WORKING_PRECISION),ALLOCATABLE :: temp_function(:)
     !Compute velocity correlation functions
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Computing microscopic charge current time autocorrelation function (CACF)."
     correlation_function(:)=0.0d0
     !$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
     !$OMP PRIVATE(temp_function)
     !$OMP SINGLE
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (overall ecaf for conductivity)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     !allocate memory for temporary functions (used for parallelisation)
     ALLOCATE(temp_function(tmax+1),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     temp_function(:)=0.0d0
     !Here, the outer loop is chosen to be the starting step, i.e. 't0'
     !$OMP DO SCHEDULE(STATIC,1)
     DO startstep=1,nsteps,sampling
      CALL iterate_timesteps_microscopic_charge_current_tacf(startstep,temp_function)
     ENDDO
     !$OMP END DO
     !CRITICAL directive to properly update the correlation_function
     !$OMP CRITICAL
     correlation_function(:)=correlation_function(:)+temp_function(:)
     !$OMP END CRITICAL
     !deallocate private memory used for parallelisation
     DEALLOCATE(temp_function,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
     !normalise by the number of startsteps / averages taken.
     DO timeline=0,tmax,1
      correlation_function(timeline+1)=correlation_function(timeline+1)/DFLOAT(1+(nsteps-(timeline+1))/sampling)
     ENDDO
     !multiply with e**2
     correlation_function(:)=correlation_function(:)*(elementary_charge**2)
    END SUBROUTINE compute_overall_conductivity

    !Calculates <j(t)*j(0)> in equation (10) in J. Chem. Phys., 2002, 116, 3018–3026.
    !Simple and quick.
    SUBROUTINE iterate_timesteps_microscopic_charge_current_tacf&
    &(startstep,temp_function)
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: initial_mcc(3),temp_mcc(3),molecular_mcc(3)
    INTEGER :: molecule_counter,local_tmax,molecule_type_counter,timeline,charge
    INTEGER,INTENT(IN) :: startstep
    REAL(WORKING_PRECISION),INTENT(INOUT) :: temp_function(:)
    REAL(WORKING_PRECISION) :: chargereal
     !first, initialise
     initial_mcc(:)=0.0d0
     DO molecule_type_counter=1,give_number_of_molecule_types()
      molecular_mcc(:)=0.0d0
      charge=give_charge_of_molecule(molecule_type_counter)
      IF (charge==0) CYCLE
      chargereal=FLOAT(charge)
      DO molecule_counter=1,give_number_of_molecules_per_step(molecule_type_counter),1
       molecular_mcc(:)=molecular_mcc(:)+give_center_of_mass(startstep,molecule_type_counter,molecule_counter)
      ENDDO
      !here, molecular_mcc contains the sum of all velocities. Now, weigh with charge and add to initial_mcc.
      initial_mcc(:)=initial_mcc(:)+molecular_mcc(:)*chargereal
     ENDDO
     !now, we have j(0)! Calculate j(t)
     !Careful: startstep+timeline can of course not exceed the number of available steps.
     IF ((startstep+tmax)>nsteps) THEN!still ok this way round, because the average positions should be allocated in the calling routine anyway.
      local_tmax=(nsteps-startstep)
     ELSE
      local_tmax=tmax
     ENDIF
     DO timeline=0,local_tmax,1
      !Calculate the j(t), which is the temp_mcc
      temp_mcc(:)=0.0d0
      DO molecule_type_counter=1,give_number_of_molecule_types()
       molecular_mcc(:)=0.0d0
       charge=give_charge_of_molecule(molecule_type_counter)
       IF (charge==0) CYCLE
       chargereal=FLOAT(charge)
       DO molecule_counter=1,give_number_of_molecules_per_step(molecule_type_counter),1
        molecular_mcc(:)=molecular_mcc(:)+give_center_of_mass(startstep+timeline,molecule_type_counter,molecule_counter)
       ENDDO
       !here, molecular_mcc contains the sum of all velocities. Now, weigh with charge and add to temp_mcc.
       temp_mcc(:)=temp_mcc(:)+molecular_mcc(:)*chargereal
      ENDDO
      !correlate, and presto.
      temp_function(timeline+1)=temp_function(timeline+1)+DOT_PRODUCT(initial_mcc(:),temp_mcc(:))
     ENDDO
    END SUBROUTINE iterate_timesteps_microscopic_charge_current_tacf

  END SUBROUTINE overall_conductivity

  !The following subroutine calculates velocity cross- and autocorrelation functions.
  !this includes distinct contributions, which are expensive to compute.
  !Handles the operation modes "vcf" and "ecaf".
  SUBROUTINE velocity_correlation()
  IMPLICIT NONE
  REAL(WORKING_PRECISION),ALLOCATABLE :: correlation_function(:,:)! first entry is time shift. second dimension of the correlation_function corresponds to the entry in vc_components
  INTEGER :: nsteps,ncomponents
  LOGICAL :: normalise
   IF (.NOT.(ALLOCATED(vc_components))) THEN
    CALL report_error(0)
    RETURN
   ENDIF
   ncomponents=SIZE(vc_components)
   CALL check_vc_components()
   !remind user of reference frames.
   IF ((ERROR_CODE==116).OR.(ERROR_CODE==117)) RETURN
   CALL initialise_velocity_correlation()
   CALL compute_velocity_correlation(sampling_interval)
   SELECT CASE (TRIM(operation_mode))
   CASE ("vcf")
    CALL report_vcfs()
   CASE ("ecaf")
    CALL report_ecaf()
   CASE DEFAULT
    CALL report_error(0)
   END SELECT
   CALL finalise_velocity_correlation

   CONTAINS

    SUBROUTINE report_vcfs()
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline,component_counter,molecule_type_counter
    REAL(KIND=WORKING_PRECISION) :: integral(ncomponents),area,conversion_factor,temperature
    CHARACTER(LEN=64) :: component_name
     !Opening output file for velocity correlation functions
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing vcf components into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"VCF_components","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"VCF_components",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the VCF components based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "See, for example, equations (A6) to (A10) in J. Phys. Chem. B, 2011, 115, 13212–13221."
     WRITE(3,FMT='(A16)',ADVANCE="NO") "timeline"
     DO component_counter=1,ncomponents,1
      IF (vc_components(component_counter)%self) THEN
       WRITE(component_name,FMT='(A,I0,A,I0,A)') "<vs",vc_components(component_counter)%reference_type,&
       &"(0)*vs",vc_components(component_counter)%observed_type,"(t)>"
       IF (DEVELOPERS_VERSION) THEN
        temperature=1.0d7*give_mass_of_molecule(vc_components(component_counter)%reference_type)*&
        &correlation_function(1,component_counter)/(3.0d0*boltzmann*avogadro)
        WRITE(*,'("  ! T (from ",A,"):",3EN10.1)') TRIM(component_name),SNGL(temperature)
       ENDIF
      ELSE
       WRITE(component_name,FMT='(A,I0,A,I0,A)') "<vd",vc_components(component_counter)%reference_type,&
       &"(0)*vd",vc_components(component_counter)%observed_type,"(t)>"
      ENDIF
      WRITE(3,FMT='(A16)',ADVANCE="NO") TRIM(component_name)
     ENDDO
     WRITE(3,*)
     integral(:)=0.0d0
     DO timeline=0,tmax-1,1
      WRITE(3,FMT='(I16)',ADVANCE="NO") timeline*TIME_SCALING_FACTOR !time variable "timeline"
      DO component_counter=1,ncomponents,1
       WRITE(3,FMT='(E16.8)',ADVANCE="NO") SNGL(correlation_function(timeline+1,component_counter)) !time correlation function
       !integrating the trapezoids:
       area=correlation_function(timeline+2,component_counter)+correlation_function(timeline+1,component_counter)
       area=area*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
       integral(component_counter)=integral(component_counter)+area
      ENDDO
      WRITE(3,*)
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
     WRITE(*,*) "Diffusion properties from integrated velocity correlation components:"
     WRITE(*,*) "(including the factor 'N' and '1/3', see J. Phys. Chem. B, 2011, 115, 13212–13221.)"
     DO component_counter=1,ncomponents,1
      IF (vc_components(component_counter)%self) THEN
       conversion_factor=1.0d0
       WRITE(component_name,FMT='(A,I0)') "Ds_",vc_components(component_counter)%reference_type
      ELSE
       conversion_factor=FLOAT(give_total_number_of_molecules_per_step())
       WRITE(component_name,FMT='(A,I0,I0)') "Dd_",&
       &vc_components(component_counter)%reference_type,vc_components(component_counter)%observed_type
      ENDIF
      conversion_factor=conversion_factor/3.0d1 !times 1/3 for dimension, and divide by 10 to convert to cm**2/s
      WRITE(*,'("   ",A," =",E17.8," cm²/s")') TRIM(component_name),integral(component_counter)*conversion_factor
     ENDDO
     WRITE(*,*) "(Assuming femtoseconds and Angströms as input units)"
    END SUBROUTINE report_vcfs

    SUBROUTINE report_ecaf()
    IMPLICIT NONE
    LOGICAL :: connected
    INTEGER :: ios,timeline,component_counter,i,j
    REAL(KIND=WORKING_PRECISION) :: integral(ncomponents),area,conversion_factor,total,self,distinct
    CHARACTER(LEN=64) :: component_name
     !Opening output file for velocity correlation functions
     INQUIRE(UNIT=3,OPENED=connected)
     IF (connected) CALL report_error(27,exit_status=3)
     IF (VERBOSE_OUTPUT) WRITE(*,*) "writing CACF (ECAF, ECCF, ...) components into file '",&
     &TRIM(ADJUSTL(OUTPUT_PREFIX))//"CACF_components","'"
     OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"CACF_components",IOSTAT=ios)
     IF (ios/=0) CALL report_error(26,exit_status=ios)
     WRITE(3,*) "This file contains the CACF components based on the input file '"&
     &,TRIM(FILENAME_AUTOCORRELATION_INPUT),"'"
     WRITE(3,*) "See also: THEORY OF SIMPLE LIQUIDS (Hansen / McDonald), fourth edition, equation (10.5.2)"
     WRITE(3,FMT='(A16)',ADVANCE="NO") "timeline"
     DO component_counter=1,ncomponents,1
      !i and j contain the molecule type index
      i=vc_components(component_counter)%reference_type
      j=vc_components(component_counter)%observed_type
      IF (vc_components(component_counter)%self) THEN
       WRITE(component_name,FMT='(A,I0,A,I0,A,I0,A,I0,A)') "z",i,"z",j,"<vs",i,"*vs",j,">"
      ELSE
       WRITE(component_name,FMT='(A,I0,A,I0,A,I0,A,I0,A)') "z",i,"z",j,"<vd",i,"*vd",j,">"
      ENDIF
      WRITE(3,FMT='(A17)',ADVANCE="NO") TRIM(component_name)
      i=give_charge_of_molecule(i)
      j=give_charge_of_molecule(j)
      !now, i and j contain the charge of the corresponding molecules / molecule types.
      correlation_function(:,component_counter)=correlation_function(:,component_counter)*FLOAT(i*j)
     ENDDO
     WRITE(3,*)
     integral(:)=0.0d0
     DO timeline=0,tmax-1,1
      WRITE(3,FMT='(I16)',ADVANCE="NO") timeline*TIME_SCALING_FACTOR !time variable "timeline"
      DO component_counter=1,ncomponents,1
       WRITE(3,FMT='(E17.8)',ADVANCE="NO") SNGL(correlation_function(timeline+1,component_counter)) !time correlation function
       !integrating the trapezoids:
       area=correlation_function(timeline+2,component_counter)+correlation_function(timeline+1,component_counter)
       area=area*(DFLOAT(TIME_SCALING_FACTOR)/2.0d0)
       integral(component_counter)=integral(component_counter)+area
      ENDDO
      WRITE(3,*)
     ENDDO
     ENDFILE 3
     CLOSE(UNIT=3)
     WRITE(*,*) "Conductivities from integrated electric current autocorrelation function components:"
     conversion_factor=(1.0d23*(elementary_charge**2))/(3.0d0*boltzmann*298.15*give_box_volume())
     total=0.0d0
     self=0.0d0
     distinct=0.0d0
     integral(:)=integral(:)*conversion_factor
     DO component_counter=1,ncomponents,1
      IF (vc_components(component_counter)%self) THEN
       self=self+integral(component_counter)
       WRITE(component_name,FMT='(A,I0)') "sigma_s_",vc_components(component_counter)%reference_type
       WRITE(*,'("   ",A,"   =",EN13.4," S/cm")') TRIM(component_name),integral(component_counter)
      ELSE
       distinct=distinct+integral(component_counter)
       WRITE(component_name,FMT='(A,I0,I0)') "sigma_d_",&
       &vc_components(component_counter)%reference_type,vc_components(component_counter)%observed_type
       WRITE(*,'("   ",A,"  =",EN13.4," S/cm")') TRIM(component_name),integral(component_counter)
      ENDIF
      total=total+integral(component_counter)
     ENDDO
     WRITE(component_name,FMT='(A10)') "sigma_self"
     WRITE(*,'("   ",A," =",EN13.4," S/cm")') TRIM(component_name),self
     WRITE(component_name,FMT='(A10)') "sigma_dist"
     WRITE(*,'("   ",A," =",EN13.4," S/cm")') TRIM(component_name),distinct
     WRITE(component_name,FMT='(A10)') "sigma_tot."
     WRITE(*,'("   ",A10," =",EN13.4," S/cm")') TRIM(component_name),total
     WRITE(*,*) "(Specific conductivity, assuming T=298.15K and femtoseconds/Angströms as input units)"
    END SUBROUTINE report_ecaf

    SUBROUTINE check_vc_components()
    IMPLICIT NONE
    INTEGER :: component_counter,ref,obs
     DO component_counter=1,ncomponents,1
      !This is a double check with extended output.
      ref=vc_components(component_counter)%reference_type
      obs=vc_components(component_counter)%observed_type
      IF ((ref<1).OR.(ref>give_number_of_molecule_types())) THEN
       CALL report_error(116,ref)
       RETURN
      ELSEIF ((obs<1).OR.(obs>give_number_of_molecule_types())) THEN
       CALL report_error(116,obs)
       RETURN
      ELSEIF ((ref/=obs).AND.(vc_components(component_counter)%self)) THEN
       CALL report_error(117)
       RETURN
      ENDIF
     ENDDO
    END SUBROUTINE check_vc_components

    SUBROUTINE initialise_velocity_correlation()
    IMPLICIT NONE
    INTEGER :: allocstatus,component_counter
     nsteps=give_number_of_timesteps()
     IF ((tmax>(nsteps-1)).OR.(tmax<1)) THEN
      tmax=(nsteps-1)
      CALL report_error(28,exit_status=INT(tmax))
     ENDIF
     !allocate memory for the correlation_function (from t=0 to t=tmax)
     ALLOCATE(correlation_function(tmax+1,ncomponents),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     correlation_function(:,:)=0.0d0
     SELECT CASE (TRIM(operation_mode))!no further output necessary here, should be covered by initialise_autocorrelation
     CASE ("vcf")
      normalise=.TRUE.
     CASE ("ecaf")
      normalise=.FALSE.
     CASE DEFAULT
      CALL report_error(0)
     END SELECT
    END SUBROUTINE initialise_velocity_correlation

    SUBROUTINE finalise_velocity_correlation()
    IMPLICIT NONE
    INTEGER :: deallocstatus
     DEALLOCATE(correlation_function,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
    END SUBROUTINE finalise_velocity_correlation

    !This subroutine computes <vi,vj>, and relevant sums thereof, as given in vc_components.
    !'sampling' is the sampling interval.
    SUBROUTINE compute_velocity_correlation(sampling)
    IMPLICIT NONE
    !$ INTERFACE
    !$  FUNCTION OMP_get_num_threads()
    !$  INTEGER :: OMP_get_num_threads
    !$  END FUNCTION OMP_get_num_threads
    !$ END INTERFACE
    INTEGER,INTENT(IN) :: sampling
    INTEGER :: startstep,timeline,component_counter,nmolecules_ref
    INTEGER :: allocstatus,deallocstatus
    REAL(WORKING_PRECISION),ALLOCATABLE :: temp_function(:)
     !Compute velocity correlation functions
     IF (VERBOSE_OUTPUT) WRITE(*,*) "Computing velocity correlation functions"
     correlation_function(:,:)=0.0d0
     !$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
     !$OMP PRIVATE(temp_function,component_counter,nmolecules_ref)
     !$OMP SINGLE
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (velocity correlation function)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     !allocate memory for temporary functions (used for parallelisation)
     ALLOCATE(temp_function(tmax+1),STAT=allocstatus)
     IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
     !Here, the outer loop is chosen to be the starting step, i.e. 't0'
     !$OMP DO SCHEDULE(STATIC,1)
     DO startstep=1,nsteps,sampling
      DO component_counter=1,ncomponents,1
       nmolecules_ref=give_number_of_molecules_per_step(vc_components(component_counter)%reference_type)
       IF (vc_components(component_counter)%self) THEN
        !self contributions
        CALL iterate_timesteps_velocity_self_correlation&
        &(component_counter,startstep,temp_function,nmolecules_ref)
       ELSE
        !distinct contributions...
        IF (vc_components(component_counter)%reference_type&
        &==vc_components(component_counter)%observed_type) THEN
         CALL iterate_timesteps_velocity_distinct_correlation_sameref&
         &(component_counter,startstep,temp_function,nmolecules_ref)
        ELSE
         CALL iterate_timesteps_velocity_distinct_correlation_differentref&
         &(component_counter,startstep,temp_function,nmolecules_ref)
        ENDIF
       ENDIF
       !CRITICAL directive to properly update the correlation_function
       !ATOMIC could also work, not sure if they allow arrays
       !$OMP CRITICAL
       correlation_function(:,component_counter)=correlation_function(:,component_counter)+temp_function(:)
       !$OMP END CRITICAL
      ENDDO
     ENDDO
     !$OMP END DO
     !deallocate private memory used for parallelisation
     DEALLOCATE(temp_function,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
     !normalise by the number of startsteps / averages taken.
     DO timeline=0,tmax,1
      correlation_function(timeline+1,:)=correlation_function(timeline+1,:)/DFLOAT(1+(nsteps-(timeline+1))/sampling)
     ENDDO
    END SUBROUTINE compute_velocity_correlation

    !Calculates <vs(t)*vs(0)> in equation A6/A7 in J. Phys. Chem. B, 2011, 115, 13212–13221.
    !Does NOT include the factor 1/3, and does NOT include the normalisation by number of molecules!!
    SUBROUTINE iterate_timesteps_velocity_self_correlation&
    &(vc_component_number,startstep,temp_function,nmolecules_ref)
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: temp_value
    INTEGER :: molecule_counter,local_tmax
    INTEGER :: molecule_type_index
    INTEGER :: timeline
    INTEGER,INTENT(IN) :: startstep,vc_component_number,nmolecules_ref
    REAL(WORKING_PRECISION),INTENT(INOUT) :: temp_function(:)
    REAL(WORKING_PRECISION) :: initial_velocities(nmolecules_ref,3),norm_factor
     !first, initialise
     norm_factor=FLOAT(nmolecules_ref)
     temp_function(:)=0.0d0
     molecule_type_index=vc_components(vc_component_number)%reference_type
     DO molecule_counter=1,nmolecules_ref,1
      initial_velocities(molecule_counter,:)=give_center_of_mass(startstep,molecule_type_index,molecule_counter)
     ENDDO
     !Careful: startstep+timeline can of course not exceed the number of available steps.
     IF ((startstep+tmax)>nsteps) THEN!still ok this way round, because the average positions should be allocated in the calling routine anyway.
      local_tmax=(nsteps-startstep)
     ELSE
      local_tmax=tmax
     ENDIF
     DO timeline=0,local_tmax,1
      temp_value=0.0d0
      DO molecule_counter=1,nmolecules_ref,1
       temp_value=temp_value+DOT_PRODUCT(&
       &give_center_of_mass(startstep+timeline,molecule_type_index,molecule_counter),&
       &initial_velocities(molecule_counter,:))
      ENDDO
      IF (normalise) THEN
       temp_function(timeline+1)=temp_function(timeline+1)+(temp_value/norm_factor)
      ELSE
       temp_function(timeline+1)=temp_function(timeline+1)+temp_value
      ENDIF
     ENDDO
    END SUBROUTINE iterate_timesteps_velocity_self_correlation
    
    !Calculates <vd(t)*vd(0)> in equation A8/A9 in J. Phys. Chem. B, 2011, 115, 13212–13221.
    !Does NOT include the factor 1/3, and does NOT include the normalisation by number of molecules!!
    SUBROUTINE iterate_timesteps_velocity_distinct_correlation_sameref&
    &(vc_component_number,startstep,temp_function,nmolecules_ref)
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: temp_value
    INTEGER :: molecule_counter_ref,molecule_counter_obs,local_tmax
    INTEGER :: molecule_type_index
    INTEGER :: timeline
    INTEGER,INTENT(IN) :: startstep,vc_component_number,nmolecules_ref
    REAL(WORKING_PRECISION),INTENT(INOUT) :: temp_function(:)
    REAL(WORKING_PRECISION) :: initial_velocities(nmolecules_ref,3),norm_factor,temp_vector(3)
     !first, initialise
     norm_factor=FLOAT(nmolecules_ref*(nmolecules_ref-1))
     temp_function(:)=0.0d0
     molecule_type_index=vc_components(vc_component_number)%reference_type
     DO molecule_counter_ref=1,nmolecules_ref,1
      initial_velocities(molecule_counter_ref,:)=&
      &give_center_of_mass(startstep,molecule_type_index,molecule_counter_ref)
     ENDDO
     !Careful: startstep+timeline can of course not exceed the number of available steps.
     IF ((startstep+tmax)>nsteps) THEN!still ok this way round, because the average positions should be allocated in the calling routine anyway.
      local_tmax=(nsteps-startstep)
     ELSE
      local_tmax=tmax
     ENDIF
     DO timeline=0,local_tmax,1
      temp_value=0.0d0
      DO molecule_counter_ref=1,nmolecules_ref,1
       temp_vector(:)=0.0d0
       DO molecule_counter_obs=1,nmolecules_ref,1
        !it can happen that the two molecules are identical.
        !in that case, the self contribution needs to be removed.
        IF (molecule_counter_ref==molecule_counter_obs) CYCLE
        temp_vector(:)=temp_vector(:)+give_center_of_mass(startstep+timeline,molecule_type_index,molecule_counter_obs)
       ENDDO
       temp_value=temp_value+DOT_PRODUCT(temp_vector(:),initial_velocities(molecule_counter_ref,:))
      ENDDO
      IF (normalise) THEN
       temp_function(timeline+1)=temp_function(timeline+1)+(temp_value/norm_factor)
      ELSE
       temp_function(timeline+1)=temp_function(timeline+1)+temp_value
      ENDIF
     ENDDO
    END SUBROUTINE iterate_timesteps_velocity_distinct_correlation_sameref

    !Calculates <vd(t)*vd(0)> in equation A8/A9/A10 in J. Phys. Chem. B, 2011, 115, 13212–13221.
    !Does NOT include the factor 1/3, and does NOT include the normalisation by number of molecules!!
    SUBROUTINE iterate_timesteps_velocity_distinct_correlation_differentref&
    &(vc_component_number,startstep,temp_function,nmolecules_ref)
    IMPLICIT NONE
    REAL(KIND=WORKING_PRECISION) :: temp_value
    INTEGER :: molecule_counter_ref,molecule_counter_obs,local_tmax
    INTEGER :: molecule_type_ref,molecule_type_obs
    INTEGER :: timeline,nmolecules_obs
    INTEGER,INTENT(IN) :: startstep,vc_component_number,nmolecules_ref
    REAL(WORKING_PRECISION),INTENT(INOUT) :: temp_function(:)
    REAL(WORKING_PRECISION) :: initial_velocities(nmolecules_ref,3),norm_factor
     !first, initialise
     temp_function(:)=0.0d0
     molecule_type_ref=vc_components(vc_component_number)%reference_type
     molecule_type_obs=vc_components(vc_component_number)%observed_type
     nmolecules_obs=give_number_of_molecules_per_step(molecule_type_obs)
     norm_factor=FLOAT(nmolecules_ref*nmolecules_obs)
     DO molecule_counter_ref=1,nmolecules_ref,1
      initial_velocities(molecule_counter_ref,:)=&
      &give_center_of_mass(startstep,molecule_type_ref,molecule_counter_ref)
     ENDDO
     !Careful: startstep+timeline can of course not exceed the number of available steps.
     IF ((startstep+tmax)>nsteps) THEN!still ok this way round, because the average positions should be allocated in the calling routine anyway.
      local_tmax=(nsteps-startstep)
     ELSE
      local_tmax=tmax
     ENDIF
     DO timeline=0,local_tmax,1
      temp_value=0.0d0
      DO molecule_counter_ref=1,nmolecules_ref,1
       DO molecule_counter_obs=1,nmolecules_obs,1
        temp_value=temp_value+DOT_PRODUCT(&
        &give_center_of_mass(startstep+timeline,molecule_type_obs,molecule_counter_obs),&
        &initial_velocities(molecule_counter_ref,:))
       ENDDO
      ENDDO
      IF (normalise) THEN
       temp_function(timeline+1)=temp_function(timeline+1)+(temp_value/norm_factor)
      ELSE
       temp_function(timeline+1)=temp_function(timeline+1)+temp_value
      ENDIF
     ENDDO
    END SUBROUTINE iterate_timesteps_velocity_distinct_correlation_differentref

  END SUBROUTINE velocity_correlation

  SUBROUTINE check_boost()
  IMPLICIT NONE
   IF (.NOT.(VERBOSE_OUTPUT)) RETURN
   IF (.NOT.((give_comboost()).AND.(PARALLEL_OPERATION))) THEN
    PRINT *,"These calculations are very (like, really) involved."
   ENDIF
   IF (.NOT.(give_comboost())) THEN
    PRINT *,"They speed up significantly when performed on a centre-of-mass trajectory."
    PRINT *,"(See also the keyword 'convert_simple')"
   ENDIF
   IF (.NOT.(PARALLEL_OPERATION)) THEN
    PRINT *,"It is *NOT* advisable to use serial operation."
   ENDIF
  END SUBROUTINE check_boost

  !WRITING input file to unit 8, which shouldn't be open.
  !has to be compliant with 'initialise_autocorrelation'
  SUBROUTINE write_simple_conductivity()
  IMPLICIT NONE
  LOGICAL :: file_exists,connected
  INTEGER :: n,ios,tstep_local
   FILENAME_AUTOCORRELATION_INPUT="prealpha_simple.inp"
   INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_AUTOCORRELATION_INPUT),EXIST=file_exists)
   IF (file_exists) CALL report_error(114)
   INQUIRE(UNIT=8,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=8)
   OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_AUTOCORRELATION_INPUT),IOSTAT=ios)!input path is added for the MSD file!
   IF (ios/=0) CALL report_error(46,exit_status=ios)
   tstep_local=(give_number_of_timesteps()-1)/10000
   IF (tstep_local<5) tstep_local=1
   WRITE(8,'("conductivity")')
   WRITE(8,'("tmax ",I0)') give_number_of_timesteps()-1
   WRITE(8,'("sampling_interval ",I0)') tstep_local
   WRITE(8,'("quit")')
   CLOSE(UNIT=8)
   IF (VERBOSE_OUTPUT) PRINT *,"File 'prealpha_simple.inp' written."
  END SUBROUTINE write_simple_conductivity

  SUBROUTINE perform_autocorrelation()
  IMPLICIT NONE
  CALL initialise_autocorrelation()
  IF ((ERROR_CODE/=21).AND.(ERROR_CODE/=39).AND.(ERROR_CODE/=33).AND.(ERROR_CODE/=83)&
  &.AND.(ERROR_CODE/=118).AND.(ERROR_CODE/=121)) THEN
   !do the actual analysis:
   SELECT CASE (TRIM(operation_mode))!no further output necessary here, should be covered by initialise_autocorrelation
   CASE ("conductivity")
    CALL check_boost()
    CALL overall_conductivity()
   CASE ("vcf","ecaf")
    CALL check_boost()
    CALL velocity_correlation()
   CASE ("dihedral")
    IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
    CALL dihedral_autocorrelation()!fill the array
    IF (.NOT.(skip_autocorr)) CALL calculate_autocorrelation_function_from_binary_array() !calculate the autocorrelation function from the array
   CASE ("reorientation")
    IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
    CALL reorientational_autocorrelation()
   CASE ("rmm-vcf")
    IF (INFORMATION_IN_TRAJECTORY=="POS") CALL report_error(56)
    CALL check_boost()
    IF (WRAP_TRAJECTORY) THEN
     CALL report_error(72)
    ELSE
     CALL cross_correlation()
    ENDIF
   CASE DEFAULT
    CALL report_error(0)
   END SELECT
   CALL finalise_autocorrelation()
  ELSE
   ERROR_CODE=ERROR_CODE_DEFAULT
   !resetting ERROR_CODE to avoid problems with a rare condition
   !(i.e. invoking perform_autocorrelation after error 21 or 39 or 33 without intermediate problems)
  ENDIF
  END SUBROUTINE perform_autocorrelation

END MODULE AUTOCORRELATION
!--------------------------------------------------------------------------------------------------------------------------------!

!This Module calculates (drift corrected) mean squared displacements.
!The diffusion implementation is shit. Use TRAVIS.
MODULE DIFFUSION ! Copyright (C) 2020 Frederik Philippi
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
 INTEGER,DIMENSION(:,:),ALLOCATABLE :: projections ! x y z molecule_type_index
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
 PUBLIC :: perform_diffusion_analysis,print_helpful_info,user_msd_input,write_simple_diffusion

 CONTAINS

  !WRITING input file to unit 8, which shouldn't be open.
  !has to be compliant with 'read_input_for_distribution'
  SUBROUTINE write_simple_diffusion()
  IMPLICIT NONE
  LOGICAL :: file_exists,connected
  INTEGER :: n,ios,tstep_local
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
 ! INTEGER,DIMENSION(:,:),ALLOCATABLE :: projections ! x y z molecule_type_index, just as in module DIFFUSION
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
   WRITE(*,'(" the time scaling factor is ",I0)') TIME_SCALING_FACTOR
   WRITE(*,'(" averages taken: max = ",I0,", min = ",I0)') (give_number_of_timesteps()-tstep),(give_number_of_timesteps()-tmax)
  END SUBROUTINE print_helpful_info

  !This subroutine generates the required diffusion functions.
  SUBROUTINE make_diffusion_functions(projection_number)
  IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
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
   !$  WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (mean squared displacement)"
   !$  CALL timing_parallel_sections(.TRUE.)
   !$ ENDIF
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
   !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
   !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
   !$  CALL timing_parallel_sections(.FALSE.)
   !$ ENDIF
   !both x_squared and x_unsquared are normalised by the number of molecules. now, account for the averaging process:
   x_squared(:)=x_squared(:)/DFLOAT(x_num(:))
   DO n=1,3,1
    x_unsquared(:,n)=x_unsquared(:,n)/DFLOAT(x_num(:))
   ENDDO
  END SUBROUTINE make_diffusion_functions

  !This subroutine is responsible for writing the output into a file.
  SUBROUTINE write_diffusion_functions(projection_number)
  IMPLICIT NONE
  !Output formats
  2 FORMAT (I15,2E15.6) !Module DIFFUSION - diffusion output file, no verbose_print
  6 FORMAT (I15,6E15.6,I15) !Module DIFFUSION - diffusion output file, verbose_print
  INTEGER,INTENT(IN) :: projection_number
  INTEGER :: array_pos,ios,i
  LOGICAL :: connected
  CHARACTER(LEN=1024) :: filename_diffusion_output
  CHARACTER(LEN=32) :: msd_exponent_str
   WRITE(msd_exponent_str,'(I0)') msd_exponent
   INQUIRE(UNIT=3,OPENED=connected)
   IF (connected) CALL report_error(27,exit_status=3)
   !generate filename, depending on operation mode
   IF (TRIM(operation_mode)=="default") THEN
    WRITE(filename_diffusion_output,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
    &//"type_",projection_number,"_diffusion"
   ELSE
    WRITE(filename_diffusion_output,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
    &//"projection_",projection_number,"_diffusion"
   ENDIF
   IF (VERBOSE_OUTPUT) WRITE(*,*) "writing diffusion data (msd) into file: ",TRIM(filename_diffusion_output)
   OPEN(UNIT=3,FILE=TRIM(filename_diffusion_output),IOSTAT=ios)
   IF (ios/=0) CALL report_error(26,exit_status=ios)
   !Write header
   WRITE(3,*) "This file contains diffusion information based on the input file '"&
   &,TRIM(FILENAME_DIFFUSION_INPUT),"'"
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
   DO array_pos=1,tmax/tstep,1 !zeroth entry is already there.
    IF (verbose_print) THEN
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
!This Module performs low-level command line tasks.
MODULE RECOGNITION ! Copyright (C) 2020 Frederik Philippi
 USE SETTINGS
 IMPLICIT NONE
 PRIVATE
 !variables
 INTEGER :: safety_shift_default=200 !how many atoms to check in the molecule recognition. Ideally the maximum number of atoms in a molecule.
 !PRIVATE/PUBLIC declarations
 PUBLIC :: molecule_recognition
 CONTAINS

  SUBROUTINE molecule_recognition(trajectory_command_line)! This subroutine builds the rotation matrix. Not accessible globally.
  IMPLICIT NONE
   !$ INTERFACE
   !$  FUNCTION OMP_get_max_threads()
   !$  INTEGER :: OMP_get_max_threads
   !$  END FUNCTION OMP_get_max_threads
   !$  SUBROUTINE OMP_set_num_threads(number_of_threads)
   !$  INTEGER,INTENT(IN) :: number_of_threads
   !$  END SUBROUTINE OMP_set_num_threads
   !$ END INTERFACE
  REAL :: box_dimensions(2,3),box_size(3)!low and high for x,y,z and difference between high and low
  REAL :: maximum_distance_squared
  CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE :: list_of_elements !--> Turned on support for  2-letter elements!
  CHARACTER(LEN=1) :: drude_symbol="0"
  REAL,DIMENSION(:,:),ALLOCATABLE :: coordinates !first dimension nlines_total, second dimension the three coordinates
  TYPE :: single_molecule
   CHARACTER(LEN=1024) :: sum_formula
   INTEGER :: number_of_atoms=0 ! = number of atoms PER SINGLE MOLECULE, not total!
   INTEGER :: total_molecule_count=0 ! = number of molecules of this type in the box
   LOGICAL :: ignore=.FALSE. ! = TRUE, when this one has been identified as duplicate, and merged to the previous one.
   LOGICAL,DIMENSION(:),ALLOCATABLE :: member(:) ! = collects all the members belonging to this molecule group
  END TYPE single_molecule
  LOGICAL,DIMENSION(:),ALLOCATABLE :: atom_assigned! = collects all the atoms that have been assigned so far
  TYPE(single_molecule),DIMENSION(:),ALLOCATABLE :: molecule_list(:) !list of molecules. There is a maximum of (nlines_total) different molecule types, each of them having a number_of_atoms and total_molecule_count
  CHARACTER(LEN=1024) :: trajectory_command_line
  LOGICAL :: file_exists,connected
  LOGICAL(KIND=1),DIMENSION(:,:),ALLOCATABLE :: connectivity(:,:) !I know, it's still a 16-fold waste of RAM, which is why it's only in the developers version.
  INTEGER :: ios,lines,nlines_total,molecule_types,number_of_drude_particles,threadnum,i
    !$ IF (DEVELOPERS_VERSION) THEN
    !$ threadnum=OMP_get_max_threads()
    !$ WRITE(*,'("  ! number of threads set to ",I0)') threadnum
    !$ CALL OMP_set_num_threads(threadnum)
    !$ ENDIF
   PRINT *,"Trying to perform molecule recognition on trajectory file:"
   trajectory_command_line=ADJUSTL(trajectory_command_line)
   WRITE(*,'(A,A,A)') ' "',TRIM(trajectory_command_line),'"'
   PRINT *,"Expecting a sorted, unwrapped lammps trajectory with cartesian coordinates!"
   PRINT *,"(Specify 'element xu yu zu' and 'sort ID' in lammps)"
   INQUIRE(FILE=TRIM(trajectory_command_line),EXIST=file_exists)
   IF (file_exists) THEN
    INQUIRE(UNIT=3,OPENED=connected)
    IF (connected) CALL report_error(27,exit_status=3)
    OPEN(UNIT=3,FILE=TRIM(trajectory_command_line),ACTION='READ',IOSTAT=ios)
    IF (ios/=0) THEN
     CALL report_error(95,exit_status=ios)
     RETURN
    ENDIF
    !trajectory file has been specified correctly. Start reading it.
    REWIND 3
    DO lines=1,9,1
     SELECT CASE (lines)
     CASE (4)
      READ(3,IOSTAT=ios,FMT=*) nlines_total
     CASE (6)
      READ(3,IOSTAT=ios,FMT=*) box_dimensions(:,1)
     CASE (7)
      READ(3,IOSTAT=ios,FMT=*) box_dimensions(:,2)
     CASE (8)
      READ(3,IOSTAT=ios,FMT=*) box_dimensions(:,3)
     CASE DEFAULT
      READ(3,IOSTAT=ios,FMT=*)
     END SELECT
     IF (ios/=0) THEN
      CALL report_error(95)
      RETURN
     ENDIF
    ENDDO
    !initialise box size
    box_size(:)=box_dimensions(2,:)-box_dimensions(1,:)
    maximum_distance_squared=box_size(2)**2+SQRT(box_size(1)**2+box_size(3)**2)
    WRITE(*,ADVANCE="NO",FMT='(" Expecting ",I0," atoms.")') nlines_total
    CALL initialise_molecule_recognition()
    IF (number_of_drude_particles/=0) &
    &WRITE(*,'(" Found ",I0," drude particles - be aware of correct sorting!")') number_of_drude_particles
    IF (ERROR_CODE==95) RETURN
    CALL recognise_molecules()
    CALL finalise_molecule_recognition()
   ELSE
    CALL report_error(95)
   ENDIF
   
   CONTAINS

    SUBROUTINE initialise_molecule_recognition()
    IMPLICIT NONE
    INTEGER :: allocstatus
     PRINT *,"Initialising."
     ALLOCATE(list_of_elements(nlines_total),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(95,exit_status=allocstatus)
      RETURN
     ENDIF
     ALLOCATE(atom_assigned(nlines_total),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(95,exit_status=allocstatus)
      RETURN
     ENDIF
     atom_assigned(:)=.FALSE.
     ALLOCATE(molecule_list(nlines_total),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(95,exit_status=allocstatus)
      RETURN
     ENDIF
     ALLOCATE(coordinates(nlines_total,3),STAT=allocstatus)
     IF (allocstatus/=0) THEN
      CALL report_error(95,exit_status=allocstatus)
      RETURN
     ENDIF
     IF (DEVELOPERS_VERSION) THEN
      ALLOCATE(connectivity(nlines_total,nlines_total),STAT=allocstatus)
      IF (allocstatus/=0) THEN
       CALL report_error(95,exit_status=allocstatus)
       RETURN
      ENDIF
     ENDIF
     number_of_drude_particles=0
     DO lines=1,nlines_total,1
      ALLOCATE(molecule_list(lines)%member(nlines_total),STAT=allocstatus)
      IF (allocstatus/=0) THEN
       CALL report_error(95,exit_status=allocstatus)
       CLOSE(UNIT=3)
       RETURN
      ENDIF
      molecule_list(lines)%member(:)=.FALSE.
      READ(3,IOSTAT=ios,FMT=*) list_of_elements(lines),coordinates(lines,:)
      CALL wrap_vector(coordinates(lines,:))
      IF ((ADJUSTL(TRIM(list_of_elements(lines)))=="X").OR.(ADJUSTL(TRIM(list_of_elements(lines)))=="D")) THEN
       !a drude particle has been found
       number_of_drude_particles=number_of_drude_particles+1
       IF (ADJUSTL(TRIM(list_of_elements(lines)))/=drude_symbol) THEN
        !drude symbol doesn't match - yet?
        IF ((drude_symbol=="X").OR.(drude_symbol=="D")) THEN
         !has already been initialised - Both X and D present!
         drude_symbol="B"
        ELSE
         drude_symbol=ADJUSTL(TRIM(list_of_elements(lines)))
        ENDIF
       ENDIF
      ENDIF
      IF (ios/=0) THEN
       CALL report_error(95)
       CLOSE(UNIT=3)
       RETURN
      ENDIF
     ENDDO
     !UNIT 3 no longer required!
     CLOSE(UNIT=3)
    END SUBROUTINE initialise_molecule_recognition

    !This FUNCTION returns the smallest squared distance of 2 atoms considering all PBCs.
    REAL FUNCTION give_smallest_atom_distance_squared(pos_1,pos_2)
    IMPLICIT NONE
    INTEGER :: a,b,c
    REAL :: pos_1(3),pos_2(3),shift(3),distance_clip
     !two atoms can be no further apart than the diagonale of the box... that's what I initialise to
     give_smallest_atom_distance_squared=maximum_distance_squared
     !Now, check all mirror images
     DO a=-1,1,1! a takes the values (-1, 0, 1)
      DO b=-1,1,1! b takes the values (-1, 0, 1)
       DO c=-1,1,1! c takes the values (-1, 0, 1)
        shift(1)=FLOAT(a)
        shift(2)=FLOAT(b)
        shift(3)=FLOAT(c)
        shift=shift*box_size
        !shift is now the translation vector to the mirror image.
        distance_clip=SUM(((pos_2+shift)-pos_1)**2)
        IF (distance_clip<give_smallest_atom_distance_squared) THEN
         !a distance has been found that's closer than the current best - amend that.
         give_smallest_atom_distance_squared=distance_clip
        ENDIF
         ENDDO
      ENDDO
     ENDDO
    END FUNCTION give_smallest_atom_distance_squared

    SUBROUTINE wrap_vector(input_vector)
    IMPLICIT NONE
    REAL,INTENT(INOUT) :: input_vector(3)
    INTEGER :: xyzcounter
     DO xyzcounter=1,3,1
      DO
       IF (input_vector(xyzcounter)<box_dimensions(1,xyzcounter)) THEN
        !input vector is outside of box (smaller)
        input_vector(xyzcounter)=input_vector(xyzcounter)+box_size(xyzcounter)
        CYCLE
       ELSE
        IF (input_vector(xyzcounter)>box_dimensions(2,xyzcounter)) THEN
         !input vector is outside of box (bigger)
         input_vector(xyzcounter)=input_vector(xyzcounter)-box_size(xyzcounter)
         CYCLE
        ELSE
         !input vector is inside box!
         EXIT
        ENDIF
       ENDIF
      ENDDO
     ENDDO
    END SUBROUTINE wrap_vector

    !The following set of functions provides the values of important variables to other routines. This serves the purpose of keeping variables local.
    SUBROUTINE assign_sum_formula(molecule_type_index,first_element)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index,first_element
    INTEGER :: outer,inner,n
    LOGICAL :: element_unused(molecule_list(molecule_type_index)%number_of_atoms)!has this atom been used up yet?
    CHARACTER(LEN=2) :: current_element
     element_unused(:)=.TRUE.
     molecule_list(molecule_type_index)%sum_formula=""
     DO outer=1,molecule_list(molecule_type_index)%number_of_atoms,1
      current_element=TRIM(list_of_elements(first_element+outer-1))
      !Print the element in outer, if not yet done:
      IF (element_unused(outer)) THEN
       !append the new element
       molecule_list(molecule_type_index)%sum_formula=&
       &TRIM(molecule_list(molecule_type_index)%sum_formula)//TRIM(current_element)
       !count how many are there, and label them as used
       n=1
       DO inner=(outer+1),molecule_list(molecule_type_index)%number_of_atoms,1
        IF (TRIM(current_element)==TRIM(list_of_elements(first_element+inner-1))) THEN
         element_unused(inner)=.FALSE.
         n=n+1
        ENDIF
       ENDDO
       !append the number
       IF (n>1) THEN
        WRITE(molecule_list(molecule_type_index)%sum_formula,'(A,I0)')&
        &TRIM(molecule_list(molecule_type_index)%sum_formula),n
       ENDIF
      ENDIF
     ENDDO
    END SUBROUTINE assign_sum_formula

    SUBROUTINE recognise_molecules()
    USE DEBUG
    IMPLICIT NONE
    !$ INTERFACE
    !$  FUNCTION OMP_get_num_threads()
    !$  INTEGER :: OMP_get_num_threads
    !$  END FUNCTION OMP_get_num_threads
    !$ END INTERFACE
    INTEGER :: molecule_type_counter,linecounter1,linecounter2,linecounter3,merged_molecule_types,newest_molecule_type
    INTEGER :: written_atoms
    LOGICAL :: file_exists,write_xyz_files
    CHARACTER(LEN=1024) :: working_directory,xyz_filename,trajectory_filename_only
     WRITE(*,'(" Starting cutoff-based molecule recognition.")')
     IF (DEVELOPERS_VERSION) THEN
      CALL initialise_connectivity()
      !$OMP SINGLE
      !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
      !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (brute force recognition)"
      !$  CALL timing_parallel_sections(.TRUE.)
      !$ ENDIF
      !$OMP END SINGLE
     ENDIF
     molecule_types=0
     DO linecounter1=1,nlines_total,1
      !look for an atom which is *not* part of a molecule yet!
      IF (.NOT.(atom_assigned(linecounter1))) THEN
       molecule_types=molecule_types+1
       molecule_list(molecule_types)%total_molecule_count=1
       molecule_list(molecule_types)%member(linecounter1)=.TRUE.
       atom_assigned(linecounter1)=.TRUE.
       !for each atom, get all those that belong to the same molecule
       IF (DEVELOPERS_VERSION) THEN
        CALL assign_atoms_to_molecules_parallelised(molecule_types,linecounter1)
       ELSE
        CALL assign_atoms_to_molecules(molecule_types,linecounter1,linecounter1)
       ENDIF
      ENDIF
     ENDDO
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION).AND.(DEVELOPERS_VERSION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='("     ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
     WRITE(*,'(" Found ",I0," molecules. Checking for consistency.")') molecule_types
     DO molecule_type_counter=1,molecule_types,1
      DO linecounter1=1,nlines_total,1
       IF (molecule_list(molecule_type_counter)%member(linecounter1)) EXIT
      ENDDO
      molecule_list(molecule_type_counter)%number_of_atoms=1
      DO linecounter2=linecounter1+1,nlines_total,1
       IF (.NOT.(molecule_list(molecule_type_counter)%member(linecounter2))) EXIT
       molecule_list(molecule_type_counter)%number_of_atoms=molecule_list(molecule_type_counter)%number_of_atoms+1
      ENDDO
      IF (.NOT.(molecule_list(molecule_type_counter)%member(linecounter2))) linecounter2=linecounter2-1
      DO linecounter3=linecounter2+1,nlines_total,1
       IF (molecule_list(molecule_type_counter)%member(linecounter3)) THEN
        CALL report_error(96,exit_status=molecule_type_counter)
        CALL report_error(95)
        RETURN
       ENDIF
      ENDDO
      CALL assign_sum_formula(molecule_type_counter,linecounter1)
     ENDDO
     WRITE(*,'(" Organising molecules into types based on sum formula and order of trajectory file.")')
     merged_molecule_types=1
     newest_molecule_type=1
     molecule_list(1)%ignore=.FALSE.
     DO molecule_type_counter=2,molecule_types,1
      IF (molecule_list(molecule_type_counter)%sum_formula==molecule_list(newest_molecule_type)%sum_formula) THEN
       !merge both!
       molecule_list(newest_molecule_type)%total_molecule_count=molecule_list(newest_molecule_type)%total_molecule_count+1
       molecule_list(molecule_type_counter)%ignore=.TRUE.
      ELSE
       !new molecule type found!
       molecule_list(molecule_type_counter)%ignore=.FALSE.
       newest_molecule_type=molecule_type_counter
       merged_molecule_types=merged_molecule_types+1
      ENDIF
     ENDDO
     WRITE(*,'(" Merged molecules into ",I0," types.")') merged_molecule_types
     working_directory=""
     trajectory_filename_only=TRIM(trajectory_command_line)
     DO i=LEN(TRIM(trajectory_command_line))-1,1,-1
      IF (IACHAR("/")==IACHAR(trajectory_command_line(i:i)))THEN
       working_directory=TRIM(trajectory_command_line(1:i))
       trajectory_filename_only=TRIM(trajectory_command_line(i+1:LEN(TRIM(trajectory_command_line))))
       WRITE(*,'(A,A,A)') ' Write into directory "',TRIM(working_directory),'"'
       EXIT
      ENDIF
      IF (i==1) working_directory=""
     ENDDO
     write_xyz_files=.TRUE.
     DO molecule_type_counter=1,merged_molecule_types,1
      WRITE(xyz_filename,'(A,I0,A)') TRIM(working_directory)//"MolRec_Type_",molecule_type_counter,".xyz"
      INQUIRE(FILE=TRIM(xyz_filename),EXIST=file_exists)
      IF (file_exists) THEN
       write_xyz_files=.FALSE.
       EXIT
      ENDIF
     ENDDO
     IF (write_xyz_files) THEN
      newest_molecule_type=1
      PRINT *,"Writing example files 'MolRec_Type_N.xyz' for each molecule type."
      DO molecule_type_counter=1,molecule_types,1
       IF (.NOT.(molecule_list(molecule_type_counter)%ignore)) THEN
        WRITE(xyz_filename,'(A,I0,A)') TRIM(working_directory)//"MolRec_Type_",newest_molecule_type,".xyz"
        IF (connected) CALL report_error(27,exit_status=3)
        OPEN(UNIT=3,FILE=TRIM(xyz_filename))
        !write temporary xyz file in SCRATCH unit
        REWIND 3
        WRITE(3,*) molecule_list(molecule_type_counter)%number_of_atoms
        WRITE(3,*)
        written_atoms=0
        DO linecounter1=1,nlines_total,1
         IF (molecule_list(molecule_type_counter)%member(linecounter1)) THEN
          written_atoms=written_atoms+1
          WRITE(3,*) list_of_elements(linecounter1),coordinates(linecounter1,:)
          IF (written_atoms==molecule_list(molecule_type_counter)%number_of_atoms) EXIT
         ENDIF
        ENDDO
        CALL center_xyz(3,.TRUE.,custom_header=&
        &"Example molecule '"//TRIM(molecule_list(molecule_type_counter)%sum_formula)//"'"&
        &,geometric_center=.TRUE.)
        CLOSE(UNIT=3)
        newest_molecule_type=newest_molecule_type+1
       ENDIF
      ENDDO
     ELSE
      PRINT *,"A file of the type 'MolRec_Type_N.xyz' already exists - no structures will be written."
     ENDIF
     !here begins the part that outputs a molecular input file.
     INQUIRE(FILE=TRIM(working_directory)//TRIM(FILENAME_MOLECULAR_INPUT),EXIST=file_exists)
     IF (file_exists) THEN
      PRINT *,"Molecular input file with name '"//TRIM(FILENAME_MOLECULAR_INPUT)//"' already exists."
      PRINT *,"Please use the following lines until the 'quit' statement as your molecular input file:"
      WRITE(*,*) "  1 ### number of timesteps in trajectory - please adjust!"
      WRITE(*,'("   ",I0," ### number of different types of molecules.")') merged_molecule_types
      DO molecule_type_counter=1,molecule_types,1
       IF (.NOT.(molecule_list(molecule_type_counter)%ignore)) THEN
        WRITE(*,'("   0 ",I0," ",I0," ### ",I0,A,I0," atoms each.")')&
        &molecule_list(molecule_type_counter)%number_of_atoms,molecule_list(molecule_type_counter)%total_molecule_count,&
        &molecule_list(molecule_type_counter)%total_molecule_count,&
        &" molecules with the formula '"//TRIM(molecule_list(molecule_type_counter)%sum_formula)//"' per step with ",&
        &molecule_list(molecule_type_counter)%number_of_atoms
       ENDIF
      ENDDO
      IF (number_of_drude_particles/=0) THEN
       IF (drude_symbol=="B") THEN
        WRITE(*,*) "  masses 2 ### Both 'D' and 'X' have been found..."
        WRITE(*,*) "  X  0.400 ### By that, the support for drude particles is turned on."
        WRITE(*,*) "  D  0.400 ### By that, the support for drude particles is turned on."
       ELSE
        WRITE(*,*) "  masses 1 ### The following line specifies a custom mass."
        WRITE(*,*) "  "//drude_symbol//"  0.400 ### By that, the support for drude particles is turned on."
       ENDIF
      ENDIF
      WRITE(*,*) "  quit"
     ELSE
      !write molecular input file
      PRINT *,"Writing molecular input file '"//TRIM(FILENAME_MOLECULAR_INPUT)//"'."
      INQUIRE(UNIT=8,OPENED=connected)
      IF (connected) CALL report_error(27,exit_status=8)
      OPEN(UNIT=8,FILE=TRIM(working_directory)//TRIM(FILENAME_MOLECULAR_INPUT),IOSTAT=ios)
      IF (ios/=0) THEN
       CALL report_error(95,exit_status=ios)
      ELSE
       WRITE(8,*) "1 ### number of timesteps in trajectory - please adjust!"
       WRITE(8,'(" ",I0," ### number of different types of molecules.")') merged_molecule_types
       DO molecule_type_counter=1,molecule_types,1
        IF (.NOT.(molecule_list(molecule_type_counter)%ignore)) THEN
         WRITE(8,'(" 0 ",I0," ",I0," ### ",I0,A,I0," atoms each.")')&
         &molecule_list(molecule_type_counter)%number_of_atoms,molecule_list(molecule_type_counter)%total_molecule_count,&
         &molecule_list(molecule_type_counter)%total_molecule_count,&
         &" molecules with the formula '"//TRIM(molecule_list(molecule_type_counter)%sum_formula)//"' per step with ",&
         &molecule_list(molecule_type_counter)%number_of_atoms
        ENDIF
       ENDDO
       IF (number_of_drude_particles/=0) THEN
        IF (drude_symbol=="B") THEN
         WRITE(8,*) "  masses 2 ### Both 'D' and 'X' have been found..."
         WRITE(8,*) "  X  0.400 ### By that, the support for drude particles is turned on."
         WRITE(8,*) "  D  0.400 ### By that, the support for drude particles is turned on."
        ELSE
         WRITE(8,*) "  masses 1 ### The following line specifies a custom mass."
         WRITE(8,*) "  "//drude_symbol//"  0.400 ### By that, the support for drude particles is turned on."
        ENDIF
       ENDIF
       WRITE(8,*) "quit"
       CLOSE(UNIT=8)
      ENDIF
     ENDIF
     PRINT *,"Charges and number of timesteps need to be adjusted manually."
     !here begins the part that outputs a general input file
     INQUIRE(FILE=TRIM(working_directory)//TRIM(FILENAME_GENERAL_INPUT),EXIST=file_exists)
     IF (file_exists) THEN
      PRINT *,"General input file with name '"//TRIM(FILENAME_GENERAL_INPUT)//"' already exists."
      PRINT *,"Please use the following lines until the 'quit' statement as your general input file:"
      WRITE(*,*) '  "'//TRIM(trajectory_filename_only)//'" ### trajectory filename'
      WRITE(*,*) '  "./molecular.inp" ### inputfile for module MOLECULAR'
      WRITE(*,*) '  "./" ### path to trajectory'
      WRITE(*,*) '  "./" ### path to other input files'
      WRITE(*,*) '  "./" ### output folder'
      IF (number_of_drude_particles/=0) &
      &WRITE(*,*) '  show_drude'
      WRITE(*,*) '  show_settings'
      WRITE(*,*) '  print_atomic_masses'
      WRITE(*,*) '  dump_example'
      WRITE(*,*) '  sequential_read T'
      WRITE(*,*) '  #parallel_operation T'
      WRITE(*,*) '  #set_threads_simple T'
      WRITE(*,*) '  #trajectory_type lmp'
      WRITE(*,*) '  quit'
     ELSE
      !write molecular input file
      PRINT *,"Writing general input file '"//TRIM(FILENAME_GENERAL_INPUT)//"'."
      INQUIRE(UNIT=8,OPENED=connected)
      IF (connected) CALL report_error(27,exit_status=8)
      OPEN(UNIT=8,FILE=TRIM(working_directory)//TRIM(FILENAME_GENERAL_INPUT),IOSTAT=ios)
      IF (ios/=0) THEN
       CALL report_error(95,exit_status=ios)
      ELSE
       WRITE(8,*) '"'//TRIM(trajectory_filename_only)//'" ### trajectory filename'
       WRITE(8,*) '"./molecular.inp" ### inputfile for module MOLECULAR'
       WRITE(8,*) '"./" ### path to trajectory'
       WRITE(8,*) '"./" ### path to other input files'
       WRITE(8,*) '"./" ### output folder'
       IF (number_of_drude_particles/=0) &
       &WRITE(8,*) 'show_drude'
       WRITE(8,*) 'show_settings'
       WRITE(8,*) 'print_atomic_masses'
       WRITE(8,*) 'dump_example'
       WRITE(8,*) 'sequential_read T'
       WRITE(8,*) '#parallel_operation T'
       WRITE(8,*) '#set_threads_simple T'
       WRITE(8,*) '#trajectory_type lmp'
       WRITE(8,*) 'quit'
       WRITE(8,*)
       WRITE(8,*) ''
       WRITE(8,*) ''
       WRITE(8,*) ''
       WRITE(8,*) 'diffusion_simple'
       WRITE(8,*) 'distribution_simple'
       WRITE(8,*) 'This is the general input file.'
       WRITE(8,*) 'It controls the behaviour of the trajectory analyser.'
       CLOSE(UNIT=8)
      ENDIF
     ENDIF
    END SUBROUTINE recognise_molecules

    !find all atoms connected to each other by close contacts.
    !firstline is the very first atom - everything before that must be included somewhere!
    !currentline contains atom whose neighbours are to be included.
    RECURSIVE SUBROUTINE assign_atoms_to_molecules(molecule_type_index,firstline,currentline)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index,firstline,currentline
    INTEGER :: new_members,linecounter1,maximum
     maximum=firstline+safety_shift_default
     IF (maximum>nlines_total) maximum=nlines_total
     new_members=0
     DO linecounter1=firstline+1,maximum,1
      IF (.NOT.(atom_assigned(linecounter1))) THEN
       IF (give_smallest_atom_distance_squared&
       &(coordinates(linecounter1,:),coordinates(currentline,:))<squared_cutoff(linecounter1,currentline)) THEN
        !The atom in linecounter1 is connected to currentline!
        new_members=new_members+1
        molecule_list(molecule_type_index)%member(linecounter1)=.TRUE.
        atom_assigned(linecounter1)=.TRUE.
        !Advance recursion
        CALL assign_atoms_to_molecules(molecule_type_index,firstline,linecounter1)
       ENDIF
      ENDIF
     ENDDO
    END SUBROUTINE assign_atoms_to_molecules

    REAL FUNCTION squared_cutoff(line1,line2)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: line1,line2
     squared_cutoff=(covalence_radius(list_of_elements(line1))+covalence_radius(list_of_elements(line2)))&
     &*VDW_RATIO_INTERMOLECULAR
    END FUNCTION squared_cutoff

    !find all atoms connected to each other by close contacts.
    !Brute force version for all contacts - parallelised.
    SUBROUTINE initialise_connectivity()
    IMPLICIT NONE
    !$ INTERFACE
    !$  FUNCTION OMP_get_num_threads()
    !$  INTEGER :: OMP_get_num_threads
    !$  END FUNCTION OMP_get_num_threads
    !$ END INTERFACE
    INTEGER :: linecounter1,linecounter2
     !$OMP SINGLE
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (brute force connectivity)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     !$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(linecounter2)
     !$OMP DO
     DO linecounter1=1,nlines_total,1
      DO linecounter2=linecounter1+1,nlines_total,1
       IF (give_smallest_atom_distance_squared&
       &(coordinates(linecounter1,:),coordinates(linecounter2,:))<squared_cutoff(linecounter1,linecounter2)) THEN
        connectivity(linecounter1,linecounter2)=.TRUE.
        connectivity(linecounter2,linecounter1)=.TRUE.
       ELSE
        connectivity(linecounter1,linecounter2)=.FALSE.
        connectivity(linecounter2,linecounter1)=.FALSE.
       ENDIF
      ENDDO
     ENDDO
     !$OMP END DO
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
     !$ ENDIF
    END SUBROUTINE initialise_connectivity

    !find all atoms connected to each other by close contacts.
    !Brute force version.
    RECURSIVE SUBROUTINE assign_atoms_to_molecules_parallelised(molecule_type_index,firstline)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: molecule_type_index,firstline
    INTEGER :: new_members,linecounter1,linecounter2
     new_members=0
     !$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(linecounter2)
     !$OMP DO
     DO linecounter1=firstline,nlines_total,1
      IF (molecule_list(molecule_type_index)%member(linecounter1)) THEN
       !the atom/line belonging to linecounter1 belongs to this molecule. Let's check everything else!
       DO linecounter2=1,nlines_total,1
        IF (.NOT.(atom_assigned(linecounter2))) THEN
         IF (connectivity(linecounter1,linecounter2)) THEN
          !New member found - linecounter2!
          !$OMP CRITICAL
          molecule_list(molecule_type_index)%member(linecounter2)=.TRUE.
          atom_assigned(linecounter2)=.TRUE.
          !$OMP END CRITICAL
          !$OMP ATOMIC
          new_members=new_members+1
         ENDIF
        ENDIF
       ENDDO
      ENDIF
     ENDDO
     !$OMP END DO
     !$OMP END PARALLEL
     IF (new_members/=0) CALL assign_atoms_to_molecules_parallelised(molecule_type_index,firstline)
    END SUBROUTINE assign_atoms_to_molecules_parallelised

    SUBROUTINE finalise_molecule_recognition()
    IMPLICIT NONE
    INTEGER :: deallocstatus
     PRINT *,"finalising molecule recognition."
     IF (DEVELOPERS_VERSION) THEN
      DEALLOCATE(connectivity,STAT=deallocstatus)
      IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     ENDIF
     DO lines=1,nlines_total,1
      DEALLOCATE(molecule_list(lines)%member,STAT=deallocstatus)
      IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     ENDDO
     DEALLOCATE(list_of_elements,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(atom_assigned,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(molecule_list,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     DEALLOCATE(coordinates,STAT=deallocstatus)
     IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
     WRITE(*,*)
    END SUBROUTINE finalise_molecule_recognition

  END SUBROUTINE molecule_recognition

END MODULE RECOGNITION
!--------------------------------------------------------------------------------------------------------------------------------!
!This Module calculates special combined distribution functions - such as cylindrical or polar.
MODULE DISTRIBUTION ! Copyright (C) 2020 Frederik Philippi
 USE SETTINGS
 USE MOLECULAR
 IMPLICIT NONE
 PRIVATE
 !default values
 INTEGER,PARAMETER :: sampling_interval_default=10
 INTEGER,PARAMETER :: bin_count_default=100
 LOGICAL,PARAMETER :: subtract_uniform_default=.FALSE.
 LOGICAL,PARAMETER :: weigh_charge_default=.FALSE.
 REAL,PARAMETER :: maxdist_default=10.0
 !variables
 CHARACTER (LEN=8) :: operation_mode="NONE"!operation mode of the distribution module.
 INTEGER :: number_of_references !number of different references. '0 0 1' uses the z-axis, etc.
 INTEGER :: bin_count_b=bin_count_default!number of bins in the final 2D combined distribution functions (in "a" direction!)
 INTEGER :: bin_count_a=bin_count_default!number of bins in the final 2D combined distribution functions (in "b" direction!)
 INTEGER :: sampling_interval=sampling_interval_default!every this many steps are used for the distribution functions
 REAL :: maxdist=maxdist_default!only molecules with a distance < maxdist are considered. Should be half the box size.
 REAL :: foolsproof_ratio!for the cdf, consider rare sampling of corners.
 LOGICAL :: subtract_uniform=subtract_uniform_default!subtract the uniform density
 LOGICAL :: weigh_charge=weigh_charge_default!weigh the distribution functions by charges.
 INTEGER,DIMENSION(:,:),ALLOCATABLE :: references ! (x y z molecule_type_index_ref molecule_type_index_obs)=1st dim, (nreference)=2nd dim
 REAL,DIMENSION(:,:),ALLOCATABLE :: distribution_function
 REAL :: step_a,step_b !step sizes in angströms or radians
 REAL :: ideal_density
 !PRIVATE/PUBLIC declarations
 PRIVATE :: operation_mode,number_of_references,sampling_interval,references,weigh_charge,subtract_uniform,maxdist
 PRIVATE :: distribution_function
 PUBLIC :: user_distribution_input,perform_distribution_analysis,write_simple_sumrules

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
   PRINT *,"Thus, the distribution functions are calculated with respect to a specified reference vector."
   PRINT *,"In cylindrical coordinates, the independent variables are"
   PRINT *,"  a) the distance perpendicular to the reference vector and"
   PRINT *,"  b) the distance collinear to the reference vector."
   PRINT *,"In polar coordinates, the independent variables are"
   PRINT *,"  a) the distance between the two molecules in the considered pair and"
   PRINT *,"  b) the angle between reference vector and the vector connecting the molecules."
   PRINT *,"Would you like to request the distribution function in polar (y) or cylindrical (n) coordinates?"
   IF (user_input_logical()) THEN
    operation_mode="pdf"
   ELSE
    operation_mode="cdf"
   ENDIF
   PRINT *,"You now need to specify the combinations of reference vector and the two molecule types for the distribution."
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
    PRINT *,"It is possible to consider *all* molecule types around this reference molecules."
    PRINT *,"Would you like to use all (y) molecules or only those of a specific type (n)?"
    IF (user_input_logical()) THEN
     references(5,n)=-1
    ELSE
     PRINT *,"Which molecule type would you like to observe / use for the distribution function?"
     references(5,n)=user_input_integer(1,maxmol)
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
    PRINT *,"You need to choose a maximum distance - molecule pairs beyond that threshold will be disregarded."
    PRINT *,"It is recommended to use half the box size where available. Would you like to use that shortcut? (y/n)"
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
   PRINT *,"Should the entries for the histogram be weighed by their charge?"
   PRINT *,"(Useful only in combination with 'all surrounding molecules')"
   weigh_charge=user_input_logical()
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
   CASE DEFAULT
    CALL report_error(0)
   END SELECT
   DO n=1,number_of_references,1
    WRITE(8,ADVANCE="NO",FMT='(" ",SP,I0," ",I0," ",I0,SS," ",I0," ",I0," ### ")') references(:,n)
    IF (references(5,n)<1) THEN
     WRITE(8,ADVANCE="NO",FMT='(" *all* molecules around ")')
    ELSE
     WRITE(8,ADVANCE="NO",FMT='(" molecules of type ",I0," around ")') references(5,n)
    ENDIF
    WRITE(8,ADVANCE="NO",FMT='(" those of type ",I0)') references(4,n)
    IF (ALL(references(1:3,n)==0)) THEN
     WRITE(8,'(", randomised reference vector.")')
    ELSE
     WRITE(8,'(".")')
    ENDIF
   ENDDO
   IF (optimize_distance) THEN
    WRITE(8,'(" maxdist_optimize ### set maximum distance to half the box size")')
   ELSE
    WRITE(8,'(" maxdist ",F0.3," ### set maximum distance to half the box size")')maxdist
   ENDIF
   IF (subtract_uniform) THEN
    WRITE(8,'(" subtract_uniform T ### subtract anisotropic density - only for pdf")')
   ELSE
    WRITE(8,'(" subtract_uniform F ### do not subtract anisotropic density.")')
   ENDIF
   IF (weigh_charge) THEN
    WRITE(8,'(" weigh_charge T ### weigh with charge - only sensible when *all* molecules are used!")')
   ELSE
    WRITE(8,'(" weigh_charge F ### no charge weighing")')
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
    CASE DEFAULT
     CALL report_error(99)
    END SELECT
    CLOSE(UNIT=3)
   ELSE
    CALL report_error(31)!No input - no output. easy as that.
   ENDIF
   CONTAINS

    SUBROUTINE read_references()!This subroutine is responsible for reading the body of the distribution input file, connected as unit 3.
    IMPLICIT NONE
    INTEGER :: n,inputbin
    CHARACTER(LEN=32) :: inputstring
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
       maxdist=give_box_limit()/2.0
       IF (maxdist<0.0d0) THEN
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3)')&
        &"  setting 'maxdist' to default (=",maxdist_default,")"
        maxdist=maxdist_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.3,A)') "   setting 'maxdist' to ",maxdist," (half the box length)"
       ENDIF
      CASE ("subtract_uniform")
       IF (TRIM(operation_mode)=="pdf") THEN
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
        WRITE(*,*) "  subtract_uniform is only available for the polar distribution function."
       ENDIF
      CASE ("weigh_charge")
       BACKSPACE 3
       READ(3,IOSTAT=ios,FMT=*) inputstring,weigh_charge
       IF (ios/=0) THEN
        CALL report_error(100,exit_status=ios)
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') "   setting 'weigh_charge' to default (=",weigh_charge_default,")"
        weigh_charge=weigh_charge_default
       ELSE
        IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'weigh_charge' to ",weigh_charge
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
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$ END INTERFACE
  INTEGER,INTENT(IN) :: reference_number
  INTEGER(KIND=DP) :: distribution_histogram(bin_count_a,bin_count_b)
  INTEGER(KIND=DP) :: within_bin!how many pairs are in the bin?
  INTEGER(KIND=DP) :: out_of_bounds!how many pairs have exceeded the bounds? should be zero for pdf.
  INTEGER :: molecule_type_index_ref,molecule_type_index_obs,timestep
  LOGICAL :: all_surrounding_molecules!IF TRUE, THEN all surrounding molecules are considered, of any type.
  LOGICAL :: randomise!IF TRUE, THEN the reference vector is randomised each step
  LOGICAL :: aligned !IF TRUE, THEN reference vector is the z direction already
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
   IF (references(5,reference_number)<1) THEN
    !invalid molecule type index, i.e. look at all available!
    all_surrounding_molecules=.TRUE.
   ELSE
    molecule_type_index_obs=references(5,reference_number)
    all_surrounding_molecules=.FALSE.
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
     !$  WRITE(*,'(A,I0,A)') "   ### Parallel execution on ",OMP_get_num_threads()," threads (distribution function)"
     !$  CALL timing_parallel_sections(.TRUE.)
     !$ ENDIF
     !$OMP END SINGLE
     distribution_histogram_local(:,:)=0
     observed_type=molecule_type_index_obs
     !$OMP DO SCHEDULE(STATIC,1)
     DO timestep=1,give_number_of_timesteps(),sampling_interval
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
     ENDDO
     !$OMP END DO
     !$OMP CRITICAL
     distribution_histogram(:,:)=distribution_histogram(:,:)+distribution_histogram_local(:,:)
     !$OMP END CRITICAL
     !$OMP END PARALLEL
     !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
     !$  WRITE(*,ADVANCE="NO",FMT='("   ### End of parallelised section, took ")')
     !$  CALL timing_parallel_sections(.FALSE.)
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
        link_vector(:)=link_vector(:)+&
        &give_center_of_mass(timestep_in,observed_type,molecule_index_obs)-&
        &give_center_of_mass(timestep_in,molecule_type_index_ref,molecule_index_ref)
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
     CASE ("pdf")
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
      WRITE(*,'("   ideal density: ",ES9.3," particles/Angström^3")')ideal_density
      IF (DEVELOPERS_VERSION) WRITE(*,'("  ! total volume:  ",ES9.3," cubic Angströms")') total_volume
     ENDIF
     DO binpos_a=1,bin_count_a,1
      !Compute chunk area / partial integral
      SELECT CASE (TRIM(operation_mode))
      CASE ("cdf")
       ! a is the distance in the 'xy' plane, i.e. perpendicular to the reference vector.
       a=(FLOAT(binpos_a)*step_a) !outer limit of a distance
       chunk_area=pi*((a**2)-((a-step_a)**2))!radial and azimuthal part of volume
      CASE ("pdf")
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
       CASE DEFAULT
        CALL report_error(0)
       END SELECT
       !correct distribution function for chunk volume
       distribution_function(binpos_a,binpos_b)=distribution_function(binpos_a,binpos_b)/chunk_volume
      ENDDO
     ENDDO
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
  LOGICAL :: connected
   reference_charge=give_charge_of_molecule(references(4,reference_number))
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
   IF (references(5,reference_number)<1) THEN
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
     IF (reference_charge/=0) WRITE(*,*) "  writing (Coulomb) energy integral into file: ",&
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
    IF (reference_charge/=0) THEN
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
     IF (reference_charge/=0) THEN
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
    IF (reference_charge/=0) CLOSE(UNIT=4)
    CLOSE(UNIT=10)
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
     IF (TRIM(operation_mode)=="pdf") b=b*degrees
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
      IF (references(5,n)<1) THEN
       name_obs="all"
      ELSE
       name_obs=TRIM(give_sum_formula(references(5,n)))
      ENDIF
      WRITE(*,'(" Calculating distribution for reference vector ",I0," out of ",I0,".")') n,number_of_references
      WRITE(*,'("   (Vector ",SP,I0," ",I0," ",I0,SS,", ",A," for ",A," around ",A,")")')&
      &references(1:3,n),TRIM(operation_mode),TRIM(name_obs),TRIM(give_sum_formula(references(4,n)))
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
CHARACTER(LEN=1024) :: filename_rmmvcf="rmmvcf.inp",filename_dihedral="dihedral.inp",filename_reorient="reorient.inp" !correlation module standard filenames
CHARACTER(LEN=1024) :: filename_vc_components="vc_components.inp",filename_cacf_components="cacf_components.inp" !correlation module standard filenames
CHARACTER(LEN=1024) :: filename_conductivity="conductivity.inp"!correlation module standard filenames
CHARACTER(LEN=1024) :: filename_msd="diffusion.inp" !diffusion module standard filename
CHARACTER(LEN=1024) :: filename_distribution="distribution.inp" !distribution module standard filename
INTEGER :: i,number_of_molecules!number_of_molecules is required for some safety checks, and initialised in generate_molecular_input()
INTEGER :: nsteps!nsteps is required again for checks (tmax...), and is initialised in generate_molecular_input()
 CALL set_default_masses()
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
 PRINT *, "   Copyright (C) 2020 Frederik Philippi (Tom Welton Group)"
 PRINT *, "   Please report any bugs."
 PRINT *,"    Suggestions and questions are also welcome. Thanks."
 PRINT *, "   Date of Release: 14_May_2020"
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
       IF (TRIM(inputstring)=="sequential_read") THEN
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
   PRINT *,"   Mean Squared Displacement: (PARALLELISED)"
   PRINT *,"   Calculates the mean squared displacement including a drift correction."
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
   PRINT *,"   Distribution Functions: (PARALLELISED)"
   PRINT *,"   Calculates distribution functions referenced to a fixed vector in space."
   PRINT *,"   cylindrical and polar coordinate systems are supported."
   PRINT *,"   Different references vectors can be chosen, including a randomised reference vector (debug purposes)."
   PRINT *,"   This could be e.g. '0 0 1' (z-direction is the reference)."
   PRINT *,"   In both coordinate systems, values are averaged over all azimuthal angles."
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
   PRINT *,"    Note that drude particles are added to N,O,C,S,P,Li - but not H and F."
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
   PRINT *,"      3) Any positive drude mass 'X', if present, is subtracted from N,O,C,S,P,Li."
   PRINT *,"      4) After that, masses of particular atoms are overwritten by 'atomic_masses'."
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
   PRINT *,"diffusion input file:"
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
   PRINT *,"distribution input file:"
   PRINT *,"The first line contains the expression 'cdf' or 'pdf', followed by the number of references N."
   PRINT *,"'cdf' requests the cylindrical distribution function, 'pdf' requests the polar distribution function."
   PRINT *,"The references are read from the N lines following the first line. The format of each line is:"
   PRINT *,"x - y - z - number of reference molecule type - number of observed molecule type."
   PRINT *,"For molecule type 1 around 2 relative to z-direction, the line would thus contain '0 0 1 2 1'."
   PRINT *,"a 'zero' reference vector, i.e. '0 0 0', triggers the randomisation of the reference vector."
   PRINT *,"If the molecule type of the observed molecule is -1, then *all* molecule types are used."
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
   PRINT *," - 'subtract_uniform'"
   PRINT *,"    If yes (T), then the uniform density / radial distribution function is subtracted."
   PRINT *,"    Only available for the polar distribution function."
   PRINT *," - 'weigh_charge'"
   PRINT *,"    The charges of observed molecules are added to the distribution histogram, rather than unity."
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
   PRINT *," - 'quit'"
   PRINT *,"    Terminates the analysis. Lines after this switch are ignored."
   PRINT *
   PRINT *,"reorientation input file:"
   PRINT *,"The molecule type index (= the molecule to observe) is given in the first line."
   PRINT *,"The expession 'reorientation' to request the appropriate analysis"
   PRINT *,"A fragment is defined by the expression 'base' or 'tip', followed by the number of atoms in this fragment."
   PRINT *,"The immediately following line must contain a list of the atom indices in this fragment."
   PRINT *,"For example, these lines define atom 16 as the base fragment and atoms 1, 3 and 4 as the tip fragment:"
   PRINT *,"  base 1"
   PRINT *,"  16"
   PRINT *,"  tip 3"
   PRINT *,"  3 1 4"
   PRINT *,"The two fragments *must* be defined, and they must appear before the quit statement (if applicable)."
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
    SELECT CASE (user_input_integer(0,14))
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
    PRINT *," 2 - Calculate mean squared displacements"
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
    PRINT *," 19 - Calculate close contact distances (inter- and intramolecular)."
    PRINT *," 20 - Calculate polar or cylindrical distribution function."
    PRINT *," 21 - Jump analysis / jump velocity distribution"
    PRINT *," 22 - Print atomic masses (in format suitable for a molecular input file)"
    PRINT *," 23 - Write trajectory with neighbouring molecules around one selected molecule"
    PRINT *," 24 - Dump close contact dimers for one step"
    SELECT CASE (user_input_integer(0,24))
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
      PRINT *,"There is a default mode available for this analysis that doesn't require additional input."
      PRINT *,"Would you like to take this shortcut? (y/n)"
      IF (user_input_logical()) THEN
       CALL append_string("diffusion_simple ### mean squared displacement with default values.")
       CYCLE
      ENDIF
      CALL user_msd_input(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps,filename_msd)
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
42 FORMAT ("    ########   #####")
44 FORMAT ("        ########")
43 FORMAT ("        #      #")
28 FORMAT ("    #   ##  #      #")
45 FORMAT ("    #      #       #")
13 FORMAT ("    ##### #### #####")
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
IMPLICIT NONE
  !$ INTERFACE
  !$  FUNCTION OMP_get_num_threads()
  !$  INTEGER :: OMP_get_num_threads
  !$  END FUNCTION OMP_get_num_threads
  !$  FUNCTION OMP_get_max_threads()
  !$  INTEGER :: OMP_get_max_threads
  !$  END FUNCTION OMP_get_max_threads
  !$  FUNCTION OMP_get_num_procs()
  !$  INTEGER :: OMP_get_num_procs
  !$  END FUNCTION OMP_get_num_procs
  !$  SUBROUTINE OMP_set_num_threads(number_of_threads)
  !$  INTEGER,INTENT(IN) :: number_of_threads
  !$  END SUBROUTINE OMP_set_num_threads
  !$ END INTERFACE
CHARACTER(LEN=32) :: inputstring,dummy
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
    IF (TRIM(inputstring)=="exit") inputstring="quit"!support for synonyms
    IF (TRIM(inputstring)=="read_sequential") inputstring="sequential_read"!support for synonyms
    IF (TRIM(inputstring)=="time_scaling_factor") inputstring="time_scaling"!support for synonyms
    IF (TRIM(inputstring)=="parallel_execution") inputstring="parallel_operation"!support for synonyms
    IF (TRIM(inputstring)=="print_atom_masses") inputstring="print_atomic_masses"!support for synonyms
    IF (TRIM(inputstring)=="print_atomic_weights") inputstring="print_atomic_masses"!support for synonyms
    IF (TRIM(inputstring)=="print_atom_weights") inputstring="print_atomic_masses"!support for synonyms
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
     CALL dump_split(startstep,endstep)
    CASE ("dump_split_simple") !Module DEBUG
     IF (INFORMATION_IN_TRAJECTORY=="VEL") CALL report_error(56)
     WRITE(*,'(A,I0,A)') " Trajectory will be split into ",&
     &give_number_of_molecule_types()," molecule types - simple mode."
     startstep=1
     endstep=give_number_of_timesteps()
     WRITE(*,'(A,I0,A,I0,A)') " (For timesteps ",startstep," to ",endstep,")"
     CALL dump_split(startstep,endstep)
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
      CALL remove_drudes(startstep,endstep)
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
      CALL remove_drudes(startstep,endstep)
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
    CASE ("cubic_box_edge")
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
     !$  CALL report_error(19,exit_status=ios)
     !$  EXIT
     !$ ENDIF
     !$ IF (inputinteger<1) THEN
     !$  inputinteger=OMP_get_max_threads()
     !$ ELSEIF (inputinteger>OMP_get_num_procs()) THEN
     !$  CALL report_error(42,exit_status=inputinteger)
     !$  inputinteger=OMP_get_num_procs()
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
     CALL dump_neighbour_traj(.FALSE.,1,10,2,1,1,2)
     !update_com,startstep_in,endstep_in,molecule_type_index_1,molecule_index_1,molecule_type_index_2,neighbour_num
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
 15 FORMAT ("    ",A," ",L1)
 16 FORMAT ("    ",A," ",I0)
   WRITE(*,*) "Printing current global settings."
   WRITE(*,15) "VERBOSE_OUTPUT      ",VERBOSE_OUTPUT
   WRITE(*,15) "TIME_OUTPUT         ",TIME_OUTPUT
   WRITE(*,15) "DEVELOPERS_VERSION  ",DEVELOPERS_VERSION
   WRITE(*,15) "ERROR_OUTPUT        ",ERROR_OUTPUT
   WRITE(*,15) "READ_SEQUENTIAL     ",READ_SEQUENTIAL
   WRITE(*,15) "BOX_VOLUME_GIVEN    ",BOX_VOLUME_GIVEN
   WRITE(*,15) "WRAP_TRAJECTORY     ",WRAP_TRAJECTORY
   WRITE(*,15) "DISCONNECTED        ",DISCONNECTED
   WRITE(*,16) "GLOBAL_ITERATIONS   ",GLOBAL_ITERATIONS
   WRITE(*,16) "TIME_SCALING_FACTOR ",TIME_SCALING_FACTOR
   WRITE(*,16) "HEADER_LINES_GINPUT ",HEADER_LINES
   WRITE(*,16) "CURRENT_ERROR_CODE  ",ERROR_CODE
   WRITE(*,16) "ERROR_COUNT         ",give_error_count()
   WRITE(*,*) '   OUTPUT_PREFIX        "',TRIM(OUTPUT_PREFIX),'"'
   WRITE(*,*) '   TRAJECTORY_TYPE      "',TRIM(TRAJECTORY_TYPE),'"'
   WRITE(*,*) '   INFO_IN_TRAJECTORY   "',TRIM(INFORMATION_IN_TRAJECTORY),'"'
   WRITE(*,15) "PARALLEL_OPERATION  ",PARALLEL_OPERATION
   !$ IF(.FALSE.) THEN
   WRITE(*,*) "-fopenmp flag not set! PARALLEL_OPERATION has no effect."
   !$ ENDIF
   !$ WRITE(*,*) "OMP flag set!"
   !$ WRITE(*,16) "NUMBER_OF_THREADS     ",OMP_get_num_threads()
   !$ WRITE(*,16) "NUMBER_OF_PROCS       ",OMP_get_num_procs()
   !$ WRITE(*,16) "MAX_NUMBER_OF_THREADS ",OMP_get_max_threads()
   CALL show_molecular_settings()
   WRITE(*,*) "(input) Filenames:"
   WRITE(*,*) '   TRAJECTORY      "',TRIM(FILENAME_TRAJECTORY),'"'
   WRITE(*,*) '   GENERAL         "',TRIM(FILENAME_GENERAL_INPUT),'"'
   WRITE(*,*) '   MOLECULAR       "',TRIM(FILENAME_MOLECULAR_INPUT),'"'
   WRITE(*,*) '   AUTOCORRELATION "',TRIM(FILENAME_AUTOCORRELATION_INPUT),'"'
   WRITE(*,*) '   DIFFUSION       "',TRIM(FILENAME_DIFFUSION_INPUT),'"'
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

PROGRAM PREALPHA ! Copyright (C) 2020 Frederik Philippi
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