
!This Module calculates a speciation distribution for donor-acceptor interactions.
!For example, considering the Li atom as acceptor and the O atoms as donors in a mixture
!of Li[PO2F2] and dimethoxyethane, these are some possible species:
!  -  Li coordinated by one DME and one anion, each with two oxygen atoms
!  -  Li coordinated by two DMEs with one oxygen each and one anion with two oxygens
!  -  Li coordinated by four anions, each with one oxygen
!  -  Li coordinated by two DMEs, with two oxygens each
! ... and so on, here I only listed those that lead to exactly four neighbours.
!This module identifies these different species and keeps track of them and how often they occur.
MODULE SPECIATION ! Copyright (C) !RELEASEYEAR! Frederik Philippi
	USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	INTEGER,PARAMETER :: nsteps_default=100 !how many steps to use from trajectory
	INTEGER,PARAMETER :: sampling_interval_default=1
	INTEGER,PARAMETER :: maximum_number_of_species_default=11!how many different species do we want to allow?
	INTEGER,PARAMETER :: maximum_number_of_neighbour_molecules_default=10
	INTEGER,PARAMETER :: maximum_number_of_connections_default=10
	INTEGER,PARAMETER :: tmax_default=10000
	LOGICAL,PARAMETER :: print_connection_beads_default=.FALSE.
	LOGICAL,PARAMETER :: calculate_autocorrelation_default=.FALSE.
	LOGICAL,PARAMETER :: use_logarithmic_spacing_default=.FALSE.
	!variables
	INTEGER :: nsteps=nsteps_default !how many steps to use from trajectory
	INTEGER :: sampling_interval=sampling_interval_default
	INTEGER :: maximum_number_of_species=maximum_number_of_species_default!how many different species do we want to allow?
	INTEGER :: maximum_number_of_connections=maximum_number_of_connections_default
	INTEGER :: N_acceptor_types
	INTEGER :: N_donor_types
	INTEGER :: bytes_needed
	INTEGER :: tmax
	INTEGER :: neighbour_atom_overflow
	LOGICAL :: atoms_grouped,grouped_by_elements
	LOGICAL :: print_connection_beads=print_connection_beads_default
	LOGICAL :: calculate_autocorrelation=calculate_autocorrelation_default
	LOGICAL :: use_logarithmic_spacing=use_logarithmic_spacing_default
	REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: squared_cutoff_list !the squared cutoffs. The dimensions are: entry of acceptor_list -- entry of acceptor atom_index_list -- entry of donor_list -- entry of donor atom_index_list.
	TYPE,PRIVATE :: connections
		!This one saves a specific pair of atoms as being connected.
		!The reference molecule_type_index does not need to be recorded since connection is owned by species which is owned by the acceptor_list or the donor_list
		INTEGER :: indices(4) !First dimension is 1/2/3/4 for:
		! 1)  observed (donor) type, i.e. n_donor
		! 2)  observed (donor) molecule_index
		! 3)  observed (donor) atom (index in atom_index_list)
		! 4)  reference (=acceptor) atom (index in atom_index_list)
		! AT THE MOMENT, the list is sorted naturally by 1)-2)-3)-4).
		!HOWEVER, note that after collapsing, the entries of acceptor/donor atom_index_list are mangled (i.e. replaced with the atom groups). first_occurrence and last_occurrence should be kept the same maybe.
	END TYPE connections
	TYPE,PRIVATE :: neighbour_molecule_status
		LOGICAL :: unmatched_entry
		INTEGER :: first_index_in_connection_list
		INTEGER :: extra_atom_entries_in_connection_list
	END TYPE neighbour_molecule_status
	TYPE,PRIVATE :: species
		INTEGER :: total_number_of_neighbour_atoms ! For example, "4" for Li[DME]2
		INTEGER :: total_number_of_neighbour_molecules ! For example, "2" for Li[DME]2
		INTEGER(KIND=WORKING_PRECISION) :: occurrences ! For example, "I have found this species 34857 times"
		INTEGER :: rank !1 for the most probable, 2 for the next, 3 for the third most probable, and so one.
		INTEGER :: timestep_first_occurrence, timestep_last_occurrence
		TYPE(neighbour_molecule_status),DIMENSION(:),ALLOCATABLE :: neighbour_molecule_starting_indices !allocated up to maximum_number_of_connections but filled only up to total_number_of_neighbour_molecules with the starting indices of said molecules.
		TYPE(connections),DIMENSION(:),ALLOCATABLE :: first_occurrence
		TYPE(connections),DIMENSION(:),ALLOCATABLE :: last_occurrence
		INTEGER :: first_occurrence_acceptorindices(2) !acceptor_type/acceptor_molecule_index
		INTEGER :: last_occurrence_acceptorindices(2) !acceptor_type/acceptor_molecule_index
		INTEGER :: number_of_logged_connections_in_species_list !this is equal to the total number of connections and serves as one of the three descriptors
		TYPE(connections),DIMENSION(:),ALLOCATABLE :: connection !First dimension goes from 1 to number_of_logged_connections_in_species_list and is allocated up to maximum_number_of_connections.
	END TYPE species
	TYPE,PRIVATE :: donors_and_acceptors
		INTEGER :: N_atoms!array size of atom_index_list and cutoff_list
		INTEGER :: N_groups!how many atom groups
		INTEGER :: molecule_type_index
		INTEGER :: species_overflow
		INTEGER(KIND=WORKING_PRECISION) :: no_neighbour_occurrences!how many times NO species was found
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: atom_index_list !First dimension is listing all relevant atom_index up to N_atoms.
		!the second dimension is 1/2 for the actual atom_index/atom_group.
		!All atoms with the same atom_group number are treated as identical.
		REAL,DIMENSION(:),ALLOCATABLE :: cutoff_list ! CONTRIBUTIONS to the cutoff, i.e. just the radius, counting up to N_atoms
		INTEGER :: number_of_logged_species_in_acceptor_list
		TYPE(species),DIMENSION(:),ALLOCATABLE :: list_of_all_species !First dimension goes from 1 to number_of_logged_species_in_acceptor_list and is allocated up to maximum_number_of_species
	END TYPE donors_and_acceptors
	TYPE(donors_and_acceptors),DIMENSION(:),ALLOCATABLE :: acceptor_list !First dimension counts the acceptor molecule types
	TYPE(donors_and_acceptors),DIMENSION(:),ALLOCATABLE :: donor_list    !First dimension counts the donor molecule types
	REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: average_h!the occurrences, including "species 0" for the overflows
	LOGICAL :: dumpfirst=.TRUE.
	LOGICAL :: dumplast=.TRUE.
	!PRIVATE/PUBLIC declarations
	PUBLIC :: perform_speciation_analysis,user_speciation_input

	CONTAINS

!WRITING input file to unit 8, which shouldn't be open.
		!has to be compliant with 'read_velocity_correlation_body' in 'AUTOCORRELATION' module
		SUBROUTINE user_speciation_input&
		&(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps_in,filename_speciation)
		IMPLICIT NONE
		CHARACTER (LEN=*) :: filename_speciation
		LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
		INTEGER,INTENT(IN) :: nsteps_in,number_of_molecules
		INTEGER :: ios,atom_counter,howmanyatoms,inputinteger1,inputinteger2,maxmol
		REAL :: inputreal
		LOGICAL :: connected
			PRINT *,"Generating input for speciation analysis."
			PRINT *,"A good range of analysis covers usually 1000 steps in total, distributed over the whole trajectory."
			PRINT *,"HOWEVER, if you want the lifetimes, you should cover the whole trajectory!"
			PRINT *,"Up to which step number of the trajectory do you want the analysis to run?"
			WRITE(*,'(" The default is currently set to ",I0,".")') nsteps_default
			nsteps=user_input_integer(1,(nsteps_in-1))
			PRINT *,"Every how many steps would you like to use?"
			WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
			sampling_interval=user_input_integer(1,nsteps_in)
			parallelisation_possible=.TRUE.
			maxmol=number_of_molecules
			IF (number_of_molecules==-1) maxmol=10000!unknown molecule number... expect the worst.
			WRITE(*,FMT='(A)',ADVANCE="NO") " opening speciation input file..."
			INQUIRE(UNIT=8,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=8)
			OPEN(UNIT=8,FILE=TRIM(PATH_INPUT)//TRIM(OUTPUT_PREFIX)//TRIM(filename_speciation),IOSTAT=ios)
			IF (ios/=0) CALL report_error(46,exit_status=ios)
			PRINT *,"How many acceptor atoms do you want to specify?"
			howmanyatoms=user_input_integer(1,20)
			WRITE(8,'(" ",I0," acceptors")') howmanyatoms
			DO atom_counter=1,howmanyatoms
				WRITE(*,'(" Reading information for acceptor atom ",I0," / ",I0,":")') atom_counter,howmanyatoms
				PRINT *,"Please enter the molecule_type_index of the molecule this atom belongs to:"
				inputinteger1=user_input_integer(1,maxmol)
				PRINT *,"Please enter the atom_index of this atom:"
				inputinteger2=user_input_integer(1,10000)
				PRINT *,"Please enter the corresponding cutoff distance of this atom:"
				inputreal=user_input_real(-2.0,15.0)
				WRITE(8,'(" ",I0," ",I0," ",F0.3)') inputinteger1,inputinteger2,inputreal
			ENDDO
			PRINT *,"How many donor atoms do you want to specify?"
			howmanyatoms=user_input_integer(1,40)
			WRITE(8,'(" ",I0," donors")') howmanyatoms
			DO atom_counter=1,howmanyatoms
				WRITE(*,'(" Reading information for donor atom ",I0," / ",I0,":")') atom_counter,howmanyatoms
				PRINT *,"Please enter the molecule_type_index of the molecule this atom belongs to:"
				inputinteger1=user_input_integer(1,maxmol)
				PRINT *,"Please enter the atom_index of this atom:"
				inputinteger2=user_input_integer(1,10000)
				PRINT *,"Please enter the corresponding cutoff distance of this atom:"
				inputreal=user_input_real(-2.0,15.0)
				WRITE(8,'(" ",I0," ",I0," ",F0.3)') inputinteger1,inputinteger2,inputreal
			ENDDO
			PRINT *,"Do you want to group atoms (y) or use every atom_index as its own group?"
			IF (user_input_logical()) THEN
				PRINT *,"Do you want to group atoms by element symbol (y) or specify groups manually (n)?"
				IF (user_input_logical()) THEN
					WRITE(8,*) "group_elements"
				ELSE
					PRINT *,"Do you want to create a new acceptor group (y/n)?"
					DO WHILE (user_input_logical())
						PRINT *,"For which molecule_type_index do you want to create this group?"
						inputinteger1=user_input_integer(1,maxmol)
						WRITE(8,'(" new_acceptor_group",I0)') inputinteger1
						PRINT *,"How many atoms do you want to group together?"
						howmanyatoms=user_input_integer(1,40)
						DO atom_counter=1,howmanyatoms
							WRITE(*,'(" Please enter the atom_index for atom ",I0," / ",I0," of this group:")')&
							&atom_counter,howmanyatoms
							inputinteger2=user_input_integer(1,maxmol)
							WRITE(8,'(" assign_to_acceptor_group",I0)') inputinteger2
						ENDDO
						PRINT *,"Do you want to add *another* new acceptor group (y/n)?"
					ENDDO
					PRINT *,"Do you want to create a new donor group (y/n)?"
					DO WHILE (user_input_logical())
						PRINT *,"For which molecule_type_index do you want to create this group?"
						inputinteger1=user_input_integer(1,maxmol)
						WRITE(8,'(" new_donor_group",I0)') inputinteger1
						PRINT *,"How many atoms do you want to group together?"
						howmanyatoms=user_input_integer(1,40)
						DO atom_counter=1,howmanyatoms
							WRITE(*,'(" Please enter the atom_index for atom ",I0," / ",I0," of this group:")')&
							&atom_counter,howmanyatoms
							inputinteger2=user_input_integer(1,maxmol)
							WRITE(8,'(" assign_to_donor_group",I0)') inputinteger2
						ENDDO
						PRINT *,"Do you want to add *another* new donor group (y/n)?"
					ENDDO
				ENDIF
			ENDIF
			PRINT *,"How many species do you want to allow per acceptor molecule type?"
			inputinteger1=user_input_integer(1,100)
			WRITE(8,'(" N_species ",I0," ### maximum number of species per acceptor molecule type")') inputinteger1
			PRINT *,"How many neighbour connections do you want to allow per acceptor molecule?"
			inputinteger1=user_input_integer(1,200)
			WRITE(8,'(" N_neighbours ",I0," ### maximum number of neighbour connections per acceptor molecule")') inputinteger1
			PRINT *,"Now let's talk about species lifetimes."
			PRINT *,"Do you want to calculate the intermittent binary autocorrelation function (y/n)?"
			calculate_autocorrelation=user_input_logical()
			IF (calculate_autocorrelation) THEN
				WRITE(8,'(" autocorrelation T ### calculate species time correlation")')
				PRINT *,"Do you want to use logarithmic spacing of timesteps (ca 8-10x faster) (y/n)?"
				IF (user_input_logical()) THEN
					WRITE(8,'(" use_log T ### logarithmic spacing of timesteps for autocorrelation")')
				ELSE
					WRITE(8,'(" use_log F ### linear spacing of timesteps for autocorrelation, maximum resolution")')
				ENDIF
				PRINT *,"Please enter the maximum time shift of the correlation function."
				tmax=user_input_integer(1,nsteps-1)
				WRITE(8,'(" tmax ",I0," ### maximum time shift of the correlation function")') tmax
			ELSE
				WRITE(8,'(" autocorrelation F ### skip the autocorrelation")')
			ENDIF
			PRINT *,"Do you want to print beads indicating the connections in the example xyz files (y/n)?"
			IF (user_input_logical()) THEN
				WRITE(8,'(" print_beads T ### show the connections as beads in the output .xyz files")')
			ELSE
				WRITE(8,'(" print_beads F ### do not indicate connections in the output .xyz files")')
			ENDIF
			IF (.NOT.(parallelisation_requested)) THEN!... but hasn't been requested so far. Thus, ask for it.
				PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
				parallelisation_requested=user_input_logical()
			ENDIF
			WRITE(8,*) "quit"
			WRITE(8,*)
			WRITE(8,*) "This is an input file to calculate speciation statistics, lifetimes, and dump structures."
			WRITE(8,*) "To actually perform the implied calculation, it has to be referenced in 'general.inp'."
			ENDFILE 8
			CLOSE(UNIT=8)
			WRITE(*,*) "done"
		END SUBROUTINE user_speciation_input

		SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
		IMPLICIT NONE
			nsteps=nsteps_default
			sampling_interval=sampling_interval_default
			maximum_number_of_species=maximum_number_of_species_default
			dumpfirst=.TRUE.!legacy
			dumplast=.TRUE.
			N_acceptor_types=0
			N_donor_types=0
			maximum_number_of_connections=maximum_number_of_connections_default
			atoms_grouped=.FALSE.
			grouped_by_elements=.FALSE.
			print_connection_beads=print_connection_beads_default
			bytes_needed=0
			tmax=tmax_default
			calculate_autocorrelation=calculate_autocorrelation_default
			use_logarithmic_spacing=use_logarithmic_spacing_default
		END SUBROUTINE set_defaults

		!initialises the speciation module by reading the specified input file.
		SUBROUTINE initialise_speciation()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,allocstatus,n,m,acceptor_maxatoms,donor_maxatoms,a,b
		REAL :: real_space_distance
			! first, check if file exists.
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_SPECIATION_INPUT),EXIST=file_exists)
			IF (file_exists) THEN
				!setting defaults to start with.
				CALL set_defaults()
				IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(PATH_INPUT)//TRIM(FILENAME_SPECIATION_INPUT),"'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_SPECIATION_INPUT),&
				&ACTION='READ',IOSTAT=ios)
				IF (ios/=0) CALL report_error(155,exit_status=ios)
				CALL read_acceptors_and_donors()
				!read rest of input file.
				CALL read_body()
				CLOSE(UNIT=3)
				acceptor_maxatoms=0
				donor_maxatoms=0
				!preparing the data structures. We only consider the "acceptors" as reference.
				DO n=1,N_acceptor_types,1
					IF (ALLOCATED(acceptor_list(n)%list_of_all_species)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%list_of_all_species(maximum_number_of_species),STAT=allocstatus)
					!these are 4 bytes, actually
					bytes_needed=bytes_needed+13*maximum_number_of_species
					bytes_needed=bytes_needed+6*maximum_number_of_species*maximum_number_of_connections
					IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
					!Initialise the list_of_all_species
					DO m=1,maximum_number_of_species,1
						IF (ALLOCATED(acceptor_list(n)%list_of_all_species(m)%&
						&connection)) CALL report_error(0)
						ALLOCATE(acceptor_list(n)%list_of_all_species(m)%&
						&connection(maximum_number_of_connections),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
						IF (ALLOCATED(acceptor_list(n)%list_of_all_species(m)%&
						&first_occurrence)) CALL report_error(0)
						ALLOCATE(acceptor_list(n)%list_of_all_species(m)%&
						&first_occurrence(maximum_number_of_connections),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
						IF (ALLOCATED(acceptor_list(n)%list_of_all_species(m)%&
						&last_occurrence)) CALL report_error(0)
						ALLOCATE(acceptor_list(n)%list_of_all_species(m)%&
						&last_occurrence(maximum_number_of_connections),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
						acceptor_list(n)%list_of_all_species(m)%first_occurrence_acceptorindices(:)=0
						acceptor_list(n)%list_of_all_species(m)%last_occurrence_acceptorindices(:)=0
						IF (ALLOCATED(acceptor_list(n)%list_of_all_species(m)%&
						&neighbour_molecule_starting_indices)) CALL report_error(0)
						ALLOCATE(acceptor_list(n)%list_of_all_species(m)%&
						&neighbour_molecule_starting_indices(maximum_number_of_connections),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
					ENDDO
					acceptor_list(n)%list_of_all_species(:)%occurrences=0
					acceptor_list(n)%list_of_all_species(:)%timestep_last_occurrence=0
					acceptor_list(n)%list_of_all_species(:)%timestep_first_occurrence=0
					acceptor_list(n)%list_of_all_species(:)%total_number_of_neighbour_atoms=0
					acceptor_list(n)%list_of_all_species(:)%total_number_of_neighbour_molecules=0
					acceptor_list(n)%list_of_all_species(:)%number_of_logged_connections_in_species_list=0
					!while we are at it, search for the most atoms.
					!this will determine the upper limit of the acceptor part of the cutoff list.
					IF (acceptor_maxatoms<acceptor_list(n)%N_atoms) THEN
						acceptor_maxatoms=acceptor_list(n)%N_atoms
					ENDIF
				ENDDO
				!similarly, search for the most atoms among the donor molecule types.
				DO n=1,N_donor_types,1
					!this will determine the upper limit of the donor part of the cutoff list.
					IF (donor_maxatoms<donor_list(n)%N_atoms) THEN
						donor_maxatoms=donor_list(n)%N_atoms
					ENDIF
				ENDDO
				!preparing the squared cutoff list.
				IF (DEVELOPERS_VERSION) THEN
					WRITE(*,'("  ! Acceptor atoms upper limit = ",I0)') acceptor_maxatoms
					WRITE(*,'("  ! Donor atoms upper limit = ",I0)') donor_maxatoms
					WRITE(*,'("  ! squared_cutoff_list number of entries = ",I0)')&
					&N_acceptor_types*acceptor_maxatoms*N_donor_types*donor_maxatoms
				ENDIF
				IF (ALLOCATED(squared_cutoff_list)) CALL report_error(0)
				ALLOCATE(squared_cutoff_list&
				&(N_acceptor_types,acceptor_maxatoms,N_donor_types,donor_maxatoms),STAT=allocstatus)
				bytes_needed=bytes_needed+N_acceptor_types*acceptor_maxatoms*N_donor_types*donor_maxatoms
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				DO n=1,N_acceptor_types,1
					DO m=1,N_donor_types,1
						DO a=1,acceptor_list(n)%N_atoms
							DO b=1,donor_list(m)%N_atoms
								real_space_distance=&
								acceptor_list(n)%cutoff_list(a)+& !acceptor contribution
								&donor_list(m)%cutoff_list(b) ! Donor contribution
								squared_cutoff_list(n,a,m,b)=real_space_distance**2
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ELSE
				CALL report_error(158)!No input - no output. easy as that.
			ENDIF

			CONTAINS

				SUBROUTINE read_body()
				IMPLICIT NONE
				CHARACTER(LEN=33) :: inputstring
				LOGICAL :: successful_group_assignment
				INTEGER :: n_acceptor,molecule_type_index,atom_index,n_donor,n_acceptor_last,n_donor_last,atom_counter
					molecule_type_index=0
					n_acceptor_last=0
					n_donor_last=0
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  End-of-file condition in ",TRIM(FILENAME_SPECIATION_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_SPECIATION_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("maximum_number_of_connections","neighbours","N_neighbours")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,maximum_number_of_connections
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'maximum_number_of_connections' to default (=",&
								&maximum_number_of_connections_default,")"
								maximum_number_of_connections=maximum_number_of_connections_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'maximum_number_of_connections' to ",&
								&maximum_number_of_connections
							ENDIF
						CASE ("group_elements","group_element_names","group_element_symbols")
							IF (atoms_grouped) THEN
								CALL report_error(162)
							ENDIF
							IF (VERBOSE_OUTPUT) WRITE(*,'("   grouping atoms together based on their element symbols.")')
							CALL group_element_names()
							atoms_grouped=.TRUE.
							grouped_by_elements=.TRUE.
						CASE ("assign_to_acceptor_group","assign_atom_to_acceptor_group",&
						&"add_to_acceptor_group","add_atom_to_acceptor_group")
							IF ((.NOT.(atoms_grouped)).OR.(molecule_type_index==0)) THEN
								CALL report_error(164)
							ELSE
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,atom_index
								IF (ios/=0) THEN
									CALL report_error(158,exit_status=ios)
								ELSE
									successful_group_assignment=.FALSE.
									DO n_acceptor=1,N_acceptor_types,1
										IF (acceptor_list(n_acceptor)%molecule_type_index==molecule_type_index) THEN
											IF (successful_group_assignment) CALL report_error(0)
											DO atom_counter=1,acceptor_list(n_acceptor)%N_atoms,1
												IF (acceptor_list(n_acceptor)%atom_index_list(atom_counter,1)==atom_index) THEN
													IF (successful_group_assignment) CALL report_error(165)
													successful_group_assignment=.TRUE.
													acceptor_list(n_acceptor)%atom_index_list(atom_counter,2)=&
													&-acceptor_list(n_acceptor)%N_groups
												ENDIF
											ENDDO
											IF ((VERBOSE_OUTPUT).AND.(successful_group_assignment)) &
											&WRITE(*,'("   Added atom_index ",I0," to acceptor group #",I0," (molecule_type_index ",I0,")")')&
											&atom_index,acceptor_list(n_acceptor)%N_groups,molecule_type_index
										ENDIF
									ENDDO
									IF (.NOT.(successful_group_assignment)) CALL report_error(111,atom_index)
								ENDIF
							ENDIF
						CASE ("assign_to_donor_group","assign_atom_to_donor_group",&
						&"add_to_donor_group","add_atom_to_donor_group")
							IF ((.NOT.(atoms_grouped)).OR.(molecule_type_index==0)) THEN
								CALL report_error(164)
							ELSE
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,atom_index
								IF (ios/=0) THEN
									CALL report_error(158,exit_status=ios)
								ELSE
									successful_group_assignment=.FALSE.
									DO n_donor=1,N_donor_types,1
										IF (donor_list(n_donor)%molecule_type_index==molecule_type_index) THEN
											IF (successful_group_assignment) CALL report_error(0)
											DO atom_counter=1,donor_list(n_donor)%N_atoms,1
												IF (donor_list(n_donor)%atom_index_list(atom_counter,1)==atom_index) THEN
													IF (successful_group_assignment) CALL report_error(165)
													successful_group_assignment=.TRUE.
													donor_list(n_donor)%atom_index_list(atom_counter,2)=&
													&-donor_list(n_donor)%N_groups
												ENDIF
											ENDDO
											IF ((VERBOSE_OUTPUT).AND.(successful_group_assignment)) &
											&WRITE(*,'("   Added atom_index ",I0," to donor group #",I0," (molecule_type_index ",I0,")")')&
											&atom_index,donor_list(n_donor)%N_groups,molecule_type_index
										ENDIF
									ENDDO
									IF (.NOT.(successful_group_assignment)) CALL report_error(111,atom_index)
								ENDIF
							ENDIF
						CASE ("new_acceptor_group","acceptor_new_group","acceptor_group","acceptors_new_group")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,molecule_type_index
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A)')&
								&"   no changes to atom groups."
							ELSE
								IF (grouped_by_elements) THEN
									CALL report_error(163)
									grouped_by_elements=.FALSE.
								ELSE
									IF (.NOT.(atoms_grouped)) THEN
										WRITE(*,'("   turned on atom groups.")')
										DO n_acceptor=1,N_acceptor_types,1
											acceptor_list(n_acceptor)%N_groups=0
										ENDDO
									ENDIF
								ENDIF
								successful_group_assignment=.FALSE.
								DO n_acceptor=1,N_acceptor_types,1
									IF (acceptor_list(n_acceptor)%molecule_type_index==molecule_type_index) THEN
										IF (successful_group_assignment) CALL report_error(0)
										successful_group_assignment=.TRUE.
										acceptor_list(n_acceptor)%N_groups=acceptor_list(n_acceptor)%N_groups+1
										IF (VERBOSE_OUTPUT) &
										&WRITE(*,'("   Added group ",I0," to acceptor #",I0," (molecule_type_index ",I0,")")')&
										&acceptor_list(n_acceptor)%N_groups,&
										&n_acceptor,molecule_type_index
										n_acceptor_last=n_acceptor
										atoms_grouped=.TRUE.
									ENDIF
								ENDDO
								IF (.NOT.(successful_group_assignment)) CALL report_error(110,molecule_type_index)
							ENDIF
						CASE ("new_donor_group","donor_new_group","donor_group","donors_new_group")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,molecule_type_index
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A)')&
								&"   no changes to atom groups."
							ELSE
								IF (grouped_by_elements) THEN
									CALL report_error(163)
									grouped_by_elements=.FALSE.
								ELSE
									IF (.NOT.(atoms_grouped)) THEN
										WRITE(*,'("   turned on atom groups.")')
										DO n_donor=1,N_donor_types,1
											donor_list(n_donor)%N_groups=0
										ENDDO
									ENDIF
								ENDIF
								successful_group_assignment=.FALSE.
								DO n_donor=1,N_donor_types,1
									IF (donor_list(n_donor)%molecule_type_index==molecule_type_index) THEN
										IF (successful_group_assignment) CALL report_error(0)
										successful_group_assignment=.TRUE.
										donor_list(n_donor)%N_groups=donor_list(n_donor)%N_groups+1
										IF (VERBOSE_OUTPUT) &
										&WRITE(*,'("   Added group ",I0," to donor #",I0," (molecule_type_index ",I0,")")')&
										&donor_list(n_donor)%N_groups,&
										&n_donor,molecule_type_index
										n_donor_last=n_donor
										atoms_grouped=.TRUE.
									ENDIF
								ENDDO
								IF (.NOT.(successful_group_assignment)) CALL report_error(110,molecule_type_index)
							ENDIF
						CASE ("print_connection_beads","print_beads","print_connections")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,print_connection_beads
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&"   setting 'print_connection_beads' to default (=",print_connection_beads_default,")"
								print_connection_beads=print_connection_beads_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'print_connection_beads' to ",&
								&print_connection_beads
							ENDIF
						CASE ("logarithmic_spacing","logarithmic","log_spacing","use_logarithmic_spacing","use_log")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,use_logarithmic_spacing
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&"   setting 'use_logarithmic_spacing' to default (=",use_logarithmic_spacing_default,")"
								use_logarithmic_spacing=use_logarithmic_spacing_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'use_logarithmic_spacing' to ",&
								&use_logarithmic_spacing
							ENDIF
						CASE ("correlate","autocorrelation","time_correlation","calculate_autocorrelation")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,calculate_autocorrelation
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&"   setting 'calculate_autocorrelation' to default (=",calculate_autocorrelation_default,")"
								calculate_autocorrelation=calculate_autocorrelation_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'calculate_autocorrelation' to ",&
								&calculate_autocorrelation
							ENDIF
						CASE ("maximum_number_of_species","species","N_species")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,maximum_number_of_species
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'maximum_number_of_species' to default (=",maximum_number_of_species_default,")"
								maximum_number_of_species=maximum_number_of_species_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'maximum_number_of_species' to ",maximum_number_of_species
							ENDIF
						CASE ("tmax")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,tmax
							IF (ios/=0) THEN
								CALL report_error(24,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') "   setting 'tmax' to default (=",tmax_default,")"
								tmax=tmax_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'tmax' to ",tmax
							ENDIF
						CASE ("sampling_interval")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'sampling_interval' to default (=",sampling_interval_default,")"
								sampling_interval=sampling_interval_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'sampling_interval' to ",sampling_interval
							ENDIF
						CASE ("nsteps")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,nsteps
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'nsteps' to default (=",nsteps_default,")"
								nsteps=nsteps_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'nsteps' to ",nsteps
							ENDIF
						CASE ("quit")
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_SPECIATION_INPUT)
							EXIT
						CASE DEFAULT
							IF (VERBOSE_OUTPUT) WRITE(*,*) "  can't interpret line - continue streaming"
						END SELECT
					ENDDO
					!first, check for sensible input.
					IF (nsteps<1) THEN
						CALL report_error(57,exit_status=nsteps)
						nsteps=1
					ELSEIF (nsteps>give_number_of_timesteps()) THEN
						CALL report_error(57,exit_status=nsteps)
						nsteps=give_number_of_timesteps()
					ENDIF
					IF ((tmax>MAX((nsteps-1+sampling_interval)/sampling_interval,0)-1).OR.(tmax<1)) THEN
						tmax=MAX((nsteps-1+sampling_interval)/sampling_interval,0)-1
						CALL report_error(28,exit_status=INT(tmax))
					ENDIF
					CALL showgroups()
				END SUBROUTINE read_body

				SUBROUTINE group_element_names()
				IMPLICIT NONE
				INTEGER :: group_number
				CHARACTER(LEN=2) :: current_element
				INTEGER :: acceptor_indexcounter,donor_indexcounter,n_acceptor,n_donor
					current_element="NO"
					DO n_acceptor=1,N_acceptor_types,1
						group_number=0
						IF (DEVELOPERS_VERSION) THEN
							WRITE(*,ADVANCE="NO",FMT='("  ! groups for acceptor #",I0,":")') n_acceptor
						ENDIF
						!set all group labels to be positive
						acceptor_list(n_acceptor)%atom_index_list(:,2)=1
						acceptor_indexcounter=1
						DO WHILE (ANY(acceptor_list(n_acceptor)%atom_index_list(:,2)>0))
							IF (acceptor_list(n_acceptor)%atom_index_list(acceptor_indexcounter,2)>0) THEN
								!This one is unassigned yet!
								IF (current_element=="NO") THEN
									!we have found a new victim
									group_number=group_number-1
									current_element=TRIM(give_element_symbol(acceptor_list(n_acceptor)%&
									&molecule_type_index,acceptor_list(n_acceptor)%&
									&atom_index_list(acceptor_indexcounter,1)))
									acceptor_list(n_acceptor)%atom_index_list(acceptor_indexcounter,2)=group_number
									IF (DEVELOPERS_VERSION) THEN
										WRITE(*,*)
										WRITE(*,ADVANCE="NO",FMT='("  !   group #",I0," (",A,") has atom indices ",I0)') -group_number,&
										&TRIM(current_element),acceptor_list(n_acceptor)%&
										&atom_index_list(acceptor_indexcounter,1)
									ENDIF
								ELSE
									!we already have a victim
									IF (TRIM(current_element)==TRIM(give_element_symbol(acceptor_list(n_acceptor)%&
										&molecule_type_index,acceptor_list(n_acceptor)%&
										&atom_index_list(acceptor_indexcounter,1)))) THEN
										acceptor_list(n_acceptor)%atom_index_list(acceptor_indexcounter,2)=group_number
										IF (DEVELOPERS_VERSION) THEN
											WRITE(*,ADVANCE="NO",FMT='(",",I0)') acceptor_list(n_acceptor)%&
											&atom_index_list(acceptor_indexcounter,1)
										ENDIF
									ENDIF
								ENDIF
							ENDIF
							IF (acceptor_indexcounter==acceptor_list(n_acceptor)%N_atoms) THEN
								acceptor_indexcounter=1
								current_element="NO"
							ELSE
								acceptor_indexcounter=acceptor_indexcounter+1
							ENDIF
							!the atom_index is:
						ENDDO
						acceptor_list(n_acceptor)%N_groups=-group_number
						IF (DEVELOPERS_VERSION) WRITE(*,*)
					ENDDO
					current_element="NO"
					DO n_donor=1,N_donor_types,1
						group_number=0
						IF (DEVELOPERS_VERSION) THEN
							WRITE(*,ADVANCE="NO",FMT='("  ! groups for donor #",I0,":")') n_donor
						ENDIF
						!set all group labels to be positive
						donor_list(n_donor)%atom_index_list(:,2)=1
						donor_indexcounter=1
						DO WHILE (ANY(donor_list(n_donor)%atom_index_list(:,2)>0))
							IF (donor_list(n_donor)%atom_index_list(donor_indexcounter,2)>0) THEN
								!This one is unassigned yet!
								IF (current_element=="NO") THEN
									!we have found a new victim
									group_number=group_number-1
									current_element=TRIM(give_element_symbol(donor_list(n_donor)%&
									&molecule_type_index,donor_list(n_donor)%&
									&atom_index_list(donor_indexcounter,1)))
									donor_list(n_donor)%atom_index_list(donor_indexcounter,2)=group_number
									IF (DEVELOPERS_VERSION) THEN
										WRITE(*,*)
										WRITE(*,ADVANCE="NO",FMT='("  !   group #",I0," (",A,") has atom indices ",I0)')&
										&-group_number,&
										&TRIM(current_element),donor_list(n_donor)%&
										&atom_index_list(donor_indexcounter,1)
									ENDIF
								ELSE
									!we already have a victim
									IF (TRIM(current_element)==TRIM(give_element_symbol(donor_list(n_donor)%&
										&molecule_type_index,donor_list(n_donor)%&
										&atom_index_list(donor_indexcounter,1)))) THEN
										donor_list(n_donor)%atom_index_list(donor_indexcounter,2)=group_number
										IF (DEVELOPERS_VERSION) THEN
											WRITE(*,ADVANCE="NO",FMT='(",",I0)') donor_list(n_donor)%&
											&atom_index_list(donor_indexcounter,1)
										ENDIF
									ENDIF
								ENDIF
							ENDIF
							IF (donor_indexcounter==donor_list(n_donor)%N_atoms) THEN
								donor_indexcounter=1
								current_element="NO"
							ELSE
								donor_indexcounter=donor_indexcounter+1
							ENDIF
							!the atom_index is:
						ENDDO
						donor_list(n_donor)%N_groups=-group_number
						IF (DEVELOPERS_VERSION) WRITE(*,*)
					ENDDO
				END SUBROUTINE group_element_names

				SUBROUTINE showgroups()
				IMPLICIT NONE
				INTEGER :: atom_counter,n_acceptor,n_donor,group_counter
				LOGICAL :: firsttime,no_ungrouped_atoms
					WRITE(*, FMT='(" For this analysis, the atoms will be grouped as follows.")')
					DO n_acceptor=1,N_acceptor_types,1
						IF (acceptor_list(n_acceptor)%N_atoms==1) THEN
							IF (N_acceptor_types==1) THEN
								WRITE(*,'("   Atom ",I0," (",A,") of molecule type ",I0," (",A,") was the sole acceptor.")')&
								&acceptor_list(n_acceptor)%atom_index_list(1,1),&
								&give_element_symbol(acceptor_list(n_acceptor)%molecule_type_index,&
								&acceptor_list(n_acceptor)%atom_index_list(1,1)),&
								&acceptor_list(n_acceptor)%molecule_type_index,&
								&TRIM(give_sum_formula(acceptor_list(n_acceptor)%molecule_type_index))
							ELSE
								WRITE(*,'("   Acceptor #",I0," was atom ",I0," (",A,") of molecule type ",I0," (",A,").")')&
								&n_acceptor,&
								&acceptor_list(n_acceptor)%atom_index_list(1,1),&
								&give_element_symbol(acceptor_list(n_acceptor)%molecule_type_index,&
								&acceptor_list(n_acceptor)%atom_index_list(1,1)),&
								&acceptor_list(n_acceptor)%molecule_type_index,&
								&TRIM(give_sum_formula(acceptor_list(n_acceptor)%molecule_type_index))
							ENDIF
							CYCLE
						ENDIF
						IF (acceptor_list(n_acceptor)%N_groups==0) THEN
							WRITE(*,FMT='("   Acceptor #",I0," has no custom groups defined.")')n_acceptor
						ELSE
							WRITE(*, FMT='("   Groups for acceptor #",I0,":")') n_acceptor
						ENDIF
						DO group_counter=1,acceptor_list(n_acceptor)%N_groups
							firsttime=.TRUE.
							DO atom_counter=1,acceptor_list(n_acceptor)%N_atoms,1
								IF (acceptor_list(n_acceptor)%atom_index_list(atom_counter,2)==-group_counter) THEN
									IF (firsttime) THEN
										firsttime=.FALSE.
										IF (grouped_by_elements) THEN
											WRITE(*,ADVANCE="NO",FMT='("     Group #",I0," (",A,") has atom indices ",I0)')&
											&group_counter,&
											&TRIM(give_element_symbol(acceptor_list(n_acceptor)%&
											&molecule_type_index,acceptor_list(n_acceptor)%&
											&atom_index_list(atom_counter,1))),&
											&acceptor_list(n_acceptor)%atom_index_list(atom_counter,1)
										ELSE
											WRITE(*,ADVANCE="NO",FMT='("     Group #",I0," has atom indices ",I0," (",A,")")') &
											&group_counter,acceptor_list(n_acceptor)%atom_index_list(atom_counter,1),&
											&TRIM(give_element_symbol(acceptor_list(n_acceptor)%&
											&molecule_type_index,acceptor_list(n_acceptor)%&
											&atom_index_list(atom_counter,1)))
										ENDIF
									ELSE
										IF (grouped_by_elements) THEN
											WRITE(*,ADVANCE="NO",FMT='(", ",I0)') acceptor_list(n_acceptor)%&
											&atom_index_list(atom_counter,1)
										ELSE
											WRITE(*,ADVANCE="NO",FMT='(", ",I0," (",A,")")')&
											&acceptor_list(n_acceptor)%atom_index_list(atom_counter,1),&
											&TRIM(give_element_symbol(acceptor_list(n_acceptor)%&
											&molecule_type_index,acceptor_list(n_acceptor)%&
											&atom_index_list(atom_counter,1)))
										ENDIF
									ENDIF
								ENDIF
							ENDDO
							IF (firsttime) THEN
								WRITE(*,FMT='("     Group #",I0," was empty (check input!).")')group_counter
							ELSE
								WRITE(*,'(".")')
							ENDIF
						ENDDO
						IF (.NOT.(grouped_by_elements)) THEN
							no_ungrouped_atoms=.TRUE.
							firsttime=.TRUE.
							DO atom_counter=1,acceptor_list(n_acceptor)%N_atoms,1
								IF (acceptor_list(n_acceptor)%atom_index_list(atom_counter,2)>0) THEN
									!sanity check
									IF (acceptor_list(n_acceptor)%atom_index_list(atom_counter,2)/=&
									&acceptor_list(n_acceptor)%atom_index_list(atom_counter,1)) CALL report_error(0)
									no_ungrouped_atoms=.FALSE.
									IF (firsttime) THEN
										WRITE(*,FMT='("     Atoms not assigned to a group, treated as their own group:")')
										WRITE(*,ADVANCE="NO",FMT='("       Atom indices ",I0)')&
										&acceptor_list(n_acceptor)%atom_index_list(atom_counter,1)
										firsttime=.FALSE.
									ELSE
										WRITE(*,ADVANCE="NO",FMT='(", ",I0)')&
										&acceptor_list(n_acceptor)%atom_index_list(atom_counter,1)
									ENDIF
								ENDIF
							ENDDO
							IF (no_ungrouped_atoms) THEN
								WRITE(*,'("     All atoms of this acceptor are assigned to a group.")')
							ELSE
								WRITE(*,'(".")')
							ENDIF
						ENDIF
					ENDDO
					DO n_donor=1,N_donor_types,1
						IF (donor_list(n_donor)%N_atoms==1) THEN
							IF (N_donor_types==1) THEN
								WRITE(*,'(" Atom ",I0," (",A,") of molecule type ",I0," (",A,") was the sole donor.")')&
								&donor_list(n_donor)%atom_index_list(1,1),&
								&give_element_symbol(donor_list(n_donor)%molecule_type_index,&
								&donor_list(n_donor)%atom_index_list(1,1)),&
								&donor_list(n_donor)%molecule_type_index,&
								&TRIM(give_sum_formula(donor_list(n_donor)%molecule_type_index))
							ELSE
								WRITE(*,'("   Donor #",I0," was atom ",I0," (",A,") of molecule type ",I0," (",A,").")')&
								&n_donor,&
								&donor_list(n_donor)%atom_index_list(1,1),&
								&give_element_symbol(donor_list(n_donor)%molecule_type_index,&
								&donor_list(n_donor)%atom_index_list(1,1)),&
								&donor_list(n_donor)%molecule_type_index,&
								&TRIM(give_sum_formula(donor_list(n_donor)%molecule_type_index))
							ENDIF
							CYCLE
						ENDIF
						IF (donor_list(n_donor)%N_groups==0) THEN
							WRITE(*,FMT='("   Donor #",I0," has no custom groups defined.")')n_donor
						ELSE
							WRITE(*, FMT='("   Groups for Donor #",I0,":")') n_donor
						ENDIF
						DO group_counter=1,donor_list(n_donor)%N_groups
							firsttime=.TRUE.
							DO atom_counter=1,donor_list(n_donor)%N_atoms,1
								IF (donor_list(n_donor)%atom_index_list(atom_counter,2)==-group_counter) THEN
									IF (firsttime) THEN
										firsttime=.FALSE.
										IF (grouped_by_elements) THEN
											WRITE(*,ADVANCE="NO",FMT='("     Group #",I0," (",A,") has atom indices ",I0)')&
											&group_counter,&
											&TRIM(give_element_symbol(donor_list(n_donor)%&
											&molecule_type_index,donor_list(n_donor)%&
											&atom_index_list(atom_counter,1))),&
											&donor_list(n_donor)%atom_index_list(atom_counter,1)
										ELSE
											WRITE(*,ADVANCE="NO",FMT='("     Group #",I0," has atom indices ",I0," (",A,")")') &
											&group_counter,donor_list(n_donor)%atom_index_list(atom_counter,1),&
											&TRIM(give_element_symbol(donor_list(n_donor)%&
											&molecule_type_index,donor_list(n_donor)%&
											&atom_index_list(atom_counter,1)))
										ENDIF
									ELSE
										IF (grouped_by_elements) THEN
											WRITE(*,ADVANCE="NO",FMT='(", ",I0)') donor_list(n_donor)%&
											&atom_index_list(atom_counter,1)
										ELSE
											WRITE(*,ADVANCE="NO",FMT='(", ",I0," (",A,")")')&
											&donor_list(n_donor)%atom_index_list(atom_counter,1),&
											&TRIM(give_element_symbol(donor_list(n_donor)%&
											&molecule_type_index,donor_list(n_donor)%&
											&atom_index_list(atom_counter,1)))
										ENDIF
									ENDIF
								ENDIF
							ENDDO
							IF (firsttime) THEN
								WRITE(*,FMT='("     Group #",I0," was empty (check input!).")')group_counter
							ELSE
								WRITE(*,'(".")')
							ENDIF
						ENDDO
						IF (.NOT.(grouped_by_elements)) THEN
							no_ungrouped_atoms=.TRUE.
							firsttime=.TRUE.
							DO atom_counter=1,donor_list(n_donor)%N_atoms,1
								IF (donor_list(n_donor)%atom_index_list(atom_counter,2)>0) THEN
									!sanity check
									IF (donor_list(n_donor)%atom_index_list(atom_counter,2)/=&
									&donor_list(n_donor)%atom_index_list(atom_counter,1)) CALL report_error(0)
									no_ungrouped_atoms=.FALSE.
									IF (firsttime) THEN
										WRITE(*,FMT='("     Atoms not assigned to a group, treated as their own group:")')
										WRITE(*,ADVANCE="NO",FMT='("       Atom indices ",I0)')&
										&donor_list(n_donor)%atom_index_list(atom_counter,1)
										firsttime=.FALSE.
									ELSE
										WRITE(*,ADVANCE="NO",FMT='(", ",I0)')&
										&donor_list(n_donor)%atom_index_list(atom_counter,1)
									ENDIF
								ENDIF
							ENDDO
							IF (no_ungrouped_atoms) THEN
								WRITE(*,'("     All atoms of this donor are assigned to a group.")')
							ELSE
								WRITE(*,'(".")')
							ENDIF
						ENDIF
					ENDDO
				END SUBROUTINE showgroups

				SUBROUTINE read_acceptors_and_donors()
				IMPLICIT NONE
				INTEGER :: ios,molecule_type_index,atom_index,N_molecule_types,counter
				INTEGER :: atom_counter,deallocstatus,failed,N_acceptors,N_donors
				REAL :: cutoff_clip
				LOGICAL :: molecule_type_exists
				TYPE(donors_and_acceptors),DIMENSION(:),ALLOCATABLE :: DA_clipboard ! here I collect the donors and acceptors while reading the file and keep track of numbers to later build the proper data structures
				!open scratch file. Structure of the scratch file is:
				!molecule_type_index - atom_index - cutoff - T/F for acceptor/donor
				INQUIRE(UNIT=10,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=10)
				OPEN(UNIT=10,STATUS="SCRATCH")
				!allocate clipboard memory
				IF (ALLOCATED(DA_clipboard)) CALL report_error(0)
				ALLOCATE(DA_clipboard(give_number_of_molecule_types()),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
				!rewind scratch file
				REWIND 10
				!initialise the counters to keep track of everything
				DA_clipboard(:)%molecule_type_index=0
				DA_clipboard(:)%N_atoms=0
				N_molecule_types=0
				READ(3,IOSTAT=ios,FMT=*) N_acceptors!read the number of acceptor atoms.
				IF (ios/=0) CALL report_error(155,exit_status=ios)
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,'(" reading ",I0," user-specified acceptor atoms...")') N_acceptors
				ENDIF
				failed=0
				!Try to read all the acceptors from the input file.
				DO n=1,N_acceptors,1
					READ(3,IOSTAT=ios,FMT=*) molecule_type_index,atom_index,cutoff_clip
					IF (ios/=0) CALL report_error(155,exit_status=ios)!ERROR 155: incorrect format in speciation.inp
					IF ((molecule_type_index>give_number_of_molecule_types()).OR.(molecule_type_index<1)) THEN
						!the specified molecule type doesn't exist. Stop execution.
						CALL report_error(33,exit_status=molecule_type_index)
						CLOSE(UNIT=3)
						CLOSE(UNIT=10)
						RETURN
					ENDIF
					IF ((atom_index>give_number_of_atoms_per_molecule(molecule_type_index))&
					&.OR.(atom_index<1)) THEN
						!the specified atom index doesn't exist. ignore line.
						CALL report_error(111,exit_status=atom_index)
						failed=failed+1
					ELSE
						!fill scratch file with neighbour infos - if the last atom_index was ok
						WRITE(10,*) molecule_type_index,atom_index,cutoff_clip
					ENDIF
					N_acceptors=N_acceptors-failed
					!check if molecule_type already there, keep track of number of atoms etc
					molecule_type_exists=.FALSE.
					!I am a little bit afraid of someone using the -onetrip but it should be fine, it will just run once with counter=1 and find a new molecule type still.
					DO counter=1,N_molecule_types,1
						IF (DA_clipboard(counter)%molecule_type_index==molecule_type_index) THEN
							!found existing molecule type
							DA_clipboard(counter)%N_atoms=DA_clipboard(counter)%N_atoms+1
							molecule_type_exists=.TRUE.
							EXIT
						ENDIF
					ENDDO
					IF (.NOT.(molecule_type_exists)) THEN
						!we have a new one or the first one!
						N_molecule_types=N_molecule_types+1
						DA_clipboard(counter)%molecule_type_index=molecule_type_index
						DA_clipboard(counter)%N_atoms=1
					ENDIF
				ENDDO
				!now we know the number of acceptor molecule types
				N_acceptor_types=N_molecule_types
				IF (N_acceptors<1) THEN
					!NO valid acceptors - abort analysis
					CALL report_error(159,exit_status=N_acceptors)
					CLOSE(UNIT=3)
					CLOSE(UNIT=10)
					DEALLOCATE(DA_clipboard,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
					RETURN
				ENDIF
				!now we now how many atoms in each molecule type. Thus, Allocate the memory.
				IF (ALLOCATED(acceptor_list)) CALL report_error(0)
				ALLOCATE(acceptor_list(N_molecule_types),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				bytes_needed=bytes_needed+N_molecule_types*4
				DO n=1,N_molecule_types,1
					acceptor_list(n)%molecule_type_index=DA_clipboard(n)%molecule_type_index
					acceptor_list(n)%N_atoms=DA_clipboard(n)%N_atoms
					IF (ALLOCATED(acceptor_list(n)%atom_index_list)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%atom_index_list(acceptor_list(n)%N_atoms,2),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					IF (ALLOCATED(acceptor_list(n)%cutoff_list)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%cutoff_list(acceptor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					bytes_needed=bytes_needed+acceptor_list(n)%N_atoms+acceptor_list(n)%N_atoms*2
					!initialise
					acceptor_list(n)%atom_index_list(:,:)=0
					acceptor_list(n)%cutoff_list(:)=0.0
					acceptor_list(n)%N_groups=0
					acceptor_list(n)%no_neighbour_occurrences=0
					acceptor_list(n)%species_overflow=0
				ENDDO
				REWIND 10
				DO n=1,N_acceptors,1
					!read the atom information from scratch file
					READ(10,*) molecule_type_index,atom_index,cutoff_clip
					!sort the atom information into the corresponding molecule type
					DO counter=1,N_molecule_types,1
						IF (acceptor_list(counter)%molecule_type_index==molecule_type_index) THEN
							!found the corresponding molecule type
							DO atom_counter=1,acceptor_list(counter)%N_atoms,1
								IF (acceptor_list(counter)%atom_index_list(atom_counter,1)==0) THEN
									!found an empty spot for the atom index
									acceptor_list(counter)%atom_index_list(atom_counter,:)=atom_index
									acceptor_list(counter)%cutoff_list(atom_counter)=cutoff_clip
									atom_index=0
								ENDIF
							ENDDO
						ENDIF
					ENDDO
				ENDDO
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,'(" Successfully read ",I0," acceptor atoms:")')N_acceptors
					DO n=1,N_molecule_types,1
						WRITE(*,FMT='("    Molecule type ",I0," (",A,"):")')&
						&acceptor_list(n)%molecule_type_index,TRIM(give_sum_formula(acceptor_list(n)%molecule_type_index))
						DO atom_counter=1,acceptor_list(n)%N_atoms,1
							WRITE(*,ADVANCE="NO",FMT='("      Atom index ",I0," (",A,"), cutoff = ")')&
							&acceptor_list(n)%atom_index_list(atom_counter,1),&
							&TRIM(give_element_symbol(acceptor_list(n)%molecule_type_index,acceptor_list(n)%&
							&atom_index_list(atom_counter,1)))
							IF (acceptor_list(n)%cutoff_list(atom_counter)<-0.1) THEN
								acceptor_list(n)%cutoff_list(atom_counter)=&
								&covalence_radius(&
								&TRIM(give_element_symbol(acceptor_list(n)%molecule_type_index,acceptor_list(n)%&
								&atom_index_list(atom_counter,1))))
								WRITE(*,'(F0.2," (covalence radius)")')acceptor_list(n)%cutoff_list(atom_counter)
							ELSE
								WRITE(*,'(F0.2)')acceptor_list(n)%cutoff_list(atom_counter)
							ENDIF
						ENDDO
					ENDDO
				ENDIF
				!rewind scratch file
				REWIND 10
				!initialise the counters to keep track of everything
				DA_clipboard(:)%molecule_type_index=0
				DA_clipboard(:)%N_atoms=0
				N_molecule_types=0
				READ(3,IOSTAT=ios,FMT=*) N_donors!read the number of donor atoms.
				IF (ios/=0) CALL report_error(155,exit_status=ios)
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,'(" reading ",I0," user-specified donor atoms...")') N_donors
				ENDIF
				failed=0
				!Try to read all the donors from the input file.
				DO n=1,N_donors,1
					READ(3,IOSTAT=ios,FMT=*) molecule_type_index,atom_index,cutoff_clip
					IF (ios/=0) CALL report_error(155,exit_status=ios)!ERROR 155: incorrect format in speciation.inp
					IF ((molecule_type_index>give_number_of_molecule_types()).OR.(molecule_type_index<1)) THEN
						!the specified molecule type doesn't exist. Stop execution.
						CALL report_error(33,exit_status=molecule_type_index)
						CLOSE(UNIT=3)
						CLOSE(UNIT=10)
						RETURN
					ENDIF
					IF ((atom_index>give_number_of_atoms_per_molecule(molecule_type_index))&
					&.OR.(atom_index<1)) THEN
						!the specified atom index doesn't exist. ignore line.
						CALL report_error(111,exit_status=atom_index)
						failed=failed+1
					ELSE
						!fill scratch file with neighbour infos - if the last atom_index was ok
						WRITE(10,*) molecule_type_index,atom_index,cutoff_clip
					ENDIF
					N_donors=N_donors-failed
					!check if molecule_type already there, keep track of number of atoms etc
					molecule_type_exists=.FALSE.
					!I am a little bit afraid of someone using the -onetrip but it should be fine, it will just run once with counter=1 and find a new molecule type still.
					DO counter=1,N_molecule_types,1
						IF (DA_clipboard(counter)%molecule_type_index==molecule_type_index) THEN
							!found existing molecule type
							DA_clipboard(counter)%N_atoms=DA_clipboard(counter)%N_atoms+1
							molecule_type_exists=.TRUE.
							EXIT
						ENDIF
					ENDDO
					IF (.NOT.(molecule_type_exists)) THEN
						!we have a new one or the first one!
						N_molecule_types=N_molecule_types+1
						DA_clipboard(counter)%molecule_type_index=molecule_type_index
						DA_clipboard(counter)%N_atoms=1
					ENDIF
				ENDDO
				!now we know the number of donor molecule types
				N_donor_types=N_molecule_types
				IF (N_donors<1) THEN
					!NO valid donors - abort analysis
					CALL report_error(159,exit_status=N_donors)
					CLOSE(UNIT=3)
					CLOSE(UNIT=10)
					DEALLOCATE(DA_clipboard,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
					RETURN
				ENDIF
				!now we now how many atoms in each molecule type. Thus, Allocate the memory.
				IF (ALLOCATED(donor_list)) CALL report_error(0)
				ALLOCATE(donor_list(N_molecule_types),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				bytes_needed=bytes_needed+N_molecule_types*4
				DO n=1,N_molecule_types,1
					donor_list(n)%molecule_type_index=DA_clipboard(n)%molecule_type_index
					donor_list(n)%N_atoms=DA_clipboard(n)%N_atoms
					IF (ALLOCATED(donor_list(n)%atom_index_list)) CALL report_error(0)
					ALLOCATE(donor_list(n)%atom_index_list(donor_list(n)%N_atoms,2),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					IF (ALLOCATED(donor_list(n)%cutoff_list)) CALL report_error(0)
					ALLOCATE(donor_list(n)%cutoff_list(donor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					bytes_needed=bytes_needed+donor_list(n)%N_atoms+donor_list(n)%N_atoms*2
					!initialise
					donor_list(n)%atom_index_list(:,:)=0
					donor_list(n)%cutoff_list(:)=0.0
					donor_list(n)%N_groups=0
				ENDDO
				REWIND 10
				DO n=1,N_donors,1
					!read the atom information from scratch file
					READ(10,*) molecule_type_index,atom_index,cutoff_clip
					!sort the atom information into the corresponding molecule type
					DO counter=1,N_molecule_types,1
						IF (donor_list(counter)%molecule_type_index==molecule_type_index) THEN
							!found the corresponding molecule type
							DO atom_counter=1,donor_list(counter)%N_atoms,1
								IF (donor_list(counter)%atom_index_list(atom_counter,1)==0) THEN
									!found an empty spot for the atom index
									donor_list(counter)%atom_index_list(atom_counter,:)=atom_index
									donor_list(counter)%cutoff_list(atom_counter)=cutoff_clip
									atom_index=0
								ENDIF
							ENDDO
						ENDIF
					ENDDO
				ENDDO
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,'(" Successfully read ",I0," donor atoms:")')N_donors
					DO n=1,N_molecule_types,1
						WRITE(*,FMT='("    Molecule type ",I0," (",A,"):")')&
						&donor_list(n)%molecule_type_index,TRIM(give_sum_formula(donor_list(n)%molecule_type_index))
						DO atom_counter=1,donor_list(n)%N_atoms,1
							WRITE(*,ADVANCE="NO",FMT='("      Atom index ",I0," (",A,"), cutoff = ")')&
							&donor_list(n)%atom_index_list(atom_counter,1),&
							&TRIM(give_element_symbol(donor_list(n)%molecule_type_index,donor_list(n)%&
							&atom_index_list(atom_counter,1)))
							IF (donor_list(n)%cutoff_list(atom_counter)<-0.1) THEN
								donor_list(n)%cutoff_list(atom_counter)=&
								&covalence_radius(&
								&TRIM(give_element_symbol(donor_list(n)%molecule_type_index,donor_list(n)%&
								&atom_index_list(atom_counter,1))))
								WRITE(*,'(F0.2," (covalence radius)")')donor_list(n)%cutoff_list(atom_counter)
							ELSE
								WRITE(*,'(F0.2)')donor_list(n)%cutoff_list(atom_counter)
							ENDIF
						ENDDO
					ENDDO
				ENDIF
				!deallocate the clipboard memory, which was used for both donors and acceptors
				DEALLOCATE(DA_clipboard,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				CLOSE(UNIT=10)
				END SUBROUTINE read_acceptors_and_donors

		END SUBROUTINE initialise_speciation

		!finalises the speciation module.
		SUBROUTINE finalise_speciation()
		IMPLICIT NONE
		INTEGER :: deallocstatus
			IF (ALLOCATED(donor_list)) THEN
				DEALLOCATE(donor_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(acceptor_list)) THEN
				DEALLOCATE(acceptor_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(squared_cutoff_list)) THEN
				DEALLOCATE(squared_cutoff_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
		END SUBROUTINE finalise_speciation

		SUBROUTINE trajectory_speciation_analysis_parallel()
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: stepcounter,acceptor_molecule_index,n_acceptor,n_donor,donor_indexcounter,ios,maxmol
		INTEGER :: donor_molecule_index,acceptor_indexcounter,allocstatus,deallocstatus,observed_species_number
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: all_observed_species_for_timestep
		LOGICAL :: existing_species,problematic_occurrence,neighbour_atom,neighbour_molecule,connected
		LOGICAL,DIMENSION(:),ALLOCATABLE :: potential_species_match !first dimension goes over all the species
		TYPE(species) :: species_clipboard
		!reset the number of species logged to zero just in case.
		acceptor_list(:)%number_of_logged_species_in_acceptor_list=0
		neighbour_atom_overflow=0
		!opening the units for the time series
		IF (calculate_autocorrelation) THEN
			maxmol=0
			DO n_acceptor=1,N_acceptor_types,1
				INQUIRE(UNIT=10+n_acceptor,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=10+n_acceptor)
				! WRITE(filename_time_series,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				! &//"acceptor",n_acceptor,"_time_series.dat"
				OPEN(UNIT=10+n_acceptor,STATUS="SCRATCH",IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				! WRITE(10+n_acceptor,'(" This file contains the species NUMBER (ordering is different to output!).")')
				! WRITE(10+n_acceptor,'(" The leftmost column is the time, after that one column per molecule index.")')
				! WRITE(10+n_acceptor,ADVANCE="NO",FMT='("timeline")')
				! DO acceptor_molecule_index=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor)%molecule_type_index),1
					! WRITE(10+n_acceptor,ADVANCE="NO",FMT='(" ",I0)')acceptor_molecule_index
				! ENDDO
				WRITE(10+n_acceptor,*)
				IF (maxmol<give_number_of_molecules_per_step(acceptor_list(n_acceptor)%molecule_type_index))&
				&maxmol=give_number_of_molecules_per_step(acceptor_list(n_acceptor)%molecule_type_index)
			ENDDO
		ENDIF
		!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL)))&
		!$OMP PRIVATE(species_clipboard,potential_species_match,stepcounter,n_acceptor,acceptor_molecule_index)&
		!$OMP PRIVATE(n_donor,donor_molecule_index,donor_indexcounter,acceptor_indexcounter,neighbour_atom,observed_species_number)&
		!$OMP PRIVATE(existing_species,problematic_occurrence,neighbour_molecule,all_observed_species_for_timestep)
		!$OMP SINGLE
	 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
	 !$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (speciation)"
	 !$ 	CALL timing_parallel_sections(.TRUE.)
	 !$ ENDIF
		IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
		!$OMP END SINGLE
		IF (calculate_autocorrelation) THEN
			IF (ALLOCATED(all_observed_species_for_timestep)) CALL report_error(0)
			ALLOCATE(all_observed_species_for_timestep(N_acceptor_types,maxmol),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		ENDIF
		!allocate memory for species clipboard
		IF (ALLOCATED(species_clipboard%connection)) CALL report_error(0)
		ALLOCATE(species_clipboard%connection(maximum_number_of_connections),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		!first_occurrence not needed!
		!allocate memory for the last occurrence
		IF (ALLOCATED(species_clipboard%last_occurrence)) CALL report_error(0)
		ALLOCATE(species_clipboard%last_occurrence(maximum_number_of_connections),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		!allocate memory for the starting indices of molecules.
		IF (ALLOCATED(species_clipboard%neighbour_molecule_starting_indices)) CALL report_error(0)
		ALLOCATE(species_clipboard%neighbour_molecule_starting_indices(maximum_number_of_connections),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		!allocate memory to keep track of possible matches
		IF (ALLOCATED(potential_species_match)) CALL report_error(0)
		ALLOCATE(potential_species_match(maximum_number_of_species),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
		!Iterate over all acceptors.
		!average over timesteps, a few should be enough.
		!$OMP DO SCHEDULE(STATIC,1) ORDERED
		DO stepcounter=1,nsteps,sampling_interval
			DO n_acceptor=1,N_acceptor_types,1
				!per ACCEPTOR molecule index, go over all its atoms and check if there are neighbour contributions.
				DO acceptor_molecule_index=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor)%molecule_type_index),1
					!everything inside THIS loop counts as neighbours towards one species!
					!Thus, at the end, we need to check if this species exist and update accordingly.
					!Furthermore, at this stage, we need to initialise the species_clipboard
					species_clipboard%total_number_of_neighbour_molecules=0
					species_clipboard%total_number_of_neighbour_atoms=0
					species_clipboard%number_of_logged_connections_in_species_list=0 !this corresponds to "species_clipboard%connection(:)%connection" being an empty list
					species_clipboard%total_number_of_neighbour_atoms=0 !every neighbour atom is counted once regardless of the number of connections.
					species_clipboard%total_number_of_neighbour_molecules=0 !every neighbour molecule is counted once regardless of the number of connections.
loop_donortypes:	DO n_donor=1,N_donor_types,1
						!Iterate over all possible neighbour (donor) molecules...
						DO donor_molecule_index=1,&
						&give_number_of_molecules_per_step(donor_list(n_donor)%molecule_type_index),1
							!now we are down to one specific donor molecule.
							!We need to check that we do not accidentally include intramolecular contributions.
							IF ((donor_list(n_donor)%molecule_type_index==acceptor_list(n_acceptor)%molecule_type_index)&
							&.AND.(donor_molecule_index==acceptor_molecule_index)) CYCLE
							!the inner loops from here on only run over atoms.
							neighbour_molecule=.FALSE.
							DO donor_indexcounter=1,donor_list(n_donor)%N_atoms,1
								neighbour_atom=.FALSE.
								DO acceptor_indexcounter=1,acceptor_list(n_acceptor)%N_atoms,1
									!Now we have a specific atom of the reference atom...
									!	... molecule type index is: acceptor_list(n_acceptor)%molecule_type_index
									!	... molecule index is: acceptor_molecule_index
									!	... atom index is: acceptor_list(n_acceptor)%atom_index_list(acceptor_indexcounter,1)
									!most donor molecules will not be neighbouring. Thus, for those non-neighbours, make sure the loop is exited with as little evaluations as possible
									!The inner loop from here on goes over different donor atoms.
									!calculate squared distance!
									!This is the part that really hurts...
									!Check if it is smaller than the cutoff
									IF (give_smallest_atom_distance_squared&
									&(stepcounter,stepcounter,&
									&acceptor_list(n_acceptor)%molecule_type_index,&
									&donor_list(n_donor)%molecule_type_index,&
									&acceptor_molecule_index,donor_molecule_index,&
									&acceptor_list(n_acceptor)%atom_index_list(acceptor_indexcounter,1),&
									&donor_list(n_donor)%atom_index_list(donor_indexcounter,1))<&
									&(squared_cutoff_list(&
									&n_acceptor,acceptor_indexcounter,n_donor,donor_indexcounter))) THEN
										!smaller than cutoff!!!
										IF (species_clipboard%number_of_logged_connections_in_species_list==maximum_number_of_connections) THEN
										!overflow
										!$OMP CRITICAL(overflow_updates)
										neighbour_atom_overflow=neighbour_atom_overflow+1
										!$OMP END CRITICAL(overflow_updates)
										ELSE
											neighbour_atom=.TRUE.
											species_clipboard%number_of_logged_connections_in_species_list=&
											species_clipboard%number_of_logged_connections_in_species_list+1
											!add this specific pair/connection to the list:
											species_clipboard%connection(species_clipboard%&
											&number_of_logged_connections_in_species_list)%indices(1)&
											&=n_donor
											species_clipboard%connection(species_clipboard%&
											&number_of_logged_connections_in_species_list)%indices(2)&
											&=donor_molecule_index
											species_clipboard%connection(species_clipboard%&
											&number_of_logged_connections_in_species_list)%indices(3)&
											&=donor_indexcounter
											species_clipboard%connection(species_clipboard%&
											&number_of_logged_connections_in_species_list)%indices(4)&
											&=acceptor_indexcounter !not directly the atom_index, but the entry in atom_index_list(x,1)
										ENDIF
									ENDIF
								ENDDO
								IF (neighbour_atom) THEN
									species_clipboard%total_number_of_neighbour_atoms=&
									&species_clipboard%total_number_of_neighbour_atoms+1
									neighbour_molecule=.TRUE.
								ENDIF
							ENDDO
							IF (neighbour_molecule) THEN
								species_clipboard%total_number_of_neighbour_molecules=&
								&species_clipboard%total_number_of_neighbour_molecules+1
							ENDIF
						ENDDO
					ENDDO loop_donortypes
					!The species clipboard is now filled with all the connections observed for this specific acceptor molecule.
					!Check if there was any neighbour at all
					IF (species_clipboard%number_of_logged_connections_in_species_list>0) THEN
						!There was at least one connection.
						!Remember the occurrence - because we will collapse indices
						species_clipboard%last_occurrence=species_clipboard%connection
						species_clipboard%last_occurrence_acceptorindices(1)=n_acceptor
						species_clipboard%last_occurrence_acceptorindices(2)=acceptor_molecule_index
						!$OMP CRITICAL(acceptor_list_access)
						CALL add_clipboard_to_list(species_clipboard,potential_species_match,n_acceptor,stepcounter,observed_species_number)
						!$OMP END CRITICAL(acceptor_list_access)
					ELSE
						observed_species_number=0
						!$OMP CRITICAL(acceptor_list_access)
						acceptor_list(n_acceptor)%no_neighbour_occurrences=&
						&acceptor_list(n_acceptor)%no_neighbour_occurrences+1
						!$OMP END CRITICAL(acceptor_list_access)
					ENDIF
					!This would be a good point to append the species corresponding to acceptor_molecule_index to the time series
					IF (calculate_autocorrelation) all_observed_species_for_timestep(n_acceptor,acceptor_molecule_index)=observed_species_number
				ENDDO
			ENDDO
			IF (VERBOSE_OUTPUT) CALL print_progress()
			IF (calculate_autocorrelation) THEN
				!$OMP ORDERED
				DO n_acceptor=1,N_acceptor_types,1
					WRITE(10+n_acceptor,ADVANCE="NO",FMT='(I0)')&
					&(stepcounter-1)*TIME_SCALING_FACTOR
					DO acceptor_molecule_index=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor)%molecule_type_index),1
						WRITE(10+n_acceptor,ADVANCE="NO",FMT='(" ",I0)')&
						&all_observed_species_for_timestep(n_acceptor,acceptor_molecule_index)
					ENDDO
					WRITE(10+n_acceptor,*)
				ENDDO
				!$OMP END ORDERED
			ENDIF
		ENDDO
		!$OMP END DO
		IF (calculate_autocorrelation) THEN
			IF (ALLOCATED(all_observed_species_for_timestep)) THEN
				DEALLOCATE(all_observed_species_for_timestep,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
		ENDIF
		IF (ALLOCATED(species_clipboard%connection)) THEN
			DEALLOCATE(species_clipboard%connection,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(potential_species_match)) THEN
			DEALLOCATE(potential_species_match,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(species_clipboard%first_occurrence)) CALL report_error(0)
		IF (ALLOCATED(species_clipboard%last_occurrence)) THEN
			DEALLOCATE(species_clipboard%last_occurrence,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(species_clipboard%neighbour_molecule_starting_indices)) THEN
			DEALLOCATE(species_clipboard%neighbour_molecule_starting_indices,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		!$OMP END PARALLEL
		IF (((MAX((nsteps-1+sampling_interval)/sampling_interval,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
	 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
	 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
	 !$ 	CALL timing_parallel_sections(.FALSE.)
	 !$ ENDIF
		END SUBROUTINE trajectory_speciation_analysis_parallel

		SUBROUTINE add_clipboard_to_list(species_clipboard,potential_species_match,n_acceptor_in,stepcounter_in,observed_species_number)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: n_acceptor_in,stepcounter_in
		INTEGER,INTENT(OUT) :: observed_species_number
		TYPE(species),INTENT(INOUT) :: species_clipboard
		LOGICAL,DIMENSION(*),INTENT(INOUT) :: potential_species_match
		LOGICAL :: new_species_found
			!First we need to make sure that atoms that are considered identical sites are treated as such.
			!This needs to be called here already, it is necessary in ANY case, regardless whether new species or not.
			!collapse_and_sort_atom_indices also assigns neighbour_molecule_starting_indices
			CALL collapse_and_sort_atom_indices(species_clipboard,n_acceptor_in)
			!Now to find out if this species already exists in the list.
			!There are many possible different cases and conditions, and I will use a jump to skip the ones which are not necessary.
			!I think in this case, the jump command is better and more easily understood than 1000 nestled if-THEN-else
			!Very first case: species list is empty. need to add.
			IF (acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list==0) THEN
				new_species_found=.TRUE.
				GOTO 60
			ELSE
				!There is already at least one entry - we now need to find out if there is an identical species already
				new_species_found=.FALSE.
			ENDIF
			!The rest is really not trivial. We will use a list of potential matches, and eliminate them one after one.
			potential_species_match(1:acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list)=.TRUE.
			!This is the first rough check we can do - check if species have the same number of connections, number of neighbour atoms, number of neighbour molecules
			CALL eliminate_number_mismatches()
			!If there are no potential species matches left, THEN this is a new one!
			IF (.NOT.(ANY(&
			&potential_species_match(1:acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list)&
			&))) THEN
				new_species_found=.TRUE.
				GOTO 60
			ENDIF
			!We now need to check all the molecules one by one,
			!and also check if the atoms which are involved in these connections are the same.
			!Molecule indices do not matter, so we need to check for permutations...
			CALL eliminate_connection_mismatches()
			!If there are no potential species matches left, THEN this is a new one!
			IF (.NOT.(ANY(&
			&potential_species_match(1:acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list)&
			&))) THEN
				new_species_found=.TRUE.
				GOTO 60
			ENDIF
		60	IF (new_species_found) THEN
				CALL append_species()
			ELSE
				IF (COUNT(potential_species_match(1:acceptor_list(n_acceptor_in)%&
				&number_of_logged_species_in_acceptor_list))/=1) CALL report_error(0)
				CALL add_to_species()
			ENDIF

		CONTAINS

			SUBROUTINE eliminate_connection_mismatches()
			IMPLICIT NONE
			INTEGER :: species_counter,clipboard_molecules,reference_species_molecules
			INTEGER :: connections_firstindex_ref,extra_atom_entries_ref
			INTEGER :: connections_firstindex_clip,extra_atom_entries_clip
			INTEGER :: connection_counter
			LOGICAL :: molecule_mismatch,successful_species_match
loop_species: 	DO species_counter=1,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list,1
					successful_species_match=.FALSE.
					!only check potential matches, of course.
					IF (.NOT.(potential_species_match(species_counter))) CYCLE loop_species
					!Now, we need to check the possible permutations of clipboard and the current entry in the list_of_all_species
					! initialise the "lists to keep track" of pairs which are ok
					species_clipboard%neighbour_molecule_starting_indices(1:species_clipboard%&
					&total_number_of_neighbour_molecules)%unmatched_entry=.TRUE.
					acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
					&neighbour_molecule_starting_indices(1:acceptor_list(n_acceptor_in)%&
					&list_of_all_species(species_counter)%total_number_of_neighbour_molecules)%unmatched_entry=.TRUE.
loop_clipboard:		DO clipboard_molecules=1,species_clipboard%total_number_of_neighbour_molecules
						IF (.NOT.(species_clipboard%neighbour_molecule_starting_indices(clipboard_molecules)%unmatched_entry)) CYCLE loop_clipboard
						!the number of molecules should be the same-This has been checked before by eliminate_number_mismatches
loop_reference:			DO reference_species_molecules=1,species_clipboard%total_number_of_neighbour_molecules
							IF ((.NOT.(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
							&neighbour_molecule_starting_indices(reference_species_molecules)%unmatched_entry)).OR.&
							&(.NOT.(species_clipboard%neighbour_molecule_starting_indices(clipboard_molecules)%unmatched_entry))) CYCLE loop_reference
							!At this point, we have selected two pairs of molecules, one from the clipboard and one from the reference.
							!The idea is the following:
							! 1) we assume that they are the same, i.e. that there is no mismatche.
							! 2) we check all the connections.
							! 3) if there is ANY connection which differs in ANY of the indices 1/3/4, THEN they are NOT the same.
							! 4) if they are still the same after checking all connections (in other words, no mismatches were found),
							! 5) THEN we set the unmatched_entry of this pair to .FALSE. in both "lists to keept track"
							molecule_mismatch=.FALSE. ! step 1)
							!get the first indices of the respective molecules
							connections_firstindex_ref=acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)&
							&%neighbour_molecule_starting_indices(reference_species_molecules)%first_index_in_connection_list
							connections_firstindex_clip=species_clipboard%neighbour_molecule_starting_indices&
							&(clipboard_molecules)%first_index_in_connection_list
							extra_atom_entries_ref=acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
							&neighbour_molecule_starting_indices(reference_species_molecules)%&
							&extra_atom_entries_in_connection_list
							extra_atom_entries_clip=species_clipboard%neighbour_molecule_starting_indices&
							&(clipboard_molecules)%extra_atom_entries_in_connection_list
							!can only be the same if they have the same number of connections
							IF (extra_atom_entries_clip/=extra_atom_entries_ref) THEN
								molecule_mismatch=.TRUE.
							ELSE
loop_connections:				DO connection_counter=0,extra_atom_entries_ref ! step 2) 
									!The two molecules we want to compare have the same number of connection entries,
									!But they differ in their starting position. connection_counter is the shift from that starting position.
									IF ((acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
										&connection(connection_counter+connections_firstindex_ref)%indices(1)&
										&/=&
										&species_clipboard%connection(connection_counter+connections_firstindex_clip)&
										&%indices(1))&
										&.OR.&
										&(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
										&connection(connection_counter+connections_firstindex_ref)%indices(3)&
										&/=&
										&species_clipboard%connection(connection_counter+connections_firstindex_clip)&
										&%indices(3))&
										&.OR.&
										&(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
										&connection(connection_counter+connections_firstindex_ref)%indices(4)&
										&/=&
										&species_clipboard%connection(connection_counter+connections_firstindex_clip)&
										&%indices(4))) THEN
											molecule_mismatch=.TRUE.
											EXIT loop_connections
										ENDIF
								ENDDO loop_connections
							ENDIF
							! step 4)
							IF (.NOT.(molecule_mismatch)) THEN
								! step 5) The pair was a match! remove this pair from the two lists.
								species_clipboard%neighbour_molecule_starting_indices(clipboard_molecules)%unmatched_entry=.FALSE.
								acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
								&neighbour_molecule_starting_indices(reference_species_molecules)%unmatched_entry=.FALSE.
								IF ((.NOT.(ANY(species_clipboard%neighbour_molecule_starting_indices(1:species_clipboard%&
								&total_number_of_neighbour_molecules)%unmatched_entry))).AND.&
								&(.NOT.(ANY(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
								&neighbour_molecule_starting_indices(1:acceptor_list(n_acceptor_in)%&
								&list_of_all_species(species_counter)%total_number_of_neighbour_molecules)%unmatched_entry)))) THEN
									!All have been matched!
									successful_species_match=.TRUE.
									EXIT loop_clipboard
								ENDIF
							ENDIF
						ENDDO loop_reference
					ENDDO loop_clipboard
					potential_species_match(species_counter)=successful_species_match
				ENDDO loop_species
			END SUBROUTINE eliminate_connection_mismatches

			SUBROUTINE eliminate_number_mismatches()
			IMPLICIT NONE
			INTEGER :: species_counter
!This can be replaced by a fancy "WHERE" statement! do this last when I have a working routine with all parts so I can compare performance...
				!We do not need an initial check of potential_species_match because this subroutine is called first, after the initialisation of potential_species_match
				DO species_counter=1,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list,1
					IF ((species_clipboard%number_of_logged_connections_in_species_list)/=&
					&(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
					&number_of_logged_connections_in_species_list)) THEN
						potential_species_match(species_counter)=.FALSE.
					ELSE
						IF ((species_clipboard%total_number_of_neighbour_atoms)/=&
						&(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
						&total_number_of_neighbour_atoms)) THEN
							potential_species_match(species_counter)=.FALSE.
						ELSE
							IF ((species_clipboard%total_number_of_neighbour_molecules)/=&
							&(acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%&
							&total_number_of_neighbour_molecules)) THEN
								potential_species_match(species_counter)=.FALSE.
							ENDIF
						ENDIF
					ENDIF
				ENDDO
			END SUBROUTINE eliminate_number_mismatches

			SUBROUTINE append_species()!appends the contents of the clipboard to the actual list
			IMPLICIT NONE
				!check for overflow
				IF (acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list<maximum_number_of_species) THEN
					acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list=&
					&acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list+1
					observed_species_number=acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
					!append the new species to the list!
					!transfer the connections and their count
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%number_of_logged_connections_in_species_list=&
					&species_clipboard%number_of_logged_connections_in_species_list
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%connection=&
					&species_clipboard%connection
					!initialise occurrences
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%occurrences=1
					!initialise the last/first occurrence - stepcounter
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%timestep_last_occurrence=stepcounter_in
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%timestep_first_occurrence=stepcounter_in
					!initialise the last/first occurrence - indices of the acceptor
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%last_occurrence_acceptorindices(:)=&
					&species_clipboard%last_occurrence_acceptorindices
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%first_occurrence_acceptorindices(:)=&
					&species_clipboard%last_occurrence_acceptorindices
					!transfer the last_occurrence
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%first_occurrence=species_clipboard%last_occurrence
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%last_occurrence=species_clipboard%last_occurrence
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%total_number_of_neighbour_atoms=&
					&species_clipboard%total_number_of_neighbour_atoms
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%total_number_of_neighbour_molecules&
					&=species_clipboard%total_number_of_neighbour_molecules
					!update the starting indices
					acceptor_list(n_acceptor_in)%list_of_all_species(acceptor_list(n_acceptor_in)%&
					&number_of_logged_species_in_acceptor_list)%neighbour_molecule_starting_indices=&
					&species_clipboard%neighbour_molecule_starting_indices
				ELSE
					acceptor_list(n_acceptor_in)%species_overflow=acceptor_list(n_acceptor_in)%species_overflow+1
					observed_species_number=0
				ENDIF
			END SUBROUTINE append_species

			SUBROUTINE add_to_species()!appends the contents of the clipboard to the actual list
			IMPLICIT NONE
			INTEGER :: species_counter,connection_counter
				DO species_counter=1,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list,1
					IF (potential_species_match(species_counter)) THEN
						observed_species_number=species_counter
						acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%occurrences=&
						&acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%occurrences+1
						!update the last occurrence - indices of the acceptor
						acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%timestep_last_occurrence=stepcounter_in
						acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%last_occurrence_acceptorindices(:)=&
						&species_clipboard%last_occurrence_acceptorindices
						acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%last_occurrence=species_clipboard%last_occurrence
					ENDIF
				ENDDO
			END SUBROUTINE add_to_species

			!The following subroutine sorts connections_inout in ascending order from firstindex to lastindex.
			!Here, the entry in position sorting_level is used for comparison.
			SUBROUTINE sort_connections(connections_inout,sorting_level,firstindex,lastindex)
			IMPLICIT NONE
			TYPE(connections),DIMENSION(*),INTENT(INOUT) :: connections_inout
			INTEGER,INTENT(IN) :: sorting_level,firstindex,lastindex
			TYPE(connections) :: sorting_clipboard
			INTEGER :: i,j
				!Now, we need to sort the elements from "firstindex" to "lastindex"
				!I chose insertionsort over quicksort for stability
				DO i=firstindex+1,lastindex,1
					!might be worth using just the indices here?
					sorting_clipboard=connections_inout(i)
					j=i
					DO WHILE ((j>firstindex).AND.&
					&(connections_inout(j-1)%indices(sorting_level)>&
					&sorting_clipboard%indices(sorting_level)))
						connections_inout(j)=connections_inout(j-1)
						j=j-1
					ENDDO
					connections_inout(j)=sorting_clipboard
				ENDDO
			END SUBROUTINE sort_connections

			SUBROUTINE collapse_and_sort_atom_indices(species_clipboard,n_acceptor_in)
			IMPLICIT NONE
			INTEGER :: current_donor_type,current_donor_molecule_index,connections_firstindex,connections_lastindex
			INTEGER :: current_donor_atom,atoms_firstindex,atoms_lastindex,molecule_counter
			INTEGER :: connection_counter,connection_counter_atoms
			INTEGER,INTENT(IN) :: n_acceptor_in
			TYPE(species),INTENT(INOUT) :: species_clipboard
				!initialise neighbour_molecule_starting_indices
				species_clipboard%neighbour_molecule_starting_indices(:)%first_index_in_connection_list=0
				species_clipboard%neighbour_molecule_starting_indices(:)%extra_atom_entries_in_connection_list=0
				species_clipboard%neighbour_molecule_starting_indices(:)%unmatched_entry=.FALSE.
				molecule_counter=1
				!At this point the natural sorting will be intact.
				!Convert atom indices to atom groups - this does not affect the sorting of donor type and donor molecule index.
				DO connection_counter=1,species_clipboard%number_of_logged_connections_in_species_list
					species_clipboard%connection(connection_counter)%indices(3)=&
					&donor_list(species_clipboard%connection(connection_counter)%indices(1))%&
					&atom_index_list(species_clipboard%connection(connection_counter)%indices(3),2)
					species_clipboard%connection(connection_counter)%indices(4)=&
					&acceptor_list(n_acceptor_in)%&
					&atom_index_list(species_clipboard%connection(connection_counter)%indices(4),2)
				ENDDO
				!Now, the atom indices should be converted to atom groups. let's sort them.
				!loop over all entries
				current_donor_type=species_clipboard%connection(1)%indices(1)
				current_donor_molecule_index=species_clipboard%connection(1)%indices(2)
				connections_firstindex=1
				DO connection_counter=1,species_clipboard%number_of_logged_connections_in_species_list-1
					IF ((species_clipboard%connection(connection_counter+1)%indices(1)/=current_donor_type).OR.&
					&(species_clipboard%connection(connection_counter+1)%indices(2)/=current_donor_molecule_index)) THEN !new molecule entry condition
						!the next entry will be a new molecule. Hence, this entry is the last index of the current one.
						connections_lastindex=connection_counter
						!update the index position
						species_clipboard%neighbour_molecule_starting_indices(molecule_counter)%first_index_in_connection_list=&
						&connections_firstindex
						species_clipboard%neighbour_molecule_starting_indices(molecule_counter)%extra_atom_entries_in_connection_list=&
						&connections_lastindex-connections_firstindex
						molecule_counter=molecule_counter+1
						!sort donor atoms from connections_firstindex to connections_lastindex.
						CALL sort_connections(species_clipboard%connection(:),3,connections_firstindex,connections_lastindex)
						!THEN, sort acceptor atoms - to be able to do so, we need to find the first and last donor atoms...
						current_donor_atom=species_clipboard%connection(connections_firstindex)%indices(3)
						atoms_firstindex=connections_firstindex
						DO connection_counter_atoms=connections_firstindex,connections_lastindex-1
							IF (species_clipboard%connection(connection_counter_atoms+1)%indices(3)/=current_donor_atom) THEN !new atom entry condition
								atoms_lastindex=connection_counter_atoms
								CALL sort_connections(species_clipboard%connection(:),4,atoms_firstindex,atoms_lastindex)
								!update the first index etc of the next atom
								atoms_firstindex=connection_counter_atoms+1
								current_donor_atom=species_clipboard%connection(connection_counter_atoms+1)%indices(3)
							ENDIF
						ENDDO
						!the very last atom will be missing!
						!hence, do the part inside the "new atom entry condition" again
						atoms_lastindex=connections_lastindex
						CALL sort_connections(species_clipboard%connection(:),4,atoms_firstindex,atoms_lastindex)
						!finally, update the first index etc of the next molecule.
						connections_firstindex=connection_counter+1
						current_donor_type=species_clipboard%connection(connection_counter+1)%indices(1)
						current_donor_molecule_index=species_clipboard%connection(connection_counter+1)%indices(2)
					ENDIF
				ENDDO
				!the very last molecule will be missing!
				!Hence, do the part inside the "new molecule entry condition" above
				connections_lastindex=species_clipboard%number_of_logged_connections_in_species_list
				species_clipboard%neighbour_molecule_starting_indices(molecule_counter)%first_index_in_connection_list=&
				&connections_firstindex
				IF (molecule_counter/=species_clipboard%total_number_of_neighbour_molecules) CALL report_error(0)
				species_clipboard%neighbour_molecule_starting_indices(molecule_counter)%extra_atom_entries_in_connection_list=&
				&connections_lastindex-connections_firstindex
				!sort donor atoms from connections_firstindex to connections_lastindex.
				CALL sort_connections(species_clipboard%connection(:),3,connections_firstindex,connections_lastindex)	
				!THEN, sort acceptor atoms - to be able to do so, we need to find the first and last donor atoms...
				current_donor_atom=species_clipboard%connection(connections_firstindex)%indices(3)
				atoms_firstindex=connections_firstindex
				DO connection_counter_atoms=connections_firstindex,connections_lastindex-1
					IF (species_clipboard%connection(connection_counter_atoms+1)%indices(3)/=current_donor_atom)THEN !new atom entry condition
						atoms_lastindex=connection_counter_atoms
						CALL sort_connections(species_clipboard%connection(:),4,atoms_firstindex,atoms_lastindex)
						!update the first index etc of the next atom
						atoms_firstindex=connection_counter_atoms+1
						current_donor_atom=species_clipboard%connection(connection_counter_atoms+1)%indices(3)
					ENDIF
				ENDDO
				!the very last atom will be missing!
				!hence, do the part inside the "new atom entry condition" again
				atoms_lastindex=connections_lastindex
				CALL sort_connections(species_clipboard%connection(:),4,atoms_firstindex,atoms_lastindex)
			END SUBROUTINE collapse_and_sort_atom_indices

		END SUBROUTINE add_clipboard_to_list

		SUBROUTINE print_acceptor_summary(n_acceptor_in)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: n_acceptor_in
		INTEGER :: entry_counter,allocstatus,deallocstatus,best_index,species_output_counter
		INTEGER :: i,previous_donor_molecule,donor_type_counter,m,double_connections_counter
		INTEGER :: species_atomcount,ios,species_charge,group_number,connection_counter,N_beads,N_beads_total
		INTEGER :: acceptorgroup,donorgroup,connections_firstindex,extra_atom_entries,species_molecules
		INTEGER :: species_molecules_doubles,extra_atom_entries_2,donorgroup_2,acceptorgroup_2,connections_firstindex_2
		INTEGER :: donortype_2,donortype
		LOGICAL :: connected,no_ungrouped_atoms,firsttime
		REAL :: total_occurrence,best_value,distance,base_atom(3),tip_atom(3),connection_vector(3),beadspan
		REAL,PARAMETER :: bead_clearance=0.5,bead_distance=0.4
		REAL(KIND=WORKING_PRECISION) :: acceptor_centre(3),shift(3)
		REAL,DIMENSION(:),ALLOCATABLE :: occurrences
		CHARACTER(LEN=1024) :: filename_speciation_output
		CHARACTER(LEN=1024) :: filename_speciation_statistics
		CHARACTER(LEN=1024) :: species_output_header
		CHARACTER(LEN=2) :: element_name
			WRITE(filename_speciation_statistics,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
			&//"acceptor",n_acceptor_in,"_speciation_statistics.dat"
			OPEN(UNIT=4,FILE=TRIM(filename_speciation_statistics),IOSTAT=ios)
			IF (ios/=0) CALL report_error(26,exit_status=ios)
			WRITE(4,'(" This file contains speciation statistics for acceptors coordinated by donors.")')
			IF (acceptor_list(n_acceptor_in)%N_atoms==1) THEN
				WRITE(4,'(" This acceptor was atom ",I0," (",A,") of molecule type ",I0," (",A,").")')&
				&acceptor_list(n_acceptor_in)%atom_index_list(1,1),&
				&give_element_symbol(acceptor_list(n_acceptor_in)%molecule_type_index,&
				&acceptor_list(n_acceptor_in)%atom_index_list(1,1)),&
				&acceptor_list(n_acceptor_in)%molecule_type_index,&
				&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
			ELSE
				IF (atoms_grouped) THEN
					WRITE(4,'(" Molecule type ",I0," (",A,") was the acceptor.")')&
					&acceptor_list(n_acceptor_in)%molecule_type_index,&
					&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
					IF (acceptor_list(n_acceptor_in)%N_groups==0) THEN
						WRITE(4,FMT='("   This acceptor (#",I0,") has no custom groups defined.")')n_acceptor_in
					ELSE
						WRITE(4, FMT='("   Groups for this acceptor (#",I0,"):")') n_acceptor_in
					ENDIF
					DO group_number=1,acceptor_list(n_acceptor_in)%N_groups,1
						no_ungrouped_atoms=.TRUE.
						firsttime=.TRUE.
						DO i=1,acceptor_list(n_acceptor_in)%N_atoms,1
							IF (group_number==-acceptor_list(n_acceptor_in)%atom_index_list(i,2)) THEN
								IF (firsttime) THEN
									firsttime=.FALSE.
									WRITE(4,ADVANCE="NO",FMT='("     Group #",I0,", atom indices ",I0,"(",A,")")')&
									&group_number,acceptor_list(n_acceptor_in)%atom_index_list(i,1),&
									&TRIM(give_element_symbol(acceptor_list(n_acceptor_in)%molecule_type_index,&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1)))
								ELSE
									WRITE(4,ADVANCE="NO",FMT='(", ",I0,"(",A,")")')&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1),&
									&TRIM(give_element_symbol(acceptor_list(n_acceptor_in)%molecule_type_index,&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1)))
								ENDIF
							ENDIF
						ENDDO
						IF (firsttime) THEN
							WRITE(4,FMT='("     Group #",I0," was empty.")')group_number
						ELSE
							WRITE(4,'(".")')
						ENDIF
					ENDDO
					IF (.NOT.(grouped_by_elements)) THEN
						no_ungrouped_atoms=.TRUE.
						firsttime=.TRUE.
						DO i=1,acceptor_list(n_acceptor_in)%N_atoms,1
							IF (acceptor_list(n_acceptor_in)%atom_index_list(i,2)>0) THEN
								!sanity check
								IF (acceptor_list(n_acceptor_in)%atom_index_list(i,2)/=&
								&acceptor_list(n_acceptor_in)%atom_index_list(i,1)) CALL report_error(0)
								no_ungrouped_atoms=.FALSE.
								IF (firsttime) THEN
									WRITE(4,FMT='("     Atoms not assigned to a group, treated as their own group:")')
									WRITE(4,ADVANCE="NO",FMT='("       Atom indices ",I0)')&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1)
									firsttime=.FALSE.
								ELSE
									WRITE(4,ADVANCE="NO",FMT='(", ",I0)')&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1)
								ENDIF
							ENDIF
						ENDDO
						IF (no_ungrouped_atoms) THEN
							WRITE(4,'("     All atoms of this acceptor were assigned to a group.")')
						ELSE
							WRITE(4,'(".")')
						ENDIF
					ENDIF
				ELSE
					WRITE(4,'(" Molecule type ",I0," (",A,") was the acceptor, with the following atom indices:")')&
					&acceptor_list(n_acceptor_in)%molecule_type_index,&
					&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
					DO i=1,acceptor_list(n_acceptor_in)%N_atoms,1
						WRITE(4,FMT='("      Atom index ",I0," (",A,")")')&
						&acceptor_list(n_acceptor_in)%atom_index_list(i,1),&
						&TRIM(give_element_symbol(acceptor_list(n_acceptor_in)%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%atom_index_list(i,1)))
					ENDDO
				ENDIF
			ENDIF
			WRITE(4,'(" The donor molecule types and atom types are given at the bottom.")')
			WRITE(4,'(" For this acceptor, ",I0," different species were observed,")')&
			&acceptor_list(n_acceptor_in)%Number_of_logged_species_in_acceptor_list
			WRITE(4,'(" Which are now given in order of decreasing probability:")')
			IF (ALLOCATED(occurrences)) CALL report_error(0)
			ALLOCATE(occurrences(acceptor_list(n_acceptor_in)%&
			&Number_of_logged_species_in_acceptor_list),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
			WRITE(*,'("molecule type ",I0," (",A,").")')&
			&acceptor_list(n_acceptor_in)%molecule_type_index,&
			&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
			WRITE(*,'("   For this acceptor, ",I0," different species were observed:")')&
			&acceptor_list(n_acceptor_in)%Number_of_logged_species_in_acceptor_list
			total_occurrence=0.0
			DO entry_counter=1,acceptor_list(n_acceptor_in)%Number_of_logged_species_in_acceptor_list,1
				occurrences(entry_counter)=FLOAT(&
				&acceptor_list(n_acceptor_in)%list_of_all_species(entry_counter)%occurrences)
				total_occurrence=total_occurrence+occurrences(entry_counter)
			ENDDO
			total_occurrence=total_occurrence+FLOAT(acceptor_list(n_acceptor_in)%species_overflow)
			total_occurrence=total_occurrence+FLOAT(acceptor_list(n_acceptor_in)%no_neighbour_occurrences)
			IF (calculate_autocorrelation) THEN
				!save the average_h for later. basically the same as the occurrences
				IF (ALLOCATED(average_h)) CALL report_error(0)
				ALLOCATE(average_h(0:acceptor_list(n_acceptor_in)%&
				&Number_of_logged_species_in_acceptor_list),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
				average_h(0)=DFLOAT(acceptor_list(n_acceptor_in)%species_overflow+&
				&acceptor_list(n_acceptor_in)%no_neighbour_occurrences)/total_occurrence
				DO entry_counter=1,acceptor_list(n_acceptor_in)%Number_of_logged_species_in_acceptor_list,1
					average_h(entry_counter)=&
					&DFLOAT(acceptor_list(n_acceptor_in)%list_of_all_species(entry_counter)%occurrences)/total_occurrence
				ENDDO
			ENDIF
			!Now, print species in order of decreasing occurrence
			!I am going to brute force this since I do not want to sort them.
			species_output_counter=0
			DO WHILE (ANY(occurrences(:)>0.0))
				!There is an unused entry left!
				species_output_counter=species_output_counter+1
				!laboriously pick out the best one and report it
				best_index=1
				best_value=0.0
				DO entry_counter=1,acceptor_list(n_acceptor_in)%&
				&Number_of_logged_species_in_acceptor_list,1
					IF (occurrences(entry_counter)>best_value) THEN
						best_value=occurrences(entry_counter)
						best_index=entry_counter
					ENDIF
				ENDDO
				acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%rank=species_output_counter
				!Remove the best_index
				occurrences(best_index)=-1.0
				!Report the best index
				WRITE(4,'("   -------")')
				IF (FLOAT(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&occurrences)/total_occurrence>0.01) THEN
					WRITE(4,'("   Species #(",I0,") in ",F0.1,"% of the cases.")')&
					&species_output_counter,100.0*FLOAT(&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)/total_occurrence
				ELSE
					WRITE(4,'("     Species #(",I0,") in ",I0," cases (less than 1%).")')&
					&species_output_counter,&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences
				ENDIF
				IF (species_output_counter<6) THEN
					IF (FLOAT(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&occurrences)/total_occurrence>0.01) THEN
						WRITE(*,'("     Species #(",I0,") in ",F0.1,"% of the cases.")')&
						&species_output_counter,100.0*FLOAT(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)/total_occurrence
					ELSE
						WRITE(*,'("     Species #(",I0,") in ",I0," cases (less than 1%).")')&
						&species_output_counter,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences
					ENDIF
				ELSE
					IF (species_output_counter==6) WRITE(*,'("     (Only the first 5 printed here)")')
				ENDIF
				WRITE(4,'("   This species has a total of ",I0," connections,")')&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&number_of_logged_connections_in_species_list
				WRITE(4,'("   with a total of ",I0," neighbour atoms in ",I0," neighbour molecules:")')&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&total_number_of_neighbour_atoms,&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&total_number_of_neighbour_molecules
				!calculate the formal total charge
				!Initialise with the acceptor charge
				species_charge=give_charge_of_molecule(&
				&acceptor_list(n_acceptor_in)%molecule_type_index)
				!Make sure we report similar entries together
				acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
				&%neighbour_molecule_starting_indices(:)%unmatched_entry=.TRUE.
				DO species_molecules=1,acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
				&%total_number_of_neighbour_molecules
					connections_firstindex=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
					&%neighbour_molecule_starting_indices(species_molecules)%first_index_in_connection_list
					extra_atom_entries=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&neighbour_molecule_starting_indices(species_molecules)%&
					&extra_atom_entries_in_connection_list
					!update charge with this neighbour molecule
					species_charge=species_charge+&
					&give_charge_of_molecule(donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
					&%connection(connections_firstindex)%indices(1))%molecule_type_index)
					IF (acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
					&%neighbour_molecule_starting_indices(species_molecules)%unmatched_entry) THEN
						double_connections_counter=1
doublespecies:			DO species_molecules_doubles=species_molecules+1,acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%total_number_of_neighbour_molecules,1
							!Compare the two molecules in species_molecules and species_molecules_doubles
							connections_firstindex_2=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
							&%neighbour_molecule_starting_indices(species_molecules_doubles)%first_index_in_connection_list
							extra_atom_entries_2=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&neighbour_molecule_starting_indices(species_molecules_doubles)%&
							&extra_atom_entries_in_connection_list
							IF (extra_atom_entries/=extra_atom_entries_2) CYCLE doublespecies
							!Only get to this point for the same number of extra_atom_entries!
							DO connection_counter=0,extra_atom_entries
								!First, get the acceptor and donor group of the "reference" molecule in species_molecules
								donortype=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
								&%connection(connection_counter+connections_firstindex)%indices(1)
								acceptorgroup=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
								&%connection(connection_counter+connections_firstindex)%indices(4)
								donorgroup=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
								&%connection(connection_counter+connections_firstindex)%indices(3)
								!THEN, compare with the acceptor and donor group of the "observed" molecule in species_molecules_doubles
								donortype_2=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
								&%connection(connection_counter+connections_firstindex_2)%indices(1)
								acceptorgroup_2=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
								&%connection(connection_counter+connections_firstindex_2)%indices(4)
								donorgroup_2=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
								&%connection(connection_counter+connections_firstindex_2)%indices(3)
								IF ((donorgroup/=donorgroup_2)&
								&.OR.(acceptorgroup/=acceptorgroup_2)&
								&.OR.(donortype/=donortype_2)) CYCLE doublespecies
							ENDDO
							!Only get to this point for the same connections in the two compared molecules!
							double_connections_counter=double_connections_counter+1
							acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&neighbour_molecule_starting_indices(species_molecules_doubles)%unmatched_entry=.FALSE.
						ENDDO doublespecies
						!Still in "unmatched_entry", but now we know how many double entries. Need to cycle through again...
						IF (double_connections_counter>1) THEN
							WRITE(4,ADVANCE="NO",FMT='("     ",I0," Molecules of type ")') double_connections_counter
						ELSE
							WRITE(4,ADVANCE="NO",FMT='("     One molecule of type ")')
						ENDIF
						WRITE(4,ADVANCE="NO",FMT='(I0,"(",A,") with ")')&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
						&%connection(connections_firstindex)%indices(1))%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
						&%connection(connections_firstindex)%indices(1))%molecule_type_index))
						IF (extra_atom_entries==0) THEN
							WRITE(4,ADVANCE="NO",FMT='("one connection")')
						ELSE
							WRITE(4,ADVANCE="NO",FMT='(I0," connections")')extra_atom_entries+1
						ENDIF
						IF (double_connections_counter>1) THEN
							WRITE(4,'(" each:")')
						ELSE
							WRITE(4,'(":")')
						ENDIF
						DO connection_counter=0,extra_atom_entries
							acceptorgroup=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
							&%connection(connection_counter+connections_firstindex)%indices(4)
							donorgroup=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
							&%connection(connection_counter+connections_firstindex)%indices(3)
							IF (donorgroup<0) THEN
								WRITE(4,ADVANCE="NO",FMT='("       atom group ",I0," (donor) ")')-donorgroup
							ELSE
								WRITE(4,ADVANCE="NO",FMT='("       atom index ",I0," (donor) ")')donorgroup
							ENDIF
							IF (acceptorgroup<0) THEN
								WRITE(4,'("to atom group ",I0," (acceptor)")')-acceptorgroup
							ELSE
								WRITE(4,'("to atom index ",I0," (acceptor)")')acceptorgroup
							ENDIF
						ENDDO
					ENDIF
				ENDDO
				WRITE(4,'("   The total formal charge including the acceptor was ",SP,I0,SS,".")')&
				&species_charge
				!print the first and last occurrence
				!Note that we are still in the "best_index" loop
				WRITE(filename_speciation_output,'(A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"acceptor",n_acceptor_in,"_species",species_output_counter,"_first.xyz"
				OPEN(UNIT=3,FILE=TRIM(filename_speciation_output),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(species_output_header,'(A,I0," ",I0," ",I0," ",I0," xx.xx",A)')&
				&"Try 'dump_cut T ",acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_first_occurrence,acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_first_occurrence,&
				&acceptor_list(n_acceptor_in)%molecule_type_index,acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%first_occurrence_acceptorindices(2)&
				&,"' (choose a good cutoff distance xx.xx)"
				!Make sure we do not print doubles, but include the acceptor
				species_atomcount=give_number_of_atoms_per_molecule(acceptor_list(n_acceptor_in)%molecule_type_index)
				previous_donor_molecule=0
				!SANITY CHECK
				IF (acceptor_list(acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%first_occurrence_acceptorindices(1))&
				&%molecule_type_index/=acceptor_list(n_acceptor_in)%molecule_type_index) CALL report_error(0)
				DO connection_counter=1,acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%&
				&number_of_logged_connections_in_species_list,1
					!Is this a new molecule? if so, we will need to reset the atom index, too
					IF (acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&first_occurrence(connection_counter)%indices(2)/=previous_donor_molecule) THEN
						species_atomcount=species_atomcount+&
						&give_number_of_atoms_per_molecule(donor_list&
						&(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(1))%molecule_type_index)
						previous_donor_molecule=&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(2)
					ENDIF
				ENDDO
				!put the acceptor in the origin
				acceptor_centre(:)=give_center_of_mass(&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_first_occurrence,acceptor_list(n_acceptor_in)%molecule_type_index,&
				&acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%first_occurrence_acceptorindices(2))
				!Write a scratch file with the beads and count their number. needs to be added to species_atomcount
				N_beads_total=0 !This needs to be outside the following IF-THEN-ELSE, because we add it to the header
				INQUIRE(UNIT=10,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=10)
				OPEN(UNIT=10,STATUS="SCRATCH")
				!rewind scratch file
				REWIND 10
				!Append all the molecules that are part of this species
				!Make sure we do not print doubles, but include the acceptor
				CALL write_molecule(10,&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_first_occurrence,acceptor_list(n_acceptor_in)%molecule_type_index,&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%first_occurrence_acceptorindices(2),&
				&include_header=.FALSE.,translate_by=-acceptor_centre(:))
				previous_donor_molecule=0
				DO connection_counter=1,acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%&
				&number_of_logged_connections_in_species_list,1
					!Is this a new molecule? if so, we will need to reset the atom index, too
					IF (acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&first_occurrence(connection_counter)%indices(2)/=previous_donor_molecule) THEN
						previous_donor_molecule=&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(2)
						!get the shift vector... needed for PBC
						distance=give_smallest_distance_squared(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,&
						&acceptor_list(n_acceptor_in)%molecule_type_index,&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(1))%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%first_occurrence_acceptorindices(2),&
						&previous_donor_molecule,&
						&shift)
						!finally, write output
						CALL write_molecule(10,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,donor_list(acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(1))%molecule_type_index,&
						&previous_donor_molecule,&
						&include_header=.FALSE.,translate_by=shift(:)-acceptor_centre(:))
					ENDIF
					IF (print_connection_beads) THEN
						!Add beads connecting the two atom vectors
						base_atom=give_atom_position(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,&
						&acceptor_list(n_acceptor_in)%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence_acceptorindices(2),&
						&acceptor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence_acceptorindices(1))%atom_index_list(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(4),1))
						tip_atom=shift(:)+give_atom_position(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(1))%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(2),&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(1))%atom_index_list(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&first_occurrence(connection_counter)%indices(3),1))
						connection_vector=tip_atom-base_atom
						distance=SQRT(SUM(connection_vector(:)**2))
						connection_vector=connection_vector/distance !this should now be a unit vector
						beadspan=distance-2.0*bead_clearance
						IF (beadspan>0.0) THEN
							!How many beads should we add?
							N_beads=INT(beadspan/bead_distance)
							N_beads_total=N_beads_total+N_beads+1
							!Convert it back to a REAL so it is nice...
							beadspan=FLOAT(N_beads)*bead_distance
							!... the reason we need this is to update the bead_clearance to (distance-beadspan)/2.0
							base_atom=base_atom+connection_vector*((distance-beadspan)/2.0)
							!Adding beads
							DO m=0,N_beads
								WRITE(10,*) "Z ",base_atom+connection_vector*FLOAT(m)*bead_distance-acceptor_centre(:)
							ENDDO
						ENDIF
					ENDIF
				ENDDO
				!Transfer the scratch file
				REWIND 10
				WRITE(3,'(I0)')species_atomcount+N_beads_total
				WRITE(3,'(A)')TRIM(species_output_header)
				DO m=1,species_atomcount+N_beads_total
					READ(10,*) element_name,base_atom
					WRITE(3,*) TRIM(element_name),base_atom
				ENDDO
				CLOSE(UNIT=10)
				WRITE(3,*)
				WRITE(3,*)
				ENDFILE 3
				CLOSE(UNIT=3)
				!very inefficient in terms of amount of code:
				!I now double everything from the "first occurrence"
				!Note that we are still in the "best_index" loop
				WRITE(filename_speciation_output,'(A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
			&//"acceptor",n_acceptor_in,"_species",species_output_counter,"_last.xyz"
				OPEN(UNIT=3,FILE=TRIM(filename_speciation_output),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(species_output_header,'(A,I0," ",I0," ",I0," ",I0," xx.xx",A)')&
				&"Try 'dump_cut T ",acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_last_occurrence,acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_last_occurrence,&
				&acceptor_list(n_acceptor_in)%molecule_type_index,acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%last_occurrence_acceptorindices(2)&
				&,"' (choose a good cutoff distance xx.xx)"
				! WRITE(species_output_header,'("Timestep=",I0,", Occurrence=",E9.3)')&
				! &acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				! &timestep_last_occurrence,FLOAT(&
				! &acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
				!Make sure we do not print doubles, but include the acceptor
				species_atomcount=give_number_of_atoms_per_molecule(acceptor_list(n_acceptor_in)%molecule_type_index)
				previous_donor_molecule=0
				!SANITY CHECK
				IF (acceptor_list(acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%last_occurrence_acceptorindices(1))&
				&%molecule_type_index/=acceptor_list(n_acceptor_in)%molecule_type_index) CALL report_error(0)
				DO connection_counter=1,acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%&
				&number_of_logged_connections_in_species_list,1
					!Is this a new molecule? if so, we will need to reset the atom index, too
					IF (acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&last_occurrence(connection_counter)%indices(2)/=previous_donor_molecule) THEN
						species_atomcount=species_atomcount+&
						&give_number_of_atoms_per_molecule(donor_list&
						&(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(1))%molecule_type_index)
						previous_donor_molecule=&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(2)
					ENDIF
				ENDDO
				!put the acceptor in the origin
				acceptor_centre(:)=give_center_of_mass(&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_last_occurrence,acceptor_list(n_acceptor_in)%molecule_type_index,&
				&acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%last_occurrence_acceptorindices(2))
				!Write a scratch file with the beads and count their number. needs to be added to species_atomcount
				N_beads_total=0 !This needs to be outside the following IF-THEN-ELSE, because we add it to the header
				INQUIRE(UNIT=10,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=10)
				OPEN(UNIT=10,STATUS="SCRATCH")
				!rewind scratch file
				REWIND 10
				!Append all the molecules that are part of this species
				!Make sure we do not print doubles, but include the acceptor
				CALL write_molecule(10,&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
				&timestep_last_occurrence,acceptor_list(n_acceptor_in)%molecule_type_index,&
				&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%last_occurrence_acceptorindices(2),&
				&include_header=.FALSE.,translate_by=-acceptor_centre(:))
				previous_donor_molecule=0
				DO connection_counter=1,acceptor_list(n_acceptor_in)%&
				&list_of_all_species(best_index)%&
				&number_of_logged_connections_in_species_list,1
					!Is this a new molecule? if so, we will need to reset the atom index, too
					IF (acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&last_occurrence(connection_counter)%indices(2)/=previous_donor_molecule) THEN
						previous_donor_molecule=&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(2)
						!get the shift vector... needed for PBC
						distance=give_smallest_distance_squared(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,&
						&acceptor_list(n_acceptor_in)%molecule_type_index,&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(1))%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%last_occurrence_acceptorindices(2),&
						&previous_donor_molecule,&
						&shift)
						!finally, write output
						CALL write_molecule(10,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,donor_list(acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(1))%molecule_type_index,&
						&previous_donor_molecule,&
						&include_header=.FALSE.,translate_by=shift(:)-acceptor_centre(:))
					ENDIF
					IF (print_connection_beads) THEN
						!Add beads connecting the two atom vectors
						base_atom=give_atom_position(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,&
						&acceptor_list(n_acceptor_in)%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence_acceptorindices(2),&
						&acceptor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence_acceptorindices(1))%atom_index_list(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(4),1))
						tip_atom=shift(:)+give_atom_position(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(1))%molecule_type_index,&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(2),&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(1))%atom_index_list(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&last_occurrence(connection_counter)%indices(3),1))
						connection_vector=tip_atom-base_atom
						distance=SQRT(SUM(connection_vector(:)**2))
						connection_vector=connection_vector/distance !this should now be a unit vector
						beadspan=distance-2.0*bead_clearance
						IF (beadspan>0.0) THEN
							!How many beads should we add?
							N_beads=INT(beadspan/bead_distance)
							N_beads_total=N_beads_total+N_beads+1
							!Convert it back to a REAL so it is nice...
							beadspan=FLOAT(N_beads)*bead_distance
							!... the reason we need this is to update the bead_clearance to (distance-beadspan)/2.0
							base_atom=base_atom+connection_vector*((distance-beadspan)/2.0)
							!Adding beads
							DO m=0,N_beads
								WRITE(10,*) "Z ",base_atom+connection_vector*FLOAT(m)*bead_distance-acceptor_centre(:)
							ENDDO
						ENDIF
					ENDIF
				ENDDO
				!Transfer the scratch file
				REWIND 10
				WRITE(3,'(I0)')species_atomcount+N_beads_total
				WRITE(3,'(A)')TRIM(species_output_header)
				DO m=1,species_atomcount+N_beads_total
					READ(10,*) element_name,base_atom
					WRITE(3,*) TRIM(element_name),base_atom
				ENDDO
				CLOSE(UNIT=10)
				WRITE(3,*)
				WRITE(3,*)
				ENDFILE 3
				CLOSE(UNIT=3)
			ENDDO
			IF ((acceptor_list(n_acceptor_in)%species_overflow)>0) THEN
				IF (FLOAT((acceptor_list(n_acceptor_in)%species_overflow))/total_occurrence>0.01) THEN
					WRITE(*,'("     There was species overflow amounting to ",F0.1,"%")')&
					&100.0*FLOAT(acceptor_list(n_acceptor_in)%species_overflow)/total_occurrence
				ELSE
					WRITE(*,'("     There was species overflow amounting to ",E9.3)')&
					&FLOAT(acceptor_list(n_acceptor_in)%species_overflow)/total_occurrence
				ENDIF
			ENDIF
			IF ((acceptor_list(n_acceptor_in)%no_neighbour_occurrences)>0) THEN
				IF (FLOAT((acceptor_list(n_acceptor_in)%no_neighbour_occurrences))/total_occurrence>0.01) THEN
					WRITE(*,'("     There were ",F0.1,"% acceptor molecules without any neighbouring donors.")')&
					&100.0*FLOAT(acceptor_list(n_acceptor_in)%no_neighbour_occurrences)/total_occurrence
				ELSE
					WRITE(*,'("     There were ",I0," acceptor molecules without neighbouring donors (less than 1%).")')&
					&acceptor_list(n_acceptor_in)%no_neighbour_occurrences
				ENDIF
			ELSE
				WRITE(*,'("     There were no acceptor molecules without neighbouring donors.")')
			ENDIF
			WRITE(*,'(A,I0,A)') "   Detailed statistics written to '"&
			&//TRIM(ADJUSTL(OUTPUT_PREFIX))//"acceptor",n_acceptor_in,"_speciation_statistics.dat'"
			WRITE(*,'(A,I0,A)')&
			&"   First observed structure per species written to '"&
			&//TRIM(ADJUSTL(OUTPUT_PREFIX))//"acceptor"&
			&,n_acceptor_in,"_speciesX_first.xyz'"
			WRITE(*,'(A,I0,A)')&
			&"   Last observed structure per species written to '"&
			&//TRIM(ADJUSTL(OUTPUT_PREFIX))//"acceptor"&
			&,n_acceptor_in,"_speciesX_last.xyz'"
			IF (ALLOCATED(occurrences)) THEN
				DEALLOCATE(occurrences,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			WRITE(4,*)
			WRITE(4,'(" The donors of this analysis were:")')
			DO donor_type_counter=1,N_donor_types,1
				IF (donor_list(donor_type_counter)%N_atoms==1) THEN
					IF (N_donor_types==1) THEN
						WRITE(4,'(" Atom ",I0," (",A,") of molecule type ",I0," (",A,") was the sole donor.")')&
						&donor_list(donor_type_counter)%atom_index_list(1,1),&
						&give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
						&donor_list(donor_type_counter)%atom_index_list(1,1)),&
						&donor_list(donor_type_counter)%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
					ELSE
						WRITE(4,'("   Donor #",I0," was atom ",I0," (",A,") of molecule type ",I0," (",A,").")')&
						&donor_type_counter,&
						&donor_list(donor_type_counter)%atom_index_list(1,1),&
						&give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
						&donor_list(donor_type_counter)%atom_index_list(1,1)),&
						&donor_list(donor_type_counter)%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
					ENDIF
				ELSE
					IF (atoms_grouped) THEN
						WRITE(4,'("   Donor #",I0," was molecule type ",I0," (",A,"), with the following atom groups:")')&
						&donor_type_counter,&
						&donor_list(donor_type_counter)%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
						DO group_number=1,donor_list(donor_type_counter)%N_groups,1
							firsttime=.TRUE.
							DO i=1,donor_list(donor_type_counter)%N_atoms,1
								IF (group_number==-donor_list(donor_type_counter)%atom_index_list(i,2)) THEN
									IF (firsttime) THEN
										firsttime=.FALSE.
										WRITE(4,ADVANCE="NO",FMT='("     Group #",I0,", atom indices ",I0,"(",A,")")')&
										&group_number,donor_list(donor_type_counter)%atom_index_list(i,1),&
										&TRIM(give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
										&donor_list(donor_type_counter)%atom_index_list(i,1)))
									ELSE
										WRITE(4,ADVANCE="NO",FMT='(", ",I0,"(",A,")")')&
										&donor_list(donor_type_counter)%atom_index_list(i,1),&
										&TRIM(give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
										&donor_list(donor_type_counter)%atom_index_list(i,1)))
									ENDIF
								ENDIF
							ENDDO
							IF (firsttime) THEN
								WRITE(4,FMT='("     Group #",I0," was empty.")')group_number
							ELSE
								WRITE(4,'(".")')
							ENDIF
						ENDDO
						IF (.NOT.(grouped_by_elements)) THEN
							no_ungrouped_atoms=.TRUE.
							firsttime=.TRUE.
							DO i=1,donor_list(donor_type_counter)%N_atoms,1
								IF (donor_list(donor_type_counter)%atom_index_list(i,2)>0) THEN
									!sanity check
									IF (donor_list(donor_type_counter)%atom_index_list(i,2)/=&
									&donor_list(donor_type_counter)%atom_index_list(i,1)) CALL report_error(0)
									no_ungrouped_atoms=.FALSE.
									IF (firsttime) THEN
										WRITE(4,FMT='("     Atoms not assigned to a group, treated as their own group:")')
										WRITE(4,ADVANCE="NO",FMT='("       Atom indices ",I0)')&
										&donor_list(donor_type_counter)%atom_index_list(i,1)
										firsttime=.FALSE.
									ELSE
										WRITE(4,ADVANCE="NO",FMT='(", ",I0)')&
										&donor_list(donor_type_counter)%atom_index_list(i,1)
									ENDIF
								ENDIF
							ENDDO
							IF (no_ungrouped_atoms) THEN
								WRITE(4,'("     All atoms of this donor were assigned to a group.")')
							ELSE
								WRITE(4,'(".")')
							ENDIF
						ENDIF
					ELSE
						WRITE(4,'("   Donor #",I0," was molecule type ",I0," (",A,"), with the following atom indices:")')&
						&donor_type_counter,&
						&donor_list(donor_type_counter)%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
						DO i=1,donor_list(donor_type_counter)%N_atoms,1
							WRITE(4,FMT='("     Atom index ",I0," (",A,")")')&
							&donor_list(donor_type_counter)%atom_index_list(i,1),&
							&TRIM(give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
							&donor_list(donor_type_counter)%atom_index_list(i,1)))
						ENDDO
					ENDIF
				ENDIF
			ENDDO
			ENDFILE 4
			CLOSE(UNIT=4)
		END SUBROUTINE print_acceptor_summary

		!The following subroutine was borrowed from MODULE AUTOCORRELATION - needed quite a few adjustments so I put it here
		SUBROUTINE calculate_autocorrelation_function_from_master_array_LOG(n_acceptor_in)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: n_acceptor_in
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: timeline,logarithmic_entry_counter,entries
		REAL(WORKING_PRECISION),DIMENSION(:,:),ALLOCATABLE :: autocorrelation_function
		!the quantity C(t), which is calculated as uncorrected and then afterwards corrected.
		!here, the first dimension are the timesteps, and the second dimension are the possible species.
		REAL(WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: temp_function !running over the possible species.
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: molecule_speciesarray !This array essentially contains the information in the time_series file. first dimension is timestep, second dimension are the molecule indices.
		INTEGER,DIMENSION(:),ALLOCATABLE :: timesteps_array !which timestep is this entry in autocorrelation_function?
		INTEGER :: allocstatus,deallocstatus,ios,stepcounter,species_counter,dummy,molecule_counter
		LOGICAL :: connected
			WRITE(*,'("   Calculating species autocorrelation function for this acceptor:")')
			!allocate memory for the molecule_speciesarray and fill from file.
			IF (ALLOCATED(molecule_speciesarray)) CALL report_error(0)
			ALLOCATE(molecule_speciesarray(MAX((nsteps-1+sampling_interval)/sampling_interval,0),&
			&give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index)),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			! INQUIRE(UNIT=3,OPENED=connected)
			! IF (connected) CALL report_error(27,exit_status=3)
			! WRITE(filename_time_series,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
			! &//"acceptor",n_acceptor_in,"_time_series.dat"
			! OPEN(UNIT=3,FILE=TRIM(filename_time_series),IOSTAT=ios)
			IF (ios/=0) CALL report_error(26,exit_status=ios)
			!the file was filled in the loop "DO stepcounter=1,nsteps,sampling_interval".
			REWIND 10+n_acceptor_in
			! WRITE(*,'("     Reading file ",A,"...")')"'"//TRIM(filename_time_series)//"'"
			DO stepcounter=1,MAX((nsteps-1+sampling_interval)/sampling_interval,0),1
				!the actual steps are stepcounter*sampling_interval*TIME_SCALING_FACTOR
				READ(10+n_acceptor_in,IOSTAT=ios,FMT=*) dummy,molecule_speciesarray(stepcounter,:)
				!sanity check
				IF (dummy/=(stepcounter-1)*TIME_SCALING_FACTOR*sampling_interval) CALL report_error(0)
				IF (ios/=0) THEN
					CALL report_error(167,exit_status=ios)
					DEALLOCATE(molecule_speciesarray,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ENDIF
			ENDDO
			CLOSE(UNIT=10+n_acceptor_in)
			!Count how many entries we will need. safer than me trying to get that number with maths.
			timeline=0
			entries=0
			DO WHILE (timeline<tmax)
				IF (use_logarithmic_spacing) THEN
					timeline=timeline+MAX(INT(LOG(FLOAT(timeline))),1)
				ELSE
					timeline=timeline+1
				ENDIF
				IF (timeline>tmax) timeline=tmax
				entries=entries+1
			ENDDO
			IF (use_logarithmic_spacing) THEN
				WRITE(*,'("     Logarithmic spacing, autocorrelation function will have ",I0," entries up to ",I0,".")')&
				&entries,tmax
			ELSE
				WRITE(*,'("     Linear spacing, autocorrelation function will have ",I0," entries up to ",I0,".")')&
				&entries,tmax
			ENDIF
			!allocate memory for the timesteps (from t=0 to t=tmax)
			IF (ALLOCATED(timesteps_array)) CALL report_error(0)
			ALLOCATE(timesteps_array(0:entries),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			timeline=0
			entries=0
			timesteps_array(0)=0
			!get the entries again - this is needed, because I do not want to work with pointers and dynamic stuff.
			DO WHILE (timeline<tmax)
				IF (use_logarithmic_spacing) THEN
					timeline=timeline+MAX(INT(LOG(FLOAT(timeline))),1)
				ELSE
					timeline=timeline+1
				ENDIF
				IF (timeline>tmax) timeline=tmax
				entries=entries+1
				timesteps_array(entries)=timeline
			ENDDO
			!allocate memory for the autocorrelation_function (from t=0 to t=tmax)
			IF (ALLOCATED(autocorrelation_function)) CALL report_error(0)
			ALLOCATE(autocorrelation_function(0:entries,0:acceptor_list(n_acceptor_in)%&
			&number_of_logged_species_in_acceptor_list),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			!initialise autocorrelation_function.
			autocorrelation_function(:,:)=0.0_DP
			!$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(temp_function,timeline,logarithmic_entry_counter,stepcounter)&
			!$OMP PRIVATE (species_counter,allocstatus)
			!$OMP SINGLE
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE (*,ADVANCE="NO",FMT='(A,I0)') "     ### Parallel execution on ",OMP_get_num_threads()
		 !$ 	WRITE (*,'(A)') " threads (intermittent autocorrelation function)"
		 !$ 	CALL timing_parallel_sections(.TRUE.)
		 !$ ENDIF
			IF (VERBOSE_OUTPUT) THEN
				IF ((give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index))>100)&
				&WRITE(*,ADVANCE="NO",FMT='("    ")')
				CALL print_progress(give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index))
				IF ((give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index))>100)&
				&WRITE(*,ADVANCE="NO",FMT='("    ")')
			ENDIF
			!$OMP END SINGLE
			!allocate the temporary autocorrelation function
			IF (ALLOCATED(temp_function)) CALL report_error(0)
			ALLOCATE(temp_function(0:acceptor_list(n_acceptor_in)%&
			&number_of_logged_species_in_acceptor_list),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			!iterate over molecules and pass chunks to the subroutine that iterates over stepcounter
			!$OMP DO SCHEDULE(STATIC,1)
			DO molecule_counter=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index),1
				!'timeline' is the argument of the autocorrelation_function, i.e. the time shift of the original function.
				!for example, if the current shift is timeline=1000, and there are 10000 timesteps in total,
				!then the argument has to be evaluated from (h(0+0)-<h>)(h(0+1000)-<h>) up to (h(9000+0)-<h>)(h(9000+1000)-<h>)
				!timeline=0 is 1.0 for all
				DO logarithmic_entry_counter=0,entries,1
					timeline=timesteps_array(logarithmic_entry_counter)
					!inner loop iterates over the whole chunk, i.e. the subset of the autocorr_array for the molecule in question.
					temp_function(:)=0.0d0
					DO stepcounter=1,tmax-timeline+1,1
						!the actual steps are stepcounter*sampling_interval*TIME_SCALING_FACTOR
						!this is the central part of the whole autocorrelation process.
						!because we have the species, it should look something like this.....
						DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list,1
							IF (molecule_speciesarray(stepcounter,molecule_counter)==species_counter)THEN
								!h(t0) is one
								IF (molecule_speciesarray(stepcounter+timeline,molecule_counter)==species_counter)THEN
									!h(t0+t) is one
									temp_function(species_counter)=temp_function(species_counter)+&
									&(1.0d0-average_h(species_counter))*(1.0d0-average_h(species_counter))
								ELSE
									!h(t0+t) is zero
									temp_function(species_counter)=temp_function(species_counter)+&
									&(1.0d0-average_h(species_counter))*(0.0d0-average_h(species_counter))
								ENDIF
							ELSE
								!h(t0) is zero
								IF (molecule_speciesarray(stepcounter+timeline,molecule_counter)==species_counter)THEN
									!h(t0+t) is one
									temp_function(species_counter)=temp_function(species_counter)+&
									&(0.0d0-average_h(species_counter))*(1.0d0-average_h(species_counter))
								ELSE
									!h(t0+t) is zero
									temp_function(species_counter)=temp_function(species_counter)+&
									&(0.0d0-average_h(species_counter))*(0.0d0-average_h(species_counter))
								ENDIF
							ENDIF
						ENDDO
					ENDDO
					!Normalise result. For the above example, this would be division by 9000
					!$OMP CRITICAL(autocorrelation_updates)
					autocorrelation_function(logarithmic_entry_counter,:)=&
					&autocorrelation_function(logarithmic_entry_counter,:)+temp_function(:)/DFLOAT(tmax-timeline+1)
					!$OMP END CRITICAL(autocorrelation_updates)
				ENDDO! End of outer loop over time shifts
				!$OMP CRITICAL(progress)
				IF (VERBOSE_OUTPUT) CALL print_progress()
				!$OMP END CRITICAL(progress)
			ENDDO
			!$OMP END DO
			DEALLOCATE(temp_function,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			!$OMP END PARALLEL
			IF (((give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index))>100)&
			&.AND.(VERBOSE_OUTPUT)) WRITE(*,*)
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='("     ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
		 !$ ENDIF
			!normalise autocorrelation_function by number of molecules
			autocorrelation_function(:,:)=autocorrelation_function(:,:)/&
			&DFLOAT(give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index))
			!print / report autocorrelation function
			CALL report_autocorrelation_function()
			DEALLOCATE(timesteps_array,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(autocorrelation_function,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(molecule_speciesarray,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			IF (ALLOCATED(average_h)) THEN
				DEALLOCATE(average_h,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF

			CONTAINS

				SUBROUTINE report_autocorrelation_function
				!Output formats - AUTOCORRELATION module
				IMPLICIT NONE
				CHARACTER(LEN=12) :: species_number_output
				CHARACTER(LEN=1024) :: filename_autocorrelationlog
				REAL :: h_low,h_high,h_stdev,h_local
					WRITE(filename_autocorrelationlog,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
					&//"acceptor",n_acceptor_in,"_species_autocorrelation.dat"
					IF (VERBOSE_OUTPUT) WRITE(*,'(A,A,A)') "     writing autocorrelation function into file '",&
					&TRIM(ADJUSTL(filename_autocorrelationlog)),"'"
					INQUIRE(UNIT=3,OPENED=connected)
					IF (connected) CALL report_error(27,exit_status=3)
					OPEN(UNIT=3,FILE=TRIM(ADJUSTL(filename_autocorrelationlog)),IOSTAT=ios)
					IF (ios/=0) CALL report_error(26,exit_status=ios)
					WRITE(3,*) "This file contains the intermittent autocorrelation function based on the input file '"&
					&,TRIM(FILENAME_SPECIATION_INPUT),"'"
			!		IF (VERBOSE_OUTPUT) WRITE(*,'(" normalise autocorrelation function: dividing by ",F5.3)') 0.0
					WRITE(3,*) "One column per species, apart from the first column. First, some statistics across molecules are printed:"
					WRITE(3,*) "Standard deviation, average/lowest/highest <h>. This is followed by normalisation and species numbers."
					WRITE(3,*) "After that, first column = timeline, other columns are C(t)=<(h(t0+t)-<h>)(h(t0)-<h>)>/<(h(t0)-<h>)**2>"!up until here, autocorrelation_function is not normalised
					!sanity check
					IF (MAX((nsteps-1+sampling_interval)/sampling_interval,0)/=&
					&SIZE(molecule_speciesarray(:,1))) CALL report_error(0)
					!calculate the standard deviation, lowest, and highest average per molecule.
					WRITE(3,ADVANCE="NO",FMT='(A15)') "<h>_stdev"
					DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
						h_stdev=0.0
						DO molecule_counter=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index)
							!how many times was this particular molecule of kind "species_counter"?
							h_local=FLOAT(COUNT(molecule_speciesarray(:,molecule_counter)==species_counter))/&
							&FLOAT(SIZE(molecule_speciesarray(:,molecule_counter)))
							h_stdev=h_stdev+(h_local-average_h(species_counter))**2
						ENDDO
						h_stdev=h_stdev/(FLOAT(give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index)-1))
						h_stdev=SQRT(h_stdev)
						WRITE(3,ADVANCE="NO",FMT='(F15.10)') h_stdev
					ENDDO
					WRITE(3,*)
					WRITE(3,ADVANCE="NO",FMT='(A15)') "<h>"
					DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
						WRITE(3,ADVANCE="NO",FMT='(F15.10)') average_h(species_counter)
					ENDDO
					WRITE(3,*)
					WRITE(3,ADVANCE="NO",FMT='(A15)') "<h>_low"
					DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
						h_low=1.0
						DO molecule_counter=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index)
							!how many times was this particular molecule of kind "species_counter"?
							h_local=FLOAT(COUNT(molecule_speciesarray(:,molecule_counter)==species_counter))/&
							&FLOAT(SIZE(molecule_speciesarray(:,molecule_counter)))
							IF (h_local<h_low) h_low=h_local
						ENDDO
						WRITE(3,ADVANCE="NO",FMT='(F15.10)') h_low
					ENDDO
					WRITE(3,*)
					WRITE(3,ADVANCE="NO",FMT='(A15)') "<h>_high"
					DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
						h_high=0.0
						DO molecule_counter=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor_in)%molecule_type_index)
							!how many times was this particular molecule of kind "species_counter"?
							h_local=FLOAT(COUNT(molecule_speciesarray(:,molecule_counter)==species_counter))/&
							&FLOAT(SIZE(molecule_speciesarray(:,molecule_counter)))
							IF (h_local>h_high) h_high=h_local
						ENDDO
						WRITE(3,ADVANCE="NO",FMT='(F15.10)') h_high
					ENDDO
					WRITE(3,*)
					WRITE(3,ADVANCE="NO",FMT='(A15)') "<(h(t0)-<h>)^2>"
					DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
						WRITE(3,ADVANCE="NO",FMT='(F15.10)')autocorrelation_function(0,species_counter)
					ENDDO
					WRITE(3,*)
					WRITE(3,ADVANCE="NO",FMT='(A15)') "timeline"
					WRITE(species_number_output,'("#",I0)') 0
					WRITE(3,ADVANCE="NO",FMT='(A15)') species_number_output
					DO species_counter=1,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
						WRITE(species_number_output,'("#",I0)')&
						&acceptor_list(n_acceptor_in)%list_of_all_species(species_counter)%rank
						WRITE(3,ADVANCE="NO",FMT='(A15)') species_number_output
					ENDDO
					WRITE(3,*)
					DO logarithmic_entry_counter=0,entries,1
						timeline=timesteps_array(logarithmic_entry_counter)
						WRITE(3,ADVANCE="NO",FMT='(I15)') timeline*TIME_SCALING_FACTOR*sampling_interval
						DO species_counter=0,acceptor_list(n_acceptor_in)%number_of_logged_species_in_acceptor_list
							WRITE(3,ADVANCE="NO",FMT='(F15.10)')&
							&SNGL(autocorrelation_function(logarithmic_entry_counter,species_counter)&
							&/autocorrelation_function(0,species_counter))
						ENDDO
						WRITE(3,*)
					ENDDO
					ENDFILE 3
					CLOSE(UNIT=3)
				END SUBROUTINE report_autocorrelation_function

		END SUBROUTINE calculate_autocorrelation_function_from_master_array_LOG

		SUBROUTINE perform_speciation_analysis()
		IMPLICIT NONE
		INTEGER :: n_acceptor
			CALL initialise_speciation()
			IF ((ERROR_CODE/=33).AND.(ERROR_CODE/=159).AND.(ERROR_CODE/=157)) THEN
				WRITE(*,'(" Using ",I0," timesteps in intervals of ",I0," for averaging.")')&
				&MAX((nsteps-1+sampling_interval)/sampling_interval,0),sampling_interval
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,ADVANCE="NO",FMT='(" The main data structures occupy approximately ")')
					CALL print_memory_requirement(FLOAT(bytes_needed)*(4.0/1024.0d0))
					WRITE(*,'(".")')
				ENDIF
				WRITE(*,'(" Starting speciation analysis.")')
				CALL refresh_IO()
				CALL trajectory_speciation_analysis_parallel() !opens UNIT=10+n_acceptor
				IF (neighbour_atom_overflow>0) CALL report_error(160,exit_status=neighbour_atom_overflow)
				IF (ANY(acceptor_list(:)%species_overflow>0))&
				&CALL report_error(161,exit_status=SUM(acceptor_list(:)%species_overflow))
				!Iterate over all acceptors. Usually, this is probably one.
				IF (N_acceptor_types==1) THEN
					WRITE(*,ADVANCE="NO",FMT='(" There was one acceptor: ")')
					CALL print_acceptor_summary(1) !allocates average_h
					IF (calculate_autocorrelation) THEN
						CALL calculate_autocorrelation_function_from_master_array_LOG(1) !deallocates average_h, closes UNIT=10+n_acceptor
					ENDIF
				ELSE
					WRITE(*,'(" There were ",I0," acceptor molecules:")') N_acceptor_types
					DO n_acceptor=1,N_acceptor_types,1
						WRITE(*,ADVANCE="NO",FMT='("  Acceptor ",I0,"/",I0," was ")')&
						&n_acceptor,N_acceptor_types
						CALL print_acceptor_summary(n_acceptor)
						IF (calculate_autocorrelation) THEN
							CALL calculate_autocorrelation_function_from_master_array_LOG(n_acceptor) !deallocates average_h, closes UNIT=10+n_acceptor
						ENDIF
					ENDDO
				ENDIF
				CALL finalise_speciation()
			ELSE
				ERROR_CODE=ERROR_CODE_DEFAULT
			ENDIF
		END SUBROUTINE perform_speciation_analysis

END MODULE SPECIATION
!--------------------------------------------------------------------------------------------------------------------------------!