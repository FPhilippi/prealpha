
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
	INTEGER,PARAMETER :: nsteps_default=1 !how many steps to use from trajectory
	INTEGER,PARAMETER :: sampling_interval_default=1
	INTEGER,PARAMETER :: maximum_number_of_species_default=11!how many different species do we want to allow?
	INTEGER,PARAMETER :: maximum_number_of_neighbour_molecules_default=10
	INTEGER,PARAMETER :: maximum_number_of_connections_default=10
	LOGICAL,PARAMETER :: print_connection_beads_default=.FALSE.
	!variables
	INTEGER :: nsteps=nsteps_default !how many steps to use from trajectory
	INTEGER :: sampling_interval=sampling_interval_default
	INTEGER :: maximum_number_of_species=maximum_number_of_species_default!how many different species do we want to allow?
	INTEGER :: maximum_number_of_connections=maximum_number_of_connections_default
	INTEGER :: N_acceptor_types
	INTEGER :: N_donor_types
	INTEGER :: species_overflow
	INTEGER :: neighbour_atom_overflow
	LOGICAL :: atoms_grouped
	LOGICAL :: print_connection_beads=print_connection_beads_default
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
		INTEGER(KIND=WORKING_PRECISION) :: problematic_occurrences ! plus one if not the same as "last_occurrence", minus one otherwise
		INTEGER :: total_number_of_neighbour_atoms ! For example, "4" for Li[DME]2
		INTEGER :: total_number_of_neighbour_molecules ! For example, "2" for Li[DME]2
		INTEGER(KIND=WORKING_PRECISION) :: occurrences ! For example, "I have found this species 34857 times"
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
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: atom_index_list !First dimension is listing all relevant atom_index up to N_atoms.
		!the second dimension is 1/2 for the actual atom_index/atom_group.
		!All atoms with the same atom_group number are treated as identical.
		REAL,DIMENSION(:),ALLOCATABLE :: cutoff_list ! CONTRIBUTIONS to the cutoff, i.e. just the radius, counting up to N_atoms
		INTEGER :: number_of_logged_species_in_acceptor_list
		TYPE(species),DIMENSION(:),ALLOCATABLE :: list_of_all_species !First dimension goes from 1 to number_of_logged_species_in_acceptor_list and is allocated up to maximum_number_of_species
	END TYPE donors_and_acceptors
	TYPE(donors_and_acceptors),DIMENSION(:),ALLOCATABLE :: acceptor_list !First dimension counts the acceptor molecule types
	TYPE(donors_and_acceptors),DIMENSION(:),ALLOCATABLE :: donor_list    !First dimension counts the donor molecule types
	LOGICAL :: dumpfirst=.TRUE.
	LOGICAL :: dumplast=.TRUE.
	!PRIVATE/PUBLIC declarations
	PUBLIC :: perform_speciation_analysis

	CONTAINS

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
			print_connection_beads=print_connection_beads_default
		END SUBROUTINE set_defaults

		!initialises the speciation module by reading the specified input file.
		SUBROUTINE initialise_speciation()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,allocstatus,n,m,acceptor_maxatoms,donor_maxatoms,a,b
		REAL :: real_space_distance
		CHARACTER(LEN=32) :: inputstring
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
					acceptor_list(n)%list_of_all_species(:)%problematic_occurrences=0
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
							IF (VERBOSE_OUTPUT) WRITE(*,'("  grouping atoms together based on their element symbols.")')
							CALL group_element_names()
							atoms_grouped=.TRUE.
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
					IF (nsteps<1) THEN
						CALL report_error(57,exit_status=nsteps)
						nsteps=1
					ELSEIF (nsteps>give_number_of_timesteps()) THEN
						CALL report_error(57,exit_status=nsteps)
						nsteps=give_number_of_timesteps()
					ENDIF
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
							WRITE(*,ADVANCE="NO",FMT='("  ! groups for acceptor # ",I0,":")') n_acceptor
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
							WRITE(*,ADVANCE="NO",FMT='("  ! groups for donor # ",I0,":")') n_donor
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
										WRITE(*,ADVANCE="NO",FMT='("  !   group #",I0," (",A,") has atom indices ",I0)') -group_number,&
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
				DA_clipboard(:)%N_groups=0
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
				DO n=1,N_molecule_types,1
					acceptor_list(n)%molecule_type_index=DA_clipboard(n)%molecule_type_index
					acceptor_list(n)%N_atoms=DA_clipboard(n)%N_atoms
					IF (ALLOCATED(acceptor_list(n)%atom_index_list)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%atom_index_list(acceptor_list(n)%N_atoms,2),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					IF (ALLOCATED(acceptor_list(n)%cutoff_list)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%cutoff_list(acceptor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					!initialise
					acceptor_list(n)%atom_index_list(:,:)=0
					acceptor_list(n)%cutoff_list(:)=0.0
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
				DO n=1,N_molecule_types,1
					donor_list(n)%molecule_type_index=DA_clipboard(n)%molecule_type_index
					donor_list(n)%N_atoms=DA_clipboard(n)%N_atoms
					IF (ALLOCATED(donor_list(n)%atom_index_list)) CALL report_error(0)
					ALLOCATE(donor_list(n)%atom_index_list(donor_list(n)%N_atoms,2),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					IF (ALLOCATED(donor_list(n)%cutoff_list)) CALL report_error(0)
					ALLOCATE(donor_list(n)%cutoff_list(donor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					!initialise
					donor_list(n)%atom_index_list(:,:)=0
					donor_list(n)%cutoff_list(:)=0.0
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
									donor_list(counter)%atom_index_list(atom_counter,1)=atom_index
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

		SUBROUTINE trajectory_speciation_analysis()
		IMPLICIT NONE
		INTEGER :: stepcounter,acceptor_molecule_index,n_acceptor,n_donor,donor_indexcounter
		INTEGER :: donor_molecule_index,acceptor_indexcounter,allocstatus,deallocstatus
		LOGICAL :: new_species_found,existing_species,problematic_occurrence,neighbour_atom,neighbour_molecule
		LOGICAL,DIMENSION(:),ALLOCATABLE :: potential_species_match !first dimension goes over all the species
		TYPE(species) :: species_clipboard
		!reset the number of species logged to zero just in case.
		acceptor_list(:)%number_of_logged_species_in_acceptor_list=0
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
		neighbour_atom_overflow=0
		species_overflow=0
		!Iterate over all acceptors.
		DO n_acceptor=1,N_acceptor_types,1
			!average over timesteps, a few should be enough.
			DO stepcounter=1,nsteps,sampling_interval
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
										neighbour_atom_overflow=neighbour_atom_overflow+1
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
						!First we need to make sure that atoms that are considered identical sites are treated as such.
						!This needs to be called here already, it is necessary in ANY case, regardless whether new species or not.
						!collapse_and_sort_atom_indices also assigns neighbour_molecule_starting_indices
						CALL collapse_and_sort_atom_indices()
						!Now to find out if this species already exists in the list.
						!There are many possible different cases and conditions, and I will use a jump to skip the ones which are not necessary.
						!I think in this case, the jump command is better and more easily understood than 1000 nestled if-then-else
						!Very first case: species list is empty. need to add.
						IF (acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list==0) THEN
							new_species_found=.TRUE.
							GOTO 60
						ELSE
							!There is already at least one entry - we now need to find out if there is an identical species already
							new_species_found=.FALSE.
						ENDIF
						!The rest is really not trivial. We will use a list of potential matches, and eliminate them one after one.
						potential_species_match(1:acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list)=.TRUE.
						!This is the first rough check we can do - check if species have the same number of connections, number of neighbour atoms, number of neighbour molecules
						CALL eliminate_number_mismatches()
						!If there are no potential species matches left, then this is a new one!
						IF (.NOT.(ANY(&
						&potential_species_match(1:acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list)&
						&))) THEN
							new_species_found=.TRUE.
							GOTO 60
						ENDIF
						!We now need to check all the molecules one by one,
						!and also check if the atoms which are involved in these connections are the same.
						!Molecule indices do not matter, so we need to check for permutations...
						CALL eliminate_connection_mismatches()
						!If there are no potential species matches left, then this is a new one!
						IF (.NOT.(ANY(&
						&potential_species_match(1:acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list)&
						&))) THEN
							new_species_found=.TRUE.
							GOTO 60
						ENDIF
					60	IF (new_species_found) THEN
							CALL append_species()
						ELSE
							IF (COUNT(potential_species_match(1:acceptor_list(n_acceptor)%&
							&number_of_logged_species_in_acceptor_list))>1) CALL report_error(0)
							CALL add_to_species()
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		ENDDO
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

		CONTAINS

			SUBROUTINE collapse_and_sort_atom_indices()
			IMPLICIT NONE
			INTEGER :: current_donor_type,current_donor_molecule_index,connections_firstindex,connections_lastindex
			INTEGER :: current_donor_atom,atoms_firstindex,atoms_lastindex,molecule_counter
			INTEGER :: connection_counter,connection_counter_atoms,species_counter
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
					&acceptor_list(n_acceptor)%&
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
						!then, sort acceptor atoms - to be able to do so, we need to find the first and last donor atoms...
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
				!then, sort acceptor atoms - to be able to do so, we need to find the first and last donor atoms...
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

			SUBROUTINE eliminate_connection_mismatches()
			IMPLICIT NONE
			INTEGER :: species_counter,clipboard_molecules,reference_species_molecules
			INTEGER :: connections_firstindex_ref,extra_atom_entries_ref
			INTEGER :: connections_firstindex_clip,extra_atom_entries_clip
			INTEGER :: connection_counter
			LOGICAL :: molecule_mismatch,successful_species_match
loop_species: 	DO species_counter=1,acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list,1
					successful_species_match=.FALSE.
					!only check potential matches, of course.
					IF (.NOT.(potential_species_match(species_counter))) CYCLE loop_species
					!Now, we need to check the possible permutations of clipboard and the current entry in the list_of_all_species
					! initialise the "lists to keep track" of pairs which are ok
					species_clipboard%neighbour_molecule_starting_indices(1:species_clipboard%&
					&total_number_of_neighbour_molecules)%unmatched_entry=.TRUE.
					acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
					&neighbour_molecule_starting_indices(1:acceptor_list(n_acceptor)%&
					&list_of_all_species(species_counter)%total_number_of_neighbour_molecules)%unmatched_entry=.TRUE.
loop_clipboard:		DO clipboard_molecules=1,species_clipboard%total_number_of_neighbour_molecules
						IF (.NOT.(species_clipboard%neighbour_molecule_starting_indices(clipboard_molecules)%unmatched_entry)) CYCLE
						!the number of molecules should be the same-This has been checked before by eliminate_number_mismatches
						DO reference_species_molecules=1,species_clipboard%total_number_of_neighbour_molecules
							IF (.NOT.(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
							&neighbour_molecule_starting_indices(reference_species_molecules)%unmatched_entry)) CYCLE
							!At this point, we have selected two pairs of molecules, one from the clipboard and one from the reference.
							!The idea is the following:
							! 1) we assume that they are the same, i.e. that there is no mismatche.
							! 2) we check all the connections.
							! 3) if there is ANY connection which differs in ANY of the indices 1/3/4, then they are NOT the same.
							! 4) if they are still the same after checking all connections (in other words, no mismatches were found),
							! 5) then we set the unmatched_entry of this pair to .FALSE. in both "lists to keept track"
							molecule_mismatch=.FALSE. ! step 1)
							!get the first indices of the respective molecules
							connections_firstindex_ref=acceptor_list(n_acceptor)%list_of_all_species(species_counter)&
							&%neighbour_molecule_starting_indices(reference_species_molecules)%first_index_in_connection_list
							connections_firstindex_clip=species_clipboard%neighbour_molecule_starting_indices&
							&(clipboard_molecules)%first_index_in_connection_list
							extra_atom_entries_ref=acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
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
									IF ((acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
										&connection(connection_counter+connections_firstindex_ref)%indices(1)&
										&/=&
										&species_clipboard%connection(connection_counter+connections_firstindex_clip)&
										&%indices(1))&
										&.OR.&
										&(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
										&connection(connection_counter+connections_firstindex_ref)%indices(3)&
										&/=&
										&species_clipboard%connection(connection_counter+connections_firstindex_clip)&
										&%indices(3))&
										&.OR.&
										&(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
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
								acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
								&neighbour_molecule_starting_indices(reference_species_molecules)%unmatched_entry=.FALSE.
								IF ((.NOT.(ANY(species_clipboard%neighbour_molecule_starting_indices(1:species_clipboard%&
								&total_number_of_neighbour_molecules)%unmatched_entry))).AND.&
								&(.NOT.(ANY(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
								&neighbour_molecule_starting_indices(1:acceptor_list(n_acceptor)%&
								&list_of_all_species(species_counter)%total_number_of_neighbour_molecules)%unmatched_entry)))) THEN
									!All have been matched!
									successful_species_match=.TRUE.
									EXIT loop_clipboard
								ENDIF
							ENDIF
						ENDDO
					ENDDO loop_clipboard
					potential_species_match(species_counter)=successful_species_match
				ENDDO loop_species
			END SUBROUTINE eliminate_connection_mismatches


			SUBROUTINE eliminate_number_mismatches()
			IMPLICIT NONE
			INTEGER :: species_counter
!This can be replaced by a fancy "WHERE" statement! do this last when I have a working routine with all parts so I can compare performance...
				!We do not need an initial check of potential_species_match because this subroutine is called first, after the initialisation of potential_species_match
				DO species_counter=1,acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list,1
					IF ((species_clipboard%number_of_logged_connections_in_species_list)/=&
					&(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
					&number_of_logged_connections_in_species_list)) THEN
						potential_species_match(species_counter)=.FALSE.
					ELSE
						IF ((species_clipboard%total_number_of_neighbour_atoms)/=&
						&(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
						&total_number_of_neighbour_atoms)) THEN
							potential_species_match(species_counter)=.FALSE.
						ELSE
							IF ((species_clipboard%total_number_of_neighbour_molecules)/=&
							&(acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
							&total_number_of_neighbour_molecules)) THEN
								potential_species_match(species_counter)=.FALSE.
							ENDIF
						ENDIF
					ENDIF
				ENDDO
			END SUBROUTINE eliminate_number_mismatches

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

			SUBROUTINE append_species()!appends the contents of the clipboard to the actual list
			IMPLICIT NONE
				!check for overflow
				IF (acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list<maximum_number_of_species) THEN
					acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list=&
					&acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list+1
					!append the new species to the list!
					!transfer the connections and their count
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%number_of_logged_connections_in_species_list=&
					&species_clipboard%number_of_logged_connections_in_species_list
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%connection=&
					&species_clipboard%connection
					!initialise (problematic) occurrences
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%problematic_occurrences=0
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%occurrences=1
					!initialise the last/first occurrence - stepcounter
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%timestep_last_occurrence=stepcounter
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%timestep_first_occurrence=stepcounter
					!initialise the last/first occurrence - indices of the acceptor
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%last_occurrence_acceptorindices(:)=&
					&species_clipboard%last_occurrence_acceptorindices
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%first_occurrence_acceptorindices(:)=&
					&species_clipboard%last_occurrence_acceptorindices
					!transfer the last_occurrence
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%first_occurrence=species_clipboard%last_occurrence
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%last_occurrence=species_clipboard%last_occurrence
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%total_number_of_neighbour_atoms=&
					&species_clipboard%total_number_of_neighbour_atoms
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%total_number_of_neighbour_molecules&
					&=species_clipboard%total_number_of_neighbour_molecules
					!update the starting indices
					acceptor_list(n_acceptor)%list_of_all_species(acceptor_list(n_acceptor)%&
					&number_of_logged_species_in_acceptor_list)%neighbour_molecule_starting_indices=&
					&species_clipboard%neighbour_molecule_starting_indices
				ELSE
					species_overflow=species_overflow+1
				ENDIF
			END SUBROUTINE append_species

			SUBROUTINE add_to_species()!appends the contents of the clipboard to the actual list
			IMPLICIT NONE
			INTEGER :: species_counter,connection_counter,problematic
				DO species_counter=1,acceptor_list(n_acceptor)%number_of_logged_species_in_acceptor_list,1
					IF (potential_species_match(species_counter)) THEN
						acceptor_list(n_acceptor)%list_of_all_species(species_counter)%occurrences=&
						&acceptor_list(n_acceptor)%list_of_all_species(species_counter)%occurrences+1
						!a problematic_occurrence is one where the molecule indices are the same, including the acceptor molecule index.
						problematic=+1
						IF (acceptor_list(n_acceptor)%list_of_all_species(species_counter)%last_occurrence_acceptorindices(2)==&
						&species_clipboard%last_occurrence_acceptorindices(2)) THEN
							inner: DO connection_counter=1,acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
								&number_of_logged_connections_in_species_list,1
								IF(species_clipboard%connection(connection_counter)%indices(2)/=&
								&acceptor_list(n_acceptor)%list_of_all_species(species_counter)%&
								&last_occurrence(connection_counter)%indices(2)) THEN
									problematic=-1
									EXIT inner
								ENDIF
							ENDDO inner
						ENDIF
						acceptor_list(n_acceptor)%list_of_all_species(species_counter)%problematic_occurrences=&
						acceptor_list(n_acceptor)%list_of_all_species(species_counter)%problematic_occurrences+problematic
						!update the last occurrence - indices of the acceptor
						acceptor_list(n_acceptor)%list_of_all_species(species_counter)%timestep_last_occurrence=stepcounter
						acceptor_list(n_acceptor)%list_of_all_species(species_counter)%last_occurrence_acceptorindices(:)=&
						&species_clipboard%last_occurrence_acceptorindices
						acceptor_list(n_acceptor)%list_of_all_species(species_counter)%last_occurrence=species_clipboard%last_occurrence
					ENDIF
				ENDDO
			END SUBROUTINE add_to_species

		END SUBROUTINE trajectory_speciation_analysis

		SUBROUTINE report_species()
		IMPLICIT NONE
		INTEGER :: n_acceptor
			!Iterate over all acceptors. Usually, this is probably one.
			IF (N_acceptor_types==1) THEN
				WRITE(*,ADVANCE="NO",FMT='(" There was one acceptor: ")')
				CALL print_acceptor_summary(1)
			ELSE
				WRITE(*,'(" There were ",I0," acceptor molecules:")') N_acceptor_types
				DO n_acceptor=1,N_acceptor_types,1
					WRITE(*,ADVANCE="NO",FMT='("  Acceptor ",I0,"/",I0," was ")')&
					&n_acceptor,N_acceptor_types
					CALL print_acceptor_summary(n_acceptor)
				ENDDO
			ENDIF

		CONTAINS

			SUBROUTINE print_acceptor_summary(n_acceptor_in)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: n_acceptor_in
			INTEGER :: entry_counter,allocstatus,deallocstatus,best_index,species_output_counter
			INTEGER :: n_donor,neighbour_entry_count,i,previous_donor_molecule,donor_type_counter,m
			INTEGER :: species_atomcount,ios,atom_counter,species_charge,group_number,connection_counter,N_beads,N_beads_total
			INTEGER :: acceptorgroup,donorgroup,connections_firstindex,extra_atom_entries,species_molecules
			LOGICAL :: connected
			REAL :: total_occurrence,best_value,distance,species_realcharge,base_atom(3),tip_atom(3),connection_vector(3),beadspan
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
					WRITE(4,'(" Atom ",I0," (",A,") of molecule type ",I0," (",A,") was the acceptor.")')&
					&acceptor_list(n_acceptor_in)%atom_index_list(1,1),&
					&give_element_symbol(acceptor_list(n_acceptor_in)%molecule_type_index,&
					&acceptor_list(n_acceptor_in)%atom_index_list(1,1)),&
					&acceptor_list(n_acceptor_in)%molecule_type_index,&
					&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
				ELSE
					IF (atoms_grouped) THEN
						WRITE(4,'(" Molecule type ",I0," (",A,") was the acceptor, with the following atom groups:")')&
						&acceptor_list(n_acceptor_in)%molecule_type_index,&
						&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
						DO group_number=1,acceptor_list(n_acceptor_in)%N_groups,1
							WRITE(4,ADVANCE="NO",FMT='("      Group #",I0,", atom indices")')group_number
							DO i=1,acceptor_list(n_acceptor_in)%N_atoms,1
								IF (group_number==-acceptor_list(n_acceptor_in)%atom_index_list(i,2)) THEN
									WRITE(4,ADVANCE="NO",FMT='(" ",I0,"(",A,")")')&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1),&
									&TRIM(give_element_symbol(acceptor_list(n_acceptor_in)%molecule_type_index,&
									&acceptor_list(n_acceptor_in)%atom_index_list(i,1)))
								ENDIF
							ENDDO
							WRITE(4,*)
						ENDDO
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
				total_occurrence=total_occurrence+FLOAT(species_overflow)
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
					DO species_molecules=1,acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
					&%total_number_of_neighbour_molecules
						connections_firstindex=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
						&%neighbour_molecule_starting_indices(species_molecules)%first_index_in_connection_list
						extra_atom_entries=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&neighbour_molecule_starting_indices(species_molecules)%&
						&extra_atom_entries_in_connection_list
						WRITE(4,'("     Molecule of type ",I0,"(",A,") with ",I0," connections:")')&
						&donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
						&%connection(connections_firstindex)%indices(1))%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
						&%connection(connections_firstindex)%indices(1))%molecule_type_index)),&
						&extra_atom_entries+1
						DO connection_counter=0,extra_atom_entries
							acceptorgroup=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
							&%connection(connection_counter+connections_firstindex)%indices(4)
							donorgroup=acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
							&%connection(connection_counter+connections_firstindex)%indices(3)
							WRITE(4,'("       atom group ",I0,"(donor) to atom group ",I0," (acceptor)")')&
							&-donorgroup,-acceptorgroup
						ENDDO
						!update charge with this neighbour molecule
						species_charge=species_charge+&
						&give_charge_of_molecule(donor_list(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)&
						&%connection(connections_firstindex)%indices(1))%molecule_type_index)
					ENDDO
					WRITE(4,'("   The similarity across all was ",SP,F5.2,SS,".")')&
					&-1.0*FLOAT(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&problematic_occurrences)/FLOAT(&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
					WRITE(4,'("   The total formal charge including the acceptor was ",SP,I0,SS,".")')&
					&species_charge
					!print the first and last occurrence
					!Note that we are still in the "best_index" loop
					WRITE(filename_speciation_output,'(A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
					&//"acceptor",n_acceptor_in,"_species",species_output_counter,"_first.xyz"
					OPEN(UNIT=3,FILE=TRIM(filename_speciation_output),IOSTAT=ios)
					IF (ios/=0) CALL report_error(26,exit_status=ios)
					WRITE(species_output_header,'("Timestep=",I0,", Occurrence=",E9.3)')&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&timestep_first_occurrence,FLOAT(&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
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
					IF (print_connection_beads) THEN
						INQUIRE(UNIT=10,OPENED=connected)
						IF (connected) CALL report_error(27,exit_status=10)
						OPEN(UNIT=10,STATUS="SCRATCH")
						!rewind scratch file
						REWIND 10
					ENDIF
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
					WRITE(species_output_header,'("Timestep=",I0,", Occurrence=",E9.3)')&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&timestep_last_occurrence,FLOAT(&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
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
					IF (print_connection_beads) THEN
						INQUIRE(UNIT=10,OPENED=connected)
						IF (connected) CALL report_error(27,exit_status=10)
						OPEN(UNIT=10,STATUS="SCRATCH")
						!rewind scratch file
						REWIND 10
					ENDIF
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
				IF (species_overflow>0) THEN
					IF (FLOAT(species_overflow)/total_occurrence>0.01) THEN
						WRITE(*,'("     There was species overflow amounting to ",F0.1,"%")')&
						&100.0*FLOAT(species_overflow)/total_occurrence
					ELSE
						WRITE(*,'("     There was species overflow amounting to ",E9.3)')&
						&FLOAT(species_overflow)/total_occurrence
					ENDIF
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
						WRITE(4,'(" Atom ",I0," (",A,") of molecule type ",I0," (",A,") was the donor.")')&
						&donor_list(donor_type_counter)%atom_index_list(1,1),&
						&give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
						&donor_list(donor_type_counter)%atom_index_list(1,1)),&
						&donor_list(donor_type_counter)%molecule_type_index,&
						&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
					ELSE
						IF (atoms_grouped) THEN
							WRITE(4,'(" Molecule type ",I0," (",A,") was the donor, with the following atom groups:")')&
							&donor_list(donor_type_counter)%molecule_type_index,&
							&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
							DO group_number=1,donor_list(donor_type_counter)%N_groups,1
								WRITE(4,ADVANCE="NO",FMT='("      Group #",I0,", atom indices")')group_number
								DO i=1,donor_list(donor_type_counter)%N_atoms,1
									IF (group_number==-donor_list(donor_type_counter)%atom_index_list(i,2)) THEN
										WRITE(4,ADVANCE="NO",FMT='(" ",I0,"(",A,")")')&
										&donor_list(donor_type_counter)%atom_index_list(i,1),&
										&TRIM(give_element_symbol(donor_list(donor_type_counter)%molecule_type_index,&
										&donor_list(donor_type_counter)%atom_index_list(i,1)))
									ENDIF
								ENDDO
								WRITE(4,*)
							ENDDO
						ELSE
							WRITE(4,'(" Molecule type ",I0," (",A,") was the donor, with the following atom indices:")')&
							&donor_list(donor_type_counter)%molecule_type_index,&
							&TRIM(give_sum_formula(donor_list(donor_type_counter)%molecule_type_index))
							DO i=1,donor_list(donor_type_counter)%N_atoms,1
								WRITE(4,FMT='("      Atom index ",I0," (",A,")")')&
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

		END SUBROUTINE report_species

		SUBROUTINE perform_speciation_analysis()
		IMPLICIT NONE
			CALL initialise_speciation()
			IF ((ERROR_CODE/=33).AND.(ERROR_CODE/=159).AND.(ERROR_CODE/=157)) THEN
				WRITE(*,'(" Using ",I0," timesteps in intervals of ",I0," for averaging.")')&
				&MAX((nsteps-1+sampling_interval)/sampling_interval,0),sampling_interval
				WRITE(*,'(" Starting speciation analysis.")')
				CALL refresh_IO()
				CALL trajectory_speciation_analysis
				IF (neighbour_atom_overflow>0) CALL report_error(160,exit_status=neighbour_atom_overflow)
				IF (species_overflow>0) CALL report_error(161,exit_status=species_overflow)
				CALL report_species()
				CALL finalise_speciation()
			ELSE
				ERROR_CODE=ERROR_CODE_DEFAULT
			ENDIF
		END SUBROUTINE perform_speciation_analysis

END MODULE SPECIATION
!--------------------------------------------------------------------------------------------------------------------------------!