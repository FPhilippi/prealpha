
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
	LOGICAL,PARAMETER :: dumpfirst_default=.TRUE.
	LOGICAL,PARAMETER :: dumplast_default=.TRUE.
	!variables
	INTEGER :: nsteps !how many steps to use from trajectory
	INTEGER :: sampling_interval=1
	INTEGER :: maximum_number_of_species=10!how many different species do we want to allow?
	INTEGER :: maximum_number_of_neighbour_molecules=10
	INTEGER :: N_acceptor_types
	INTEGER :: N_donor_types
	INTEGER :: species_overflow
	REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: squared_cutoff_list !the squared cutoffs. The dimensions are: entry of acceptor_list -- entry of acceptor atom_index_list -- entry of donor_list -- entry of donor atom_index_list
	TYPE,PRIVATE :: species
		LOGICAL :: entry_used=.FALSE. !False for empty entries. occurrences etc are only allocated for non-empty entries.
		INTEGER(KIND=WORKING_PRECISION) :: problematic_occurrences ! plus one if not the same as "last_occurrence", minus one otherwise
		INTEGER :: total_number_of_neighbour_atoms ! For example, "4" for Li[DME]2
		INTEGER(KIND=WORKING_PRECISION) :: occurrences ! For example, "I have found this species 34857 times"
		INTEGER :: timestep_first_occurrence, timestep_last_occurrence
		INTEGER :: occurrence_entry_counter
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: last_occurrence !first dimension counts the molecules from 1 (acceptor) up to (maximum_number_of_neighbour_molecules*N_donor_types)+1, The second dimension INDEXES molecule_type_index and molecule_index, which are then saved as the integer value of the array.
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: first_occurrence !first dimension counts the molecules from 1 (acceptor) up to (maximum_number_of_neighbour_molecules*N_donor_types)+1. The second dimension INDEXES molecule_type_index and molecule_index, which are then saved as the integer value of the array.
		INTEGER,DIMENSION(:),ALLOCATABLE :: neighbour_entry_count !first dimension INDEXES the donor molecule types. this is a list of how many entries have been used in neighbouring_atoms.
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: neighbouring_atoms ! first dimension INDEXES the donor molecule types, second dimension INDEXES the distinct neighbours of this type.
		!for example, if there is a species where one molecule of type 3 donates 1 atom and another one donates 4 atoms, then the entry here would be (3,1)=4, (3,2)=1, and (3,3:maximum_number_of_neighbour_molecules)=0.
		!If there are two molecules of type 1 and three molecules of type 4, each with one atom, then the entry would become:
		!(1,1)=1, (1,2)=1, (4,3)=1, (4,4)=1, (4,5)=1, and all the remaining ones are zero. 
	END TYPE species
	TYPE,PRIVATE :: donors_and_acceptors
		INTEGER :: N_atoms!array size of atom_index_list and cutoff_list
		INTEGER :: molecule_type_index
		INTEGER,DIMENSION(:),ALLOCATABLE :: atom_index_list !counting all relevant atom_index up to N_atoms
		REAL,DIMENSION(:),ALLOCATABLE :: cutoff_list ! CONTRIBUTIONS to the cutoff, i.e. just the radius, counting up to N_atoms
		INTEGER :: Number_of_logged_species_in_acceptor_list
		TYPE(species),DIMENSION(:),ALLOCATABLE :: list_of_all_species !First dimension goes from 1 to Number_of_logged_species_in_acceptor_list and is allocated up to maximum_number_of_species
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
			maximum_number_of_neighbour_molecules=maximum_number_of_neighbour_molecules_default
			dumpfirst=dumpfirst_default
			dumplast=dumplast_default
			N_acceptor_types=0
			N_donor_types=0
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
					acceptor_list(n)%list_of_all_species(:)%entry_used=.FALSE.
					acceptor_list(n)%list_of_all_species(:)%problematic_occurrences=0
					acceptor_list(n)%list_of_all_species(:)%occurrences=0
					acceptor_list(n)%list_of_all_species(:)%timestep_last_occurrence=0
					acceptor_list(n)%list_of_all_species(:)%timestep_first_occurrence=0
					acceptor_list(n)%list_of_all_species(:)%total_number_of_neighbour_atoms=0
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
				ENDIF
				IF (ALLOCATED(squared_cutoff_list)) CALL report_error(0)
				ALLOCATE(squared_cutoff_list&
				&(N_acceptor_types,acceptor_maxatoms,N_donor_types,donor_maxatoms),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
				DO n=1,N_acceptor_types,1
					DO m=1,N_donor_types,1
						DO a=1,acceptor_list(n)%N_atoms
							DO b=1,donor_list(n)%N_atoms
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
				CHARACTER(LEN=32) :: inputstring
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  End-of-file condition in ",TRIM(FILENAME_SPECIATION_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_SPECIATION_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("maximum_number_of_neighbours","neighbours","N_neighbours")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,maximum_number_of_neighbour_molecules
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
								&"  setting 'maximum_number_of_neighbour_molecules' to default (=",maximum_number_of_neighbour_molecules_default,")"
								maximum_number_of_neighbour_molecules=maximum_number_of_neighbour_molecules_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'maximum_number_of_neighbour_molecules' to ",&
								&maximum_number_of_neighbour_molecules
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
						CASE ("dumpfirst","dump_first")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,dumpfirst
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								dumpfirst=dumpfirst_default
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&"   setting 'dumpfirst' to default (",dumpfirst,")"
							ELSE
								IF (dumpfirst) THEN
									WRITE(*,'("   write the first occurrence of each species.")')
								ELSE
									WRITE(*,'("   do not write the first occurrence of each species.")')
								ENDIF
							ENDIF
						CASE ("dumplast","dump_last")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,dumplast
							IF (ios/=0) THEN
								CALL report_error(158,exit_status=ios)
								dumplast=dumplast_default
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&"   setting 'dumplast' to default (",dumplast,")"
							ELSE
								IF (dumplast) THEN
									WRITE(*,'("   write the last occurrence of each species.")')
								ELSE
									WRITE(*,'("   do not write the last occurrence of each species.")')
								ENDIF
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
				DO n=1,N_molecule_types,1
					acceptor_list(n)%molecule_type_index=DA_clipboard(n)%molecule_type_index
					acceptor_list(n)%N_atoms=DA_clipboard(n)%N_atoms
					IF (ALLOCATED(acceptor_list(n)%atom_index_list)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%atom_index_list(acceptor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					IF (ALLOCATED(acceptor_list(n)%cutoff_list)) CALL report_error(0)
					ALLOCATE(acceptor_list(n)%cutoff_list(acceptor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					!initialise
					acceptor_list(n)%atom_index_list(:)=0
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
								IF (acceptor_list(counter)%atom_index_list(atom_counter)==0) THEN
									!found an empty spot for the atom index
									acceptor_list(counter)%atom_index_list(atom_counter)=atom_index
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
							&acceptor_list(n)%atom_index_list(atom_counter),&
							&TRIM(give_element_symbol(acceptor_list(n)%molecule_type_index,acceptor_list(n)%atom_index_list(atom_counter)))
							IF (acceptor_list(n)%cutoff_list(atom_counter)<-0.1) THEN
								acceptor_list(n)%cutoff_list(atom_counter)=&
								&covalence_radius(&
								&TRIM(give_element_symbol(acceptor_list(n)%molecule_type_index,acceptor_list(n)%atom_index_list(atom_counter))))
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
					ALLOCATE(donor_list(n)%atom_index_list(donor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					IF (ALLOCATED(donor_list(n)%cutoff_list)) CALL report_error(0)
					ALLOCATE(donor_list(n)%cutoff_list(donor_list(n)%N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					!initialise
					donor_list(n)%atom_index_list(:)=0
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
								IF (donor_list(counter)%atom_index_list(atom_counter)==0) THEN
									!found an empty spot for the atom index
									donor_list(counter)%atom_index_list(atom_counter)=atom_index
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
							&donor_list(n)%atom_index_list(atom_counter),&
							&TRIM(give_element_symbol(donor_list(n)%molecule_type_index,donor_list(n)%atom_index_list(atom_counter)))
							IF (donor_list(n)%cutoff_list(atom_counter)<-0.1) THEN
								donor_list(n)%cutoff_list(atom_counter)=&
								&covalence_radius(&
								&TRIM(give_element_symbol(donor_list(n)%molecule_type_index,donor_list(n)%atom_index_list(atom_counter))))
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
		END SUBROUTINE finalise_speciation

		SUBROUTINE trajectory_speciation_analysis()
		IMPLICIT NONE
		INTEGER :: stepcounter,acceptor_molecule_index,indexcounter,n_acceptor,n_donor,neighbour_counter
		INTEGER :: donor_molecule_index,donor_indexcounter,counter,allocstatus,deallocstatus,i,j,sorting_clipboard
		INTEGER :: neighbour_overflow,first_occurrence_entry_counter,entry_counter,previous_entry
		LOGICAL :: new_species_found,existing_species,problematic_occurrence
		TYPE(species) :: species_clipboard
		!reset the number of species logged to zero just in case.
		acceptor_list(:)%Number_of_logged_species_in_acceptor_list=0
		!allocate memory for species clipboard
		IF (ALLOCATED(species_clipboard%neighbouring_atoms)) CALL report_error(0)
		ALLOCATE(species_clipboard%neighbouring_atoms(&
		&N_donor_types,maximum_number_of_neighbour_molecules),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		!allocate memory for the number of neighbour entries
		IF (ALLOCATED(species_clipboard%neighbour_entry_count)) CALL report_error(0)
		ALLOCATE(species_clipboard%neighbour_entry_count(&
		&N_donor_types),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		neighbour_overflow=0
		species_overflow=0
		!also reserve some space for the potential "first" and "last"
		IF ((dumpfirst).OR.(dumplast)) THEN
			IF (ALLOCATED(species_clipboard%first_occurrence)) CALL report_error(0)
			ALLOCATE(species_clipboard%first_occurrence(&
			&(maximum_number_of_neighbour_molecules*N_donor_types)+1,3),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
		ENDIF
		!Iterate over all acceptors. Usually, this is probably one.
		DO n_acceptor=1,N_acceptor_types,1
			!average over timesteps, a few should be enough.
			!Parallelisation will be horrible probably.
			DO stepcounter=1,nsteps,sampling_interval
				!per ACCEPTOR molecule index, go over all its atoms and check if there are neighbour contributions.
				DO acceptor_molecule_index=1,give_number_of_molecules_per_step(acceptor_list(n_acceptor)%molecule_type_index),1
					!everything inside THIS loop counts as neighbours towards one species!
					!Thus, at the end, we need to check if this species exist and update accordingly.
					!Furthermore, at this stage, we need to initialise the species_clipboard
					IF ((dumpfirst).OR.(dumplast)) THEN
						species_clipboard%neighbouring_atoms(:,:)=0
					ENDIF
					!The species clipboard should not be updated inside the atom loop - because I want each acceptor MOLECULE to be considered as a whole.
					species_clipboard%entry_used=.FALSE.
					first_occurrence_entry_counter=0
					species_clipboard%total_number_of_neighbour_atoms=0
					DO indexcounter=1,acceptor_list(n_acceptor)%N_atoms,1
						!previous_entry after each loop through n_donor, previous_entry will be updated.
						!Thus, to recognise the first set of donors as such, we need to initialise this to zero:
						previous_entry=0
						!Now we have a specific atom of the reference atom...
						!	... molecule type index is: acceptor_list(n_acceptor)%molecule_type_index
						!	... molecule index is: acceptor_molecule_index
						!	... atom index is: acceptor_list(n_acceptor)%atom_index_list(indexcounter)
						!THUS, iterate over all possible neighbour (donor) molecules...
						DO n_donor=1,N_donor_types,1
							species_clipboard%neighbour_entry_count(n_donor)=0
							DO donor_molecule_index=1,&
							&give_number_of_molecules_per_step(donor_list(n_donor)%molecule_type_index),1
								!now we are down to one specific donor molecule.
								!We need to check that we do not accidentally include intramolecular contributions.
								IF ((donor_list(n_donor)%molecule_type_index==acceptor_list(n_acceptor)%molecule_type_index)&
								&.AND.(donor_molecule_index==acceptor_molecule_index)) CYCLE
								!most donor molecules will not be neighbouring. Thus, for those non-neighbours, make sure the loop is exited with as little evaluations as possible
								neighbour_counter=0
								DO donor_indexcounter=1,donor_list(n_donor)%N_atoms,1
									!calculate squared distance!
									!This is the part that really hurts...
									!Check if it is smaller than the cutoff
									IF (give_smallest_atom_distance_squared&
									&(stepcounter,stepcounter,&
									&acceptor_list(n_acceptor)%molecule_type_index,&
									&donor_list(n_donor)%molecule_type_index,&
									&acceptor_molecule_index,donor_molecule_index,&
									&acceptor_list(n_acceptor)%atom_index_list(indexcounter),&
									&donor_list(n_donor)%atom_index_list(donor_indexcounter))<&
									&(squared_cutoff_list(&
									&n_acceptor,indexcounter,n_donor,donor_indexcounter))) THEN
										neighbour_counter=neighbour_counter+1
									ENDIF
								ENDDO
								!If there was a neighbour atom, then neighbour_counter is >0.
								!We now need to append an entry to neighbouring_atoms
								IF (neighbour_counter>0) THEN
									species_clipboard%entry_used=.TRUE.
									!find a free entry
									species_clipboard%neighbour_entry_count(n_donor)=&
									&species_clipboard%neighbour_entry_count(n_donor)+1
									IF (species_clipboard%neighbour_entry_count(n_donor)<maximum_number_of_neighbour_molecules) THEN
										species_clipboard%total_number_of_neighbour_atoms=species_clipboard%total_number_of_neighbour_atoms+neighbour_counter
										species_clipboard%neighbouring_atoms(n_donor,species_clipboard%neighbour_entry_count(n_donor))=neighbour_counter
										!keep track of the identities of the involved molecules if required.
										IF ((dumpfirst).OR.(dumplast)) THEN
											IF (first_occurrence_entry_counter<1) THEN
												first_occurrence_entry_counter=1
												species_clipboard%timestep_first_occurrence=stepcounter
												!save the acceptor
												species_clipboard%first_occurrence(1,1)=&
												&acceptor_list(n_acceptor)%molecule_type_index
												species_clipboard%first_occurrence(1,2)=&
												&acceptor_molecule_index
											ENDIF
											!save the donor
											first_occurrence_entry_counter=first_occurrence_entry_counter+1
											species_clipboard%first_occurrence(first_occurrence_entry_counter,1)=&
											&donor_list(n_donor)%molecule_type_index
											species_clipboard%first_occurrence(first_occurrence_entry_counter,2)=&
											&donor_molecule_index
										ENDIF
									ELSE
										!keep track of overflows
										neighbour_overflow=neighbour_overflow+1
									ENDIF
								ENDIF
							ENDDO
							!This is a good place to sort...
							!specifically, if species_clipboard%neighbour_entry_count(n_donor) is bigger than previous_entry,
							!then some entries have been added. If only one has been added, then
							!species_clipboard%neighbour_entry_count(n_donor)=previous_entry+1. We only need to sort if there are at least 2 elements.
							!Thus, we only need to sort (or check for the need of sorting), if:
							IF ((species_clipboard%neighbour_entry_count(n_donor)-previous_entry)>1) THEN
								!Now, we need to sort the elements from "previous_entry+1" to "species_clipboard%neighbour_entry_count(n_donor)"
								!I chose insertionsort over quicksort for stability
								DO i=previous_entry+2,species_clipboard%neighbour_entry_count(n_donor),1
									sorting_clipboard=species_clipboard%neighbouring_atoms(n_donor,i)
									j=i
									DO WHILE ((j>previous_entry+1).AND.(species_clipboard%neighbouring_atoms(n_donor,j-1)<sorting_clipboard))
										!we only need to swap the neighbour atom counts, since n_donor is already sorted by design, and the counter variable is meaningless anyways
										species_clipboard%neighbouring_atoms(n_donor,j)=species_clipboard%neighbouring_atoms(n_donor,j-1)
										j=j-1
									ENDDO
									species_clipboard%neighbouring_atoms(n_donor,j)=sorting_clipboard
								ENDDO
							ENDIF
							!Remember to update for the next round! This is necessary - even if we did not sort.
							previous_entry=species_clipboard%neighbour_entry_count(n_donor)
						ENDDO
					ENDDO
					!Check if there was any neighbour at all
					IF (species_clipboard%entry_used) THEN
						!at this point, the species_clipboard is filled.
						! The sorting should already have brought "species_clipboard%neighbouring_atoms" into a consistent standard form.
						!Specifically, the entries are 1) sorted by donor types and 2) their number of neighbour atoms in descending order.
						!see if we have this species already in the "acceptor_list(n_acceptor)%list_of_all_species(????),"
						!and accordingly either update the "occurrences" count or make a new entry.
						new_species_found=.TRUE.
						i=0
						!The following loop continues running under the assumption that species_clipboard contains a previously unobserved species.
						DO WHILE ((new_species_found).AND.&
						&(i<acceptor_list(n_acceptor)%Number_of_logged_species_in_acceptor_list))
							i=i+1
							!The list of neighbour entry counts serves as a checksum of sorts here,
							!Because most entries will differ in it and do not need a detailed comparison.
							IF (ALL(acceptor_list(n_acceptor)%list_of_all_species(i)%&
							&neighbour_entry_count(:)==species_clipboard%neighbour_entry_count(:)))THEN
								!The neighbour entries are the same - we should thus compare all the entries in neighbouring_atoms.
								!We need two loops. We will run them under the assumption that there is NO difference between the entries.
								!the reason for this is that, even if this is an existing species,
								!it will be only one entry among others. Thus, switching to existing_species
								!allows us to exit the inner loops as fast as possible.
								!Thus, we initialise:
								existing_species=.TRUE.
								!The first loop runs from 1 up to N_donor_types
								n_donor=0
								DO WHILE ((existing_species).AND.(n_donor<N_donor_types))
									n_donor=n_donor+1
									!The second loop runs, for each donor type, from 1 up to neighbour_entry_count.
									!we can use species_clipboard here because we already checked that its neighbouring_atoms list has the same dimensions as list_of_all_species.
									entry_counter=0
									DO WHILE ((existing_species).AND.&
									&(entry_counter<species_clipboard%neighbour_entry_count(n_donor)))
										entry_counter=entry_counter+1
										!now, make the comparison!
										IF (.NOT.(species_clipboard%neighbouring_atoms(n_donor,entry_counter)==&
										&acceptor_list(n_acceptor)%list_of_all_species(i)%&
										&neighbouring_atoms(n_donor,entry_counter))) THEN
											!They are not the same, contrary to our assumption - exit the loops.
											existing_species=.FALSE.
										ENDIF
									ENDDO
								ENDDO
								!At this point, we need to check if a difference has been found in the preceding two while loops or not.
								IF (existing_species) THEN
									new_species_found=.FALSE.
								ENDIF
							ENDIF
						ENDDO
						IF (new_species_found) THEN
							!sanity check
							IF (acceptor_list(n_acceptor)%Number_of_logged_species_in_acceptor_list/=i) CALL report_error(0)
							!check for overflows
							IF (maximum_number_of_species==i) THEN
								species_overflow=species_overflow+1
							ELSE!NO overflow - continue adding entries
								i=i+1
								acceptor_list(n_acceptor)%Number_of_logged_species_in_acceptor_list=i
								!now we need to add the entry at acceptor_list(n_acceptor)%Number_of_logged_species_in_acceptor_list
								!(which is the same as variable i)
								acceptor_list(n_acceptor)%list_of_all_species(i)%occurrences=1
								acceptor_list(n_acceptor)%list_of_all_species(i)%&
								&total_number_of_neighbour_atoms=species_clipboard%&
								&total_number_of_neighbour_atoms
								!the entry_used is not technically needed for the acceptor_list, but I kept it for consistency with species_clipboard
								acceptor_list(n_acceptor)%list_of_all_species(i)%entry_used=.TRUE.
								!Need to allocate and initialise neighbour_entry_count
								IF (ALLOCATED(acceptor_list(n_acceptor)%list_of_all_species(i)%&
								&neighbour_entry_count)) CALL report_error(0)
								ALLOCATE(acceptor_list(n_acceptor)%list_of_all_species(i)%&
								&neighbour_entry_count(&
								&N_donor_types),STAT=allocstatus)
								IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
								!Need to allocate and initialise neighbouring_atoms
								IF (ALLOCATED(acceptor_list(n_acceptor)%list_of_all_species(i)%&
								&neighbouring_atoms)) CALL report_error(0)
								ALLOCATE(acceptor_list(n_acceptor)%list_of_all_species(i)%&
								&neighbouring_atoms(N_donor_types,maximum_number_of_neighbour_molecules&
								&),STAT=allocstatus)
								IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
								acceptor_list(n_acceptor)%list_of_all_species(i)%neighbour_entry_count(:)=&
								&species_clipboard%neighbour_entry_count(:)
								!and, most importantly, the neighbouring atoms.
								acceptor_list(n_acceptor)%list_of_all_species(i)%neighbouring_atoms(:,:)=&
								&species_clipboard%neighbouring_atoms(:,:)
								!We will allocate the memory for both first_occurrence and last_occurrence as required.
								!If we are dumping the first one, then we need to allocate and transfer.
								IF (dumpfirst) THEN
									acceptor_list(n_acceptor)%list_of_all_species(i)%&
									&timestep_first_occurrence=stepcounter
									!allocate memory for first_occurrence
									IF (ALLOCATED(acceptor_list(n_acceptor)%list_of_all_species(i)%&
									&first_occurrence)) CALL report_error(0)
									ALLOCATE(acceptor_list(n_acceptor)%&
									&list_of_all_species(i)%first_occurrence(&
									&(maximum_number_of_neighbour_molecules*N_donor_types)+1,3),STAT=allocstatus)
									IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
									!transfer...
									DO j=1,first_occurrence_entry_counter,1
										acceptor_list(n_acceptor)%list_of_all_species(i)%first_occurrence(j,:)=&
										&species_clipboard%first_occurrence(j,:)
									ENDDO
								ENDIF
								!If we are dumping the last one, then we need to allocate and transfer.
								IF (dumplast) THEN
									acceptor_list(n_acceptor)%list_of_all_species(i)%&
									&timestep_last_occurrence=stepcounter
									!allocate memory for last_occurrence
									IF (ALLOCATED(acceptor_list(n_acceptor)%list_of_all_species(i)%&
									&last_occurrence)) CALL report_error(0)
									ALLOCATE(acceptor_list(n_acceptor)%&
									&list_of_all_species(i)%last_occurrence(&
									&(maximum_number_of_neighbour_molecules*N_donor_types)+1,3),STAT=allocstatus)
									IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
									!transfer...
									DO j=1,first_occurrence_entry_counter,1
										acceptor_list(n_acceptor)%list_of_all_species(i)%last_occurrence(j,:)=&
										&species_clipboard%first_occurrence(j,:)
									ENDDO
								ENDIF
								!transfer complete. update the entry counter for the global variable
								acceptor_list(n_acceptor)%list_of_all_species(i)%&
								&occurrence_entry_counter=first_occurrence_entry_counter
							ENDIF
						ELSE !NO new species found
							!species_clipboard is the same as acceptor_list(n_acceptor)%list_of_all_species(i)
							!Thus, update acceptor_list(n_acceptor)%list_of_all_species(i)
							acceptor_list(n_acceptor)%list_of_all_species(i)%occurrences=&
							&acceptor_list(n_acceptor)%list_of_all_species(i)%occurrences+1
							!if we are dumping the last one, then transfer species_clipboard%first_occurrence
							!to acceptor_list(n_acceptor)%list_of_all_species(i)%last_occurrence.
							!While doing so, check if there is a difference. If no, increase problematic_occurrences counter. if yes, decrease problematic_occurrences counter.
							problematic_occurrence=.TRUE.
							!sanity check
							IF (acceptor_list(n_acceptor)%list_of_all_species(i)%&
							&occurrence_entry_counter/=first_occurrence_entry_counter) CALL report_error(0)
							DO j=1,first_occurrence_entry_counter,1 !at this point I am just praying that j is an unused variable
								!check for difference... only necessary if none found so far.
								IF (problematic_occurrence) THEN
									IF (.NOT.(ALL(species_clipboard%first_occurrence(j,:)==&
									&acceptor_list(n_acceptor)%list_of_all_species(i)%&
									&last_occurrence(j,:)))) THEN
										!they are different.
										problematic_occurrence=.FALSE.
									ENDIF
								ENDIF
								!transfer...
								acceptor_list(n_acceptor)%list_of_all_species(i)%last_occurrence(j,:)=&
								&species_clipboard%first_occurrence(j,:)
							ENDDO
							IF (problematic_occurrence) THEN
								acceptor_list(n_acceptor)%list_of_all_species(i)%problematic_occurrences=&
								&acceptor_list(n_acceptor)%list_of_all_species(i)%problematic_occurrences+1
							ELSE
								acceptor_list(n_acceptor)%list_of_all_species(i)%problematic_occurrences=&
								&acceptor_list(n_acceptor)%list_of_all_species(i)%problematic_occurrences-1
							ENDIF
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		IF (neighbour_overflow>0) CALL report_error(160,exit_status=neighbour_overflow)
		IF (species_overflow>0) CALL report_error(161,exit_status=species_overflow)
		IF (ALLOCATED(species_clipboard%neighbouring_atoms)) THEN
			DEALLOCATE(species_clipboard%neighbouring_atoms,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(species_clipboard%neighbour_entry_count)) THEN
			DEALLOCATE(species_clipboard%neighbour_entry_count,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(species_clipboard%first_occurrence)) THEN
			DEALLOCATE(species_clipboard%first_occurrence,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		END SUBROUTINE trajectory_speciation_analysis

		SUBROUTINE report_species()
		IMPLICIT NONE
		INTEGER :: n_acceptor
			!Iterate over all acceptors. Usually, this is probably one.
			IF (N_acceptor_types==1) THEN
				WRITE(*,ADVANCE="NO",FMT='(" There was one acceptor molecule:")')
				CALL print_acceptor_summary(1)
			ELSE
				WRITE(*,'(" There were ",I0," acceptor molecules:")')
				DO n_acceptor=1,N_acceptor_types,1
					WRITE(*,ADVANCE="NO",FMT='("  Acceptor ",I0,"/",I0," was ")')&
					&n_acceptor,N_acceptor_types
					CALL print_acceptor_summary(n_acceptor)
				ENDDO
			ENDIF

			!Furthermore, if we are dumping structures, then we need to remove molecules
			!which are double. that can happen for more than one acceptor atom.

			CONTAINS

			SUBROUTINE print_acceptor_summary(n_acceptor_in)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: n_acceptor_in
			INTEGER :: entry_counter,allocstatus,deallocstatus,best_index,species_output_counter
			INTEGER :: n_donor,neighbour_entry_count,previous_entries,previous_neighbourcount,i
			INTEGER :: species_atomcount,ios
			REAL :: total_occurrence,best_value,distance
			REAL(KIND=WORKING_PRECISION) :: acceptor_centre(3),shift(3)
			REAL,DIMENSION(:),ALLOCATABLE :: occurrences
			CHARACTER(LEN=1024) :: filename_speciation_output
			CHARACTER(LEN=1024) :: filename_speciation_statistics
			CHARACTER(LEN=1024) :: species_output_header
				WRITE(filename_speciation_statistics,'(A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
				&//"acceptor",n_acceptor_in,"_speciation_statistics.dat"
				OPEN(UNIT=4,FILE=TRIM(filename_speciation_statistics),IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(4,'(" This file contains speciation statistics for acceptors coordinated by donors.")')
				WRITE(4,'(" The acceptor molecule was molecule type ",I0," (",A,").")')&
				&acceptor_list(n_acceptor_in)%molecule_type_index,&
				&TRIM(give_sum_formula(acceptor_list(n_acceptor_in)%molecule_type_index))
				WRITE(4,'("   For this acceptor, ",I0," different species were observed,")')&
				&acceptor_list(n_acceptor_in)%Number_of_logged_species_in_acceptor_list
				WRITE(4,'("   Which are now given in order of decreasing probability:")')
				IF (ALLOCATED(occurrences)) CALL report_error(0)
				ALLOCATE(occurrences(acceptor_list(n_acceptor_in)%&
				&Number_of_logged_species_in_acceptor_list),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(156,exit_status=allocstatus)
				WRITE(*,'(" molecule type ",I0," (",A,").")')&
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
				DO WHILE (ANY(acceptor_list(n_acceptor_in)%list_of_all_species&
				&(1:acceptor_list(n_acceptor_in)%&
				&Number_of_logged_species_in_acceptor_list)%entry_used))
					!There is an unused entry left!
					species_output_counter=species_output_counter+1
					!laboriously pick out the best one and report it
					best_index=1
					best_value=0.0
					DO entry_counter=1,acceptor_list(n_acceptor_in)%&
					&Number_of_logged_species_in_acceptor_list,1
						IF (.NOT.(acceptor_list(n_acceptor_in)%&
						&list_of_all_species(entry_counter)%entry_used)) CYCLE
						IF (occurrences(entry_counter)>best_value) THEN
							best_value=occurrences(entry_counter)
							best_index=entry_counter
						ENDIF
					ENDDO
					!Remove the best_index
					acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%entry_used=.FALSE.
					!Report the best index
					WRITE(4,'("   -------")')
					WRITE(4,'("   Species #(",I0,") in ",F0.1,"% of the cases.")')&
					&species_output_counter,100.0*FLOAT(&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)/total_occurrence
					WRITE(4,'("   The similarity across all was ",SP,F5.2,SS,".")')&
					&-1.0*FLOAT(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&problematic_occurrences)/FLOAT(&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
					IF (species_output_counter<6) THEN
						WRITE(*,'("     Species #(",I0,") in ",F0.1,"% of the cases.")')&
						&species_output_counter,100.0*FLOAT(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)/total_occurrence
					ELSE
						IF (species_output_counter==6) WRITE(*,'("   (Only the first 5 printed here)")')
					ENDIF
					WRITE(4,'("     This species has a total of ",I0," neighbour atoms:")')&
					&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
					&total_number_of_neighbour_atoms
					DO n_donor=1,N_donor_types,1
						neighbour_entry_count=acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%neighbour_entry_count(n_donor)
						IF (neighbour_entry_count>0) THEN
							IF (neighbour_entry_count>1) THEN
								WRITE(4,FMT='("       ",I0," neighbour molecules of type ",I0)')&
								&neighbour_entry_count,donor_list(n_donor)%molecule_type_index
							ELSE
								!Just one...
								WRITE(4,ADVANCE="NO",FMT='("       ",I0," molecule of type ",I0,", ")')&
								&neighbour_entry_count,donor_list(n_donor)%molecule_type_index
							ENDIF
							!report how many neighbours etc
							!Here, for every molcule type, I need to go through the neighbours and list them.
							!importantly, I want to sum all those with the same number of neighbour atoms together.
							!We initialise with the first entry...
							entry_counter=1
							previous_entries=1
							previous_neighbourcount=&
							&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&neighbouring_atoms(n_donor,entry_counter)
							!... and thus start the loop at the second one.
							!every time there is a difference, we report the previous batch
							DO WHILE (entry_counter<neighbour_entry_count)
								entry_counter=entry_counter+1
								IF (acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
								&neighbouring_atoms(n_donor,entry_counter)/=previous_neighbourcount) THEN
									!There was a difference in neighbour atom counts.
									!Thus, report all the ones we have so far.
									WRITE(4,'("         ",I0," of them donating ",I0," atoms.")')&
									&previous_entries,previous_neighbourcount
									!Reset "previous entries" to the current one
									previous_entries=1
									previous_neighbourcount=&
									&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
									&neighbouring_atoms(n_donor,entry_counter)
								ELSE
									!this one is the same as the previous! increase counter:
									previous_entries=previous_entries+1
								ENDIF
							ENDDO
							!importantly, since we only reported the previous one, here the last one too:
							IF (neighbour_entry_count>1) THEN
								IF (previous_entries>1) THEN
									IF (previous_neighbourcount>1) THEN
										WRITE(4,'("         ",I0," of them donating ",I0," atoms each.")')&
										&previous_entries,previous_neighbourcount
									ELSE
										WRITE(4,'("         ",I0," of them donating ",I0," atom each.")')&
										&previous_entries,previous_neighbourcount
									ENDIF
										
								ELSE
									IF (previous_neighbourcount>1) THEN
										WRITE(4,'("         ",I0," of them donating ",I0," atoms.")')&
										&previous_entries,previous_neighbourcount
									ELSE
										WRITE(4,'("         ",I0," of them donating ",I0," atom.")')&
										&previous_entries,previous_neighbourcount
									ENDIF
								ENDIF
							ELSE
								!Just one...
								IF (previous_entries>1) THEN
									CALL report_error(0)
								ELSE
									IF (previous_neighbourcount>1) THEN
										WRITE(4,'("donating ",I0," atoms.")')&
										&previous_neighbourcount
									ELSE
										WRITE(4,'("donating ",I0," atom.")')&
										&previous_neighbourcount
									ENDIF
								ENDIF
							ENDIF
						ENDIF
					ENDDO
					!print the first and last occurrence, if required
					!Note that we are still in the "best_index" loop
					IF (dumpfirst) THEN
						WRITE(filename_speciation_output,'(A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
						&//"acceptor",n_acceptor_in,"_species",species_output_counter,"_first.xyz"
						OPEN(UNIT=3,FILE=TRIM(filename_speciation_output),IOSTAT=ios)
						IF (ios/=0) CALL report_error(26,exit_status=ios)
						WRITE(species_output_header,'("Timestep=",I0,", Occurrence=",E9.3,A)')&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,FLOAT(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
						!get the total number of atoms
						species_atomcount=0
						DO i=1,acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%&
						&occurrence_entry_counter,1
							species_atomcount=species_atomcount+&
							&give_number_of_atoms_per_molecule(acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%&
							&first_occurrence(i,1))
						ENDDO
						!put the acceptor in the origin
						acceptor_centre(:)=give_center_of_mass(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_first_occurrence,&
						&acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%first_occurrence(1,1),&
						&acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%first_occurrence(1,2))
						!First, write the header
						WRITE(3,'(I0)')species_atomcount
						WRITE(3,'(A)')TRIM(species_output_header)
						!then, append all the molecules that are part of this species
						DO i=1,acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%&
						&occurrence_entry_counter,1
							distance=give_smallest_distance_squared&
							&(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&timestep_first_occurrence,acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&timestep_first_occurrence,&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%first_occurrence(1,1),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%first_occurrence(i,1),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%first_occurrence(1,2),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%first_occurrence(i,2),&
							&shift)
							CALL write_molecule(3,&
							&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&timestep_first_occurrence,acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%first_occurrence(i,1),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%first_occurrence(i,2),&
							&include_header=.FALSE.,translate_by=shift(:)-acceptor_centre(:))
						ENDDO
						WRITE(3,*)
						WRITE(3,*)
						ENDFILE 3
						CLOSE(UNIT=3)
					ENDIF
					IF (dumplast) THEN
						WRITE(filename_speciation_output,'(A,I0,A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))&
						&//"acceptor",n_acceptor_in,"_species",species_output_counter,"_last.xyz"
						OPEN(UNIT=3,FILE=TRIM(filename_speciation_output),IOSTAT=ios)
						IF (ios/=0) CALL report_error(26,exit_status=ios)
						WRITE(species_output_header,'("Timestep=",I0,", Occurrence=",E9.3,A)')&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,FLOAT(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%occurrences)
						!get the total number of atoms
						species_atomcount=0
						DO i=1,acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%&
						&occurrence_entry_counter,1
							species_atomcount=species_atomcount+&
							&give_number_of_atoms_per_molecule(acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%&
							&last_occurrence(i,1))
						ENDDO
						!put the acceptor in the origin
						acceptor_centre(:)=give_center_of_mass(&
						&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
						&timestep_last_occurrence,&
						&acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%last_occurrence(1,1),&
						&acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%last_occurrence(1,2))
						!First, write the header
						WRITE(3,'(I0)')species_atomcount
						WRITE(3,'(A)')TRIM(species_output_header)
						!then, append all the molecules that are part of this species
						DO i=1,acceptor_list(n_acceptor_in)%&
						&list_of_all_species(best_index)%&
						&occurrence_entry_counter,1
							distance=give_smallest_distance_squared&
							&(acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&timestep_last_occurrence,acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&timestep_last_occurrence,&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%last_occurrence(1,1),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%last_occurrence(i,1),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%last_occurrence(1,2),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%last_occurrence(i,2),&
							&shift)
							CALL write_molecule(3,&
							&acceptor_list(n_acceptor_in)%list_of_all_species(best_index)%&
							&timestep_last_occurrence,acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%last_occurrence(i,1),&
							&acceptor_list(n_acceptor_in)%&
							&list_of_all_species(best_index)%last_occurrence(i,2),&
							&include_header=.FALSE.,translate_by=shift(:)-acceptor_centre(:))
						ENDDO
						WRITE(3,*)
						WRITE(3,*)
						ENDFILE 3
						CLOSE(UNIT=3)
					ENDIF
				ENDDO
				WRITE(*,'(A,I0,A)') "   Detailed statistics written to '"&
				&//TRIM(ADJUSTL(OUTPUT_PREFIX))//"acceptor",n_acceptor_in,"_speciation_statistics.dat'"
				IF (dumpfirst) WRITE(*,'(A,I0,A)')&
				&"   First observed structure per species written to '"&
				&//TRIM(ADJUSTL(OUTPUT_PREFIX))//"acceptor"&
				&,n_acceptor_in,"_speciesX_first.xyz'"
				IF (dumplast) WRITE(*,'(A,I0,A)')&
				&"   Last observed structure per species written to '"&
				&//TRIM(ADJUSTL(OUTPUT_PREFIX))//"acceptor"&
				&,n_acceptor_in,"_speciesX_last.xyz'"
				IF (species_overflow>0) THEN
					WRITE(*,'("   (There was species overflow amounting to ",F0.1,"%)")')&
					&100.0*FLOAT(species_overflow)/total_occurrence
				ENDIF
				IF (ALLOCATED(occurrences)) THEN
					DEALLOCATE(occurrences,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
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
				CALL report_species()
				CALL finalise_speciation()
			ELSE
				ERROR_CODE=ERROR_CODE_DEFAULT
			ENDIF
		END SUBROUTINE perform_speciation_analysis

END MODULE SPECIATION
!--------------------------------------------------------------------------------------------------------------------------------!