
!This Module performs a detailed cluster analysis.
MODULE CLUSTER ! Copyright (C) !RELEASEYEAR! Frederik Philippi
	USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	INTEGER,PARAMETER :: nsteps_default=1 !how many steps to use from trajectory
	INTEGER,PARAMETER :: sampling_interval_default=1
	LOGICAL,PARAMETER :: print_statistics_default=.FALSE.
	LOGICAL,PARAMETER :: print_spectators_default=.FALSE.
	!variables
	CHARACTER (LEN=10) :: operation_mode="NONE"!operation mode of the cluster module.
	!operation_mode="PAIRS": needs pairs of atom indices and their fixed cutoffs
	!operation_mode="GLOBAL": needs a global cutoff operating on a set of atoms. The latter are specified by molecule_type_index and atom_index, and the code iterates over molecule_index.
	INTEGER :: nsteps=nsteps_default !how many steps to use from trajectory
	INTEGER :: sampling_interval=sampling_interval_default
	INTEGER :: N_cutoffs,N_atoms,N_steps_to_print,N_Pairs
	INTEGER :: bytes_needed
	LOGICAL :: print_statistics=print_statistics_default
	LOGICAL :: print_spectators=print_spectators_default !if T, print also the molecules which are not members of list_of_atoms
	LOGICAL :: custom_weights=.FALSE.
	TYPE :: cluster_atom
		INTEGER :: molecule_type_index
		INTEGER :: molecule_index
		INTEGER :: atom_index
		REAL :: atom_weight !usually, this would be the charge.
	END TYPE cluster_atom
	TYPE :: cluster_pair
		INTEGER :: list_of_atoms_starting_index_A
		INTEGER :: list_of_atoms_last_index_A
		INTEGER :: molecule_type_index_A
		INTEGER :: atom_index_A
		INTEGER :: list_of_atoms_starting_index_B
		INTEGER :: list_of_atoms_last_index_B
		INTEGER :: molecule_type_index_B
		INTEGER :: atom_index_B
	END TYPE cluster_pair
	TYPE(cluster_atom),DIMENSION(:),ALLOCATABLE :: list_of_atoms
	TYPE(cluster_pair),DIMENSION(:),ALLOCATABLE :: list_of_pairs
	INTEGER(KIND=DP),DIMENSION(:,:),ALLOCATABLE :: overall_occurrence
	REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE :: average_weights
	REAL,DIMENSION(:),ALLOCATABLE :: cutoff_list
	INTEGER(KIND=DP),DIMENSION(:),ALLOCATABLE :: ClusterCountDistributionFunction
	REAL,DIMENSION(:),ALLOCATABLE :: squared_cutoff_list
	INTEGER,DIMENSION(:),ALLOCATABLE :: steps_toprint
	!PRIVATE/PUBLIC declarations
	PUBLIC :: perform_cluster_analysis,user_cluster_input

	CONTAINS

		!WRITING input file to unit 8, which shouldn't be open.
		SUBROUTINE user_cluster_input&
		&(parallelisation_possible,parallelisation_requested,number_of_molecules,nsteps_in,filename_speciation)
		IMPLICIT NONE
		CHARACTER (LEN=*) :: filename_speciation
		LOGICAL,INTENT(INOUT) :: parallelisation_possible,parallelisation_requested
		INTEGER,INTENT(IN) :: nsteps_in,number_of_molecules
		REAL :: inputreal
		LOGICAL :: connected

		END SUBROUTINE user_cluster_input

		SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
		IMPLICIT NONE
			nsteps=nsteps_default
			sampling_interval=sampling_interval_default
			print_statistics=print_statistics_default
			custom_weights=.FALSE.
			bytes_needed=0
			operation_mode="NONE"
		END SUBROUTINE set_defaults

		!initialises the cluster module by reading the specified input file.
		SUBROUTINE initialise_cluster()
		IMPLICIT NONE
		INTEGER :: i,allocstatus
		LOGICAL :: file_exists,connected
		INTEGER :: ios,current_step
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
			INQUIRE(FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_CLUSTER_INPUT),EXIST=file_exists)
			IF (file_exists) THEN
				CALL set_defaults()
				IF (VERBOSE_OUTPUT) WRITE(*,*) "reading file '",TRIM(PATH_INPUT)//TRIM(FILENAME_CLUSTER_INPUT),"'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_INPUT)//TRIM(FILENAME_CLUSTER_INPUT),&
				&ACTION='READ',IOSTAT=ios)
				IF (ios/=0) CALL report_error(172,exit_status=ios)
				CALL read_first_pass()
				IF (ERROR_CODE==177) RETURN !couldn't read operation mode...
				!N_atoms, N_cutoffs, N_steps_to_print are initialised
				IF (N_atoms==0) THEN
					CALL report_error(173)
				ELSEIF (N_cutoffs==0) THEN
					CALL report_error(174)
				ELSE
					IF (TRIM(operation_mode)=="GLOBAL") THEN
						WRITE(*,'("   Expecting ",I0," atoms and ",I0," cutoffs.")') N_atoms,N_cutoffs
					ELSE
						WRITE(*,'("   Expecting ",I0," pairs and ",I0," atoms.")') N_pairs,N_atoms
					ENDIF
					!allocate memory which is needed now - that's why I do the firstpass!
					IF (TRIM(operation_mode)=="PAIRS") THEN
						ALLOCATE(list_of_pairs(N_pairs),STAT=allocstatus)
						bytes_needed=bytes_needed+N_pairs*8 !multiply with 4 bytes later
						IF (allocstatus/=0) THEN
							CALL report_error(171,exit_status=allocstatus)
							RETURN
						ENDIF
					ENDIF
					ALLOCATE(cutoff_list(N_cutoffs),STAT=allocstatus)
					bytes_needed=bytes_needed+N_cutoffs !multiply with 4 bytes later
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=allocstatus)
						RETURN
					ENDIF
					IF (N_steps_to_print>0) THEN
						CALL initialise_print_members()
						ALLOCATE(steps_toprint(N_steps_to_print),STAT=allocstatus)
						bytes_needed=bytes_needed+N_steps_to_print !multiply with 4 bytes later
						IF (allocstatus/=0) THEN
							CALL report_error(171,exit_status=allocstatus)
							RETURN
						ENDIF
					ENDIF
					ALLOCATE(list_of_atoms(N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=allocstatus)
						RETURN
					ENDIF
					CALL read_body()
					CLOSE(UNIT=3)
					ALLOCATE(squared_cutoff_list(N_cutoffs),STAT=allocstatus)
					bytes_needed=bytes_needed+N_cutoffs !multiply with 4 bytes later
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=allocstatus)
						RETURN
					ENDIF
					bytes_needed=bytes_needed+4*N_atoms !multiply with 4 bytes later
					!list_of_clusters cannot be larger than N_atoms (=all monomers...)
					!the following lines are for the variables local to each thread!
				 !$ IF (PARALLEL_OPERATION) THEN
				 !$ 	bytes_needed=bytes_needed+2*N_atoms*&
				 !$ 	&*OMP_get_num_threads()!Here we have the list_of_clusters
				 !$ 	bytes_needed=bytes_needed+OMP_get_num_threads()*N_atoms !this is for the cluster_number list
				 !$ ELSE
					bytes_needed=bytes_needed+2*N_atoms
					bytes_needed=bytes_needed+N_atoms
				 !$ ENDIF
					IF (print_statistics) THEN
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							ALLOCATE(ClusterCountDistributionFunction(N_cutoffs),STAT=allocstatus)
							IF (allocstatus/=0) THEN
								CALL report_error(171,exit_status=allocstatus)
								RETURN
							ENDIF
							ALLOCATE(average_weights(N_atoms,N_cutoffs),STAT=allocstatus)
							IF (allocstatus/=0) THEN
								CALL report_error(171,exit_status=allocstatus)
								RETURN
							ENDIF
							ALLOCATE(overall_occurrence(N_atoms,N_cutoffs),STAT=allocstatus)
							IF (allocstatus/=0) THEN
								CALL report_error(171,exit_status=allocstatus)
								RETURN
							ENDIF
							bytes_needed=bytes_needed+2*2*N_atoms*N_cutoffs+2*N_cutoffs !multiply with 4 bytes later
						ELSE
							ALLOCATE(ClusterCountDistributionFunction(1),STAT=allocstatus)
							IF (allocstatus/=0) THEN
								CALL report_error(171,exit_status=allocstatus)
								RETURN
							ENDIF
							ALLOCATE(average_weights(N_atoms,1),STAT=allocstatus)
							IF (allocstatus/=0) THEN
								CALL report_error(171,exit_status=allocstatus)
								RETURN
							ENDIF
							ALLOCATE(overall_occurrence(N_atoms,1),STAT=allocstatus)
							IF (allocstatus/=0) THEN
								CALL report_error(171,exit_status=allocstatus)
								RETURN
							ENDIF
							bytes_needed=bytes_needed+2*2*N_atoms+2 !multiply with 4 bytes later
						ENDIF
						IF (.NOT.(custom_weights)) THEN
							WRITE(*,'(" no custom weights have been specified, using molecular charge.")')
							WRITE(*,'(" (Note: even for more than one atom, they are all assigned the full charge)")')
							DO i=1,N_atoms,1
								list_of_atoms(i)%atom_weight=&
								&FLOAT(give_charge_of_molecule(list_of_atoms(i)%molecule_type_index))
							ENDDO
						ENDIF
						ClusterCountDistributionFunction(:)=0
						average_weights(:,:)=0.0
						overall_occurrence(:,:)=0
					ENDIF
					IF (TRIM(operation_mode)=="GLOBAL") THEN
						CALL sort_cutoff_list(1,N_cutoffs,cutoff_list)
					ELSE
						WRITE(*,'(" The following pairs have been read:")')
						DO i=1,N_pairs
							WRITE(*,'("   Pair ",I0," out of ",I0,":")')i,N_pairs
							WRITE(*,'("     From atom #",I0,"(",A,") of molecule type #",I0,"(",A,")...")')&
							&list_of_pairs(i)%atom_index_A,&
							&TRIM(give_element_symbol(list_of_pairs(i)%molecule_type_index_A,&
							&list_of_pairs(i)%atom_index_A)),&
							&list_of_pairs(i)%molecule_type_index_A,&
							&TRIM(give_sum_formula(list_of_pairs(i)%molecule_type_index_A))
							WRITE(*,'("     ... to atom #",I0,"(",A,") of molecule type #",I0,"(",A,"),")')&
							&list_of_pairs(i)%atom_index_B,&
							&TRIM(give_element_symbol(list_of_pairs(i)%molecule_type_index_B,&
							&list_of_pairs(i)%atom_index_B)),&
							&list_of_pairs(i)%molecule_type_index_B,&
							&TRIM(give_sum_formula(list_of_pairs(i)%molecule_type_index_B))
							IF ((cutoff_list(i)<0.01).OR.(cutoff_list(i)>9.99)) THEN
								WRITE(*,'("     with cutoff = ",E10.3,".")')cutoff_list(i)
							ELSE
								WRITE(*,'("     with cutoff =",F5.2,".")')cutoff_list(i)
							ENDIF
						ENDDO
					ENDIF
					CALL sort_steps_to_print_list(1,N_steps_to_print,steps_toprint)
					DO i=1,N_cutoffs,1
						squared_cutoff_list(i)=cutoff_list(i)**2
					ENDDO
					IF (N_steps_to_print>1) THEN
						current_step=1
						DO	i=1,nsteps,sampling_interval
							IF (steps_toprint(current_step)==i) THEN
								current_step=current_step+1
							ELSE
								IF (i>steps_toprint(current_step)) THEN
									IF (VERBOSE_OUTPUT) THEN
										WRITE(*,'(A,I0,A,I0)') "  Using step #",i," instead of step #",steps_toprint(current_step)
									ENDIF
									steps_toprint(current_step)=i
									current_step=current_step+1
								ENDIF
							ENDIF
						ENDDO
						DO i=current_step,N_steps_to_print,1
							IF (steps_toprint(i)>nsteps) THEN
								IF (VERBOSE_OUTPUT) THEN
									WRITE(*,'(A,I0,A,I0)') "  Using step #",nsteps," instead of step #",steps_toprint(i)
								ENDIF
								steps_toprint(i)=nsteps
							ENDIF
						ENDDO
					ENDIF
				ENDIF
			ELSE
				CALL report_error(170)!No input - no output. easy as that.
			ENDIF

		CONTAINS

			SUBROUTINE read_first_pass()
			IMPLICIT NONE
			CHARACTER(LEN=33) :: inputstring
			INTEGER :: scansteps,inputinteger,inputinteger2,n
			INTEGER :: molecule_type_index,atom_index,molecule_type_index2,atom_index2
			REAL :: dummy,cutlow,cuthigh,cutoff_in
				REWIND 3
				READ(3,IOSTAT=ios,FMT=*) inputstring
				SELECT CASE (TRIM(inputstring))
				CASE ("PAIRS","Pairs","pairs","PAIR","Pair","pair")
					operation_mode="PAIRS"
					WRITE(*,*) "PAIRS operation mode:"
					WRITE(*,*) "Uses a set of pairs of atoms with their corresponding cutoffs."
				CASE ("GLOBAL","Global","global")
					operation_mode="GLOBAL"
					WRITE(*,*) "GLOBAL operation mode:"
					WRITE(*,*) "Uses a set of atoms and a global cutoff, applied to any combination of atoms."
				CASE DEFAULT
					CALL report_error(177)
					RETURN
				END SELECT
				molecule_type_index=0
				N_atoms=0
				N_Pairs=0
				N_cutoffs=0
				N_steps_to_print=0
				DO n=1,MAXITERATIONS,1
					READ(3,IOSTAT=ios,FMT=*) inputstring
					IF (ios/=0) EXIT
					SELECT CASE (TRIM(inputstring))
					CASE ("pair","Pair","add_pair")
						IF (TRIM(operation_mode)=="PAIRS") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,molecule_type_index,atom_index,&
							&molecule_type_index2,atom_index2,cutoff_in
							IF (ios==0) THEN
								IF ((valid_atom_index(molecule_type_index,atom_index))&
								&.AND.(valid_atom_index(molecule_type_index2,atom_index2))) THEN
									N_atoms=N_atoms+give_number_of_molecules_per_step(molecule_type_index)
									N_atoms=N_atoms+give_number_of_molecules_per_step(molecule_type_index2)
									N_pairs=N_pairs+1
									N_cutoffs=N_cutoffs+1
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("scan_cutoff")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,cutlow,cuthigh,scansteps
							IF (ios==0) THEN
								IF ((scansteps>0).OR.(cutlow<0.0).OR.(cutlow>cuthigh))&
								& N_cutoffs=N_cutoffs+scansteps
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("single_cutoff","cutoff")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,dummy
							IF (ios==0) THEN
								IF (dummy>0.0d0) N_cutoffs=N_cutoffs+1
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("print_step","printstep","dump","dump_step")
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
						IF (ios==0) THEN
							N_steps_to_print=N_steps_to_print+1
						ENDIF
					CASE ("add_molecule_type","add_molecule_type_index","add_molecule")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
							IF (ios==0) THEN
								IF (valid_molecule_type_index(inputinteger)) THEN
									N_atoms=N_atoms+&
									&give_number_of_atoms_per_molecule(inputinteger)*&
									&give_number_of_molecules_per_step(inputinteger)
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("add_atom_index","add_atom")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger,inputinteger2
							IF (ios==0) THEN
								IF (valid_atom_index(inputinteger,inputinteger2)) THEN
									N_atoms=N_atoms+give_number_of_molecules_per_step(inputinteger)
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("quit")
						EXIT
					END SELECT
				ENDDO
			END SUBROUTINE read_first_pass

			SUBROUTINE read_body()
			IMPLICIT NONE
			CHARACTER(LEN=33) :: inputstring
			INTEGER :: current_cutoff,current_printstep,current_atom,molecule_index
			INTEGER :: scansteps,n,inputinteger,inputinteger2,m,current_pair,pair_counter
			INTEGER :: molecule_type_index,atom_index,molecule_type_index2,atom_index2
			LOGICAL :: A_double,B_double
			REAL :: cutlow,cuthigh,interval,cutoff_in,current_weight,inputreal
				REWIND 3
				READ(3,IOSTAT=ios,FMT=*) inputstring
				current_cutoff=0
				current_pair=0
				current_printstep=0
				current_atom=0
				current_weight=0.0
				DO n=1,MAXITERATIONS,1
					READ(3,IOSTAT=ios,FMT=*) inputstring
					IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  End-of-file condition in ",TRIM(FILENAME_CLUSTER_INPUT)
					IF (ios/=0) THEN
						IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_CLUSTER_INPUT)
						EXIT
					ENDIF
					SELECT CASE (TRIM(inputstring))
					CASE ("scan_cutoff")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,cutlow,cuthigh,scansteps
							IF (ios/=0) THEN
								CALL report_error(170)
							ELSE
								IF ((scansteps>0).OR.(cutlow<0.0).OR.(cutlow>cuthigh)) THEN
									interval=(cuthigh-cutlow)/FLOAT(scansteps-1)
									IF (DEVELOPERS_VERSION) WRITE(*,ADVANCE="NO",FMT='("  ! Adding cutoffs:")')
									DO m=0,scansteps-1,1
										current_cutoff=current_cutoff+1
										cutoff_list(current_cutoff)=cutlow+FLOAT(m)*interval
										IF (DEVELOPERS_VERSION) WRITE(*,ADVANCE="NO",FMT='(" ",F0.2)') cutoff_list(current_cutoff)
									ENDDO
									IF (DEVELOPERS_VERSION) WRITE(*,*)
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A,I0,A)')&
									&"   ",scansteps,"/",N_cutoffs," cutoffs added via scan keyword."
								ELSE
									CALL report_error(175)
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("single_cutoff","cutoff")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,cutlow
							IF (ios/=0) THEN
								CALL report_error(170)
							ELSE
								IF (cutlow>0.0d0) THEN
									current_cutoff=current_cutoff+1
									cutoff_list(current_cutoff)=cutlow
								ELSE
									CALL report_error(175)
								ENDIF
							ENDIF
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A,I0,A)')&
							&"   ",1,"/",N_cutoffs," cutoff added via single_cutoff keyword."
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("print_step","printstep","dump","dump_step")
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
						IF (ios/=0) THEN
							CALL report_error(170)
						ELSE
							current_printstep=current_printstep+1
							steps_toprint(current_printstep)=inputinteger
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
							&"   Added step ",inputinteger," to the list of steps to print."
						ENDIF
					CASE ("custom_weight","custom_weights","weight","weights","set_weight","set_weights")
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,inputreal
						IF (ios/=0) THEN
							CALL report_error(170)
						ELSE
							current_weight=inputreal
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,E11.3,A)')&
							&"   Setting weight to",current_weight," for all atoms below until new 'custom_weight'."
						ENDIF
					CASE ("pair","Pair","add_pair")
						IF (TRIM(operation_mode)=="PAIRS") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,molecule_type_index,atom_index,&
							&molecule_type_index2,atom_index2,cutoff_in
							IF (ios==0) THEN
								IF ((valid_atom_index(molecule_type_index,atom_index))&
								&.AND.(valid_atom_index(molecule_type_index2,atom_index2))) THEN
									current_pair=current_pair+1
									list_of_pairs(current_pair)%molecule_type_index_A=molecule_type_index
									list_of_pairs(current_pair)%atom_index_A=atom_index
									list_of_pairs(current_pair)%molecule_type_index_B=molecule_type_index2
									list_of_pairs(current_pair)%atom_index_B=atom_index2
									cutoff_list(current_pair)=cutoff_in
									!Check for doubles present in previous pairs
									A_double=.FALSE.
									B_double=.FALSE.
									IF (current_pair>1) THEN
check_A:								DO pair_counter=1,current_pair-1,1
											IF ((list_of_pairs(current_pair)%molecule_type_index_A&
											&==list_of_pairs(pair_counter)%molecule_type_index_A).AND.&
											&(list_of_pairs(current_pair)%atom_index_A&
											&==list_of_pairs(pair_counter)%atom_index_A)) THEN
												!We found a double entry!
												N_atoms=N_atoms-give_number_of_molecules_per_step&
												&(list_of_pairs(pair_counter)%molecule_type_index_A)
												IF (VERBOSE_OUTPUT) WRITE(*,'(3A,I0)') "   Double entry for '",&
												&TRIM(give_element_symbol(list_of_pairs(current_pair)%molecule_type_index_A,&
												&list_of_pairs(current_pair)%atom_index_A)),"': reducing N_atoms to ",N_atoms
												IF (DEVELOPERS_VERSION) WRITE(*,FMT='("  ! Use indices of A for new A")')
												A_double=.TRUE.
												list_of_pairs(current_pair)%list_of_atoms_starting_index_A=&
												&list_of_pairs(pair_counter)%list_of_atoms_starting_index_A
												list_of_pairs(current_pair)%list_of_atoms_last_index_A=&
												&list_of_pairs(pair_counter)%list_of_atoms_last_index_A
												EXIT check_A
											ELSEIF ((list_of_pairs(current_pair)%molecule_type_index_A&
											&==list_of_pairs(pair_counter)%molecule_type_index_B).AND.&
											&(list_of_pairs(current_pair)%atom_index_A&
											&==list_of_pairs(pair_counter)%atom_index_B)) THEN
												!We found a double entry!
												N_atoms=N_atoms-give_number_of_molecules_per_step(&
												&list_of_pairs(pair_counter)%molecule_type_index_B)
												IF (VERBOSE_OUTPUT) WRITE(*,'(3A,I0)') "   Double entry for '",&
												&TRIM(give_element_symbol(list_of_pairs(current_pair)%molecule_type_index_A,&
												&list_of_pairs(current_pair)%atom_index_A)),"': reducing N_atoms to ",N_atoms
												IF (DEVELOPERS_VERSION) WRITE(*,FMT='("  ! Use indices of B for new A")')
												A_double=.TRUE.
												list_of_pairs(current_pair)%list_of_atoms_starting_index_A=&
												&list_of_pairs(pair_counter)%list_of_atoms_starting_index_B
												list_of_pairs(current_pair)%list_of_atoms_last_index_A=&
												&list_of_pairs(pair_counter)%list_of_atoms_last_index_B
												EXIT check_A
											ENDIF
										ENDDO check_A
check_B:								DO pair_counter=1,current_pair-1,1
											IF ((list_of_pairs(current_pair)%molecule_type_index_B&
											&==list_of_pairs(pair_counter)%molecule_type_index_A).AND.&
											&(list_of_pairs(current_pair)%atom_index_B&
											&==list_of_pairs(pair_counter)%atom_index_A)) THEN
												!We found a double entry!
												N_atoms=N_atoms-give_number_of_molecules_per_step(&
												&list_of_pairs(pair_counter)%molecule_type_index_A)
												IF (VERBOSE_OUTPUT) WRITE(*,'(3A,I0)') "   Double entry for '",&
												&TRIM(give_element_symbol(list_of_pairs(current_pair)%molecule_type_index_B,&
												&list_of_pairs(current_pair)%atom_index_B)),"': reducing N_atoms to ",N_atoms
												IF (DEVELOPERS_VERSION) WRITE(*,FMT='("  ! Use indices of A for new B")')
												B_double=.TRUE.
												list_of_pairs(current_pair)%list_of_atoms_starting_index_B=&
												&list_of_pairs(pair_counter)%list_of_atoms_starting_index_A
												list_of_pairs(current_pair)%list_of_atoms_last_index_B=&
												&list_of_pairs(pair_counter)%list_of_atoms_last_index_A
												EXIT check_B
											ELSEIF ((list_of_pairs(current_pair)%molecule_type_index_B&
											&==list_of_pairs(pair_counter)%molecule_type_index_B).AND.&
											&(list_of_pairs(current_pair)%atom_index_B&
											&==list_of_pairs(pair_counter)%atom_index_B)) THEN
												!We found a double entry!
												N_atoms=N_atoms-give_number_of_molecules_per_step(&
												&list_of_pairs(pair_counter)%molecule_type_index_B)
												IF (VERBOSE_OUTPUT) WRITE(*,'(3A,I0)') "   Double entry for '",&
												&TRIM(give_element_symbol(list_of_pairs(current_pair)%molecule_type_index_B,&
												&list_of_pairs(current_pair)%atom_index_B)),"': reducing N_atoms to ",N_atoms
												IF (DEVELOPERS_VERSION) WRITE(*,FMT='("  ! Use indices of B for new B")')
												B_double=.TRUE.
												list_of_pairs(current_pair)%list_of_atoms_starting_index_B=&
												&list_of_pairs(pair_counter)%list_of_atoms_starting_index_B
												list_of_pairs(current_pair)%list_of_atoms_last_index_B=&
												&list_of_pairs(pair_counter)%list_of_atoms_last_index_B
												EXIT check_B
											ENDIF
										ENDDO check_B
									ENDIF !check with old pairs complete, still need to check for A=B
									IF (.NOT.(A_double)) THEN
										!assign the molecules of type A to the list_of_atoms
										list_of_pairs(current_pair)%list_of_atoms_starting_index_A=current_atom+1
										list_of_pairs(current_pair)%list_of_atoms_last_index_A=current_atom+&
										&give_number_of_molecules_per_step(molecule_type_index)
										DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index)
											current_atom=current_atom+1
											list_of_atoms(current_atom)%molecule_type_index=molecule_type_index
											list_of_atoms(current_atom)%molecule_index=molecule_index
											list_of_atoms(current_atom)%atom_index=atom_index
											list_of_atoms(current_atom)%atom_weight=current_weight
										ENDDO
									ENDIF
									!Here we check for A=B
									IF (.NOT.(B_double)) THEN
										!why do I do this? B is not a double of old members, but can be a double of the new A!
										IF ((list_of_pairs(current_pair)%molecule_type_index_B&
										&==list_of_pairs(current_pair)%molecule_type_index_A).AND.&
										&(list_of_pairs(current_pair)%atom_index_B&
										&==list_of_pairs(current_pair)%atom_index_A)) THEN
											!We found a double entry!
											N_atoms=N_atoms-give_number_of_molecules_per_step(molecule_type_index)
											IF (VERBOSE_OUTPUT) WRITE(*,'(3A,I0)') "   Double entry for '",&
											&TRIM(give_element_symbol(list_of_pairs(current_pair)%molecule_type_index_A,&
											&list_of_pairs(current_pair)%atom_index_A)),"': reducing N_atoms to ",N_atoms
											IF (DEVELOPERS_VERSION) WRITE(*,FMT='("  ! Use indices of new A for new B")')
											B_double=.TRUE.
											list_of_pairs(current_pair)%list_of_atoms_starting_index_B=&
											&list_of_pairs(current_pair)%list_of_atoms_starting_index_A
											list_of_pairs(current_pair)%list_of_atoms_last_index_B=&
											&list_of_pairs(current_pair)%list_of_atoms_last_index_A
										ENDIF
									ENDIF
									IF (.NOT.(B_double)) THEN
										!assign the molecules of type B to the list_of_atoms
										list_of_pairs(current_pair)%list_of_atoms_starting_index_B=current_atom+1
										list_of_pairs(current_pair)%list_of_atoms_last_index_B=current_atom+&
										&give_number_of_molecules_per_step(molecule_type_index2)
										DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index2)
											current_atom=current_atom+1
											list_of_atoms(current_atom)%molecule_type_index=molecule_type_index2
											list_of_atoms(current_atom)%molecule_index=molecule_index
											list_of_atoms(current_atom)%atom_index=atom_index2
											list_of_atoms(current_atom)%atom_weight=current_weight
										ENDDO
									ENDIF
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("add_molecule_type","add_molecule_type_index","add_molecule")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
							IF (ios/=0) THEN
								CALL report_error(170)
							ELSE
								IF (valid_molecule_type_index(inputinteger)) THEN
									DO molecule_index=1,give_number_of_molecules_per_step(inputinteger)
										DO atom_index=1,give_number_of_atoms_per_molecule(inputinteger)
											current_atom=current_atom+1
											list_of_atoms(current_atom)%molecule_type_index=inputinteger
											list_of_atoms(current_atom)%molecule_index=molecule_index
											list_of_atoms(current_atom)%atom_index=atom_index
											list_of_atoms(current_atom)%atom_weight=current_weight
										ENDDO
									ENDDO
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
									&"   Added ",&
									give_number_of_molecules_per_step(inputinteger)*&
									&give_number_of_atoms_per_molecule(inputinteger)&
									&," atoms to be considered for the cluster analysis."
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("add_atom_index","add_atom")
						IF (TRIM(operation_mode)=="GLOBAL") THEN
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger,atom_index
							IF (ios/=0) THEN
								CALL report_error(170)
							ELSE
								IF (valid_atom_index(inputinteger,atom_index)) THEN
									DO molecule_index=1,give_number_of_molecules_per_step(inputinteger)
										current_atom=current_atom+1
										list_of_atoms(current_atom)%molecule_type_index=inputinteger
										list_of_atoms(current_atom)%molecule_index=molecule_index
										list_of_atoms(current_atom)%atom_index=atom_index
										list_of_atoms(current_atom)%atom_weight=current_weight
									ENDDO
									IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
									&"   Added ",&
									give_number_of_molecules_per_step(inputinteger)&
									&," atoms to be considered for the cluster analysis."
								ENDIF
							ENDIF
						ELSE
							CALL report_error(178)
						ENDIF
					CASE ("print_statistics","statistics")
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,print_statistics
						IF (ios/=0) THEN
							CALL report_error(170,exit_status=ios)
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
							&"   setting 'print_statistics' to default (=",print_statistics_default,")"
							print_statistics=print_statistics_default
						ELSE
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'print_statistics' to ",&
							&print_statistics
						ENDIF
					CASE ("print_spectators","spectators","print_spectator","spectator")
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,print_spectators
						IF (ios/=0) THEN
							CALL report_error(170,exit_status=ios)
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
							&"   setting 'print_spectators' to default (=",print_spectators_default,")"
							print_spectators=print_spectators_default
						ELSE
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'print_spectators' to ",&
							&print_spectators
						ENDIF
					CASE ("sampling_interval")
						BACKSPACE 3
						READ(3,IOSTAT=ios,FMT=*) inputstring,sampling_interval
						IF (ios/=0) THEN
							CALL report_error(170,exit_status=ios)
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
							CALL report_error(170,exit_status=ios)
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)')&
							&"  setting 'nsteps' to default (=",nsteps_default,")"
							nsteps=nsteps_default
						ELSE
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0)') "   setting 'nsteps' to ",nsteps
						ENDIF
					CASE ("quit")
						IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_CLUSTER_INPUT)
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
			END SUBROUTINE read_body

			!The following subroutine is a quicksort algorithm to sort the cutoff list.
			RECURSIVE SUBROUTINE sort_cutoff_list(left,right,cutoff_list_inout)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: left,right
			INTEGER :: a,b
			REAL,DIMENSION(*),INTENT(INOUT) :: cutoff_list_inout
			REAL pivot
			REAL :: cutoff_clip
				IF (left<right) THEN
					pivot=cutoff_list_inout(left)
					a=left
					b=right
					DO
						DO WHILE (cutoff_list_inout(a)<pivot)
							a=a+1
						ENDDO
						DO WHILE (cutoff_list_inout(b)>pivot)
							b=b-1
						ENDDO
						IF (a>=b) EXIT
						!swap elements, unless they're the same
						IF (cutoff_list_inout(a)==cutoff_list_inout(b)) THEN
							b=b-1
							IF (a==b) EXIT
						ELSE
							cutoff_clip=cutoff_list_inout(a)
							cutoff_list_inout(a)=cutoff_list_inout(b)
							cutoff_list_inout(b)=cutoff_clip
						ENDIF
					ENDDO
					CALL sort_cutoff_list(left,b,cutoff_list_inout)
					CALL sort_cutoff_list(b+1,right,cutoff_list_inout)
				ENDIF
			END SUBROUTINE sort_cutoff_list

			!The following subroutine is a quicksort algorithm to sort the steps_toprint list.
			RECURSIVE SUBROUTINE sort_steps_to_print_list(left,right,steps_to_print_list_inout)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: left,right
			INTEGER :: a,b
			INTEGER,DIMENSION(*),INTENT(INOUT) :: steps_to_print_list_inout
			INTEGER pivot
			REAL :: steps_to_print_clip
				IF (left<right) THEN
					pivot=steps_to_print_list_inout(left)
					a=left
					b=right
					DO
						DO WHILE (steps_to_print_list_inout(a)<pivot)
							a=a+1
						ENDDO
						DO WHILE (steps_to_print_list_inout(b)>pivot)
							b=b-1
						ENDDO
						IF (a>=b) EXIT
						!swap elements, unless they're the same
						IF (steps_to_print_list_inout(a)==steps_to_print_list_inout(b)) THEN
							b=b-1
							IF (a==b) EXIT
						ELSE
							steps_to_print_clip=steps_to_print_list_inout(a)
							steps_to_print_list_inout(a)=steps_to_print_list_inout(b)
							steps_to_print_list_inout(b)=steps_to_print_clip
						ENDIF
					ENDDO
					CALL sort_steps_to_print_list(left,b,steps_to_print_list_inout)
					CALL sort_steps_to_print_list(b+1,right,steps_to_print_list_inout)
				ENDIF
			END SUBROUTINE sort_steps_to_print_list

		END SUBROUTINE initialise_cluster

		!finalises the cluster module.
		SUBROUTINE finalise_cluster()
		IMPLICIT NONE
		INTEGER :: deallocstatus
			IF (TRIM(operation_mode)=="PAIRS") THEN
				IF (ALLOCATED(list_of_pairs)) THEN
					DEALLOCATE(list_of_pairs,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
			ENDIF
			IF (ALLOCATED(list_of_atoms)) THEN
				DEALLOCATE(list_of_atoms,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(cutoff_list)) THEN
				DEALLOCATE(cutoff_list,STAT=deallocstatus)
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
			IF (N_steps_to_print>0) THEN
				CALL finalise_print_members()
				IF (ALLOCATED(steps_toprint)) THEN
					DEALLOCATE(steps_toprint,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
			ENDIF
			IF (print_statistics) THEN
				IF (ALLOCATED(ClusterCountDistributionFunction)) THEN
					DEALLOCATE(ClusterCountDistributionFunction,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
				IF (ALLOCATED(average_weights)) THEN
					DEALLOCATE(average_weights,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
				IF (ALLOCATED(overall_occurrence)) THEN
					DEALLOCATE(overall_occurrence,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
			ENDIF
		END SUBROUTINE finalise_cluster

		!This subroutine looks at ALL connections between ALL members of list_of_atoms
		SUBROUTINE cluster_analysis_GLOBAL()
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: stepcounter,allocstatus,deallocstatus,m,n,current_cluster,i,ios,cluster_todelete,cluster_counts
		INTEGER :: cutoff_counter
		LOGICAL :: steps_printed
		REAL :: normalisation,CCF
		CHARACTER(LEN=1024) :: fstring,outstring
		TYPE :: one_cluster
			REAL :: cluster_weight
			INTEGER :: members
		END TYPE one_cluster
		TYPE(one_cluster),DIMENSION(:),ALLOCATABLE :: list_of_clusters
		INTEGER,DIMENSION(:),ALLOCATABLE :: cluster_number
		steps_printed=.FALSE.
			!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL)))&
			!$OMP PRIVATE(list_of_clusters,cluster_number,m,n,i,current_cluster,cutoff_counter,cluster_todelete,cluster_counts)
			!$OMP SINGLE
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (GLOBAL cluster analysis)"
		 !$ 	CALL timing_parallel_sections(.TRUE.)
		 !$ ENDIF
			IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
			!$OMP END SINGLE
			ALLOCATE(list_of_clusters(N_atoms),STAT=allocstatus)
			IF (allocstatus/=0) THEN
				CALL report_error(171,exit_status=allocstatus)
			ENDIF
			ALLOCATE(cluster_number(N_atoms),STAT=allocstatus)
			IF (allocstatus/=0) THEN
				CALL report_error(171,exit_status=allocstatus)
			ENDIF
			!$OMP DO SCHEDULE(STATIC,1) ORDERED
			DO stepcounter=1,nsteps,sampling_interval
				current_cluster=0
				cluster_number(:)=0
				list_of_clusters(:)%cluster_weight=0.0
				list_of_clusters(:)%members=0
				!Iterate over cutoffs...
				DO cutoff_counter=1,N_cutoffs,1
					DO m=1,N_atoms
						DO n=m,N_atoms
							IF (m==n) CYCLE
							!do not need to check pairs that belong together already
							IF (((cluster_number(m)==0).OR.(cluster_number(n)==0))&
							&.OR.(cluster_number(m)/=cluster_number(n))) THEN
								!check if closer than cutoff
								IF (give_smallest_atom_distance_squared&
								&(stepcounter,stepcounter,&
								&list_of_atoms(m)%molecule_type_index,list_of_atoms(n)%molecule_type_index,&
								&list_of_atoms(m)%molecule_index,list_of_atoms(n)%molecule_index,&
								&list_of_atoms(m)%atom_index,list_of_atoms(n)%atom_index)&
								&<squared_cutoff_list(cutoff_counter)) THEN
									IF (cluster_number(m)==0) THEN
									! m is a monomer!
										IF (cluster_number(n)==0) THEN
										! n is a monomer!
											!--> appropriate action: assign new cluster.
											current_cluster=current_cluster+1
											cluster_number(m)=current_cluster
											cluster_number(n)=current_cluster
											list_of_clusters(current_cluster)%cluster_weight=&
											&list_of_atoms(m)%atom_weight+&
											&list_of_atoms(n)%atom_weight
											list_of_clusters(current_cluster)%members=2
										ELSE
										! n is NOT a monomer!
											!--> appropriate action: add entry "m" to cluster "n"
											cluster_number(m)=cluster_number(n)
											list_of_clusters(cluster_number(n))%cluster_weight=&
											&list_of_clusters(cluster_number(n))%cluster_weight+&
											&list_of_atoms(m)%atom_weight
											list_of_clusters(cluster_number(n))%members=&
											&list_of_clusters(cluster_number(n))%members+1
										ENDIF
									ELSE
									! m is NOT a monomer!
										IF (cluster_number(n)==0) THEN
										! n is a monomer!
											!--> appropriate action: add entry "n" to cluster "m"
											cluster_number(n)=cluster_number(m)
											list_of_clusters(cluster_number(m))%cluster_weight=&
											&list_of_clusters(cluster_number(m))%cluster_weight+&
											&list_of_atoms(n)%atom_weight
											list_of_clusters(cluster_number(m))%members=&
											&list_of_clusters(cluster_number(m))%members+1
										ELSE
										! n is NOT a monomer!
											!--> appropriate action: merge clusters...
											!based on absolutely nothing, merge n into m
											cluster_todelete=cluster_number(n)
											list_of_clusters(cluster_number(m))%members=&
											&list_of_clusters(cluster_number(m))%members+&
											&list_of_clusters(cluster_todelete)%members
											list_of_clusters(cluster_number(m))%cluster_weight=&
											&list_of_clusters(cluster_number(m))%cluster_weight+&
											&list_of_clusters(cluster_todelete)%cluster_weight
											DO i=1,N_atoms,1
												IF (cluster_number(i)==cluster_todelete) THEN
													cluster_number(i)=cluster_number(m)
												ENDIF
											ENDDO
											!here, entry n is not needed anymore.
											!rename the cluster with the highest number to be n
											IF (current_cluster/=cluster_todelete) THEN
												list_of_clusters(cluster_todelete)%cluster_weight=&
												&list_of_clusters(current_cluster)%cluster_weight
												list_of_clusters(cluster_todelete)%members=&
												&list_of_clusters(current_cluster)%members
												DO i=1,N_atoms,1
													IF (cluster_number(i)==current_cluster) THEN
														cluster_number(i)=cluster_todelete
													ENDIF
												ENDDO
											ENDIF
											!delete the entry for the current_cluster. probably not necessary but safer to have clean memory.
											!list_of_clusters(current_cluster)%cluster_weight=0.0
											!list_of_clusters(current_cluster)%members=0
											current_cluster=current_cluster-1
										ENDIF
									ENDIF
								ENDIF
							ENDIF
						ENDDO
					ENDDO
					IF (print_statistics) THEN
						cluster_counts=current_cluster !we found this many non-monomeric clusters!
						DO i=1,current_cluster,1
							!$OMP CRITICAL(weight_updates)
							average_weights(list_of_clusters(i)%members,cutoff_counter)=&
							&average_weights(list_of_clusters(i)%members,cutoff_counter)+list_of_clusters(i)%cluster_weight
							!$OMP END CRITICAL(weight_updates)
							!$OMP CRITICAL(occurrence_updates)
							overall_occurrence(list_of_clusters(i)%members,cutoff_counter)=&
							&overall_occurrence(list_of_clusters(i)%members,cutoff_counter)+1
							!$OMP END CRITICAL(occurrence_updates)
						ENDDO
						DO i=1,N_atoms,1
							IF (cluster_number(i)==0) THEN
								cluster_counts=cluster_counts+1!every monomer counts as a cluster
								!$OMP CRITICAL(occurrence_updates)
								overall_occurrence(1,cutoff_counter)=overall_occurrence(1,cutoff_counter)+1
								!$OMP END CRITICAL(occurrence_updates)
								!$OMP CRITICAL(weight_updates)
								average_weights(1,cutoff_counter)=average_weights(1,cutoff_counter)+&
								&list_of_atoms(i)%atom_weight
								!$OMP END CRITICAL(weight_updates)
							ENDIF
						ENDDO
						!$OMP CRITICAL(CCDF_update)
						ClusterCountDistributionFunction(cutoff_counter)=&
						&ClusterCountDistributionFunction(cutoff_counter)+cluster_counts
						!$OMP END CRITICAL(CCDF_update)
					ENDIF
					!print clusters to xyz file if necessary
					IF (N_steps_to_print>0) THEN
						IF (ANY(steps_toprint(:)==stepcounter)) THEN
							!$OMP CRITICAL(print_members)
							steps_printed=.TRUE.
							!print this step
							!print spectators
							IF (print_spectators) THEN
								WRITE(fstring,'(2A,I0,A,E9.3,A)') &
								&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",&
								&stepcounter,"_cutoff=",cutoff_list(cutoff_counter),"_spectators.xyz"
								WRITE(outstring,'(A)') &
								&"These are the molecules which are not considered in the cluster analysis"
								CALL reset_print_members()
								DO i=1,N_atoms,1
									CALL add_print_member&
									&(list_of_atoms(i)%molecule_type_index,list_of_atoms(i)%molecule_index)
								ENDDO
								CALL invert_print_members()
								CALL print_members(stepcounter,TRIM(fstring),TRIM(outstring))
							ENDIF
							!print monomers
							WRITE(fstring,'(2A,I0,A,E9.3,A)') &
							&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",&
							&stepcounter,"_cutoff=",cutoff_list(cutoff_counter),"_monomers.xyz"
							WRITE(outstring,'(A)') &
							&"These are the monomers"
							CALL reset_print_members()
							DO i=1,N_atoms,1
								IF (cluster_number(i)==0) THEN
									CALL add_print_member&
									&(list_of_atoms(i)%molecule_type_index,list_of_atoms(i)%molecule_index)
								ENDIF
							ENDDO
							CALL print_members(stepcounter,TRIM(fstring),TRIM(outstring))
							!print oligomers
							DO i=2,N_atoms,1
								CALL reset_print_members()
								WRITE(fstring,'(2A,I0,A,E9.3,A,I0,A)') &
								&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",&
								&stepcounter,"_cutoff=",cutoff_list(cutoff_counter),"_",i,"-mers.xyz"
								WRITE(outstring,'(A,I0,A)') &
								&"These are the ",i,"-mers"
								DO m=1,current_cluster
									!go through all the clusters, find the ones which have "i" members
									IF (list_of_clusters(m)%members==i) THEN
										!here, cluster "m" has "i" members, now find all the atoms that belong to it
										DO n=1,N_atoms,1
											IF (cluster_number(n)==m) THEN
												CALL add_print_member&
												&(list_of_atoms(n)%molecule_type_index,list_of_atoms(n)%molecule_index)
											ENDIF
										ENDDO
										!here we have added all atoms belonging to that cluster to print_members
									ENDIF
								ENDDO
								!here we should have found all relevant clusters
								CALL print_members(stepcounter,TRIM(fstring),TRIM(outstring))
							ENDDO
							!$OMP END CRITICAL(print_members)
						ENDIF
					ENDIF
				ENDDO
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			!$OMP END DO
			IF (ALLOCATED(cluster_number)) THEN
				DEALLOCATE(cluster_number,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(list_of_clusters)) THEN
				DEALLOCATE(list_of_clusters,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			!$OMP END PARALLEL
			IF ((VERBOSE_OUTPUT).AND.(steps_printed)) WRITE(*,'(" various xyz files of clusters written to output folder.")')
			IF (print_statistics) THEN
				WRITE(fstring,'(2A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"cluster_statistics.dat"
				IF (VERBOSE_OUTPUT) WRITE(*,'(3A)')" Writing cluster statistics to '",TRIM(fstring),"'"
				OPEN(UNIT=4,FILE=TRIM(fstring))
				WRITE(4,*) "This file contains cluster statistics based on the input file '"&
					&,TRIM(FILENAME_CLUSTER_INPUT),"'"
				WRITE(4,*) "atom_fraction    = #(specified atoms in this cluster size) / total #(specified atoms)"
				WRITE(4,*) "cluster_fraction = #(clusters of this size) / total number of clusters"
				IF (N_cutoffs==1) THEN
					WRITE(4,'(" The cutoff was ",E9.3)') cutoff_list(1)
				ELSE
					WRITE(4,ADVANCE="NO",FMT='(A9)') "cutoff"
				ENDIF
				WRITE(4,'(A8,2A15,A17,A15)')"members","weight_average","atom_fraction","cluster_fraction"
				DO cutoff_counter=1,N_cutoffs,1
					normalisation=0.0
					DO i=1,N_atoms,1
						normalisation=normalisation+FLOAT(overall_occurrence(i,cutoff_counter)*i)
					ENDDO
					DO i=1,N_atoms,1
						IF (overall_occurrence(i,cutoff_counter)>0) THEN
							IF (N_cutoffs>1) WRITE(4,ADVANCE="NO",FMT='(E9.3)') cutoff_list(cutoff_counter)
							WRITE(4,'(I8,2E15.6,E17.6)')i,average_weights(i,cutoff_counter)/&
							&FLOAT(overall_occurrence(i,cutoff_counter)),&
							&(FLOAT(overall_occurrence(i,cutoff_counter)*i))/normalisation,&
							&FLOAT(overall_occurrence(i,cutoff_counter))/&
							&FLOAT(SUM(overall_occurrence(:,cutoff_counter)))
						ENDIF
					ENDDO
				ENDDO
				ENDFILE 4
				CLOSE(UNIT=4)
				WRITE(*,*) "The following values of the CCF were observed:"
				WRITE(*,'(2A10)') "cutoff","CCF"
				DO cutoff_counter=1,N_cutoffs,1
					IF (cutoff_list(cutoff_counter)<0.01) THEN
						WRITE(*,ADVANCE="NO",FMT='(E10.4)')cutoff_list(cutoff_counter)
					ELSEIF (cutoff_list(cutoff_counter)<99999.9) THEN
						WRITE(*,ADVANCE="NO",FMT='(F10.3)')cutoff_list(cutoff_counter)
					ELSE
						WRITE(*,ADVANCE="NO",FMT='(E10.4)')cutoff_list(cutoff_counter)
					ENDIF
					CCF=FLOAT(ClusterCountDistributionFunction(cutoff_counter))/&
					&FLOAT(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
					IF (CCF<0.01) THEN
						WRITE(*,'(E10.4)')CCF
					ELSEIF (CCF<99999.9) THEN
						WRITE(*,'(F10.3)')CCF
					ELSE
						WRITE(*,'(E10.4)')CCF
					ENDIF
				ENDDO
				WRITE(*,*) "( CCF = cluster count function, for comparison with TRAVIS:"
				WRITE(*,*) "J. Chem. Inf. Model., 2022, 62, 5634–5644, dx.doi.org/10.1021/acs.jcim.2c01244 )"
			ENDIF
		END SUBROUTINE cluster_analysis_GLOBAL

		!This is the version that looks only at pair connections. I duplicated the whole routine to not accidentally mess up the GLOBAL analysis.
		SUBROUTINE cluster_analysis_PAIRS()
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: stepcounter,allocstatus,deallocstatus,m,n,current_cluster,i,ios,cluster_todelete,cluster_counts
		INTEGER :: pair_counter
		LOGICAL :: steps_printed
		REAL :: normalisation,CCF
		CHARACTER(LEN=1024) :: fstring,outstring
		TYPE :: one_cluster
			REAL :: cluster_weight
			INTEGER :: members
		END TYPE one_cluster
		TYPE(one_cluster),DIMENSION(:),ALLOCATABLE :: list_of_clusters
		INTEGER,DIMENSION(:),ALLOCATABLE :: cluster_number
			steps_printed=.FALSE.
			!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL)))&
			!$OMP PRIVATE(list_of_clusters,cluster_number,m,n,i,current_cluster,pair_counter,cluster_todelete,cluster_counts)
			!$OMP SINGLE
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (GLOBAL cluster analysis)"
		 !$ 	CALL timing_parallel_sections(.TRUE.)
		 !$ ENDIF
			IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
			!$OMP END SINGLE
			ALLOCATE(list_of_clusters(N_atoms),STAT=allocstatus)
			IF (allocstatus/=0) THEN
				CALL report_error(171,exit_status=allocstatus)
			ENDIF
			ALLOCATE(cluster_number(N_atoms),STAT=allocstatus)
			IF (allocstatus/=0) THEN
				CALL report_error(171,exit_status=allocstatus)
			ENDIF
			!$OMP DO SCHEDULE(STATIC,1) ORDERED
			DO stepcounter=1,nsteps,sampling_interval
				current_cluster=0
				cluster_number(:)=0
				list_of_clusters(:)%cluster_weight=0.0
				list_of_clusters(:)%members=0
				!Iterate over pairs
				DO pair_counter=1,N_pairs,1
					DO m=list_of_pairs(pair_counter)%list_of_atoms_starting_index_A,&
					&list_of_pairs(pair_counter)%list_of_atoms_last_index_A,1
						DO n=list_of_pairs(pair_counter)%list_of_atoms_starting_index_B,&
						&list_of_pairs(pair_counter)%list_of_atoms_last_index_B,1
							!do not need to check pairs that belong together already
							IF (((cluster_number(m)==0).OR.(cluster_number(n)==0))&
							&.OR.(cluster_number(m)/=cluster_number(n))) THEN
								!check if closer than cutoff
								IF (give_smallest_atom_distance_squared&
								&(stepcounter,stepcounter,&
								&list_of_atoms(m)%molecule_type_index,list_of_atoms(n)%molecule_type_index,&
								&list_of_atoms(m)%molecule_index,list_of_atoms(n)%molecule_index,&
								&list_of_atoms(m)%atom_index,list_of_atoms(n)%atom_index)&
								&<squared_cutoff_list(pair_counter)) THEN
									IF (cluster_number(m)==0) THEN
									! m is a monomer!
										IF (cluster_number(n)==0) THEN
										! n is a monomer!
											!--> appropriate action: assign new cluster.
											current_cluster=current_cluster+1
											cluster_number(m)=current_cluster
											cluster_number(n)=current_cluster
											list_of_clusters(current_cluster)%cluster_weight=&
											&list_of_atoms(m)%atom_weight+&
											&list_of_atoms(n)%atom_weight
											list_of_clusters(current_cluster)%members=2
										ELSE
										! n is NOT a monomer!
											!--> appropriate action: add entry "m" to cluster "n"
											cluster_number(m)=cluster_number(n)
											list_of_clusters(cluster_number(n))%cluster_weight=&
											&list_of_clusters(cluster_number(n))%cluster_weight+&
											&list_of_atoms(m)%atom_weight
											list_of_clusters(cluster_number(n))%members=&
											&list_of_clusters(cluster_number(n))%members+1
										ENDIF
									ELSE
									! m is NOT a monomer!
										IF (cluster_number(n)==0) THEN
										! n is a monomer!
											!--> appropriate action: add entry "n" to cluster "m"
											cluster_number(n)=cluster_number(m)
											list_of_clusters(cluster_number(m))%cluster_weight=&
											&list_of_clusters(cluster_number(m))%cluster_weight+&
											&list_of_atoms(n)%atom_weight
											list_of_clusters(cluster_number(m))%members=&
											&list_of_clusters(cluster_number(m))%members+1
										ELSE
										! n is NOT a monomer!
											!--> appropriate action: merge clusters...
											!based on absolutely nothing, merge n into m
											cluster_todelete=cluster_number(n)
											list_of_clusters(cluster_number(m))%members=&
											&list_of_clusters(cluster_number(m))%members+&
											&list_of_clusters(cluster_todelete)%members
											list_of_clusters(cluster_number(m))%cluster_weight=&
											&list_of_clusters(cluster_number(m))%cluster_weight+&
											&list_of_clusters(cluster_todelete)%cluster_weight
											DO i=1,N_atoms,1
												IF (cluster_number(i)==cluster_todelete) THEN
													cluster_number(i)=cluster_number(m)
												ENDIF
											ENDDO
											!here, entry n is not needed anymore.
											!rename the cluster with the highest number to be n
											IF (current_cluster/=cluster_todelete) THEN
												list_of_clusters(cluster_todelete)%cluster_weight=&
												&list_of_clusters(current_cluster)%cluster_weight
												list_of_clusters(cluster_todelete)%members=&
												&list_of_clusters(current_cluster)%members
												DO i=1,N_atoms,1
													IF (cluster_number(i)==current_cluster) THEN
														cluster_number(i)=cluster_todelete
													ENDIF
												ENDDO
											ENDIF
											!delete the entry for the current_cluster. probably not necessary but safer to have clean memory.
											!list_of_clusters(current_cluster)%cluster_weight=0.0
											!list_of_clusters(current_cluster)%members=0
											current_cluster=current_cluster-1
										ENDIF
									ENDIF
								ENDIF
							ENDIF
						ENDDO
					ENDDO
				ENDDO !end of pair counter loop
				IF (print_statistics) THEN
					cluster_counts=current_cluster !we found this many non-monomeric clusters!
					DO i=1,current_cluster,1
						!$OMP CRITICAL(weight_updates)
						average_weights(list_of_clusters(i)%members,1)=&
						&average_weights(list_of_clusters(i)%members,1)+list_of_clusters(i)%cluster_weight
						!$OMP END CRITICAL(weight_updates)
						!$OMP CRITICAL(occurrence_updates)
						overall_occurrence(list_of_clusters(i)%members,1)=&
						&overall_occurrence(list_of_clusters(i)%members,1)+1
						!$OMP END CRITICAL(occurrence_updates)
					ENDDO
					DO i=1,N_atoms,1
						IF (cluster_number(i)==0) THEN
							cluster_counts=cluster_counts+1!every monomer counts as a cluster
							!$OMP CRITICAL(occurrence_updates)
							overall_occurrence(1,1)=overall_occurrence(1,1)+1
							!$OMP END CRITICAL(occurrence_updates)
							!$OMP CRITICAL(weight_updates)
							average_weights(1,1)=average_weights(1,1)+&
							&list_of_atoms(i)%atom_weight
							!$OMP END CRITICAL(weight_updates)
						ENDIF
					ENDDO
					!$OMP CRITICAL(CCDF_update)
					ClusterCountDistributionFunction(1)=&
					&ClusterCountDistributionFunction(1)+cluster_counts
					!$OMP END CRITICAL(CCDF_update)
				ENDIF
				!print clusters to xyz file if necessary
				IF (N_steps_to_print>0) THEN
					IF (ANY(steps_toprint(:)==stepcounter)) THEN
						!$OMP CRITICAL(print_members)
						steps_printed=.TRUE.
						!print this step
						!print spectators
						IF (print_spectators) THEN
							WRITE(fstring,'(2A,I0,A,A)') &
							&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",&
							&stepcounter,"_spectators.xyz"
							WRITE(outstring,'(A)') &
							&"These are the molecules which are not considered in the cluster analysis"
							CALL reset_print_members()
							DO i=1,N_atoms,1
								CALL add_print_member&
								&(list_of_atoms(i)%molecule_type_index,list_of_atoms(i)%molecule_index)
							ENDDO
							CALL invert_print_members()
							CALL print_members(stepcounter,TRIM(fstring),TRIM(outstring))
						ENDIF
						!print monomers
						WRITE(fstring,'(2A,I0,A)') &
						&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",&
						&stepcounter,"_monomers.xyz"
						WRITE(outstring,'(A)') &
						&"These are the monomers"
						CALL reset_print_members()
						DO i=1,N_atoms,1
							IF (cluster_number(i)==0) THEN
								CALL add_print_member&
								&(list_of_atoms(i)%molecule_type_index,list_of_atoms(i)%molecule_index)
							ENDIF
						ENDDO
						CALL print_members(stepcounter,TRIM(fstring),TRIM(outstring))
						!print oligomers
						DO i=2,N_atoms,1
							CALL reset_print_members()
							WRITE(fstring,'(2A,I0,A,I0,A)') &
							&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",&
							&stepcounter,"_",i,"-mers.xyz"
							WRITE(outstring,'(A,I0,A)') &
							&"These are the ",i,"-mers"
							DO m=1,current_cluster
								!go through all the clusters, find the ones which have "i" members
								IF (list_of_clusters(m)%members==i) THEN
									!here, cluster "m" has "i" members, now find all the atoms that belong to it
									DO n=1,N_atoms,1
										IF (cluster_number(n)==m) THEN
											CALL add_print_member&
											&(list_of_atoms(n)%molecule_type_index,list_of_atoms(n)%molecule_index)
										ENDIF
									ENDDO
									!here we have added all atoms belonging to that cluster to print_members
								ENDIF
							ENDDO
							!here we should have found all relevant clusters
							CALL print_members(stepcounter,TRIM(fstring),TRIM(outstring))
						ENDDO
						!$OMP END CRITICAL(print_members)
					ENDIF
				ENDIF
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			!$OMP END DO
			IF (ALLOCATED(cluster_number)) THEN
				DEALLOCATE(cluster_number,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(list_of_clusters)) THEN
				DEALLOCATE(list_of_clusters,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			!$OMP END PARALLEL
			IF ((VERBOSE_OUTPUT).AND.(steps_printed)) WRITE(*,'(" various xyz files of clusters written to output folder.")')
			IF (print_statistics) THEN
				WRITE(fstring,'(2A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"cluster_statistics.dat"
				IF (VERBOSE_OUTPUT) WRITE(*,'(3A)')" Writing cluster statistics to '",TRIM(fstring),"'"
				OPEN(UNIT=4,FILE=TRIM(fstring))
				WRITE(4,*) "This file contains cluster statistics based on the input file '"&
					&,TRIM(FILENAME_CLUSTER_INPUT),"'"
				WRITE(4,*) "atom_fraction    = #(specified atoms in this cluster size) / total #(specified atoms)"
				WRITE(4,*) "cluster_fraction = #(clusters of this size) / total number of clusters"
				WRITE(4,'(A8,2A15,A17,A15)')"members","weight_average","atom_fraction","cluster_fraction"
				normalisation=0.0
				DO i=1,N_atoms,1
					normalisation=normalisation+FLOAT(overall_occurrence(i,1)*i)
				ENDDO
				DO i=1,N_atoms,1
					IF (overall_occurrence(i,1)>0) THEN
						WRITE(4,'(I8,2E15.6,E17.6)')i,average_weights(i,1)/&
						&FLOAT(overall_occurrence(i,1)),&
						&(FLOAT(overall_occurrence(i,1)*i))/normalisation,&
						&FLOAT(overall_occurrence(i,1))/&
						&FLOAT(SUM(overall_occurrence(:,1)))
					ENDIF
				ENDDO
				ENDFILE 4
				CLOSE(UNIT=4)
				WRITE(*,ADVANCE="NO",FMT='(A)') " The 'CCF' for this analysis was"
				CCF=FLOAT(ClusterCountDistributionFunction(1))/&
				&FLOAT(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
				IF ((CCF<0.01).OR.(CCF>99.99)) THEN
					WRITE(*,'(E11.4)')CCF
				ELSE
					WRITE(*,'(F6.2)')CCF
				ENDIF
				WRITE(*,*) "( CCF = cluster count function = average number of separate clusters per box)"
			ENDIF
		END SUBROUTINE cluster_analysis_PAIRS

		SUBROUTINE perform_cluster_analysis()
		IMPLICIT NONE
			CALL initialise_cluster()
			IF ((ERROR_CODE/=170).AND.(ERROR_CODE/=171).AND.(ERROR_CODE/=173).AND.(ERROR_CODE/=174).AND.(ERROR_CODE/=177)) THEN
				WRITE(*,'(" Using ",I0," timesteps in intervals of ",I0," for averaging.")')&
				&MAX((nsteps-1+sampling_interval)/sampling_interval,0),sampling_interval
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,ADVANCE="NO",FMT='(" The main data structures occupy approximately ")')
					CALL print_memory_requirement(FLOAT(bytes_needed)*(4.0/1024.0d0))
					WRITE(*,'(".")')
				ENDIF
				WRITE(*,'(" Starting cluster analysis.")')
				CALL refresh_IO()
				IF (TRIM(operation_mode)=="GLOBAL") THEN
					CALL cluster_analysis_GLOBAL()
				ELSE
					CALL cluster_analysis_PAIRS()
				ENDIF
				CALL finalise_cluster()
			ELSE
				ERROR_CODE=ERROR_CODE_DEFAULT
			ENDIF
		END SUBROUTINE perform_cluster_analysis

END MODULE CLUSTER
!--------------------------------------------------------------------------------------------------------------------------------!