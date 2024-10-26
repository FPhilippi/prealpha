
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
	LOGICAL,PARAMETER :: print_CDistanceDF_default=.FALSE. !cluster distance distribution function
	LOGICAL,PARAMETER :: print_CCountF_default=.FALSE. !cluster count function
	LOGICAL,PARAMETER :: print_CSizeDF_default=.FALSE. !cluster size distribution function
	LOGICAL,PARAMETER :: print_CChargeDF_default=.FALSE. !cluster charge distribution function
	LOGICAL,PARAMETER :: complete_molecules_default=.TRUE.
	!variables
	CHARACTER (LEN=10) :: operation_mode="NONE"!operation mode of the cluster module.
	!operation_mode="PAIRS": needs pairs of atom indices and their fixed cutoffs
	!operation_mode="GLOBAL": needs a global cutoff operating on a set of atoms. The latter are specified by molecule_type_index and atom_index, and the code iterates over molecule_index.
	INTEGER :: nsteps=nsteps_default !how many steps to use from trajectory
	INTEGER :: sampling_interval=sampling_interval_default
	INTEGER :: N_cutoffs,N_atoms,N_steps_to_print
	INTEGER :: bytes_needed
	LOGICAL :: print_statistics=print_statistics_default
	LOGICAL :: print_spectators=print_spectators_default !if T, print also the molecules which are not members of list_of_atoms
	LOGICAL :: complete_molecules=complete_molecules_default !if T, print whole molecules. If T, only print the atoms from the list.
	LOGICAL :: print_CDistanceDF=print_CDistanceDF_default
	LOGICAL :: print_CCountF=print_CCountF_default
	LOGICAL :: print_CSizeDF=print_CSizeDF_default
	LOGICAL :: print_CChargeDF=print_CChargeDF_default
	LOGICAL :: custom_weights=.FALSE.
	TYPE :: cluster_atom
		INTEGER :: molecule_type_index
		INTEGER :: molecule_index
		INTEGER :: atom_index
		INTEGER :: cluster_number
		REAL :: atom_weight !usually, this would be the charge.
	END TYPE cluster_atom
	TYPE :: one_cluster
		REAL :: cluster_weight
		INTEGER :: members
	END TYPE one_cluster
	TYPE(cluster_atom),DIMENSION(:),ALLOCATABLE :: list_of_atoms
	TYPE(one_cluster),DIMENSION(:),ALLOCATABLE :: list_of_clusters
	INTEGER(KIND=DP),DIMENSION(:,:),ALLOCATABLE :: overall_occurrence
	REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE :: average_weights
	REAL,DIMENSION(:),ALLOCATABLE :: cutoff_list
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
			print_CDistanceDF=print_CDistanceDF_default
			print_CCountF=print_CCountF_default
			print_CSizeDF=print_CSizeDF_default
			print_CChargeDF=print_CChargeDF_default
			print_statistics=print_statistics_default
			custom_weights=.FALSE.
			bytes_needed=0
		END SUBROUTINE set_defaults

		!initialises the cluster module by reading the specified input file.
		SUBROUTINE initialise_cluster()
		IMPLICIT NONE
		INTEGER :: i,allocstatus
		LOGICAL :: file_exists,connected
		INTEGER :: ios
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
				!N_atoms, N_cutoffs, N_steps_to_print are initialised
				IF (N_atoms==0) THEN
					CALL report_error(173)
				ELSEIF (N_cutoffs==0) THEN
					CALL report_error(174)
				ELSE
					WRITE(*,'("   Expecting ",I0," atoms and ",I0," cutoffs.")') N_atoms,N_cutoffs
					!allocate memory which is needed now - that's why I do the firstpass!
					ALLOCATE(cutoff_list(N_cutoffs),STAT=allocstatus)
					bytes_needed=bytes_needed+N_cutoffs !multiply with 4 bytes later
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=ios)
						RETURN
					ENDIF
					IF (N_steps_to_print>0) THEN
						ALLOCATE(steps_toprint(N_steps_to_print),STAT=allocstatus)
						bytes_needed=bytes_needed+N_steps_to_print !multiply with 4 bytes later
						IF (allocstatus/=0) THEN
							CALL report_error(171,exit_status=ios)
							RETURN
						ENDIF
					ENDIF
					ALLOCATE(list_of_atoms(N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=ios)
						RETURN
					ENDIF
					CALL read_body()
					ALLOCATE(squared_cutoff_list(N_cutoffs),STAT=allocstatus)
					bytes_needed=bytes_needed+N_cutoffs !multiply with 4 bytes later
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=ios)
						RETURN
					ENDIF
					bytes_needed=bytes_needed+5*N_atoms !multiply with 4 bytes later
					!list_of_clusters cannot be larger than N_atoms (=all monomers...)
					ALLOCATE(list_of_clusters(N_atoms),STAT=allocstatus)
					IF (allocstatus/=0) THEN
						CALL report_error(171,exit_status=ios)
						RETURN
					ENDIF
					bytes_needed=bytes_needed+2*N_atoms !multiply with 4 bytes later
					IF (print_statistics) THEN
						ALLOCATE(average_weights(N_atoms,N_cutoffs),STAT=allocstatus)
						IF (allocstatus/=0) THEN
							CALL report_error(171,exit_status=ios)
							RETURN
						ENDIF
						ALLOCATE(overall_occurrence(N_atoms,N_cutoffs),STAT=allocstatus)
						IF (allocstatus/=0) THEN
							CALL report_error(171,exit_status=ios)
							RETURN
						ENDIF
						average_weights(:,:)=0.0
						overall_occurrence(:,:)=0
						bytes_needed=bytes_needed+2*2*N_atoms*N_cutoffs !multiply with 4 bytes later
					ENDIF
					CALL sort_cutoff_list(1,N_cutoffs,cutoff_list)
					CALL sort_steps_to_print_list(1,N_steps_to_print,steps_toprint)
					DO i=1,N_cutoffs,1
						squared_cutoff_list(i)=cutoff_list(i)**2
					ENDDO
				ENDIF
			ELSE
				CALL report_error(170)!No input - no output. easy as that.
			ENDIF

		CONTAINS

				SUBROUTINE read_first_pass()
				IMPLICIT NONE
				CHARACTER(LEN=33) :: inputstring
				INTEGER :: molecule_type_index,atom_index,scansteps,inputinteger,inputinteger2,n
				REAL :: dummy,cutlow,cuthigh
					REWIND 3
					molecule_type_index=0
					N_atoms=0
					N_cutoffs=0
					N_steps_to_print=0
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF (ios/=0) EXIT
						SELECT CASE (TRIM(inputstring))
						CASE ("scan_cutoff")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,cutlow,cuthigh,scansteps
							IF (ios==0) THEN
								IF ((scansteps>0).OR.(cutlow<0.0).OR.(cutlow>cuthigh))&
								& N_cutoffs=N_cutoffs+scansteps
							ENDIF
						CASE ("single_cutoff","cutoff")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,dummy
							IF (ios==0) THEN
								IF (dummy>0.0d0) N_cutoffs=N_cutoffs+1
							ENDIF
						CASE ("N_steps_to_print","N_printsteps","n_steps_to_print","n_printsteps")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
							IF (ios==0) THEN
								N_steps_to_print=N_steps_to_print+1
							ENDIF
						CASE ("add_molecule_type","add_molecule_type_index","add_molecule")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
							IF (ios==0) THEN
								IF (valid_molecule_type_index(inputinteger)) THEN
									N_atoms=N_atoms+&
									&give_number_of_atoms_per_molecule(inputinteger)*&
									&give_number_of_molecules_per_step(inputinteger)
								ENDIF
							ENDIF
						CASE ("add_atom_index","add_atom")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger,inputinteger2
							IF (ios==0) THEN
								IF (valid_atom_index(inputinteger,inputinteger2)) THEN
									N_atoms=N_atoms+give_number_of_molecules_per_step(inputinteger)
								ENDIF
							ENDIF
						CASE ("quit")
							EXIT
						END SELECT
					ENDDO
				END SUBROUTINE read_first_pass

				SUBROUTINE read_body()
				IMPLICIT NONE
				CHARACTER(LEN=33) :: inputstring
				INTEGER :: atom_index,current_cutoff,current_printstep,current_atom,molecule_index
				INTEGER :: scansteps,n,inputinteger,inputinteger2,m
				REAL :: cutlow,cuthigh,interval
					REWIND 3
					current_cutoff=0
					current_printstep=0
					current_atom=0
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "  End-of-file condition in ",TRIM(FILENAME_CLUSTER_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_CLUSTER_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("scan_cutoff")
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
						CASE ("single_cutoff","cutoff")
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
							IF (VERBOSE_OUTPUT) WRITE(*,'(A,F0.2,A)')&
							&"   single cutoff (",cutlow,") added."
						CASE ("N_steps_to_print","N_printsteps","n_steps_to_print","n_printsteps")
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
						CASE ("add_molecule_type","add_molecule_type_index","add_molecule")
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
											list_of_atoms(current_atom)%cluster_number=0
											list_of_atoms(current_atom)%atom_weight=0.0
										ENDDO
									ENDDO
								ENDIF
							ENDIF
						CASE ("add_atom_index","add_atom")
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
										list_of_atoms(current_atom)%cluster_number=0
										list_of_atoms(current_atom)%atom_weight=0.0
									ENDDO
								ENDIF
							ENDIF
						CASE ("complete_molecules")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,complete_molecules
							IF (ios/=0) THEN
								CALL report_error(170,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)')&
								&"   setting 'complete_molecules' to default (=",complete_molecules_default,")"
								complete_molecules=complete_molecules_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1)') "   setting 'complete_molecules' to ",&
								&complete_molecules
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
			INTEGER pivot
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
			IF (ALLOCATED(list_of_atoms)) THEN
				DEALLOCATE(list_of_atoms,STAT=deallocstatus)
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
				IF (ALLOCATED(steps_toprint)) THEN
					DEALLOCATE(steps_toprint,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				ELSE
					CALL report_error(0)
				ENDIF
			ENDIF
			IF (print_statistics) THEN
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
		INTEGER :: stepcounter,allocstatus,deallocstatus,m,n,current_cluster,i,ios,cluster_todelete,output_atom_count
		INTEGER :: cutoff_counter
		LOGICAL :: first_atom
		REAL :: normalisation,dummy
			IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
			DO stepcounter=1,nsteps,sampling_interval
				current_cluster=0
				list_of_atoms(:)%cluster_number=0
				list_of_clusters(:)%cluster_weight=0.0
				list_of_clusters(:)%members=0
				!Iterate over cutoffs...
				DO cutoff_counter=1,N_cutoffs,1
					DO m=1,N_atoms
						DO n=m,N_atoms
							IF (m==n) CYCLE
							!do not need to check pairs that belong together already
							IF (((list_of_atoms(m)%cluster_number==0).OR.(list_of_atoms(n)%cluster_number==0))&
							&.OR.(list_of_atoms(m)%cluster_number/=list_of_atoms(n)%cluster_number)) THEN
								!check if closer than cutoff
								IF (give_smallest_atom_distance_squared&
								&(stepcounter,stepcounter,&
								&list_of_atoms(m)%molecule_type_index,list_of_atoms(n)%molecule_type_index,&
								&list_of_atoms(m)%molecule_index,list_of_atoms(n)%molecule_index,&
								&list_of_atoms(m)%atom_index,list_of_atoms(n)%atom_index)&
								&<squared_cutoff_list(cutoff_counter)) THEN
									IF (list_of_atoms(m)%cluster_number==0) THEN
									! m is a monomer!
										IF (list_of_atoms(n)%cluster_number==0) THEN
										! n is a monomer!
											!--> appropriate action: assign new cluster.
											current_cluster=current_cluster+1
											list_of_atoms(m)%cluster_number=current_cluster
											list_of_atoms(n)%cluster_number=current_cluster
											list_of_clusters(current_cluster)%cluster_weight=&
											&list_of_atoms(m)%atom_weight+&
											&list_of_atoms(n)%atom_weight
											list_of_clusters(current_cluster)%members=2
										ELSE
										! n is NOT a monomer!
											!--> appropriate action: add entry "m" to cluster "n"
											list_of_atoms(m)%cluster_number=list_of_atoms(n)%cluster_number
											list_of_clusters(list_of_atoms(n)%cluster_number)%cluster_weight=&
											&list_of_clusters(list_of_atoms(n)%cluster_number)%cluster_weight+&
											&give_charge_of_molecule(list_of_atoms(m)%molecule_type_index)
											list_of_clusters(list_of_atoms(n)%cluster_number)%members=&
											&list_of_clusters(list_of_atoms(n)%cluster_number)%members+1
										ENDIF
									ELSE
									! m is NOT a monomer!
										IF (list_of_atoms(n)%cluster_number==0) THEN
										! n is a monomer!
											!--> appropriate action: add entry "n" to cluster "m"
											list_of_atoms(n)%cluster_number=list_of_atoms(m)%cluster_number
											list_of_clusters(list_of_atoms(m)%cluster_number)%cluster_weight=&
											&list_of_clusters(list_of_atoms(m)%cluster_number)%cluster_weight+&
											&list_of_atoms(n)%atom_weight
											list_of_clusters(list_of_atoms(m)%cluster_number)%members=&
											&list_of_clusters(list_of_atoms(m)%cluster_number)%members+1
										ELSE
										! n is NOT a monomer!
											!--> appropriate action: merge clusters...
											!based on absolutely nothing, merge n into m
											cluster_todelete=list_of_atoms(n)%cluster_number
											list_of_clusters(list_of_atoms(m)%cluster_number)%members=&
											&list_of_clusters(list_of_atoms(m)%cluster_number)%members+&
											&list_of_clusters(cluster_todelete)%members
											list_of_clusters(list_of_atoms(m)%cluster_number)%cluster_weight=&
											&list_of_clusters(list_of_atoms(m)%cluster_number)%cluster_weight+&
											&list_of_clusters(cluster_todelete)%cluster_weight
											DO i=1,N_atoms,1
												IF (list_of_atoms(i)%cluster_number==cluster_todelete) THEN
													list_of_atoms(i)%cluster_number=list_of_atoms(m)%cluster_number
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
													IF (list_of_atoms(i)%cluster_number==current_cluster) THEN
														list_of_atoms(i)%cluster_number=cluster_todelete
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
						DO i=1,current_cluster,1
							average_weights(list_of_clusters(i)%members,cutoff_counter)=&
							&average_weights(list_of_clusters(i)%members,cutoff_counter)+list_of_clusters(i)%cluster_weight
							overall_occurrence(list_of_clusters(i)%members,cutoff_counter)=&
							&overall_occurrence(list_of_clusters(i)%members,cutoff_counter)+1
						ENDDO
						DO i=1,N_atoms,1
							IF (list_of_atoms(i)%cluster_number==0) THEN
								overall_occurrence(1,cutoff_counter)=overall_occurrence(1,cutoff_counter)+1
								average_weights(1,cutoff_counter)=average_weights(1,cutoff_counter)+&
								&list_of_atoms(i)%atom_weight
							ENDIF
						ENDDO
					ENDIF
				ENDDO
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			IF (print_statistics) THEN
				DO cutoff_counter=1,N_cutoffs,1
					PRINT *,"CUTOFF #",cutoff_counter
					normalisation=0.0
					DO i=1,N_atoms,1
						WRITE(*,'(" Members= ",I0,", charge=",I0,", occurrence=",I0)')&
						&i,average_weights(i,cutoff_counter),overall_occurrence(i,cutoff_counter)
						normalisation=normalisation+FLOAT(overall_occurrence(i,cutoff_counter)*i)
					ENDDO
					DO i=1,N_atoms,1
						WRITE(*,*)"Members ","average_weights ","fraction_of_molecules ","fraction_of_clusters "
						WRITE(*,*)i,average_weights(i,cutoff_counter)/&
						&FLOAT(overall_occurrence(i,cutoff_counter)),&
						&(FLOAT(overall_occurrence(i,cutoff_counter)*i))/normalisation,&
						&FLOAT(overall_occurrence(i,cutoff_counter))/&
						&FLOAT(SUM(overall_occurrence(:,cutoff_counter)))
					ENDDO
				ENDDO
			ENDIF
		END SUBROUTINE cluster_analysis_GLOBAL

		SUBROUTINE print_clusters(logged_clusters)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: logged_clusters
		REAL(KIND=DP) :: current_centre(3)
		REAL :: dummy
		LOGICAL :: first_atom
		INTEGER :: first_entry,ip,np,mp,output_atom_count,stepcounter
		CHARACTER(LEN=1024) :: fstring
				!print the glymes
				WRITE(fstring,'(2A,I0,A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",stepcounter,"_glymes.xyz"
				OPEN(UNIT=4,FILE=TRIM(fstring))
				output_atom_count=3360
				CALL write_header(4,1,output_atom_count,"xyz")
				DO ip=1,210,1
					CALL write_molecule(4,stepcounter,2,ip,include_header=.FALSE.)
				ENDDO
				ENDFILE 4
				CLOSE(UNIT=4)
						
						
				!print the clusters
				WRITE(fstring,'(2A,I0,A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",stepcounter,"_monomers.xyz"
				OPEN(UNIT=4,FILE=TRIM(fstring))
				output_atom_count=0
				DO ip=1,N_atoms,1
					IF (list_of_atoms(ip)%cluster_number==0) THEN
						output_atom_count=output_atom_count+&
						&give_number_of_atoms_per_molecule(list_of_atoms(ip)%molecule_type_index)
					ENDIF
				ENDDO
				CALL write_header(4,1,output_atom_count,"xyz")
				DO ip=1,N_atoms,1
					IF (list_of_atoms(ip)%cluster_number==0) THEN
						CALL write_molecule(4,stepcounter,list_of_atoms(ip)%molecule_type_index,&
						&list_of_atoms(ip)%molecule_index,include_header=.FALSE.)
					ENDIF
				ENDDO
				ENDFILE 4
				CLOSE(UNIT=4)
				DO ip=2,N_atoms,1
					output_atom_count=0
					DO mp=1,logged_clusters
						!iterate over the specified timesteps
						!check if cluster has the right number of members
						IF (list_of_clusters(mp)%members==ip) THEN
							!count atoms
							DO np=1,N_atoms,1
								IF (list_of_atoms(np)%cluster_number==mp) THEN
									output_atom_count=output_atom_count+&
									&give_number_of_atoms_per_molecule(list_of_atoms(np)%molecule_type_index)
								ENDIF
							ENDDO
						ENDIF
					ENDDO
					IF (output_atom_count/=0) THEN
						WRITE(fstring,'(2A,I0,A,I0,A)') &
						&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"step_",stepcounter,"_clusters_",ip,"-mers.xyz"
						OPEN(UNIT=4,FILE=TRIM(fstring))
						CALL write_header(4,1,output_atom_count,"xyz")
						DO mp=1,logged_clusters
							!iterate over the specified timesteps
							!check if cluster has the right number of members
							IF (list_of_clusters(mp)%members==ip) THEN
								!write body
								first_atom=.TRUE.
								DO np=1,N_atoms,1
									IF (list_of_atoms(np)%cluster_number==mp) THEN
										IF (first_atom) THEN
											current_centre(:)=0.0d0
											first_atom=.FALSE.
											first_entry=np
										ELSE
											dummy=give_smallest_atom_distance_squared&
											&(stepcounter,stepcounter,&
											&list_of_atoms(first_entry)%molecule_type_index,list_of_atoms(np)%molecule_type_index,&
											&list_of_atoms(first_entry)%molecule_index,list_of_atoms(np)%molecule_index,&
											&list_of_atoms(first_entry)%atom_index,list_of_atoms(np)%atom_index,translation=current_centre)
										ENDIF
										CALL write_molecule(4,stepcounter,list_of_atoms(np)%molecule_type_index,&
										&list_of_atoms(np)%molecule_index&
										&,include_header=.FALSE.,translate_by=current_centre(:))
									ENDIF
								ENDDO
							ENDIF
						ENDDO
						ENDFILE 4
						CLOSE(UNIT=4)
					ENDIF
				ENDDO
		END SUBROUTINE print_clusters

		SUBROUTINE perform_cluster_analysis()
		IMPLICIT NONE
			CALL initialise_cluster()
			IF ((ERROR_CODE/=170).AND.(ERROR_CODE/=171).AND.(ERROR_CODE/=173).AND.(ERROR_CODE/=174)) THEN
				WRITE(*,'(" Using ",I0," timesteps in intervals of ",I0," for averaging.")')&
				&MAX((nsteps-1+sampling_interval)/sampling_interval,0),sampling_interval
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,ADVANCE="NO",FMT='(" The main data structures occupy approximately ")')
					CALL print_memory_requirement(FLOAT(bytes_needed)*(4.0/1024.0d0))
					WRITE(*,'(".")')
				ENDIF
				WRITE(*,'(" Starting cluster analysis.")')
				CALL refresh_IO()
				CALL cluster_analysis_GLOBAL()
				CALL finalise_cluster()
			ELSE
				ERROR_CODE=ERROR_CODE_DEFAULT
			ENDIF
		END SUBROUTINE perform_cluster_analysis

END MODULE CLUSTER
!--------------------------------------------------------------------------------------------------------------------------------!