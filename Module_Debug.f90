
!This module contains procedures for debugging and technical purposes.
MODULE DEBUG ! Copyright (C) !RELEASEYEAR! Frederik Philippi
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
	PUBLIC dump_neighbour_traj,remove_cores,dump_atomic_properties
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
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
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
					!$ 	WRITE(*,'(A,I0,A)') "   ### Parallel execution on ",OMP_get_num_threads()," threads (jump histogram)"
					!$ 	CALL timing_parallel_sections(.TRUE.)
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
				 !$ 	WRITE(*,ADVANCE="NO",FMT='("   ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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
			IF (VERBOSE_OUTPUT) CALL print_progress(startstep_in-endstep_in)
			CALL initialise_neighbourtraj()
			DO stepcounter=startstep_in,endstep_in,1
				IF (update_com) current_centre(:)=give_center_of_mass(stepcounter,molecule_type_index_1,molecule_index_1)
				CALL make_neighbour_list()
				CALL write_neighbours()
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			IF (((startstep_in-endstep_in)>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
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

		!This SUBROUTINE writes atomic charges and masses to two separate files.
		SUBROUTINE dump_atomic_properties()
		IMPLICIT NONE
		INTEGER :: molecule_type_index,molecule_index,atom_index
		CHARACTER(LEN=1024) :: fstring
		LOGICAL :: connected
			!First, do the fools-proof check
			INQUIRE(UNIT=4,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=4)
			WRITE(*,'(" Writing atomic properties to files:")')
			WRITE(fstring,'(2A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"atomic_charges.dat"
			WRITE(*,*) TRIM(fstring)
			OPEN(UNIT=4,FILE=TRIM(fstring))
			DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
				DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
					DO atom_index=1,give_number_of_atoms_per_molecule(molecule_type_index),1
						WRITE(4,*) give_charge_of_atom(molecule_type_index,atom_index)
					ENDDO
				ENDDO
			ENDDO
			CLOSE(UNIT=4)
			WRITE(fstring,'(2A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"atomic_masses.dat"
			WRITE(*,*) TRIM(fstring)
			OPEN(UNIT=4,FILE=TRIM(fstring))
			DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
				DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
					DO atom_index=1,give_number_of_atoms_per_molecule(molecule_type_index),1
						WRITE(4,*) give_mass_of_atom(molecule_type_index,atom_index)
					ENDDO
				ENDDO
			ENDDO
			CLOSE(UNIT=4)
		END SUBROUTINE dump_atomic_properties

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
			CHARACTER (LEN=2) :: atom_type='X'
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
	 !$ 	FUNCTION OMP_get_wtime()
	 !$ 	REAL(8) :: OMP_get_wtime
	 !$ 	END FUNCTION OMP_get_wtime
	 !$ END INTERFACE
		LOGICAL,INTENT(IN),OPTIONAL :: total
		INTEGER :: clipboard
	 !$ REAL(8) :: clipboard_real
		INTEGER,SAVE :: timeline=0
	 !$ REAL(8),SAVE :: timeline_real=0.0d0
		 !$ clipboard_real=OMP_get_wtime()
		 !$ IF (PRESENT(total).AND.(TIME_OUTPUT)) THEN
		 !$ 	IF (total) THEN
		 !$ 		WRITE(*,ADVANCE="NO",FMT='(" TOTAL elapsed time: ")')
		 !$ 		CALL user_friendly_time_output(clipboard_real-timeline_begin_real)
		 !$ 	ENDIF
		 !$ 	WRITE(*,*)
		 !$ 	RETURN
		 !$ ENDIF
		 !$ IF ((timeline_real>0.0d0).AND.(TIME_OUTPUT)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='(A)') " elapsed time: "
		 !$ 	CALL user_friendly_time_output(clipboard_real-timeline_real)
		 !$ ELSE
		 !$ 	timeline_begin_real=clipboard_real
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
			IF (VERBOSE_OUTPUT) CALL print_progress(startstep_in-endstep_in)
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
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			WRITE(4,*)
			WRITE(4,*)
			ENDFILE 4
			CLOSE(UNIT=4)
			CLOSE(UNIT=10)
			IF (((startstep_in-endstep_in)>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
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
			IF (VERBOSE_OUTPUT) CALL print_progress(startstep_in-endstep_in)
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
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			IF (((startstep_in-endstep_in)>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
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

		SUBROUTINE dump_split(startstep_in,endstep_in,output_format)
		IMPLICIT NONE
		INTEGER :: startstep,endstep
		INTEGER :: molecule_type_index,stepcounter,moleculecounter,atomcount
		INTEGER,INTENT(IN) :: startstep_in,endstep_in
		LOGICAL :: connected
		CHARACTER(LEN=1024) :: fstring
		CHARACTER(LEN=3),INTENT(IN) :: output_format
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
				IF (VERBOSE_OUTPUT) THEN
					IF ((endstep-startstep)>100) WRITE(*,*)
					CALL print_progress(endstep-startstep)
				ENDIF
				WRITE(fstring,'(2A,I0,A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),&
				&"traj_",molecule_type_index,"."//output_format
				INQUIRE(UNIT=4,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=4)
				OPEN(UNIT=4,FILE=TRIM(fstring),STATUS="REPLACE")
				atomcount=give_number_of_atoms_per_molecule(molecule_type_index)*give_number_of_molecules_per_step(molecule_type_index)
				DO stepcounter=startstep,endstep,1
					!Write head, depending on which type the trajectory has...
					CALL write_header(4,stepcounter*TIME_SCALING_FACTOR,atomcount,output_format)
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
					IF (VERBOSE_OUTPUT) CALL print_progress()
				ENDDO
				ENDFILE 4
				CLOSE(UNIT=4)
				IF (VERBOSE_OUTPUT) THEN
					IF ((give_number_of_timesteps())>100) THEN
						WRITE(*,*)
						WRITE(*,FMT='(" ")',ADVANCE="NO")
					ENDIF
				ENDIF
				WRITE(*,'("done.")')
			ENDDO
			CLOSE(UNIT=3)
		END SUBROUTINE dump_split

		SUBROUTINE convert(writemolecularinputfile,output_format,useCOC_in)
		IMPLICIT NONE
		CHARACTER(LEN=3),INTENT(IN) :: output_format
		INTEGER :: molecule_type_index,stepcounter,moleculecounter,output(3),ncentres,valid_COC
		LOGICAL,INTENT(IN) :: writemolecularinputfile
		LOGICAL,INTENT(IN),OPTIONAL :: useCOC_in
		LOGICAL :: connected,useCOC
		CHARACTER(LEN=1024) :: fstring
		CHARACTER(LEN=1) :: element
			IF (PRESENT(useCOC_in)) THEN
				useCOC=useCOC_in
			ELSE
				useCOC=.FALSE.
			ENDIF
			IF (DEVELOPERS_VERSION) THEN
				PRINT *," ! CAREFUL: 'PARALLELISED' CONVERTER USED"
				CALL convert_parallel()
				RETURN
			ELSE
				!first, get the number of lines that will be added.
				ncentres=0
				valid_COC=0
				IF (useCOC) THEN
					!count only those molecules with nonzero charge
					DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
						IF (give_charge_of_molecule(molecule_type_index)/=0) THEN
							valid_COC=valid_COC+1
							ncentres=ncentres+give_number_of_molecules_per_step(molecule_type_index)
						ENDIF
					ENDDO
				ELSE
					!count everything
					DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
						ncentres=ncentres+give_number_of_molecules_per_step(molecule_type_index)
					ENDDO
				ENDIF
				IF (ncentres==0) CALL report_error(143)
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				IF (useCOC) THEN
					WRITE(fstring,'(3A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_COC.",output_format
				ELSE
					WRITE(fstring,'(3A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_COM.",output_format
				ENDIF
				OPEN(UNIT=3,FILE=TRIM(fstring))
				IF (VERBOSE_OUTPUT) THEN
					WRITE(*,ADVANCE="NO",FMT='(" Writing new trajectory to file ",A,"...")') "'"//TRIM(fstring)//"'"
					IF (give_number_of_timesteps()>100) WRITE(*,*)
					CALL print_progress(give_number_of_timesteps())
				ENDIF
				DO stepcounter=1,give_number_of_timesteps(),1
					!Write head, depending on which type the trajectory has...
					CALL write_header(3,stepcounter*TIME_SCALING_FACTOR,ncentres,output_format)
					DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
						element=CHAR(ALPHABET_small(MOD((molecule_type_index-1),26)+1)) !assign the element names a,b,c,... to the centred molecules.
						IF (useCOC) THEN
							!use the centre of charge - but only for charged molecules.
							IF (give_charge_of_molecule(molecule_type_index)/=0) THEN
								DO moleculecounter=1,give_number_of_molecules_per_step(molecule_type_index),1
									!Sort of high accuracy should be kept here because of the way I use this routine.
									!reduced to 16.8 from 18.10
									WRITE(3,'(A1,3E16.8)') element,give_center_of_charge(stepcounter,molecule_type_index,moleculecounter)
								ENDDO
							ENDIF
						ELSE
							DO moleculecounter=1,give_number_of_molecules_per_step(molecule_type_index),1
								!Sort of high accuracy should be kept here because of the way I use this routine.
								!reduced to 16.8 from 18.10
								WRITE(3,'(A1,3E16.8)') element,give_center_of_mass(stepcounter,molecule_type_index,moleculecounter)
							ENDDO
						ENDIF
					ENDDO
					IF (VERBOSE_OUTPUT) CALL print_progress()
				ENDDO
				IF (VERBOSE_OUTPUT) THEN
					IF ((give_number_of_timesteps())>100) THEN
						WRITE(*,*)
						WRITE(*,FMT='(" ")',ADVANCE="NO")
					ENDIF
					WRITE(*,'("done.")')
				ENDIF
				CLOSE(UNIT=3)
			ENDIF
			IF (writemolecularinputfile) THEN
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				IF (useCOC) THEN
					WRITE(fstring,'(2A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"COC_molecular.inp"
				ELSE
					WRITE(fstring,'(2A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"COM_molecular.inp"
				ENDIF
				OPEN(UNIT=3,FILE=TRIM(fstring))
				WRITE(*,ADVANCE="NO",FMT='(" Writing new molecular input file in ",A,"...")') "'"//TRIM(fstring)//"'"
				WRITE(3,'(I0," ### number of timesteps")') give_number_of_timesteps()
				WRITE(3,'(I0," ### number of different types of molecules. Followed by list of molecules.")') give_number_of_molecule_types()
				output(2)=1
				DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
					output(1)=give_charge_of_molecule(molecule_type_index)
					output(3)=give_number_of_molecules_per_step(molecule_type_index)
					IF (useCOC) THEN
						!use the centre of charge - but only for charged molecules.
						IF (give_charge_of_molecule(molecule_type_index)/=0) THEN
							WRITE(3,ADVANCE="NO",FMT='(SP,I2,SS," ",I0," ",I0," ### ")') output(:) !write the crucial part
							WRITE(3,'("There are ",I0," molecules per step with charge ",SP,I2," (",F0.3,SS,"), given as centre of charge.")')& !write the comments
							& output(3),output(1),give_realcharge_of_molecule(molecule_type_index)
						ENDIF
					ELSE
						WRITE(3,ADVANCE="NO",FMT='(SP,I2,SS," ",I0," ",I0," ### ")') output(:) !write the crucial part
						WRITE(3,'("There are ",I0," molecules per step with charge ",SP,I2,SS,", given as centre of mass.")')& !write the comments
						& output(3),output(1)
					ENDIF
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
			WRITE(*,FMT='(A)')" ('maxdist' is the largest intramolecular distance of any atom from the COM)."
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
		SUBROUTINE remove_drudes(startstep_in,endstep_in,output_format,writemolecularinputfile)
		IMPLICIT NONE
		INTEGER :: stepcounter,molecule_type_index,molecule_index
		INTEGER,INTENT(IN) :: startstep_in,endstep_in
		LOGICAL,INTENT(IN) :: writemolecularinputfile
		LOGICAL :: connected
		CHARACTER(LEN=1024) :: fstring
		CHARACTER(LEN=3),INTENT(IN) :: output_format
			!open the output file
			INQUIRE(UNIT=4,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=4)
			IF (startstep_in==endstep_in) THEN
				WRITE(fstring,'(2A,I0,A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_nodrudes_step_",startstep_in,"."//output_format
			ELSE
				WRITE(fstring,'(2A,I0,A,I0,A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_nodrudes_step_",startstep_in,"-",endstep_in,"."//output_format
			ENDIF
			OPEN(UNIT=4,FILE=TRIM(fstring))
			!iterate over the specified timesteps
			IF (VERBOSE_OUTPUT) THEN
				WRITE(*,FMT='(A)',ADVANCE="NO") " writing new trajectory to file '"//TRIM(fstring)//"'..."
				IF ((endstep_in-startstep_in)>100) WRITE(*,*)
				CALL print_progress(endstep_in-startstep_in)
			ENDIF
			DO stepcounter=startstep_in,endstep_in,1
				!Write head, depending on which type the trajectory has...
				CALL write_header(4,stepcounter*TIME_SCALING_FACTOR,&
				&(give_number_of_atoms_per_step()-give_number_of_drude_particles()),output_format)
				DO molecule_type_index=1,give_number_of_molecule_types(),1
					DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
						!the following line appends one molecule to the output trajectory.
						CALL write_molecule_merged_drudes(4,stepcounter,molecule_type_index,molecule_index,include_header=.FALSE.)
					ENDDO
				ENDDO
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			IF (VERBOSE_OUTPUT) THEN
				IF ((endstep_in-startstep_in)>100) THEN
					WRITE(*,*)
					WRITE(*,FMT='(" ")',ADVANCE="NO")
				ENDIF
				WRITE(*,'("done.")')
			ENDIF
			ENDFILE 4
			CLOSE(UNIT=4)
			IF (writemolecularinputfile) CALL write_molecule_input_file_without_drudes(endstep_in-startstep_in+1)
		END SUBROUTINE

		!This SUBROUTINE writes a SUBROUTINE with removed core particles. This requires assigned drude particles!
		SUBROUTINE remove_cores(startstep_in,endstep_in,output_format)
		IMPLICIT NONE
		INTEGER :: stepcounter,molecule_type_index,molecule_index
		INTEGER,INTENT(IN) :: startstep_in,endstep_in
		LOGICAL :: connected
		CHARACTER(LEN=1024) :: fstring
		CHARACTER(LEN=3),INTENT(IN) :: output_format
			!open the output file
			INQUIRE(UNIT=4,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=4)
			IF (startstep_in==endstep_in) THEN
				WRITE(fstring,'(2A,I0,A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_nocores_step_",startstep_in,"."//output_format
			ELSE
				WRITE(fstring,'(2A,I0,A,I0,A)') &
				&TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX)),"traj_nocores_step_",startstep_in,"-",endstep_in,"."//output_format
			ENDIF
			OPEN(UNIT=4,FILE=TRIM(fstring))
			!iterate over the specified timesteps
			IF (VERBOSE_OUTPUT) THEN
				WRITE(*,FMT='(A)',ADVANCE="NO") " writing new trajectory to file '"//TRIM(fstring)//"'..."
				IF ((endstep_in-startstep_in)>100) WRITE(*,*)
				CALL print_progress(endstep_in-startstep_in)
			ENDIF
			DO stepcounter=startstep_in,endstep_in,1
				!Write head, depending on which type the trajectory has...
				CALL write_header(4,stepcounter*TIME_SCALING_FACTOR,&
				&(give_number_of_drude_particles()),output_format)
				DO molecule_type_index=1,give_number_of_molecule_types(),1
					DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
						!the following line appends the drudes of one molecule to the output trajectory.
						CALL write_only_drudes_relative_to_core(4,stepcounter,molecule_type_index,molecule_index,include_header=.FALSE.)
					ENDDO
				ENDDO
				IF (VERBOSE_OUTPUT) CALL print_progress()
			ENDDO
			IF (VERBOSE_OUTPUT) THEN
				IF ((endstep_in-startstep_in)>100) THEN
					WRITE(*,*)
					WRITE(*,FMT='(" ")',ADVANCE="NO")
				ENDIF
				WRITE(*,'("done.")')
			ENDIF
			ENDFILE 4
			CLOSE(UNIT=4)
		END SUBROUTINE remove_cores

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
				54	FORMAT (EN9.1," K")
				55	FORMAT (F5.1," K")
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
	50	FORMAT ("   T",A1,": ",EN9.1," K (",EN9.1," K)")
	51	FORMAT ("   T",A1,": ",F5.1," K (",F5.1," K)")
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
				52	FORMAT ("     T",A,": ",EN9.1," K")
				53	FORMAT ("     T",A,": ",F5.1," K")
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