
!This Module performs low-level command line tasks.
MODULE RECOGNITION ! Copyright (C) !RELEASEYEAR! Frederik Philippi
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
		 !$ 	FUNCTION OMP_get_max_threads()
		 !$ 	INTEGER :: OMP_get_max_threads
		 !$ 	END FUNCTION OMP_get_max_threads
		 !$ 	SUBROUTINE OMP_set_num_threads(number_of_threads)
		 !$ 	INTEGER,INTENT(IN) :: number_of_threads
		 !$ 	END SUBROUTINE OMP_set_num_threads
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
			 !$ 	FUNCTION OMP_get_num_threads()
			 !$ 	INTEGER :: OMP_get_num_threads
			 !$ 	END FUNCTION OMP_get_num_threads
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
						!$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (brute force recognition)"
						!$ 	CALL timing_parallel_sections(.TRUE.)
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
				 !$ 	WRITE(*,ADVANCE="NO",FMT='("     ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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
						WRITE(*,*) '  #set_threads_simple'
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
							WRITE(8,*) '#set_threads_simple'
							WRITE(8,*) '#trajectory_type lmp'
							WRITE(8,*) 'quit'
							WRITE(8,*)
							WRITE(8,*)
							IF (number_of_drude_particles/=0) THEN
								WRITE(8,*) 'remove_drudes_simple'
								WRITE(8,*) 'remove_cores_simple'
								WRITE(8,*) 'drude_temp_simple'
							ENDIF
							WRITE(8,*)
							WRITE(8,*) 'diffusion_simple'
							WRITE(8,*) 'distribution_simple'
							WRITE(8,*) 'dump_snapshot_simple'
							WRITE(8,*) 'dump_split_simple'
							WRITE(8,*) 'contact_distance_simple'
							WRITE(8,*) 'dump_neighbour_traj_simple'
							WRITE(8,*) 'convert_simple'
							WRITE(8,*) 'temperature_simple'
							WRITE(8,*) 'gyradius_simple'
							WRITE(8,*) 'charge_arm_simple'
							WRITE(8,*) 'conductivity_simple'
							WRITE(8,*) 'jump_velocity_simple'
							WRITE(8,*)
							WRITE(8,*)
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
			 !$ 	FUNCTION OMP_get_num_threads()
			 !$ 	INTEGER :: OMP_get_num_threads
			 !$ 	END FUNCTION OMP_get_num_threads
			 !$ END INTERFACE
				INTEGER :: linecounter1,linecounter2
					!$OMP SINGLE
					!$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
					!$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (brute force connectivity)"
					!$ 	CALL timing_parallel_sections(.TRUE.)
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
				 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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