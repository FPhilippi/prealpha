
!This module is responsible for handling the trajectory and passing information to other modules.
MODULE MOLECULAR ! Copyright (C) !RELEASEYEAR! Frederik Philippi
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
	REAL,PARAMETER :: default_mass_sodium=22.990

	REAL,PARAMETER :: default_charge_hydrogen=0.0
	REAL,PARAMETER :: default_charge_fluorine=0.0
	REAL,PARAMETER :: default_charge_boron=0.0
	REAL,PARAMETER :: default_charge_chlorine=0.0
	REAL,PARAMETER :: default_charge_bromine=0.0
	REAL,PARAMETER :: default_charge_sodium=0.0
	REAL,PARAMETER :: default_charge_iodine=0.0
	REAL,PARAMETER :: default_charge_nitrogen=0.0
	REAL,PARAMETER :: default_charge_oxygen=0.0
	REAL,PARAMETER :: default_charge_carbon=0.0
	REAL,PARAMETER :: default_charge_sulfur=0.0
	REAL,PARAMETER :: default_charge_phosphorus=0.0
	REAL,PARAMETER :: default_charge_lithium=0.0
	!variables
	REAL :: charge_hydrogen=default_charge_hydrogen
	REAL :: charge_fluorine=default_charge_fluorine
	REAL :: charge_boron=default_charge_boron
	REAL :: charge_chlorine=default_charge_chlorine
	REAL :: charge_bromine=default_charge_bromine
	REAL :: charge_iodine=default_charge_iodine
	REAL :: charge_nitrogen=default_charge_nitrogen
	REAL :: charge_oxygen=default_charge_oxygen
	REAL :: charge_carbon=default_charge_carbon
	REAL :: charge_sulfur=default_charge_sulfur
	REAL :: charge_phosphorus=default_charge_phosphorus
	REAL :: charge_lithium=default_charge_lithium
	REAL :: charge_sodium=default_charge_sodium

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
	REAL :: mass_sodium=default_mass_sodium
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
		REAL :: realcharge=0.0
		INTEGER :: number_of_atoms=0 ! = extent of dimension one of trajectory, number of atoms PER SINGLE MOLECULE, not total!
		INTEGER :: total_molecule_count=0 ! = extent of dimension two of trajectory, number of molecules of this type in the box
		INTEGER :: number_of_drudes_in_molecule=0 ! = number of drude particles in this molecule (counter for successfully assigned drude pairs)
        REAL(KIND=GENERAL_PRECISION) :: mass=0.0d0
		CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE :: list_of_elements !--> Turned on support for  2-letter elements!
		TYPE(drude_pair),DIMENSION(:),ALLOCATABLE :: list_of_drude_pairs ! list of pairs of drude particles / drude cores / drudes
		REAL(KIND=GENERAL_PRECISION),DIMENSION(:),ALLOCATABLE :: list_of_atom_masses !corresponding masses for the atoms
		REAL(KIND=GENERAL_PRECISION),DIMENSION(:),ALLOCATABLE :: list_of_atom_charges !corresponding charges for the atoms
		CHARACTER(LEN=5) :: residue_name !The residue name adhering to GROMACS specifications
		LOGICAL,DIMENSION(:),ALLOCATABLE :: manual_atom_mass_specified !any elements that are .TRUE. will not be initialised from the trajectory!
		LOGICAL,DIMENSION(:),ALLOCATABLE :: manual_atom_charge_specified !any elements that are .TRUE. will not be initialised from the trajectory!
		TYPE(atom),DIMENSION(:,:),ALLOCATABLE :: snapshot !like trajectory, but for one timestep only. Required for READ_SEQUENTIAL.
		TYPE(atom),DIMENSION(:,:,:),ALLOCATABLE :: trajectory
		!trajectory is now in column major order!
		!first dimension: Atom index in molecule
		!second dimension: index of the molecule (not to be confused with molecule_type_index)
		!third dimension: timestep
		TYPE(atom),DIMENSION(:,:),ALLOCATABLE :: snapshot_2 !alternative storage if position and velocity are required simultaneously
		TYPE(atom),DIMENSION(:,:,:),ALLOCATABLE :: trajectory_2 !alternative storage if position and velocity are required simultaneously
		TYPE(atom),DIMENSION(:,:,:),ALLOCATABLE :: queue!the queue for circular parallel operation. conceptually the same as 'trajectory', but only big enough to act as a buffer.
	END TYPE molecule
	REAL(KIND=STORAGE_PRECISION) :: box_dimensions(2,3)!low and high for x,y,z
	REAL(KIND=STORAGE_PRECISION) :: box_size(3) !size of the box.
	REAL(KIND=STORAGE_PRECISION) :: maximum_distance !maximum possible distance in box.
	REAL(KIND=STORAGE_PRECISION) :: maximum_distance_squared !square of maximum possible distance in box.
    REAL(KIND=SP) :: drude_mass=0.0e0
	REAL(KIND=SP) :: drude_charge=0.0e0
	REAL(KIND=SP) :: COM_mass_list(IACHAR("a"):(IACHAR("a")+25)) !A list of user-specified masses for centre-of-mass trajectories.
	LOGICAL :: custom_default_masses !if 'T', then the user has specified his own default masses for lowercase letters or elements.
	LOGICAL :: custom_atom_masses !if 'T', then the user has manually specified atom masses.
	LOGICAL :: custom_default_charges !like custom_default_masses, but for the charges - and not for lowercase letters.
	LOGICAL :: custom_atom_charges !like custom_atom_masses, but for charges.
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
	PUBLIC :: print_atomic_masses,give_comboost,switch_to_barycenter,print_atomic_charges,set_default_charges,charge_arm,check_charges
	PUBLIC :: print_dipole_statistics,write_molecule_input_file_without_drudes,write_only_drudes_relative_to_core
	PUBLIC :: give_number_of_specific_atoms_per_molecule,give_indices_of_specific_atoms_per_molecule,give_number_of_specific_atoms
	PUBLIC :: give_indices_of_specific_atoms,give_element_symbol,give_center_of_charge,give_realcharge_of_molecule,write_trajectory
	PUBLIC :: give_charge_of_atom,give_mass_of_atom,give_center_of_charge_2,give_center_of_mass_2,charge_arm_2,write_molecule_in_slab
	PUBLIC :: give_box_boundaries
	CONTAINS

		LOGICAL FUNCTION check_charges(molecule_type_index)
		IMPLICIT NONE
		INTEGER :: atom_index
		INTEGER,INTENT(IN) :: molecule_type_index
			check_charges=.FALSE.
			DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
				IF (ABS(molecule_list(molecule_type_index)%list_of_atom_charges(atom_index))>0.001) THEN
					check_charges=.TRUE.
					RETURN
				ENDIF
			ENDDO
		END FUNCTION check_charges

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
			mass_sodium=default_mass_sodium
		END SUBROUTINE set_default_masses

		SUBROUTINE set_default_charges()
		IMPLICIT NONE
			charge_hydrogen=default_charge_hydrogen
			charge_fluorine=default_charge_fluorine
			charge_boron=default_charge_boron
			charge_chlorine=default_charge_chlorine
			charge_bromine=default_charge_bromine
			charge_iodine=default_charge_iodine
			charge_nitrogen=default_charge_nitrogen
			charge_oxygen=default_charge_oxygen
			charge_carbon=default_charge_carbon
			charge_sulfur=default_charge_sulfur
			charge_phosphorus=default_charge_phosphorus
			charge_lithium=default_charge_lithium
			charge_sodium=default_charge_sodium
		END SUBROUTINE set_default_charges

		SUBROUTINE subtract_drude_masses()
		IMPLICIT NONE
			IF (VERBOSE_OUTPUT) PRINT *,"Subtracting drude masses from N,O,C,S,P,Li,F"
			mass_nitrogen=mass_nitrogen-drude_mass
			mass_oxygen=mass_oxygen-drude_mass
			mass_carbon=mass_carbon-drude_mass
			mass_sulfur=mass_sulfur-drude_mass
			mass_phosphorus=mass_phosphorus-drude_mass
			mass_lithium=mass_lithium-drude_mass
			mass_fluorine=mass_fluorine-drude_mass
		END SUBROUTINE subtract_drude_masses

		SUBROUTINE write_molecule_input_file_without_drudes(nsteps)
		LOGICAL :: connected
		CHARACTER(LEN=1024) :: fstring,chargestring
		INTEGER,INTENT(IN) :: nsteps
		REAL :: DrudeCharge,CoreCharge
		INTEGER :: output(3),atom_index,natoms,counter,drude_flag,molecule_type_index
			INQUIRE(UNIT=3,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=3)
			WRITE(fstring,'(3A)') TRIM(PATH_OUTPUT),TRIM(ADJUSTL(OUTPUT_PREFIX)),"molecular_nodrude.inp"
			OPEN(UNIT=3,FILE=TRIM(fstring))
			WRITE(*,ADVANCE="NO",FMT='(" Writing new molecular input file in ",A,"...")') "'"//TRIM(fstring)//"'"
			WRITE(3,'(I0," ### number of timesteps")') nsteps
			WRITE(3,'(I0," ### number of different types of molecules. Followed by list of molecules.")') give_number_of_molecule_types()
			natoms=0
			DO molecule_type_index=1,give_number_of_molecule_types(),1 !iterate over number of molecule types. (i.e. cation and anion, usually)
				output(1)=give_charge_of_molecule(molecule_type_index)
				output(2)=0
				DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
					!drude_flag is the atom index of the drude particle attached to this core.
					!Will be -1 if no drude particle is found, and 0 if this is a drude itself.
					drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
					IF (drude_flag/=0) output(2)=output(2)+1
					!get total number of output atoms
					IF (drude_flag/=0) natoms=natoms+1
				ENDDO
				output(3)=give_number_of_molecules_per_step(molecule_type_index)
				WRITE(3,ADVANCE="NO",FMT='(SP,I2,SS," ",I0," ",I0," ### ")')output(:)
				!write the crucial part
				WRITE(3,'("There are ",I0," molecules per step with charge ",SP,I2,SS,".")')output(3),output(1)
			ENDDO
			!write the custom charges section.
			IF (((custom_atom_charges).OR.(custom_default_charges))) THEN
				WRITE(3,'("atomic_charges ",I0," ### The following lines contain the atomic charges.")') natoms
				DO molecule_type_index=1,give_number_of_molecule_types(),1
					counter=1
					DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
						drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
						SELECT CASE (drude_flag)
						CASE (-1) !non-polarisable atom - print just as it is.
							WRITE(3,'("   ",I0," ",I0," ",F0.4," (Element ",A,", non-polarisable)")')&
							&molecule_type_index,counter,&
							&molecule_list(molecule_type_index)%list_of_atom_charges(atom_index),&
							&TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))
							counter=counter+1
						CASE (0) !this is a drude particle - skip this one by cycling.
							CYCLE
						CASE DEFAULT
							!everything else should be drude cores - merge with the drude particle.
							CoreCharge=molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
							DrudeCharge=molecule_list(molecule_type_index)%list_of_atom_charges(drude_flag)
							IF ((CoreCharge<1.0).AND.(CoreCharge>-1.0)) THEN
								IF (CoreCharge>=0.0) THEN
									WRITE(chargestring,'(F6.4)') CoreCharge
								ELSE
									WRITE(chargestring,'(F7.4)') CoreCharge
								ENDIF
							ELSE
								WRITE(chargestring,'(F0.4)') CoreCharge
							ENDIF

							IF ((DrudeCharge<1.0).AND.(DrudeCharge>-1.0)) THEN
								IF (DrudeCharge>=0.0) THEN
									WRITE(chargestring,'(A,"+",F6.4)') TRIM(chargestring),DrudeCharge
								ELSE
									WRITE(chargestring,'(A,F7.4)') TRIM(chargestring),DrudeCharge
								ENDIF
							ELSE
								IF (DrudeCharge>=0.0) THEN
									WRITE(chargestring,'(A,"+",F0.4)') TRIM(chargestring),DrudeCharge
								ELSE
									WRITE(chargestring,'(A,F0.4)') TRIM(chargestring),DrudeCharge
								ENDIF
							ENDIF
							WRITE(3,'("   ",I0," ",I0," ",F0.4," (Element ",A,", Core+Drude=",A,")")')&
							&molecule_type_index,counter,&
							&molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)+&
							&molecule_list(molecule_type_index)%list_of_atom_charges(drude_flag),&
							&TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index)),&
							&TRIM(chargestring)
							counter=counter+1
						END SELECT
					ENDDO
				ENDDO
			ENDIF
			!write the custom masses section.
			WRITE(3,'("atomic_masses ",I0," ### The following lines contain the atomic masses.")') natoms
			DO molecule_type_index=1,give_number_of_molecule_types(),1
				counter=1
				DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
					drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
					SELECT CASE (drude_flag)
					CASE (-1) !non-polarisable atom - print just as it is.
						WRITE(3,'("   ",I0," ",I0," ",F0.3)') molecule_type_index,counter,&
						&molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
						counter=counter+1
					CASE (0) !this is a drude particle - skip this one by cycling.
						CYCLE
					CASE DEFAULT
						!everything else should be drude cores - merge with the drude particle.
						WRITE(3,'("   ",I0," ",I0," ",F0.3)') molecule_type_index,counter,&
						&molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)+&
						&molecule_list(molecule_type_index)%list_of_atom_masses(drude_flag)
						counter=counter+1
					END SELECT
				ENDDO
			ENDDO
			CLOSE(UNIT=3)
			WRITE(*,*) "done"
		END SUBROUTINE write_molecule_input_file_without_drudes

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

		SUBROUTINE print_atomic_charges()
		IMPLICIT NONE
		INTEGER :: molecule_type_index,atom_index
		REAL :: native_charge,atom_charge
		CHARACTER(LEN=2) :: element_name
			WRITE(*,'(" The following lines can be used/altered to request custom atomic charges:")')
			WRITE(*,'("   atomic_charges ",I0)') SUM(molecule_list(:)%number_of_atoms)
			DO molecule_type_index=1,number_of_molecule_types,1
				DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
					atom_charge=molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
					native_charge=atomic_charge(molecule_list(molecule_type_index)%list_of_elements(atom_index))
					element_name=molecule_list(molecule_type_index)%list_of_elements(atom_index)
					WRITE(*,FMT='("   ",I0," ",I0," ",F0.4," (")',ADVANCE="NO") molecule_type_index,atom_index,atom_charge
					IF ((ABS(atom_charge-native_charge)>0.001).AND.(ABS(native_charge)>0.001)) THEN
						WRITE(*,'(A,", default charge was ",F0.4,")")') TRIM(element_name),native_charge
					ELSE
						WRITE(*,'("Element ",A,")")') TRIM(element_name)
					ENDIF
				ENDDO
			ENDDO
		END SUBROUTINE print_atomic_charges

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
		REAL(KIND=WORKING_PRECISION),INTENT(OUT) :: rgy_sq
		REAL(KIND=WORKING_PRECISION),INTENT(OUT),OPTIONAL :: maxdist
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
			IF (PRESENT(maxdist)) maxdist=SQRT(maxdist_squared)
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
	 !$ 	FUNCTION OMP_get_thread_num()
	 !$ 	INTEGER :: OMP_get_thread_num
	 !$ 	END FUNCTION OMP_get_thread_num
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER,INTENT(IN) :: timestep,molecule_type_index_1
		REAL,INTENT(OUT) :: smallest
		INTEGER :: molecule_index_1
		REAL :: smallest_local
			!two atoms can be no further apart than the diagonale of the box... that's what I initialise to
			!$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(smallest_local)
			!$OMP SINGLE
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE (*,'(A,I0,A)') "     ### Parallel execution on ",OMP_get_num_threads()," threads (intermolecular contact distance)"
		 !$ 	CALL timing_parallel_sections(.TRUE.)
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
		 !$ 	WRITE(*,ADVANCE="NO",FMT='("     ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
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

		!This FUNCTION returns the smallest squared distance of 2 atoms considering all PBCs.
		REAL(KIND=GENERAL_PRECISION) FUNCTION give_smallest_atom_distance_squared&
		&(timestep1,timestep2,molecule_type_index_1,molecule_type_index_2,&
		&molecule_index_1,molecule_index_2,atom_index_1,atom_index_2,translation)
		IMPLICIT NONE
		INTEGER :: a,b,c
		INTEGER,INTENT(IN) :: timestep1,timestep2,molecule_type_index_1,molecule_type_index_2
		INTEGER,INTENT(IN) :: molecule_index_1,molecule_index_2,atom_index_1,atom_index_2
		REAL(KIND=WORKING_PRECISION) :: pos_1(3),pos_2(3),shift(3),distance_clip,wrapshift(3)
		REAL(KIND=WORKING_PRECISION),INTENT(OUT),OPTIONAL :: translation(3)
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
			! The following is always needed, because the trajectory is wrapped molecule-wise.
			CALL wrap_vector(pos_1,shift(:))
			CALL wrap_vector(pos_2,wrapshift(:))
			!"shift" is now the vector to translate the reference atom into the box by wrapping.
			!"wrapshift" is the same for the second, observed atom.
			!now, store in wrapshift the vector to bring the second *into the same box* as the first molecule:
			wrapshift(:)=wrapshift(:)-shift(:)
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
							IF (PRESENT(translation)) translation(:)=shift(:)
						ENDIF
				   ENDDO
				ENDDO
			ENDDO
			IF (PRESENT(translation)) translation(:)=translation(:)+wrapshift(:)
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
	 !$ 	FUNCTION OMP_get_thread_num()
	 !$ 	INTEGER :: OMP_get_thread_num
	 !$ 	END FUNCTION OMP_get_thread_num
	 !$ 	FUNCTION OMP_get_max_threads()
	 !$ 	INTEGER :: OMP_get_max_threads
	 !$ 	END FUNCTION OMP_get_max_threads
	 !$ 	SUBROUTINE OMP_set_nested(enable)
	 !$ 	LOGICAL,INTENT(IN) :: enable
	 !$ 	END SUBROUTINE OMP_set_nested
	 !$ END INTERFACE
		 !$ nprocs=OMP_get_max_threads()
		 !$ IF (nprocs<3) THEN
		 !$ 	WRITE(*,*) " ! less than 3 threads available. Returning control to main module."
		 !$ 	RETURN
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
			1	IF (element_status((MOD(stepcounter-1,queuesize)+1))==2) THEN
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
				2	skip=.FALSE.
					IF (completed) THEN
				 !$ 	WRITE(*,'("  ! Thread number ",I0," caught termination signal.")') OMP_get_thread_num()
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
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
			!$OMP PARALLEL IF(PARALLEL_OPERATION) PRIVATE(molecule_type_index,molecule_index,centre_of_mass,xyzcounter)
			!$OMP SINGLE
			!start the timer for the parallel section.
		 !$ CALL timing_parallel_sections(.TRUE.)
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION))&
		 !$ &WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (wrapping)"
			IF (VERBOSE_OUTPUT) CALL print_progress(number_of_steps)
			!$OMP END SINGLE
			!$OMP DO SCHEDULE(STATIC,1)
			DO stepcounter=1,number_of_steps,1
				DO molecule_type_index=1,number_of_molecule_types,1
					DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
						centre_of_mass=give_center_of_mass(stepcounter,molecule_type_index,molecule_index)
						DO xyzcounter=1,3,1
						!See? I got rid of it! (The jump)
							DO
								IF (centre_of_mass(xyzcounter)<box_dimensions(1,xyzcounter)) THEN
									!centre of mass is outside of box (smaller)
									centre_of_mass(xyzcounter)=centre_of_mass(xyzcounter)+box_size(xyzcounter)
									molecule_list(molecule_type_index)%trajectory(:,molecule_index,stepcounter)%coordinates(xyzcounter)=&
									&molecule_list(molecule_type_index)%trajectory(:,molecule_index,stepcounter)%coordinates(xyzcounter)&
									&+box_size(xyzcounter)
									CYCLE
								ELSE
									IF (centre_of_mass(xyzcounter)>box_dimensions(2,xyzcounter)) THEN
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
				!$OMP CRITICAL
				IF (VERBOSE_OUTPUT) CALL print_progress()
				!$OMP END CRITICAL
			ENDDO
			!$OMP END DO
			!$OMP END PARALLEL
			IF ((number_of_steps>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
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
					10  IF (centre_of_mass(xyzcounter)<box_dimensions(1,xyzcounter)) THEN
							centre_of_mass(xyzcounter)=centre_of_mass(xyzcounter)+box_size(xyzcounter)
							molecule_list(molecule_type_index)%snapshot(:,molecule_index)%coordinates(xyzcounter)=&
							&molecule_list(molecule_type_index)%snapshot(:,molecule_index)%coordinates(xyzcounter)+box_size(xyzcounter)
							GOTO 10
						ELSE
							IF (centre_of_mass(xyzcounter)>box_dimensions(2,xyzcounter)) THEN
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
			WRITE(*,'("    ",A," ",A)') "dihedrals_initialised     ",TRIM(logical_to_yesno(dihedrals_initialised))
			WRITE(*,'("    ",A," ",A)') "custom_constraints        ",TRIM(logical_to_yesno(custom_constraints))
			WRITE(*,'("    ",A," ",A)') "custom_default_masses     ",TRIM(logical_to_yesno(custom_default_masses))
			WRITE(*,'("    ",A," ",A)') "custom_atom_masses        ",TRIM(logical_to_yesno(custom_atom_masses))
			WRITE(*,'("    ",A," ",A)') "custom_default_charges    ",TRIM(logical_to_yesno(custom_default_charges))
			WRITE(*,'("    ",A," ",A)') "custom_atom_charges       ",TRIM(logical_to_yesno(custom_atom_charges))
			WRITE(*,'("    ",A," ",A)') "drudes_assigned           ",TRIM(logical_to_yesno(drudes_assigned))
			WRITE(*,'("    ",A," ",A)') "drude_details             ",TRIM(logical_to_yesno(drude_details))
			WRITE(*,'("    ",A," ",A)') "drudes_allocated          ",TRIM(logical_to_yesno(drudes_allocated))
			WRITE(*,'("    ",A," ",A)') "use_firstatom_as_com      ",TRIM(logical_to_yesno(use_firstatom_as_com))
			WRITE(*,'("    ",A," ",A)') "use_barycentre            ",TRIM(logical_to_yesno(use_barycentre))
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
				IF (EXTRA_VELOCITY) required_storage=required_storage*2.0
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
				IF (EXTRA_VELOCITY) THEN
					WRITE(*,'(" 6*",I0,"*",I0,"*4Byte =",F6.1,A2," (single precision)")')& !printing memory requirement
					&total_number_of_atoms,number_of_steps,required_storage,memory_suffix
				ELSE
					WRITE(*,'(" 3*",I0,"*",I0,"*4Byte =",F6.1,A2," (single precision)")')& !printing memory requirement
					&total_number_of_atoms,number_of_steps,required_storage,memory_suffix
				ENDIF
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
				IF (EXTRA_VELOCITY) THEN
					CALL read_snapshot_body_posvel()
				ELSE
					CALL read_snapshot_body()
				ENDIF
			ELSE!timestep MUST be smaller than the file position - have to go back!
				!How far back is that new, requested timestep?
				IF (timestep>(file_position/2)) THEN
					!Use backspacing.
					!Dear future self, the following line might look weird, but it works. Trust me.
					DO counter=-1,lines_to_skip*(timestep-file_position-1),-1!
						BACKSPACE 9
					ENDDO
					IF (EXTRA_VELOCITY) THEN
						CALL read_snapshot_body_posvel()
					ELSE
						CALL read_snapshot_body()
					ENDIF
				ELSE
					!complete REWIND is advisable.
					REWIND 9
					!in principle, one could now set the file_position flag to zero and let the subroutine call itself recursively.
					!However, I prefer iteration:
					!Skip over the unnecessary (including body) lines, i.e. (header+natoms)*(steps to skip)
					DO counter=1,lines_to_skip*(timestep-1),1!'file_position' is zero here due to REWIND
						READ(9,*)
					ENDDO
					IF (EXTRA_VELOCITY) THEN
						CALL read_snapshot_body_posvel()
					ELSE
						CALL read_snapshot_body()
					ENDIF
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

				SUBROUTINE read_snapshot_body_posvel()
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
								&molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates,&
								&molecule_list(molecule_type_index)%snapshot_2(atom_index,molecule_index)%coordinates
							ENDDO
						ENDDO	
					ENDDO
					!And now, the file is positioned just after the timestep that has been read.
				END SUBROUTINE read_snapshot_body_posvel

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

		!returns lower and upper boundary in the given dimension
		FUNCTION give_box_boundaries(dimension_in)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: dimension_in
		REAL(KIND=WORKING_PRECISION) :: give_box_boundaries(2)
			IF (BOX_VOLUME_GIVEN) THEN
				give_box_boundaries(:)=box_dimensions(:,dimension_in)!low and high for x,y,z
			ELSE
				give_box_boundaries(:)=0.0d0
				CALL report_error(41)
			ENDIF
		END FUNCTION give_box_boundaries

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

		!This subroutine reports the molecule_type_index and atom_index values for the atoms in "give_number_of_specific_atoms"
		SUBROUTINE give_indices_of_specific_atoms(element_name_input,indices_array)
		IMPLICIT NONE
		CHARACTER(LEN=*),INTENT(IN) :: element_name_input
		INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) :: indices_array
		INTEGER :: atom_index,counter,molecule_type_index
			IF (SIZE(indices_array(:,2))/=give_number_of_specific_atoms(element_name_input))&
			&CALL report_error(138)
			counter=0
			DO molecule_type_index=1,number_of_molecule_types,1
				DO atom_index=1,give_number_of_atoms_per_molecule(molecule_type_index),1
					IF (TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))==TRIM(element_name_input)) THEN
						counter=counter+1
						indices_array(counter,2)=atom_index
						indices_array(counter,1)=molecule_type_index
					ENDIF
				ENDDO
			ENDDO
		END SUBROUTINE give_indices_of_specific_atoms

		CHARACTER(LEN=2) FUNCTION give_element_symbol(molecule_type_index,atom_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: atom_index,molecule_type_index
			give_element_symbol=molecule_list(molecule_type_index)%list_of_elements(atom_index)
		END FUNCTION give_element_symbol

		!This subroutine gives out how many specific atoms with label "element_name_input" are present.
		!A "specific atom" is one defined by its combination of molecule_type_index and atom_index, ignoring molecule_index.
		INTEGER FUNCTION give_number_of_specific_atoms(element_name_input)
		IMPLICIT NONE
		CHARACTER(LEN=*),INTENT(IN) :: element_name_input
		INTEGER :: atom_index,molecule_type_index
			give_number_of_specific_atoms=0
			DO molecule_type_index=1,number_of_molecule_types,1
				DO atom_index=1,give_number_of_atoms_per_molecule(molecule_type_index),1
					IF (TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))==TRIM(element_name_input))&
					&give_number_of_specific_atoms=give_number_of_specific_atoms+1
				ENDDO
			ENDDO
		END FUNCTION give_number_of_specific_atoms

		!This subroutine gives all the atom indices in "molecule_type_index" which have the label "element_name_input"
		SUBROUTINE give_indices_of_specific_atoms_per_molecule(molecule_type_index,element_name_input,indices_array)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: molecule_type_index
		CHARACTER(LEN=*),INTENT(IN) :: element_name_input
		INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) :: indices_array
		INTEGER :: atom_index,counter
			IF (SIZE(indices_array(:,2))/=give_number_of_specific_atoms_per_molecule(molecule_type_index,element_name_input))&
			&CALL report_error(138)
			counter=0
			DO atom_index=1,give_number_of_atoms_per_molecule(molecule_type_index),1
				IF (TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))==TRIM(element_name_input)) THEN
					counter=counter+1
					indices_array(counter,2)=atom_index
				ENDIF
			ENDDO
			indices_array(:,1)=molecule_type_index
		END SUBROUTINE give_indices_of_specific_atoms_per_molecule

		!This subroutine gives out how many specific atoms with label "element_name_input" are present in the given molecule type.
		!A "specific atom" is one defined by its combination of molecule_type_index and atom_index, ignoring molecule_index.
		INTEGER FUNCTION give_number_of_specific_atoms_per_molecule(molecule_type_index,element_name_input)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: molecule_type_index
		CHARACTER(LEN=*),INTENT(IN) :: element_name_input
		INTEGER :: atom_index
			give_number_of_specific_atoms_per_molecule=0
			DO atom_index=1,give_number_of_atoms_per_molecule(molecule_type_index),1
				IF (TRIM(molecule_list(molecule_type_index)%list_of_elements(atom_index))==TRIM(element_name_input))&
				&give_number_of_specific_atoms_per_molecule=give_number_of_specific_atoms_per_molecule+1
			ENDDO
		END FUNCTION give_number_of_specific_atoms_per_molecule

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

		REAL FUNCTION give_realcharge_of_molecule(molecule_type_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: molecule_type_index
			give_realcharge_of_molecule=molecule_list(molecule_type_index)%realcharge
		END FUNCTION give_realcharge_of_molecule

		REAL FUNCTION give_charge_of_atom(molecule_type_index,atom_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: molecule_type_index,atom_index
			give_charge_of_atom=molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
		END FUNCTION give_charge_of_atom

		REAL FUNCTION give_mass_of_atom(molecule_type_index,atom_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: molecule_type_index,atom_index
			give_mass_of_atom=molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
		END FUNCTION give_mass_of_atom

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
			IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
			IF (use_firstatom_as_com) THEN
				IF (READ_SEQUENTIAL) THEN
					give_center_of_mass(:)=DBLE(molecule_list(molecule_type_index)%snapshot(1,molecule_index)%coordinates(:))
				ELSE
					give_center_of_mass(:)=DBLE(molecule_list(molecule_type_index)%trajectory(1,molecule_index,timestep)%coordinates(:))
				ENDIF
				RETURN
			ENDIF
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
		!Alternative version using the extra velocity
		FUNCTION give_center_of_mass_2(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		REAL(KIND=WORKING_PRECISION) :: give_center_of_mass_2(3),weighted_pos(3)!higher precision, because intermediate result.
		INTEGER :: atom_index
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
			IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
			IF (use_firstatom_as_com) THEN
				IF (READ_SEQUENTIAL) THEN
					give_center_of_mass_2(:)=DBLE(molecule_list(molecule_type_index)%snapshot_2(1,molecule_index)%coordinates(:))
				ELSE
					give_center_of_mass_2(:)=DBLE(molecule_list(molecule_type_index)%trajectory_2(1,molecule_index,timestep)%coordinates(:))
				ENDIF
				RETURN
			ENDIF
			give_center_of_mass_2(:)=0.0d0
			DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
				!first, the current atom's position is stored in weighted_pos.
				!added support for sequential read.
				IF (READ_SEQUENTIAL) THEN
					weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%snapshot_2(atom_index,molecule_index)%coordinates(:))
				ELSE
					weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%trajectory_2(atom_index,molecule_index,timestep)%coordinates(:))
				ENDIF
				!then, this position is weighted with the atom's mass
				weighted_pos(:)=weighted_pos(:)*molecule_list(molecule_type_index)%list_of_atom_masses(atom_index)
				!this weighted position is now added to the center of mass.
				give_center_of_mass_2(:)=give_center_of_mass_2(:)+weighted_pos(:)
			ENDDO
			!finally, the center of mass has to be normalised by the total mass.
			give_center_of_mass_2(:)=give_center_of_mass_2(:)/DBLE(molecule_list(molecule_type_index)%mass)
		END FUNCTION give_center_of_mass_2


		!The following function gives what would be the dipole moment of a neutral molecule
		FUNCTION give_qd_vector(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		REAL(KIND=WORKING_PRECISION) :: give_qd_vector(3),weighted_pos(3)!higher precision, because intermediate result.
		INTEGER :: atom_index
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
			IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
			give_qd_vector(:)=0.0d0
			DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
				!first, the current atom's position is stored in weighted_pos.
				!added support for sequential read.
				IF (READ_SEQUENTIAL) THEN
					weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:))
				ELSE
					weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:))
				ENDIF
				!then, this position is weighted with the atom's charge
				weighted_pos(:)=weighted_pos(:)*molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
				!this weighted position is now added to the center of mass.
				give_qd_vector(:)=give_qd_vector(:)+weighted_pos(:)
			ENDDO
			!for neutral molecules, the result is already the dipole moment. But for charged ones, we'll need to normalise with charge.
		END FUNCTION give_qd_vector
		!Alternative version using the extra velocity
		FUNCTION give_qd_vector_2(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		REAL(KIND=WORKING_PRECISION) :: give_qd_vector_2(3),weighted_pos(3)!higher precision, because intermediate result.
		INTEGER :: atom_index
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
			IF ((READ_SEQUENTIAL).AND.((timestep/=file_position))) CALL goto_timestep(timestep)
			give_qd_vector_2(:)=0.0d0
			DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
				!first, the current atom's position is stored in weighted_pos.
				!added support for sequential read.
				IF (READ_SEQUENTIAL) THEN
					weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%snapshot_2(atom_index,molecule_index)%coordinates(:))
				ELSE
					weighted_pos(:)=DBLE(molecule_list(molecule_type_index)%trajectory_2(atom_index,molecule_index,timestep)%coordinates(:))
				ENDIF
				!then, this position is weighted with the atom's charge
				weighted_pos(:)=weighted_pos(:)*molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
				!this weighted position is now added to the center of mass.
				give_qd_vector_2(:)=give_qd_vector_2(:)+weighted_pos(:)
			ENDDO
			!for neutral molecules, the result is already the dipole moment. But for charged ones, we'll need to normalise with charge.
		END FUNCTION give_qd_vector_2

		!charge arm |Q*lq| as defined by kobrak in ECS Proc. Vol., 2004, 2004–24, 417–425.
		FUNCTION charge_arm_length(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
		REAL(KIND=WORKING_PRECISION) :: charge_arm_length
			charge_arm_length=SQRT(SUM((charge_arm(timestep,molecule_type_index,molecule_index))**2))
		END FUNCTION charge_arm_length

		!charge arm vector Q*lq as defined by kobrak in ECS Proc. Vol., 2004, 2004–24, 417–425.
		!if he molecule is not charged, then the dipole moment is given.
		FUNCTION charge_arm(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
		REAL(KIND=WORKING_PRECISION) :: charge_arm(3)
			IF (molecule_list(molecule_type_index)%charge==0) THEN
				charge_arm(:)=give_qd_vector(timestep,molecule_type_index,molecule_index)
			ELSE
				charge_arm(:)=give_qd_vector(timestep,molecule_type_index,molecule_index)/molecule_list(molecule_type_index)%realcharge
				charge_arm(:)=molecule_list(molecule_type_index)%realcharge*&
				&(charge_arm(:)-give_center_of_mass(timestep,molecule_type_index,molecule_index))
			ENDIF
		END FUNCTION charge_arm
		!Alternative version using the extra velocity
		FUNCTION charge_arm_2(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
		REAL(KIND=WORKING_PRECISION) :: charge_arm_2(3)
			IF (molecule_list(molecule_type_index)%charge==0) THEN
				charge_arm_2(:)=give_qd_vector_2(timestep,molecule_type_index,molecule_index)
			ELSE
				charge_arm_2(:)=give_qd_vector_2(timestep,molecule_type_index,molecule_index)/molecule_list(molecule_type_index)%realcharge
				charge_arm_2(:)=molecule_list(molecule_type_index)%realcharge*&
				&(charge_arm_2(:)-give_center_of_mass_2(timestep,molecule_type_index,molecule_index))
			ENDIF
		END FUNCTION charge_arm_2

		!gives back the center of charge.
		!if he molecule is not charged, then the dipole moment is given.
		FUNCTION give_center_of_charge(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
		REAL(KIND=WORKING_PRECISION) :: give_center_of_charge(3)
			IF (molecule_list(molecule_type_index)%charge==0) THEN
				give_center_of_charge(:)=give_qd_vector(timestep,molecule_type_index,molecule_index)
			ELSE
				give_center_of_charge(:)=give_qd_vector(timestep,molecule_type_index,molecule_index)/&
				&molecule_list(molecule_type_index)%realcharge
			ENDIF
		END FUNCTION give_center_of_charge
		!Alternative version using the extra velocity
		FUNCTION give_center_of_charge_2(timestep,molecule_type_index,molecule_index)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: timestep,molecule_type_index,molecule_index
		REAL(KIND=WORKING_PRECISION) :: give_center_of_charge_2(3)
			IF (molecule_list(molecule_type_index)%charge==0) THEN
				give_center_of_charge_2(:)=give_qd_vector_2(timestep,molecule_type_index,molecule_index)
			ELSE
				give_center_of_charge_2(:)=give_qd_vector_2(timestep,molecule_type_index,molecule_index)/&
				&molecule_list(molecule_type_index)%realcharge
			ENDIF
		END FUNCTION give_center_of_charge_2

		SUBROUTINE print_dipole_statistics()
		IMPLICIT NONE
		REAL :: maximum,minimum,average,stdev,chargearm
		INTEGER :: molecule_index,molecule_type_index
		CHARACTER(LEN=32) :: dipole_name
			IF (.NOT.((custom_atom_charges).OR.(custom_default_charges))) CALL report_error(128)
			WRITE(*,'(" Dipole Moment / Charge Arm statistics from first step in trajectory:")')
			DO molecule_type_index=1,number_of_molecule_types,1
				!initialise statistics for this molecule type
				maximum=0.0
				minimum=maximum_distance
				average=0.0
				stdev=0.0
				DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
					chargearm=SQRT(SUM(charge_arm(1,molecule_type_index,molecule_index)**2))
					IF (chargearm<minimum) minimum=chargearm
					IF (chargearm>maximum) maximum=chargearm
					average=average+chargearm
				ENDDO
				average=average/FLOAT(molecule_list(molecule_type_index)%total_molecule_count)
				DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1
					chargearm=SQRT(SUM(charge_arm(1,molecule_type_index,molecule_index)**2))
					stdev=stdev+(chargearm-average)**2
				ENDDO
				stdev=SQRT(stdev/FLOAT(molecule_list(molecule_type_index)%total_molecule_count-1))
				IF (give_charge_of_molecule(molecule_type_index)/=0) THEN
					dipole_name="charge arm (COC-COM)"
				ELSE
					dipole_name="dipole moment (sum of q*d)"
				ENDIF
				WRITE(*,'(5A,I0,A)') "   ",TRIM(dipole_name)," of Molecule ",&
				&TRIM(give_sum_formula(molecule_type_index))," (#",molecule_type_index,"):"
				IF (molecule_list(molecule_type_index)%total_molecule_count==1) THEN
					WRITE(*,'("     average =",ES10.3)') average
				ELSE
					WRITE(*,'("     average =",ES10.3)') average
					WRITE(*,'("     st.dev. =",ES10.3)') stdev
					WRITE(*,'("     minimum =",ES10.3)') minimum
					WRITE(*,'("     maximum =",ES10.3)') maximum
				ENDIF
			ENDDO
		END SUBROUTINE print_dipole_statistics

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
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: stepcounter,molecule_type_index,molecule_index,atom_index
		REAL :: box_COM(3)
			IF (READ_SEQUENTIAL) CALL report_error(0)
			IF (VERBOSE_OUTPUT) WRITE(*,*) "Removing box centre-of-mass from every step."
			!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL))) &
			!$OMP PRIVATE(box_COM,molecule_type_index,molecule_index,atom_index)
			!$OMP SINGLE
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (remove barycentre)"
		 !$ 	CALL timing_parallel_sections(.TRUE.)
		 !$ ENDIF
			IF (VERBOSE_OUTPUT) CALL print_progress(number_of_steps)
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
				!$OMP CRITICAL
				IF (VERBOSE_OUTPUT) CALL print_progress()
				!$OMP END CRITICAL
			ENDDO
			!$OMP END DO
			!$OMP END PARALLEL
			IF ((number_of_steps>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
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

		!writes the specified molecule in xyz format into the specified unit IF it is crossing a slab.
		!Unit is not overwritten if opened with APPEND!
		!the blanks are not added!
		SUBROUTINE write_molecule_in_slab(unit_number,timestep_in,molecule_type_index,molecule_index,&
		&lower_boundary,upper_boundary,slab_dimension,molecule_was_on_slab)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: unit_number,molecule_type_index,molecule_index,timestep_in,slab_dimension
		INTEGER :: atom_index,natoms,timestep
		REAL(KIND=WORKING_PRECISION),INTENT(IN) :: lower_boundary,upper_boundary
		REAL(KIND=WORKING_PRECISION) :: shift(3),atom_position(3)
		LOGICAL :: morethanlower,lessthanupper
		LOGICAL,INTENT(OUT) :: molecule_was_on_slab
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
			!must be centre-of-mass wise... otherwise it will catch atom-wise wrapped molecules...
			atom_position(:)=give_center_of_mass(timestep,molecule_type_index,molecule_index)
			CALL wrap_vector(atom_position(:),shift(:))
			natoms=molecule_list(molecule_type_index)%number_of_atoms
			!check if the molecule is actually on the slab
			morethanlower=.FALSE.
			lessthanupper=.FALSE.
			DO atom_index=1,natoms,1
				IF (READ_SEQUENTIAL) THEN
					atom_position=&
					&molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
				ELSE
					atom_position=&
					&molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
				ENDIF
				atom_position(:)=atom_position(:)+shift(:)
				IF (atom_position(slab_dimension)>lower_boundary) THEN
					!at least one atom is within the lower boundary
					morethanlower=.TRUE.
					EXIT
				ENDIF
			ENDDO
			DO atom_index=1,natoms,1
				IF (READ_SEQUENTIAL) THEN
					atom_position=&
					&molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)
				ELSE
					atom_position=&
					&molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)
				ENDIF
				atom_position(:)=atom_position(:)+shift(:)
				IF (atom_position(slab_dimension)<upper_boundary) THEN
					!at least one atom is within the upper boundary
					lessthanupper=.TRUE.
					EXIT
				ENDIF
			ENDDO
			IF ((morethanlower).AND.(lessthanupper)) THEN
				molecule_was_on_slab=.TRUE.
				!Print the molecule!
				DO atom_index=1,natoms,1
					IF (READ_SEQUENTIAL) THEN
						WRITE(unit_number,*) molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),&
						&SNGL(molecule_list(molecule_type_index)%snapshot(atom_index,molecule_index)%coordinates(:)+shift(:))
					ELSE
						WRITE(unit_number,*) molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1),&
						&SNGL(molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,timestep)%coordinates(:)+shift(:))
					ENDIF
				ENDDO
			ELSE
				molecule_was_on_slab=.FALSE.
			ENDIF
		END SUBROUTINE write_molecule_in_slab

		!Writes the trajectory to some format.
		SUBROUTINE write_trajectory(startstep_in,endstep_in,output_format)
		IMPLICIT NONE
		INTEGER :: startstep,endstep
		INTEGER :: stepcounter
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
			SELECT CASE (output_format)
			CASE ("lmp")
				WRITE(*,ADVANCE="NO",FMT='(" Writing trajectory in LAMMPS format (.",A,")...")') TRIM(ADJUSTL(output_format))
			CASE ("gro")
				CALL report_error(145)
				WRITE(*,ADVANCE="NO",FMT='(" Writing trajectory in GROMACS format (.",A,")...")') TRIM(ADJUSTL(output_format))
			CASE ("xyz")
				WRITE(*,ADVANCE="NO",FMT='(" Writing trajectory in xyz format (.",A,")...")') TRIM(ADJUSTL(output_format))
			CASE DEFAULT
				CALL report_error(144)!unknown trajectory output format, which should never be passed to this subroutine.
				RETURN
			END SELECT
			WRITE(fstring,'(A)') TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"traj."//TRIM(ADJUSTL(output_format))
			INQUIRE(UNIT=4,OPENED=connected)
			IF (connected) CALL report_error(27,exit_status=4)
			OPEN(UNIT=4,FILE=TRIM(fstring),STATUS="REPLACE")
			IF (VERBOSE_OUTPUT) THEN
				IF ((endstep-startstep)>100) WRITE(*,*)
				CALL print_progress(endstep-startstep)
			ENDIF
			DO stepcounter=startstep,endstep,1
				!First, write header
				CALL write_header(4,stepcounter,give_number_of_atoms_per_step(),output_format)
				CALL write_body(4,stepcounter,output_format)
				CALL print_progress()
			ENDDO
			ENDFILE 4
			CLOSE(UNIT=4)
			IF (VERBOSE_OUTPUT) THEN
				IF ((endstep-startstep)>100) THEN
					WRITE(*,*)
					WRITE(*,FMT='(" ")',ADVANCE="NO")
				ENDIF
			ENDIF
			WRITE(*,'("done.")')
			CLOSE(UNIT=3)

		END SUBROUTINE write_trajectory

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

		!Writes all drude particles into a trajectory file - with their respective core positions subtracted! good for power spectra...
		SUBROUTINE write_only_drudes_relative_to_core(unit_number,timestep_in,molecule_type_index,molecule_index,include_header)
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
					WRITE(unit_number,'(I0)') molecule_list(molecule_type_index)%number_of_drudes_in_molecule
					WRITE(unit_number,'("Drude positions - Core positions.")')
				ENDIF
			ENDIF
			DO atom_index=1,natoms,1
				!drude_flag is the atom index of the drude particle attached to this core.
				!Will be -1 if no drude particle is found, and 0 if this is a drude itself.
				drude_flag=molecule_list(molecule_type_index)%list_of_drude_pairs(atom_index)%drude_flag
				SELECT CASE (drude_flag)
				CASE (-1) !non-polarisable atom - skip by cycling.
					CYCLE
				CASE (0) !this is a drude particle - skip by cycling.
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
					!subtract corepos from drudepos - write to output.
					WRITE(unit_number,ADVANCE="NO",FMT='("X",I0," ")') MODULO(drude_flag-1,natoms)+1
					WRITE(unit_number,*) drudepos(:)-corepos(:)
				END SELECT
			ENDDO
		END SUBROUTINE write_only_drudes_relative_to_core

		!initialises the molecular module by reading the input file 'molecular.inp'
		SUBROUTINE initialise_molecular()
		IMPLICIT NONE
		LOGICAL :: file_exists,connected
		INTEGER :: ios,n,allocstatus,a,b,c,totalcharge,headerlines_molecular,m
		CHARACTER(LEN=16) :: inputstring
		REAL(KIND=SP) :: mass_input,charge_input
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
				!Read default charges!
				CALL read_molecular_input_file_default_charges()
				!Read atomic charges
				CALL read_molecular_input_file_atomic_charges()
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
				INTEGER :: nanions,ncations,nneutrals
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
					nanions=0
					ncations=0
					nneutrals=0
					!Iterate over all the molecule types - n is the molecule index here.
					DO n=1,number_of_molecule_types,1
						READ(3,IOSTAT=ios,FMT=*) a,b,c
						IF (ios/=0) CALL report_error(7,exit_status=ios)
						molecule_list(n)%charge=a
						molecule_list(n)%realcharge=FLOAT(a)
						molecule_list(n)%number_of_atoms=b
						IF (b<1) CALL report_error(80,exit_status=b)
						molecule_list(n)%total_molecule_count=c
						IF (c<1) CALL report_error(119,exit_status=c)
						total_number_of_atoms=total_number_of_atoms+b*c
						totalcharge=totalcharge+a*c
						!assign generic GROMACS residue name
						IF (a<0) THEN
							nanions=nanions+1
							WRITE(molecule_list(n)%residue_name,'("ANI",I0)') nanions
						ELSE
							IF (a>0) THEN
								ncations=ncations+1
								WRITE(molecule_list(n)%residue_name,'("CAT",I0)') ncations
							ELSE
								nneutrals=nneutrals+1
								WRITE(molecule_list(n)%residue_name,'("MOL",I0)') nneutrals
							ENDIF
						ENDIF
						ALLOCATE(molecule_list(n)%list_of_elements(b),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
						ALLOCATE(molecule_list(n)%list_of_atom_masses(b),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
						ALLOCATE(molecule_list(n)%list_of_atom_charges(b),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
						ALLOCATE(molecule_list(n)%manual_atom_mass_specified(b),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
						ALLOCATE(molecule_list(n)%manual_atom_charge_specified(b),STAT=allocstatus)
						IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
						IF (READ_SEQUENTIAL) THEN
							ALLOCATE(molecule_list(n)%snapshot(b,c),STAT=allocstatus)
							IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
							IF (EXTRA_VELOCITY) THEN
								ALLOCATE(molecule_list(n)%snapshot_2(b,c),STAT=allocstatus)
								IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
							ENDIF
						ELSE
							ALLOCATE(molecule_list(n)%trajectory(b,c,number_of_steps),STAT=allocstatus)
							IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
							IF (EXTRA_VELOCITY) THEN
								ALLOCATE(molecule_list(n)%trajectory_2(b,c,number_of_steps),STAT=allocstatus)
								IF (allocstatus/=0) CALL report_error(6,exit_status=allocstatus)
							ENDIF
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
						PRINT *,"All molecules have only 1 'atom' - turn on com boost."
						PRINT *,"(centre-of-mass position = position of this atom)"
					ENDIF
				END SUBROUTINE read_molecular_input_file_header

				SUBROUTINE read_molecular_input_file_default_charges()
				IMPLICIT NONE
				CHARACTER(LEN=2) :: shortstring
					custom_default_charges=.FALSE.
					WRITE(*,ADVANCE="NO",FMT='(" Searching for ",A," statement...")') "'default_charges'"
					!Skip over headerlines in molecular input file
					REWIND 3
					DO n=1,headerlines_molecular,1
						READ(3,*)
					ENDDO
					!search for charges section.
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF (ios<0) THEN
							!end of file encountered
							WRITE(*,'("done (none found, end of file encountered).")')
							EXIT
						ENDIF
						!support for synonyms
						IF (TRIM(inputstring)=="default_charges") inputstring="charges"
						IF (ios==0) THEN
							IF (TRIM(inputstring)=="charges") THEN
								WRITE(*,'("found in line ",I0,".")') n+headerlines_molecular
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,a
								IF (ios/=0) THEN
									!something went wrong
									custom_default_charges=.FALSE.
								ELSE
									!keyword ok - read the section.
									IF (a>0) THEN
										custom_default_charges=.TRUE.
										WRITE(*,FMT='(" Trying to read ",I0," custom default charges...")',ADVANCE="NO") a
										DO m=1,a,1
											READ(3,IOSTAT=ios,FMT=*) shortstring,charge_input
											IF (ios/=0) THEN
												!wrong format... abort.
												custom_default_charges=.FALSE.
												EXIT
											ELSE
												IF (ANY(ALPHABET_small==IACHAR(shortstring(1:1)))) THEN
													!lowercase letter found!
													IF (ERROR_CODE/=122) CALL report_error(122)
												ELSE
													CALL change_default_charge(shortstring,charge_input)
												ENDIF
											ENDIF
										ENDDO
										!test if custom_default_charges still true. If not, print error message.
										IF (custom_default_charges) THEN
											WRITE(*,'(" done.")') 
										ELSE
											WRITE(*,'(" failed.")') 
											CALL report_error(123)
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
				END SUBROUTINE read_molecular_input_file_default_charges

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

				SUBROUTINE change_default_charge(element_name_input,new_charge)
				IMPLICIT NONE
				REAL,INTENT(IN) :: new_charge
				REAL :: old_charge
				LOGICAL,SAVE :: print_blank=.TRUE.
				CHARACTER(LEN=*),INTENT(IN) :: element_name_input
				CHARACTER(LEN=32) :: element_name_full
					IF (print_blank) THEN
						IF (VERBOSE_OUTPUT) WRITE(*,*)
						print_blank=.FALSE.
					ENDIF
					SELECT CASE (TRIM(element_name_input))
					CASE ("H")
						old_charge=charge_hydrogen
						charge_hydrogen=new_charge
						element_name_full="Hydrogen"
					CASE ("F")
						old_charge=charge_fluorine
						charge_fluorine=new_charge
						element_name_full="Fluorine"
					CASE ("B")
						old_charge=charge_boron
						charge_boron=new_charge
						element_name_full="Boron"
					CASE ("Cl")
						old_charge=charge_chlorine
						charge_chlorine=new_charge
						element_name_full="Chlorine"
					CASE ("Br")
						old_charge=charge_bromine
						charge_bromine=new_charge
						element_name_full="Bromine"
					CASE ("I")
						old_charge=charge_iodine
						charge_iodine=new_charge
						element_name_full="Iodine"
					CASE ("N")
						old_charge=charge_nitrogen
						charge_nitrogen=new_charge
						element_name_full="Nitrogen"
					CASE ("O")
						old_charge=charge_oxygen
						charge_oxygen=new_charge
						element_name_full="Oxygen"
					CASE ("C")
						old_charge=charge_carbon
						charge_carbon=new_charge
						element_name_full="Carbon"
					CASE ("S")
						old_charge=charge_sulfur
						charge_sulfur=new_charge
						element_name_full="Sulfur"
					CASE ("P")
						old_charge=charge_phosphorus
						charge_phosphorus=new_charge
						element_name_full="Phosphorus"
					CASE ("Li")
						old_charge=charge_lithium
						charge_lithium=new_charge
						element_name_full="Lithium"
					CASE ("Na")
						old_charge=charge_sodium
						charge_sodium=new_charge
						element_name_full="Sodium"
					CASE ("D","X")
						old_charge=drude_charge
						drude_charge=new_charge
						element_name_full="drude particle"
					CASE DEFAULT
						CALL report_error(60,exit_status=IACHAR(element_name_input(1:1)))
						RETURN
					END SELECT
					IF (VERBOSE_OUTPUT) WRITE(*,'(" Changing charge of ",A," from ",F0.3," to ",F0.4,"")')&
					&TRIM(element_name_full),old_charge,new_charge
				END SUBROUTINE change_default_charge

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
					CASE ("Na")
						old_mass=mass_sodium
						mass_sodium=new_mass
						element_name_full="Sodium"
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
				LOGICAL :: printwarning
					!initialise logicals
					printwarning=.FALSE.
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
														IF (molecule_list(molecule_type_index)%manual_atom_mass_specified(atom_index)) THEN
															IF (.NOT.(printwarning)) WRITE(*,*)
															printwarning=.TRUE.
															CALL report_error(112)
														ENDIF
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
										IF (printwarning) WRITE(*,ADVANCE="NO",FMT='(" ")')
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

				SUBROUTINE read_molecular_input_file_atomic_charges()
				IMPLICIT NONE
				INTEGER :: molecule_type_index,atom_index
				LOGICAL :: printwarning
					!initialise logicals
					printwarning=.FALSE.
					DO molecule_type_index=1,number_of_molecule_types,1
						molecule_list(molecule_type_index)%manual_atom_charge_specified(:)=.FALSE.
					ENDDO
					custom_atom_charges=.FALSE.
					WRITE(*,ADVANCE="NO",FMT='(" Searching for ",A," statement...")') "'atomic_charges'"
					!Skip over headerlines in molecular input file
					REWIND 3
					DO n=1,headerlines_molecular,1
						READ(3,*)
					ENDDO
					!search for charges section.
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF (ios<0) THEN
							!end of file encountered
							WRITE(*,'("done (none found, end of file encountered).")')
							EXIT
						ENDIF
						IF (ios==0) THEN
							!support for synonyms
							IF (TRIM(inputstring)=="atom_charges") inputstring="atomic_charges"
							IF (TRIM(inputstring)=="atomic_charges") THEN
								WRITE(*,'("found in line ",I0,".")') n+headerlines_molecular
								BACKSPACE 3
								READ(3,IOSTAT=ios,FMT=*) inputstring,a
								IF (ios/=0) THEN
									!something went wrong
									CALL report_error(124)
									custom_atom_charges=.FALSE.
								ELSE
									!keyword ok - read the section.
									IF (a>0) THEN
										custom_atom_charges=.TRUE.
										WRITE(*,FMT='(" Trying to read ",I0," custom atomic charges...")',ADVANCE="NO") a
										DO m=1,a,1
											READ(3,IOSTAT=ios,FMT=*) molecule_type_index,atom_index,charge_input
											IF (ios/=0) THEN
												!wrong format... abort.
												custom_atom_charges=.FALSE.
												EXIT
											ELSE
												!check molecule_type_index and atom_index
												IF ((molecule_type_index>0).AND.(molecule_type_index<=number_of_molecule_types)) THEN
													!valid molecule type - check atom_index
													IF ((atom_index>0).AND.(atom_index<=molecule_list(molecule_type_index)%number_of_atoms)) THEN
														IF (molecule_list(molecule_type_index)%manual_atom_charge_specified(atom_index)) THEN
															IF (.NOT.(printwarning)) WRITE(*,*)
															printwarning=.TRUE.
															CALL report_error(125)
														ENDIF
														molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)=charge_input
														molecule_list(molecule_type_index)%manual_atom_charge_specified(atom_index)=.TRUE.
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
										!test if custom_atom_charges still true. If not, print error message.
										IF (printwarning) WRITE(*,ADVANCE="NO",FMT='(" ")')
										IF (custom_atom_charges) THEN
											WRITE(*,'("done.")') 
										ELSE
											WRITE(*,'("failed.")') 
											CALL report_error(124)
											DO molecule_type_index=1,number_of_molecule_types,1
												molecule_list(molecule_type_index)%manual_atom_charge_specified(:)=.FALSE.
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
				END SUBROUTINE read_molecular_input_file_atomic_charges

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
				DEALLOCATE(molecule_list(n)%list_of_atom_charges,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
				DEALLOCATE(molecule_list(n)%manual_atom_mass_specified,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
				DEALLOCATE(molecule_list(n)%manual_atom_charge_specified,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
				IF (READ_SEQUENTIAL) THEN
					DEALLOCATE(molecule_list(n)%snapshot,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
					IF (EXTRA_VELOCITY) THEN
						DEALLOCATE(molecule_list(n)%snapshot_2,STAT=deallocstatus)
						IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
					ENDIF
				ELSE
					DEALLOCATE(molecule_list(n)%trajectory,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
					IF (EXTRA_VELOCITY) THEN
						DEALLOCATE(molecule_list(n)%trajectory_2,STAT=deallocstatus)
						IF (deallocstatus/=0) CALL report_error(8,exit_status=deallocstatus)
					ENDIF
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
				fragments_initialised=.FALSE.
			ENDIF
			IF (READ_SEQUENTIAL) THEN
				IF (VERBOSE_OUTPUT) WRITE(*,*) "closing file '",TRIM(FILENAME_MOLECULAR_INPUT),"'"!no input path is added for the molecular file!
				CLOSE(UNIT=9)
			ENDIF
		END SUBROUTINE finalise_molecular

		!reports properties: Box size, density, molecule types and masses, formulae, charge...
		SUBROUTINE report_trajectory_properties()
		IMPLICIT NONE
		INTEGER :: molecule_type_index,atom_index
		REAL(KIND=GENERAL_PRECISION) :: volume,box_weight,sum_of_atom_charges
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
				IF ((custom_atom_charges).OR.(custom_default_charges)) THEN
					!there are charges present!
					sum_of_atom_charges=0.0
					DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
						sum_of_atom_charges=sum_of_atom_charges+molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
					ENDDO
					IF (DEVELOPERS_VERSION) WRITE(*,'("  ! sum of atomic charges: ",E9.3)') sum_of_atom_charges
					IF (ABS(sum_of_atom_charges-FLOAT(molecule_list(molecule_type_index)%charge))>0.001) THEN
						CALL report_error(126)
						molecule_list(molecule_type_index)%realcharge=sum_of_atom_charges
						IF (ABS(sum_of_atom_charges)>=1.0) THEN
							WRITE(*,'(" Formal charge = ",I0,", using real charge = ",F0.3," for charge arm.")')&
							&molecule_list(molecule_type_index)%charge,molecule_list(molecule_type_index)%realcharge
						ELSE
							WRITE(*,'(" Formal charge = ",I0,", using real charge = ",E10.3," for charge arm.")')&
							&molecule_list(molecule_type_index)%charge,molecule_list(molecule_type_index)%realcharge
						ENDIF
					ENDIF
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
				!WRAP_TRAJECTORY check here
				!error report 73 also sets WRAP_TRAJECTORY to .FALSE.
				IF ((INFORMATION_IN_TRAJECTORY=="VEL").AND.(WRAP_TRAJECTORY)) CALL report_error(73)
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
				REAL(KIND=SP) :: element_mass,element_charge
				INTEGER :: current_atom_index,n,skipped_masses,skipped_charges
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
						IF (DEVELOPERS_VERSION) WRITE(*,'("  ! load_trajectory_header_information is faulty")') skipped_charges
						CALL report_error(0)!unknown trajectory format, which should never be passed to this subroutine.
					END SELECT
					!define the number of lines to skip when advancing through the trajectory file
					lines_to_skip=headerlines_to_skip+total_number_of_atoms
					number_of_drude_particles=0
					!skipped_masses counts how many masses were already specified manually. Same for skipped_charges
					skipped_masses=0
					skipped_charges=0
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
							IF (molecule_list(molecule_type_index)%number_of_atoms==1) THEN
								element_charge=FLOAT(molecule_list(molecule_type_index)%charge)
							ELSE
								element_charge=atomic_charge(element_name)
							ENDIF
							IF (molecule_list(molecule_type_index)%manual_atom_charge_specified(atom_index)) THEN
								skipped_charges=skipped_charges+1
								element_charge=molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)
							ELSE
								molecule_list(molecule_type_index)%list_of_atom_charges(atom_index)=element_charge
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
					IF (DEVELOPERS_VERSION) WRITE(*,'("  ! skipped ",I0," previously specified atomic charges.")') skipped_charges
					IF ((skipped_masses>0).AND.(.NOT.(custom_atom_masses))) CALL report_error(0)
					IF ((skipped_charges>0).AND.(.NOT.(custom_atom_charges))) CALL report_error(0)
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
						IF (EXTRA_VELOCITY) THEN
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
										&molecule_list(molecule_type_index)%trajectory(atom_index,molecule_index,stepcounter)%coordinates,&
										&molecule_list(molecule_type_index)%trajectory_2(atom_index,molecule_index,stepcounter)%coordinates
									ENDDO
								ENDDO	
							ENDDO
						ELSE
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
						ENDIF
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
					IF ((use_scaling).AND.(step>0)) THEN
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
					WRITE(unit_number,'("Velocities:")')
				CASE("POS")
					WRITE(unit_number,'("Positions:")')
				CASE DEFAULT
					CALL report_error(0)
				END SELECT
			CASE ("gro")
				WRITE(unit_number,'("converted from prealpha - check units. t= ",I0,".0")') TIME_SCALING_FACTOR*step_number
				WRITE(unit_number,'(I0)') natoms
			CASE DEFAULT
				CALL report_error(0)!unknown trajectory output format, which should never be passed to this subroutine.
			END SELECT
		END SUBROUTINE write_header

		SUBROUTINE write_body(unit_number,step_number,output_format)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: unit_number,step_number
		CHARACTER(LEN=3),INTENT(IN) :: output_format
		CHARACTER(LEN=5) :: atom_name !The atom name adhering to GROMACS specifications
		INTEGER :: molecule_type_index,molecule_index,atom_index,atom_number,natoms,residue_number
			IF ((READ_SEQUENTIAL).AND.((step_number/=file_position))) CALL goto_timestep(step_number)
			!Write body, depending on which type the trajectory has...
			atom_number=0
			residue_number=0
			DO molecule_type_index=1,number_of_molecule_types,1
				natoms=molecule_list(molecule_type_index)%number_of_atoms
				SELECT CASE (output_format)
				CASE ("lmp","xyz")
					DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory
						DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
							WRITE(unit_number,FMT='(A)',ADVANCE="NO")&
							&TRIM(molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1))//" "
							IF (READ_SEQUENTIAL) THEN
								WRITE(unit_number,*) SNGL(molecule_list(molecule_type_index)%&
								&snapshot(atom_index,molecule_index)%coordinates(:))
							ELSE
								WRITE(unit_number,*) SNGL(molecule_list(molecule_type_index)%&
								&trajectory(atom_index,molecule_index,step_number)%coordinates(:))
							ENDIF
						ENDDO
					ENDDO
				CASE ("gro")
					DO molecule_index=1,molecule_list(molecule_type_index)%total_molecule_count,1 !gives dimension 2 of trajectory
						residue_number=residue_number+1
						DO atom_index=1,molecule_list(molecule_type_index)%number_of_atoms,1
							!First we write the columns before the coordinates. Thank you GROMACS.
							WRITE(atom_name,'(A,I0)')&
							&TRIM(molecule_list(molecule_type_index)%list_of_elements(MODULO(atom_index-1,natoms)+1)),&
							&atom_index
							atom_number=atom_number+1
							WRITE(unit_number,FMT='(I5,2A5,I5)',ADVANCE="NO")&
							&residue_number,molecule_list(molecule_type_index)%residue_name,atom_name,atom_number
							!Then the coordinates.
							IF (READ_SEQUENTIAL) THEN
								WRITE(unit_number,FMT='(3F8.3)') SNGL(molecule_list(molecule_type_index)%&
								&snapshot(atom_index,molecule_index)%coordinates(:)/10.0)
							ELSE
								WRITE(unit_number,FMT='(3F8.3)') SNGL(molecule_list(molecule_type_index)%&
								&trajectory(atom_index,molecule_index,step_number)%coordinates(:)/10.0)
							ENDIF
						ENDDO
					ENDDO
				CASE DEFAULT
					CALL report_error(0)!unknown trajectory output format, which should never be passed to this subroutine.
				END SELECT
			ENDDO
			IF (output_format=="gro") WRITE(unit_number,*) box_size(:)/10.0
		END SUBROUTINE write_body

		REAL(KIND=SP) FUNCTION atomic_weight(element_name) !this function returns the atomic weight for a given element.
		IMPLICIT NONE
		!If you change this part, then also change Module_SETTINGS and change_default_mass!
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
			CASE ("Na")
				atomic_weight=(mass_sodium)
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

		REAL(KIND=SP) FUNCTION atomic_charge(element_name) !this function returns the atomic charge for a given element.
		IMPLICIT NONE
		CHARACTER(LEN=*),INTENT(IN) :: element_name
			SELECT CASE (TRIM(element_name))
			CASE ("H")
				atomic_charge=charge_hydrogen
			CASE ("F")
				atomic_charge=charge_fluorine
			CASE ("B")
				atomic_charge=charge_boron
			CASE ("Cl")
				atomic_charge=charge_chlorine
			CASE ("Br")
				atomic_charge=charge_bromine
			CASE ("I")
				atomic_charge=charge_iodine
			CASE ("N")
				atomic_charge=charge_nitrogen
			CASE ("O")
				atomic_charge=charge_oxygen
			CASE ("C")
				atomic_charge=charge_carbon
			CASE ("S")
				atomic_charge=charge_sulfur
			CASE ("P")
				atomic_charge=charge_phosphorus
			CASE ("Li")
				atomic_charge=charge_lithium
			CASE ("Na")
				atomic_charge=charge_sodium
			CASE ("X")
				atomic_charge=drude_charge
			CASE ("D")
				atomic_charge=drude_charge
			CASE DEFAULT
				!the 'convert' keyword produces a trajectory with a,b,c,...,z as element names.
				IF (ANY(ALPHABET_small==IACHAR(element_name(1:1)))) THEN
					atomic_charge=0.0
				ELSE
					CALL report_error(4)
				ENDIF
			END SELECT
		END FUNCTION atomic_charge

END MODULE MOLECULAR
!--------------------------------------------------------------------------------------------------------------------------------!
