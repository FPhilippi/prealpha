
!This Module is capable of calculating autocorrelation functions. Currently implemented:
!the intermittent autocorrelation function for a binary operator, based on an arbitrary set of dihedral constraints.
!these constraints can be only one (e.g. for simple chain conformer analyses), two (e.g. two dihedrals plus folding for cisoid/transoid transitions), or more.
!reorientational autocorrelation functions.
!relative mean molecular velocity correlation functions.
MODULE AUTOCORRELATION ! Copyright (C) !RELEASEYEAR! Frederik Philippi
    USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	LOGICAL,PARAMETER :: fold_default=.FALSE.
	LOGICAL,PARAMETER :: dump_verbose_default=.FALSE.
	LOGICAL,PARAMETER :: skip_autocorr_default=.FALSE.
	LOGICAL,PARAMETER :: export_dihedral_default=.FALSE.
	LOGICAL,PARAMETER :: export_angle_default=.FALSE.
	LOGICAL,PARAMETER :: jump_analysis_dihedral_default=.FALSE.
	LOGICAL,PARAMETER :: use_dipole_moment_default=.FALSE.
	INTEGER,PARAMETER :: sampling_interval_default=1
	INTEGER,PARAMETER :: legendre_order_default=2
	INTEGER,PARAMETER :: maximum_number_of_jump_analyses=20
	INTEGER(KIND=GENERAL_PRECISION),PARAMETER :: tmax_default=1000
	INTEGER(KIND=GENERAL_PRECISION),PARAMETER :: bin_count_default=100
	!Variables.
	TYPE,PRIVATE :: statistics_component
		REAL :: average
		REAL :: stdev
		INTEGER :: counts
	END TYPE statistics_component
	TYPE,PRIVATE :: vc_component
        INTEGER :: reference_type=-1 ! reference molecule_type_index for velocity correlation
		INTEGER :: observed_type=-1 ! observed molecule_type_index
		LOGICAL :: self=.FALSE. ! if true, then this component is a self-correlation - the types must be identical
    END TYPE vc_component
	TYPE(vc_component),ALLOCATABLE :: vc_components(:) ! list of velocity correlation components to be calculated.
	INTEGER :: legendre_order !the order of the legendre polynomial to use, usually 2, maybe 1.
	LOGICAL,ALLOCATABLE :: autocorr_array(:,:)!first dimension: number of timesteps. second dimension: number of molecules per step.
	LOGICAL :: use_dipole_moment=use_dipole_moment_default !uses vector from centre of mass to centre of charge when true.
	LOGICAL :: fold=fold_default !when true, then on top of the values a,b,... specified in the dihedral_list, (360-b),(360-a)... will be considered, too.
	LOGICAL :: dump_verbose=dump_verbose_default!controls if additional information is dumped into separate files for not.
	LOGICAL :: export_dihedral=export_dihedral_default!Controls whether dihedrals are to be exported.
	LOGICAL :: export_angle=export_angle_default!Controls whether angles are to be exported.
	LOGICAL :: skip_autocorr=skip_autocorr_default!skips the actual autocorrelation. The preparation is still done, which is useful when only PES subsets or dihedral shares are required.
	LOGICAL :: jump_analysis_dihedral=jump_analysis_dihedral_default !Jump analysis has been requested.
	INTEGER :: molecule_type_index_b!molecule_type_indices for the second molecule, 'b'
	REAL(KIND=GENERAL_PRECISION),ALLOCATABLE :: boundaries(:,:)!first dimension: index of dihedral. second dimension: upper and lower boundary.
	INTEGER,ALLOCATABLE :: export_list(:) !the list of molecules to export.
	INTEGER,ALLOCATABLE :: jump_length_list(:) !list of jump lengths to consider
	INTEGER :: jump_length_total_number !how many entries jump_length_list has
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
		INTEGER :: maxmol,ios,allocstatus,deallocstatus,n,number_of_tip_atoms,number_of_base_atoms,inputinteger
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
			PRINT *,"This vector can be defined either between the centres of mass of two fragments (base + tip),"
			PRINT *,"or alternatively the charge arm / dipole moment vector can be used."
			PRINT *,"(For the latter, atomic charges need to be defined in some way)"
			PRINT *,"Would you like to use the charge arm (y) or the two fragments (n) to define a vector?"
			IF (user_input_logical()) THEN
				!charge arm / dipole moment
				use_dipole_moment=.TRUE.
			ELSE
				use_dipole_moment=.FALSE.
				!two fragments
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
			ENDIF
			PRINT *,"Would you like to print the orientation evolution of a particular molecule? (y/n)"
			PRINT *,"(Written into a file containing timestep, unit vector, vector length, angle and Pl[u(0)u(t)])"
			export_angle=user_input_logical()
			skip_autocorr=.FALSE.
			IF (export_angle) THEN
				IF (number_of_molecules==1) THEN
					PRINT *,"Just one molecule, which will be exported."
					export_total_number=1
				ELSE
					PRINT *,"Of how many molecules would you like to export the orientation evolution?"
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
				PRINT *,"Do you want to skip the autocorrelation? (y/n)"
				skip_autocorr=user_input_logical()
			ENDIF
			IF (.NOT.(skip_autocorr)) THEN
				PRINT *,"Every how many steps would you like to use for the time correlation function?"
				WRITE(*,'(A54,I0,A2)') " (Type '1' for full accuracy. The current default is '",sampling_interval_default,"')"
				sampling_interval=user_input_integer(1,nsteps)
				!parallelisation is possible, because autocorrelation.
				parallelisation_possible=.TRUE.
				IF (.NOT.(parallelisation_requested)) THEN! ask for parallelisation, if not yet requested.
					PRINT *,"The requested feature benefits from parallelisation. Would you like to turn on parallelisation? (y/n)"
					IF (user_input_logical()) parallelisation_requested=.TRUE.
				ENDIF
			ENDIF
			!sufficient information collected
			WRITE(*,FMT='(A)',ADVANCE="NO") " writing input file for orientation analysis..."
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
			IF (export_angle) THEN
				DO n=1,export_total_number,1
					WRITE(8,'(" export ",I0," ### print orientation evolution for molecule index ",I0," into separate file.")')&
					& export_list(n),export_list(n)
				ENDDO
				!DEALLOCATE memory again
				DEALLOCATE(export_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
				!reset variables to avoid interference
				export_angle=.FALSE.
				export_total_number=0
				WRITE(8,FMT='(" skip_autocorrelation ",L1)',ADVANCE="NO") skip_autocorr
				IF (skip_autocorr) THEN
					WRITE(8,*) " ### no autocorrelation, just orientation evolution."
				ELSE
					WRITE(8,*) " ### compute the autocorrelation function."
				ENDIF
			ENDIF
			IF (use_dipole_moment) THEN
				WRITE(8,*) "charge_arm ### use charge arm (or dipole moment for neutral molecules)"
			ELSE
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
			ENDIF
			WRITE(8,*) "quit"
			WRITE(8,*)
			WRITE(8,*) "This is an input file for the calculation of a vector reorientational time correlation function."
			WRITE(8,*) "To actually perform the implied calculations, it has to be referenced in 'general.inp'."
			ENDFILE 8
			CLOSE(UNIT=8)
			WRITE(*,*) "done"
			IF (.NOT.(use_dipole_moment)) THEN
				DEALLOCATE(fragment_list_base,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
				DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
			ENDIF
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
					IF (ERROR_CODE==121) THEN
						RETURN
						CLOSE(UNIT=3)
					ENDIF
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
						CLOSE(UNIT=3)
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
					IF (ERROR_CODE==118) THEN
						RETURN
						CLOSE(UNIT=3)
					ENDIF
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
					IF (ERROR_CODE==118) THEN
						RETURN
						CLOSE(UNIT=3)
					ENDIF
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
				CHARACTER(LEN=1024) :: inputstring
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
				CHARACTER(LEN=1024) :: inputstring
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
				CHARACTER(LEN=1024) :: inputstring
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
						CASE ("skip_autocorrelation","skip_autocorr")
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
					jump_analysis_dihedral=jump_analysis_dihedral_default
					jump_length_total_number=0
					use_dipole_moment=use_dipole_moment_default
					export_angle=export_angle_default
					export_dihedral=export_dihedral_default
				END SUBROUTINE set_defaults

				SUBROUTINE read_input_for_reorientation()
				IMPLICIT NONE
				LOGICAL :: tip_read,base_read
				INTEGER :: number_of_tip_atoms,number_of_base_atoms,counter,allocstatus,deallocstatus,n,inputinteger
				INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_base(:) !temporary list of centre-of-mass fragments (defined as atom_indices) for base atom
				INTEGER,DIMENSION(:),ALLOCATABLE :: fragment_list_tip(:) !temporary list of centre-of-mass fragments (defined as atom_indices) for tip atom
				CHARACTER(LEN=1024) :: inputstring
					tip_read=.FALSE.
					base_read=.FALSE.
					legendre_order=legendre_order_default
					!initialise variables for angle export.
					export_total_number=0
					export_angle=.FALSE.
					DO n=1,MAXITERATIONS,1
						READ(3,IOSTAT=ios,FMT=*) inputstring
						IF ((ios<0).AND.(VERBOSE_OUTPUT)) WRITE(*,*) "End-of-file condition in ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
						IF (ios/=0) THEN
							IF (VERBOSE_OUTPUT) WRITE(*,*) "Done reading ",TRIM(FILENAME_AUTOCORRELATION_INPUT)
							EXIT
						ENDIF
						SELECT CASE (TRIM(inputstring))
						CASE ("skip_autocorrelation","skip_autocorr")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,skip_autocorr
							IF (ios/=0) THEN
								CALL report_error(24,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,L1,A)') " setting 'skip_autocorr' to default (=",skip_autocorr_default,")"
								skip_autocorr=skip_autocorr_default
							ELSE
								IF (VERBOSE_OUTPUT) WRITE(*,*) "setting 'skip_autocorr' to ",skip_autocorr
							ENDIF
						CASE ("export")
							!prepare the list of molecule indices whose angles are to be exported in a separate file.
							!molecule_type_index is initialised elsewhere.
							IF (export_angle) THEN !export_angle is already active.
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
										&" Angle information for molecule number '",inputinteger,"' of type '",&
										&molecule_type_index,"' will also be exported."
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
									&" Angle information for molecule number '",inputinteger,"' of type '",&
									&molecule_type_index,"' will be exported."
									export_total_number=1
									export_angle=.TRUE.
									!first occurrence of "export". Allocate memory for list.
									ALLOCATE(export_list(give_number_of_molecules_per_step(molecule_type_index)),STAT=allocstatus)
									IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
									!add molecule_index to list
									export_list(export_total_number)=inputinteger
								ENDIF
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
									IF (use_dipole_moment) THEN
										IF (VERBOSE_OUTPUT) WRITE(*,'(" (turned off use of dipole moment.)")')
										use_dipole_moment=.FALSE.
									ENDIF
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
									IF (use_dipole_moment) THEN
										IF (VERBOSE_OUTPUT) WRITE(*,'(" (turned off use of dipole moment.)")')
										use_dipole_moment=.FALSE.
									ENDIF
								ELSE
									CALL report_error(81,exit_status=n+2)
									DEALLOCATE(fragment_list_base,STAT=deallocstatus)
									IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
								ENDIF
							ELSE
								CALL report_error(81,exit_status=n+2)
							ENDIF
						CASE ("dipole","charge_arm")
							use_dipole_moment=.TRUE.
							IF (VERBOSE_OUTPUT) WRITE(*,'(" Using dipole moment / charge arm vector for reorientation.")')
							!reset everything
							IF (base_read) THEN
								DEALLOCATE(fragment_list_base,STAT=deallocstatus)
								IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
								IF (VERBOSE_OUTPUT) WRITE(*,'(" (fragment_list_base reset)")')
							ENDIF
							IF (tip_read) THEN
								DEALLOCATE(fragment_list_tip,STAT=deallocstatus)
								IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
								IF (VERBOSE_OUTPUT) WRITE(*,'(" (fragment_list_tip reset)")')
							ENDIF
							!start with the assumption that something will go wrong:
						CASE DEFAULT
							IF (VERBOSE_OUTPUT) WRITE(*,*) "can't interpret line - continue streaming"
						END SELECT
					ENDDO
					!If necessary, check for issues with export_list
					IF (export_angle) THEN
						DO n=1,export_total_number-1,1
							IF (ANY(export_list((n+1):(export_total_number))==export_list(n))) THEN
								CALL report_error(76)
								EXIT
							ENDIF
						ENDDO
					ENDIF
					IF ((base_read).AND.(tip_read)) THEN
						!both fragment lists should be filled now, call the initialisation in MODULE MOLECULAR
						CALL initialise_fragments(fragment_list_tip,fragment_list_base,number_of_tip_atoms,number_of_base_atoms,molecule_type_index)
					ELSE
						IF (.NOT.(use_dipole_moment)) CALL report_error(83)
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
				CHARACTER(LEN=1024) :: inputstring
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
						CASE ("skip_autocorrelation","skip_autocorr")
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
				USE DEBUG
				IMPLICIT NONE
				INTEGER :: n,deallocstatus,inputinteger
				CHARACTER(LEN=1024) :: inputstring
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
						WRITE(*,'("   ",A,2F6.1)') TRIM(inputstring),boundaries(n,:)
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
										&" Dihedrals for molecule number '",inputinteger,"' of type '",molecule_type_index,"' will also be exported."
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
						CASE ("skip_autocorrelation","skip_autocorr")
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
						CASE ("jump_analysis","jump_velocity","jump")
							BACKSPACE 3
							READ(3,IOSTAT=ios,FMT=*) inputstring,inputinteger
							IF (ios/=0) THEN
								CALL report_error(24,exit_status=ios)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') " jump analysis requires the desired jump time (in number of timesteps)"
							ELSEIF (jump_length_total_number>maximum_number_of_jump_analyses) THEN
								!too many jump analyses requested
								CALL report_error(130,exit_status=maximum_number_of_jump_analyses)
							ELSE
								CALL check_timestep(inputinteger)
								IF (VERBOSE_OUTPUT) WRITE(*,'(A,I0,A)') &
								&" Will perform jump analysis for jump time of ",inputinteger," timesteps."
								IF (inputinteger<100) CALL report_error(131)
								IF (jump_analysis_dihedral) THEN !jump analysis is already active.
									jump_length_total_number=jump_length_total_number+1
									!add to list
									jump_length_list(jump_length_total_number)=inputinteger
								ELSE !first element of jump_length_list!
									jump_length_total_number=1
									jump_analysis_dihedral=.TRUE.
									!first occurrence of "jump_analysis". Allocate memory for list.
									ALLOCATE(jump_length_list(maximum_number_of_jump_analyses),STAT=allocstatus)
									IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
									!add molecule_index to list
									jump_length_list(jump_length_total_number)=inputinteger
								ENDIF
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
					!If necessary, check for issues with jump_length_list
					IF (jump_analysis_dihedral) THEN
						DO n=1,jump_length_total_number-1,1
							IF (ANY(jump_length_list((n+1):(jump_length_total_number))==jump_length_list(n))) THEN
								CALL report_error(129)
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
					export_dihedral=.FALSE.
				ENDIF
				IF (jump_analysis_dihedral) THEN
					DEALLOCATE(jump_length_list,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
					jump_analysis_dihedral=.FALSE.
				ENDIF
			ENDIF
			IF (export_angle) THEN
				DEALLOCATE(export_list,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(15,exit_status=deallocstatus)
				export_angle=.FALSE.
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
			 !$ 	FUNCTION OMP_get_num_threads()
			 !$ 	INTEGER :: OMP_get_num_threads
			 !$ 	END FUNCTION OMP_get_num_threads
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
				 !$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (autocorrelation)"
				 !$ 	CALL timing_parallel_sections(.TRUE.)
				 !$ ENDIF
					IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling)/sampling,0))
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
						!$OMP CRITICAL
						IF (VERBOSE_OUTPUT) CALL print_progress()
						!$OMP END CRITICAL
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
					IF (((MAX((nsteps-1+sampling)/sampling,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
				 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
				 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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
						IF (timeline==(tmax-1)) EXIT
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
			13	FORMAT ("   ",A5," ",E15.6)
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
			14	FORMAT ("   ",A8," ",E15.6," (",EN13.4," Scm²/mol)")!molar
			15	FORMAT ("   ",A8," ",E15.6," (",EN13.4," S/cm)")!specific
			16	FORMAT ("   ",A5," ",E15.6,A14)
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
			IF (use_dipole_moment) THEN
				IF (give_charge_of_molecule(molecule_type_index)==0) THEN
					WRITE(*,'(" Vector for reorientation is the dipole moment of ",A,".")')&
					&TRIM(give_sum_formula(molecule_type_index))
					WRITE(*,'(" (defined as SUM(q*r), q=atom charge, r=atom position)")')
				ELSE
					WRITE(*,'(" Vector for reorientation is the charge arm of ",A,".")')&
					&TRIM(give_sum_formula(molecule_type_index))
					WRITE(*,'(" ( defined as (SUM(q*r)/SUM(q))-(SUM(m*r)/SUM(m)),")')
					WRITE(*,'("   q=atom charge, r=atom position, m=atom mass )")')
				ENDIF
				IF (.NOT.(check_charges(molecule_type_index))) CALL report_error(127)
			ELSE
				!Report fragment information
				CALL give_fragment_information(tip_fragment=.FALSE.)
				CALL give_fragment_information(tip_fragment=.TRUE.)
			ENDIF
			IF (export_angle) CALL export_angles()
			IF (skip_autocorr) THEN
				IF (VERBOSE_OUTPUT) WRITE(*,*) "skip time correlation function."
			ELSE
				!Main part: compute the tcf
				CALL compute_reorientational_tcf_parallel(sampling_interval)
				!Report the tcf
				CALL report_time_correlation_function()
				!Deallocate memory for tcf and average counter
			ENDIF
			DEALLOCATE(time_correlation_function,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			DEALLOCATE(x_num,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			CONTAINS

				!This subroutine merely writes the angle information.
				SUBROUTINE export_angles()
				IMPLICIT NONE
				INTEGER :: unit_number,export_counter,stepcounter,ios
				CHARACTER(LEN=1024) :: filename_export
				CHARACTER(LEN=15) :: Pl_name
				LOGICAL :: connected
				REAL(WORKING_PRECISION) :: orientation(3),initial_orientation(3),length,u0ut_dotproduct
					IF (.NOT.(export_angle)) CALL report_error(0)
					IF (VERBOSE_OUTPUT) WRITE(*,*) "Orientation evolution is written to separate files..."
					IF (VERBOSE_OUTPUT) CALL print_progress(MAX(nsteps*export_total_number,0))
					!open files for angle export.
					!Files will be opened in unit 10,11,12,...
					DO export_counter=1,export_total_number,1
						unit_number=export_counter+9
						INQUIRE(UNIT=unit_number,OPENED=connected)
						IF (connected) CALL report_error(27,exit_status=unit_number)
						WRITE(filename_export,'("orientation_export",I0,"_",I0)') export_counter,export_list(export_counter)
						OPEN(UNIT=unit_number,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//filename_export,IOSTAT=ios)
						IF (ios/=0) CALL report_error(26,exit_status=ios)
						WRITE(unit_number,'(A,I0,A,I0)') "This file contains orientational information for the molecule number "&
						&,export_list(export_counter)," of type ",molecule_type_index
						IF (use_dipole_moment) THEN
							WRITE(unit_number,*) "'Vector' is the dipole moment (or charge arm for charged species)."
						ELSE
							WRITE(unit_number,*) "'Vector' connects base point to tip point."
						ENDIF
						WRITE(Pl_name,'("P",I0,"[u(t)u(0)]")') legendre_order
						WRITE(unit_number,FMT='(7A17)') "step","unit_vector_x","unit_vector_y","unit_vector_z",&
						&"vector_length","angle/degrees",TRIM(Pl_name)
						IF (use_dipole_moment) THEN
							initial_orientation(:)=charge_arm(1,molecule_type_index,export_list(export_counter))
						ELSE
							initial_orientation(:)=&
							&give_tip_fragment(1,export_list(export_counter))-give_base_fragment(1,export_list(export_counter))
						ENDIF
						CALL normalize3D(initial_orientation(:))
						DO stepcounter=1,nsteps,1
							WRITE(unit_number,ADVANCE="NO",FMT='(I17)') stepcounter
							IF (use_dipole_moment) THEN
								!orientation = dipole moment or "charge arm" vector
								orientation(:)=charge_arm(stepcounter,molecule_type_index,export_list(export_counter))
							ELSE
								!orientation = normalised vector from 'base atom' to 'tip atom' at the initial timestep
								orientation(:)=&
								&give_tip_fragment(stepcounter,export_list(export_counter))-&
								&give_base_fragment(stepcounter,export_list(export_counter))
							ENDIF
							length=SQRT(SUM(orientation(:)**2))
							IF (length==0.0_WORKING_PRECISION) CALL report_error(1)
							orientation=(orientation/length)
							!the dot product is equal to the cosine of the angle between the two vectors
							u0ut_dotproduct=DOT_PRODUCT(orientation(:),initial_orientation(:))
							WRITE(unit_number,*) SNGL(orientation),SNGL(length),&
							&SNGL(ACOS(u0ut_dotproduct)*degrees),& 
							&SNGL(legendre_polynomial(u0ut_dotproduct,legendre_order))!take the two normalised values and compute the legendre polynomial of the dot product.
							IF (VERBOSE_OUTPUT) CALL print_progress()
						ENDDO
						ENDFILE unit_number
						CLOSE(UNIT=unit_number)
					ENDDO
					IF (((MAX(nsteps*export_total_number,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
					IF (VERBOSE_OUTPUT) WRITE(*,*) "...done."
				END SUBROUTINE export_angles

				!This subroutine computes eq (11.11.1) from "THEORY OF SIMPLE LIQUIDS" (Hansen / McDonald)
				!'sampling' is the sampling interval.
				SUBROUTINE compute_reorientational_tcf_parallel(sampling)
				IMPLICIT NONE
			 !$ INTERFACE
			 !$ 	FUNCTION OMP_get_num_threads()
			 !$ 	INTEGER :: OMP_get_num_threads
			 !$ 	END FUNCTION OMP_get_num_threads
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
			!		!$OMP PRIVATE(iterate_timesteps_self_contributions)
					!$OMP SINGLE
				 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
				 !$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (time correlation function)"
				 !$ 	CALL timing_parallel_sections(.TRUE.)
				 !$ ENDIF
					IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling)/sampling,0))
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
							IF (use_dipole_moment) THEN
								!initial orientation = dipole moment or charge arm vector
								initial_orientation(molecule_counter,:)=charge_arm(startstep,molecule_type_index,molecule_counter)
							ELSE
								!initial orientation = normalised vector from 'base atom' to 'tip atom' at the initial timestep
								initial_orientation(molecule_counter,:)=&
								&give_tip_fragment(startstep,molecule_counter)-give_base_fragment(startstep,molecule_counter)
							ENDIF
							CALL normalize3D(initial_orientation(molecule_counter,:))
						ENDDO
						CALL iterate_timesteps_tcf(startstep,temp_function,initial_orientation,x_num_temp)
						!$OMP CRITICAL
						IF (VERBOSE_OUTPUT) CALL print_progress()
						!$OMP END CRITICAL
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
					IF (((MAX((nsteps-1+sampling)/sampling,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
				 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
				 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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
							IF (use_dipole_moment) THEN
								second_vector(:)=&
								&charge_arm(startstep+timeline,molecule_type_index,molecule_counter)
							ELSE
								second_vector(:)=&
								&give_tip_fragment(startstep+timeline,molecule_counter)-&
								&give_base_fragment(startstep+timeline,molecule_counter)
							ENDIF
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
		CHARACTER(LEN=1024) :: filename_export
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
			IF (VERBOSE_OUTPUT) WRITE(*,*) "Filling dihedral array..."
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
					IF (VERBOSE_OUTPUT) CALL print_progress(give_number_of_timesteps()-1)
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
						IF (VERBOSE_OUTPUT) CALL print_progress()
					ENDDO
					IF (((give_number_of_timesteps()-1)>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
					DEALLOCATE(dihedral_list,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				END SUBROUTINE fill_dihedral_array

				SUBROUTINE dump_PES()
				IMPLICIT NONE
				!Output formats - AUTOCORRELATION module
			3	FORMAT (2F15.3,I15) !Module AUTOCORRELATION - PES_subset_dependent
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
			5	FORMAT (2I15,F15.3) !Module AUTOCORRELATION - local_incidence
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

		SUBROUTINE calculate_jump_histograms_from_binary_array()
		IMPLICIT NONE
		INTEGER :: jump_number
		TYPE(statistics_component),DIMENSION(:),ALLOCATABLE :: jump_vs_share,jump_vs_changes
		REAL :: step_size
		!step size only required for jump_vs_share, which goes from 0 to 100%
		step_size=(1.0)/FLOAT(bin_count)
		DO jump_number=1,jump_length_total_number,1
			WRITE(*,'(" Jump analysis ",I0," out of ",I0,".")')&
			&jump_number,jump_length_total_number
			WRITE(*,'(" (jumps of ",I0," time steps, which equates to ",I0,"*",I0,"=",I0," time units)")')&
			&jump_length_list(jump_number),jump_length_list(jump_number),&
			&TIME_SCALING_FACTOR,jump_length_list(jump_number)*TIME_SCALING_FACTOR
			CALL initialise_jump_histograms(jump_length_list(jump_number))
			CALL calculate_jump_statistics_parallel(jump_length_list(jump_number))
			CALL report_jump_statistics(jump_length_list(jump_number))
			CALL finalise_jump_histograms()
		ENDDO
		CONTAINS

			SUBROUTINE calculate_jump_statistics_parallel(jump_length)
			IMPLICIT NONE
		 !$ INTERFACE
		 !$ 	FUNCTION OMP_get_num_threads()
		 !$ 	INTEGER :: OMP_get_num_threads
		 !$ 	END FUNCTION OMP_get_num_threads
		 !$ END INTERFACE
			REAL :: jump_velocity,jump_vector(3),jump_time,jump_length_real,share
			INTEGER,INTENT(IN) :: jump_length
			INTEGER :: number_of_changes,molecule_index,time_shift,out_of_bounds,allocstatus,deallocstatus,bin_position,timestep
			TYPE(statistics_component),DIMENSION(:),ALLOCATABLE :: jump_vs_share_local,jump_vs_changes_local
				out_of_bounds=0
				jump_time=FLOAT(TIME_SCALING_FACTOR*jump_length)
				jump_length_real=FLOAT(jump_length)
				!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL)))&
				!$OMP PRIVATE(jump_vs_share_local,jump_vs_changes_local,jump_velocity,number_of_changes,jump_vector,molecule_index)&
				!$OMP PRIVATE(share,time_shift,allocstatus,deallocstatus,bin_position)
				!$OMP SINGLE
			 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
			 !$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (jump analysis)"
			 !$ 	CALL timing_parallel_sections(.TRUE.)
			 !$ ENDIF
				!$OMP END SINGLE
				!initialise local variables
				ALLOCATE(jump_vs_share_local(bin_count),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				!for a jump length of n steps, there can be no more than n-1 changes!
				ALLOCATE(jump_vs_changes_local(0:jump_length-1),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				jump_vs_share_local(:)%average=0.0
				jump_vs_share_local(:)%stdev=0.0
				jump_vs_share_local(:)%counts=0
				jump_vs_changes_local(:)%average=0.0
				jump_vs_changes_local(:)%stdev=0.0
				jump_vs_changes_local(:)%counts=0
				!need to iterate over all molecule indices and possible starting timesteps.
				!$OMP DO SCHEDULE(STATIC,1)
				DO timestep=1,give_number_of_timesteps()-jump_length,1
					DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
						!calculate distance/length of jump, and from that the jump velocity
						jump_vector(:)=&
						&give_center_of_mass(timestep+jump_length,molecule_type_index,molecule_index)-&
						&give_center_of_mass(timestep,molecule_type_index,molecule_index)
						jump_velocity=SQRT(SUM(jump_vector(:)**2))/jump_time
						!now, calculate the number of changes from one state to the other for the duration of this jump.
						number_of_changes=0
						DO time_shift=0,jump_length-1,1
							IF (autocorr_array(timestep+time_shift,molecule_index)&
							&.NEQV.autocorr_array(timestep+time_shift+1,molecule_index)) number_of_changes=number_of_changes+1
						ENDDO
						!calculate the share of fulfilled criteria
						share=FLOAT(COUNT(autocorr_array(timestep:timestep+jump_length-1,molecule_index)))/jump_length_real
						!update share and number histograms
						bin_position=(INT(share/step_size)+1)
						IF ((bin_position>0).AND.(bin_position<=bin_count+1)) THEN
							IF (bin_position==bin_count+1) bin_position=bin_count
							jump_vs_share_local(bin_position)%average=jump_vs_share_local(bin_position)%average+jump_velocity
							jump_vs_share_local(bin_position)%counts=jump_vs_share_local(bin_position)%counts+1
						ELSE
							!$OMP ATOMIC
							out_of_bounds=out_of_bounds+1
						ENDIF
						jump_vs_changes_local(number_of_changes)%average=jump_vs_changes_local(number_of_changes)%average+jump_velocity
						jump_vs_changes_local(number_of_changes)%counts=jump_vs_changes_local(number_of_changes)%counts+1
					ENDDO
				ENDDO
				!$OMP END DO
				!CRITICAL directive to properly update the histograms
				!$OMP CRITICAL
				jump_vs_changes(:)%average=jump_vs_changes(:)%average+jump_vs_changes_local(:)%average
				jump_vs_changes(:)%counts=jump_vs_changes(:)%counts+jump_vs_changes_local(:)%counts
				jump_vs_share(:)%average=jump_vs_share(:)%average+jump_vs_share_local(:)%average
				jump_vs_share(:)%counts=jump_vs_share(:)%counts+jump_vs_share_local(:)%counts
				!$OMP END CRITICAL
				!The explicit synchronisation and it's position here is important - we are still in the parralel section!
				!$OMP BARRIER
				!now, calculate the final average and the standard deviations
				!$OMP SINGLE
				jump_vs_changes(:)%average=jump_vs_changes(:)%average/FLOAT(jump_vs_changes(:)%counts)
				jump_vs_share(:)%average=jump_vs_share(:)%average/FLOAT(jump_vs_share(:)%counts)
				!$OMP END SINGLE
				!Another barrier - everyone needs to wait until the averages are properly updated.
				!$OMP BARRIER
				!$OMP DO SCHEDULE(STATIC,1)
				DO timestep=1,give_number_of_timesteps()-jump_length,1
					DO molecule_index=1,give_number_of_molecules_per_step(molecule_type_index),1
						!calculate distance/length of jump, and from that the jump velocity
						jump_vector(:)=&
						&give_center_of_mass(timestep+jump_length,molecule_type_index,molecule_index)-&
						&give_center_of_mass(timestep,molecule_type_index,molecule_index)
						jump_velocity=SQRT(SUM(jump_vector(:)**2))/jump_time
						!now, calculate the number of changes from one state to the other for the duration of this jump.
						number_of_changes=0
						DO time_shift=0,jump_length-1,1
							IF (autocorr_array(timestep+time_shift,molecule_index)&
							&.NEQV.autocorr_array(timestep+time_shift+1,molecule_index)) number_of_changes=number_of_changes+1
						ENDDO
						!calculate the share of fulfilled criteria
						share=FLOAT(COUNT(autocorr_array(timestep:timestep+jump_length-1,molecule_index)))/jump_length_real
						!update share and number histograms
						bin_position=(INT(share/step_size)+1)
						IF ((bin_position>0).AND.(bin_position<=bin_count)) THEN
							jump_vs_share_local(bin_position)%stdev=jump_vs_share_local(bin_position)%stdev+&
							&(jump_velocity-jump_vs_share(bin_position)%average)**2
						ENDIF
						jump_vs_changes_local(number_of_changes)%stdev=jump_vs_changes_local(number_of_changes)%stdev+&
						&(jump_velocity-jump_vs_changes(number_of_changes)%average)**2
					ENDDO
				ENDDO
				!$OMP END DO
				!CRITICAL directive to properly update the standard deviation histograms
				!$OMP CRITICAL
				jump_vs_changes(:)%stdev=jump_vs_changes(:)%stdev+jump_vs_changes_local(:)%stdev
				jump_vs_share(:)%stdev=jump_vs_share(:)%stdev+jump_vs_share_local(:)%stdev
				!$OMP END CRITICAL
				!deallocate local variables
				DEALLOCATE(jump_vs_share_local,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(jump_vs_changes_local,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				!$OMP END PARALLEL
			 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
			 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
			 !$ 	CALL timing_parallel_sections(.FALSE.)
			 !$ ENDIF
				IF ((VERBOSE_OUTPUT).AND.(out_of_bounds>0)) WRITE (*,'(" ",I0," entries exceeded the array boundary.")') out_of_bounds
				!update the standard deviation - so far, we only have the sum of (x-xaverage)**2
				jump_vs_changes(:)%stdev=SQRT(jump_vs_changes(:)%stdev/FLOAT(jump_vs_changes(:)%counts-1))
				jump_vs_share(:)%stdev=SQRT(jump_vs_share(:)%stdev/FLOAT(jump_vs_share(:)%counts-1))
			END SUBROUTINE calculate_jump_statistics_parallel

			SUBROUTINE report_jump_statistics(jump_length)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: jump_length
			INTEGER :: bin_counter,number_counter,first,last,ios
			LOGICAL :: connected
			REAL :: share
			CHARACTER(LEN=16) :: jump_length_string
				WRITE(jump_length_string,'(I0,"fs")') jump_length*TIME_SCALING_FACTOR
				IF (VERBOSE_OUTPUT) WRITE(*,*) "writing jump velocity vs. number of changes into file '",&
				&TRIM(ADJUSTL(OUTPUT_PREFIX))//"jump_"//TRIM(jump_length_string)//"_number_of_changes","'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"jump_"//TRIM(jump_length_string)//"_number_of_changes",&
				&IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(3,*) "This file contains the jump velocity data for a jump length of "&
				&//TRIM(jump_length_string)//" time steps vs. number of changes of dihedral condition."
				WRITE(3,'(4A15)') "#changes","average","stdev","#counts"
				!find first element to print
				first=0
				last=jump_length-1
				DO bin_counter=0,jump_length-1,1
					IF (jump_vs_changes(bin_counter)%counts>0) THEN
						first=bin_counter
						EXIT
					ENDIF
				ENDDO
				!find last element to print
				DO bin_counter=jump_length-1,0,-1
					IF (jump_vs_changes(bin_counter)%counts>0) THEN
						last=bin_counter
						EXIT
					ENDIF
				ENDDO
				DO bin_counter=first,last,1
					WRITE(3,'(I15,2E15.6,I15)') bin_counter,jump_vs_changes(bin_counter)
				ENDDO
				ENDFILE 3
				CLOSE(UNIT=3)
				IF (VERBOSE_OUTPUT) WRITE(*,*) "writing jump velocity vs. share of fulfilled condition into file '",&
				&TRIM(ADJUSTL(OUTPUT_PREFIX))//"jump_"//TRIM(jump_length_string)//"_share","'"
				INQUIRE(UNIT=3,OPENED=connected)
				IF (connected) CALL report_error(27,exit_status=3)
				OPEN(UNIT=3,FILE=TRIM(PATH_OUTPUT)//TRIM(ADJUSTL(OUTPUT_PREFIX))//"jump_"//TRIM(jump_length_string)//"_share",IOSTAT=ios)
				IF (ios/=0) CALL report_error(26,exit_status=ios)
				WRITE(3,*) "This file contains the jump velocity data for a jump length of "&
				&//TRIM(jump_length_string)//" time steps vs. share of fulfilled dihedral condition."
				WRITE(3,'(4A15)') "share","average","stdev","#counts"
				DO bin_counter=1,bin_count,1
					share=(FLOAT(bin_counter)*step_size)-(step_size*0.5) !middle of the range
					WRITE(3,'(F15.6,2E15.6,I15)') share,jump_vs_share(bin_counter)
				ENDDO
				ENDFILE 3
				CLOSE(UNIT=3)
			END SUBROUTINE report_jump_statistics

			SUBROUTINE initialise_jump_histograms(jump_length)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: jump_length
			INTEGER :: allocstatus
				ALLOCATE(jump_vs_share(bin_count),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				!for a jump length of n steps, there can be no more than n-1 changes!
				ALLOCATE(jump_vs_changes(0:jump_length-1),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				jump_vs_share(:)%average=0.0
				jump_vs_share(:)%stdev=0.0
				jump_vs_share(:)%counts=0
				jump_vs_changes(:)%average=0.0
				jump_vs_changes(:)%stdev=0.0
				jump_vs_changes(:)%counts=0
			END SUBROUTINE initialise_jump_histograms

			SUBROUTINE finalise_jump_histograms()
			IMPLICIT NONE
			INTEGER :: deallocstatus
				DEALLOCATE(jump_vs_share,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
				DEALLOCATE(jump_vs_changes,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			END SUBROUTINE finalise_jump_histograms

		END SUBROUTINE calculate_jump_histograms_from_binary_array

		SUBROUTINE calculate_autocorrelation_function_from_binary_array()!This subroutine converts the binary autocorr_array to the correlation function.
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
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
		 !$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (intermittent autocorrelation function)"
		 !$ 	CALL timing_parallel_sections(.TRUE.)
		 !$ ENDIF
			IF (VERBOSE_OUTPUT) CALL print_progress(give_number_of_molecules_per_step(molecule_type_index))
			!$OMP END SINGLE
			!allocate memory and initialise temp_function for every member of the team.
			ALLOCATE(temp_function(tmax+1),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			temp_function(:)=0.0_DP
			!iterate over molecules and pass chunks to the subroutine that iterates over timesteps
			!$OMP DO SCHEDULE(STATIC,1)
			DO n=1,give_number_of_molecules_per_step(molecule_type_index),1
				temp_function(:)=temp_function(:)+iterate_timesteps(autocorr_array(:,n))
				!$OMP CRITICAL
				IF (VERBOSE_OUTPUT) CALL print_progress()
				!$OMP END CRITICAL
			ENDDO
			!$OMP END DO
			!CRITICAL directive to properly update the autocorrelation_function
			!$OMP CRITICAL
			autocorrelation_function(:)=autocorrelation_function(:)+temp_function(:)
			!$OMP END CRITICAL
			DEALLOCATE(temp_function,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			!$OMP END PARALLEL
			IF (((give_number_of_molecules_per_step(molecule_type_index))>100)&
			&.AND.(VERBOSE_OUTPUT)) WRITE(*,*)
		 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
		 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
		 !$ 	CALL timing_parallel_sections(.FALSE.)
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
			1	FORMAT (I15,3F15.10) !Module AUTOCORRELATION - autocorrelation_function file
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
			 !$ 	FUNCTION OMP_get_num_threads()
			 !$ 	INTEGER :: OMP_get_num_threads
			 !$ 	END FUNCTION OMP_get_num_threads
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
				 !$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (overall ecaf for conductivity)"
				 !$ 	CALL timing_parallel_sections(.TRUE.)
				 !$ ENDIF
					IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling)/sampling,0))
					!$OMP END SINGLE
					!allocate memory for temporary functions (used for parallelisation)
					ALLOCATE(temp_function(tmax+1),STAT=allocstatus)
					IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
					temp_function(:)=0.0d0
					!Here, the outer loop is chosen to be the starting step, i.e. 't0'
					!$OMP DO SCHEDULE(STATIC,1)
					DO startstep=1,nsteps,sampling
						CALL iterate_timesteps_microscopic_charge_current_tacf(startstep,temp_function)
						!$OMP CRITICAL
						IF (VERBOSE_OUTPUT) CALL print_progress()
						!$OMP END CRITICAL
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
					IF (((MAX((nsteps-1+sampling)/sampling,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
				 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
				 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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
			 !$ 	FUNCTION OMP_get_num_threads()
			 !$ 	INTEGER :: OMP_get_num_threads
			 !$ 	END FUNCTION OMP_get_num_threads
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
				 !$ 	WRITE (*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (velocity correlation function)"
				 !$ 	CALL timing_parallel_sections(.TRUE.)
				 !$ ENDIF
					IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling)/sampling,0))
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
						!$OMP CRITICAL
						IF (VERBOSE_OUTPUT) CALL print_progress()
						!$OMP END CRITICAL
					ENDDO
					!$OMP END DO
					!deallocate private memory used for parallelisation
					DEALLOCATE(temp_function,STAT=deallocstatus)
					IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
					!$OMP END PARALLEL
					IF (((MAX((nsteps-1+sampling)/sampling,0))>100).AND.(VERBOSE_OUTPUT)) WRITE(*,*)
				 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
				 !$ 	WRITE(*,ADVANCE="NO",FMT='(" ### End of parallelised section, took ")')
				 !$ 	CALL timing_parallel_sections(.FALSE.)
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
				IF (jump_analysis_dihedral) CALL calculate_jump_histograms_from_binary_array()
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
