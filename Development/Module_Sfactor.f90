
!This Module calculates properties that require the microscopic density fourier components.
MODULE SFACTOR ! Copyright (C) !RELEASEYEAR! Frederik Philippi
	USE SETTINGS
	USE MOLECULAR
	IMPLICIT NONE
	PRIVATE
	!default values
	INTEGER,PARAMETER :: nsteps_default=4000 !how many steps to use from trajectory
	INTEGER,PARAMETER :: sampling_interval_default=1
	REAL,PARAMETER :: maximum_q_default=3.0d0!max q in reciprocal angströms
	REAL,PARAMETER :: resolution_in_q_default=0.05d0!target resolution for smart binning
	INTEGER,PARAMETER :: tmax_default=10000 !needed for intermediate scattering function
	LOGICAL,PARAMETER :: use_logarithmic_spacing_default=.FALSE.!Keep for later - intermediate scattering function
	!variables
	INTEGER :: nsteps=nsteps_default !how many steps to use from trajectory
	INTEGER :: sampling_interval=sampling_interval_default
	REAL :: maximum_q=maximum_q_default
	REAL :: resolution_in_q=resolution_in_q_default
	INTEGER :: number_of_bins_in_q !this should be allocated dynamically!
	INTEGER :: tmax=tmax_default
	INTEGER :: number_of_wavevectors!how many explicit combinations of nx,ny,nz will be calculated
	INTEGER :: number_of_sites!how many different sites there are. for the simple structure factor, this equals the number of element types.
	LOGICAL :: use_logarithmic_spacing=use_logarithmic_spacing_default

	REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: total_structure_factor !S(q), allocated up to number_of_bins_in_q
	REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: q_modulus !|q|, allocated from 0 up to number_of_bins_in_q

	TYPE,PRIVATE :: wavevector
		!Here we keep track of nx,ny,nz for each wavevector, the vector itself, and which entry in the q bins it belongs to
		INTEGER :: n(3) !nx,ny,nz
		REAL(KIND=WORKING_PRECISION) :: q(3) !the wavevector itself. should be (2Pi/L)*(nx ny nz)
		REAL(KIND=WORKING_PRECISION) :: q_modulus !the modulus of the wavevector.
		!important: the first entry should always be (1 0 0), (0 1 0) and (0 0 1) only, i.e. the shortest possible q.
	END TYPE wavevector

	TYPE,PRIVATE :: site
		INTEGER :: N_site_members !extent of first dimension of atoms.
		INTEGER :: N_atoms !total number of atoms in the box that are part of this site
		INTEGER,DIMENSION(:,:),ALLOCATABLE :: atoms !first dimension: counts from 1 to N_site_members; second dimension: 1=molecule_type_index, 2=atom_index. 
		REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: form_factor(:) !the atomic form factors for this site, if required.
	END TYPE site

	TYPE(wavevector),DIMENSION(:),ALLOCATABLE :: wavevectors !allocated up to number_of_wavevectors
	TYPE(site),DIMENSION(:),ALLOCATABLE :: sites !allocated up to number_of_sites
	!PRIVATE/PUBLIC declarations
	PUBLIC :: perform_sfactor_analysis

	CONTAINS

		SUBROUTINE set_defaults()!setting defaults, so that there are no bad surprises between subsequent calls.
		IMPLICIT NONE
			nsteps=nsteps_default
			sampling_interval=sampling_interval_default
			maximum_q=maximum_q_default
			tmax=tmax_default
			use_logarithmic_spacing=use_logarithmic_spacing_default
		END SUBROUTINE set_defaults

		SUBROUTINE choose_wavevectors()!dynamic bins
		IMPLICIT NONE
		INTEGER :: nx,ny,nz,n_max(3)
		INTEGER :: m,n,allocstatus,deallocstatus,current_index
		INTEGER :: number_of_wavevectors_in_bin,current_qbin,lastindex,startindex
		REAL :: current_qstart
		REAL(KIND=WORKING_PRECISION) :: q(3),prefactor(3),L,current_q_modulus,average_q
		REAL(KIND=WORKING_PRECISION) :: box_boundary(2)
		LOGICAL :: isfound
			!first, find the maximum n that make sense at all.
			DO m=1,3 !loop over the three dimensions
				box_boundary=give_box_boundaries(m)
				L=ABS(box_boundary(2)-box_boundary(1))
				prefactor(m)=twopi/L
				n=0
	inner1:		DO
					n=n+1
					IF (prefactor(m)*DSQRT(DFLOAT(n**2))>maximum_q) THEN
						n_max(m)=n-1
						EXIT inner1
					ENDIF
				ENDDO inner1
				PRINT *,"Box length L_",m,"=",L," with n_max=",n_max(m)
			ENDDO
			!count the valid wavevectors
			number_of_wavevectors=0
			DO nx=0,n_max(1) !loop over n_x
				q(1)=prefactor(1)*DFLOAT(nx)
				DO ny=0,n_max(2) !loop over n_y
					q(2)=prefactor(2)*DFLOAT(ny)
					DO nz=0,n_max(3) !loop over n_z
						q(3)=prefactor(3)*DFLOAT(nz)
						IF (nx+ny+nz>0) THEN
							IF (SUM(q(:)**2)<(maximum_q)**2) THEN
								number_of_wavevectors=number_of_wavevectors+1
							ENDIF
						ENDIF
					ENDDO
				ENDDO
			ENDDO
			PRINT *,"found a total of ",number_of_wavevectors," valid wavevectors."
			!lets allocate memory for the wavevectors and fill it
			current_index=0
			IF (ALLOCATED(wavevectors)) CALL report_error(0)
			ALLOCATE(wavevectors(number_of_wavevectors),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
xloop:		DO nx=0,n_max(1) !loop over n_x
				q(1)=prefactor(1)*DFLOAT(nx)
				DO ny=0,n_max(2) !loop over n_y
					q(2)=prefactor(2)*DFLOAT(ny)
					DO nz=0,n_max(3) !loop over n_z
						q(3)=prefactor(3)*DFLOAT(nz)
						IF (nx+ny+nz>0) THEN
							current_q_modulus=DSQRT(SUM(q(:)**2))
							IF (current_q_modulus<maximum_q) THEN
								current_index=current_index+1
								IF (current_index>number_of_wavevectors) EXIT xloop
								wavevectors(current_index)%n(1)=nx
								wavevectors(current_index)%n(2)=ny
								wavevectors(current_index)%n(3)=nz
								wavevectors(current_index)%q(:)=q(:)
								wavevectors(current_index)%q_modulus=current_q_modulus
							ENDIF
						ENDIF
					ENDDO
				ENDDO
			ENDDO xloop
			IF (current_index<number_of_wavevectors) THEN
				number_of_wavevectors=current_index
			ENDIF
			PRINT *,"wavevectors has ",number_of_wavevectors," entries"
			!make them nice
			PRINT *,"sorting wavevectors by ascending modulus"
			CALL sort_wavevectors(1,number_of_wavevectors,wavevectors)
			PRINT *,"done sorting. lowest 3 wavevectors:"
			PRINT *,wavevectors(1:3)%q_modulus
			PRINT *,"setting up equitable binning."
			number_of_bins_in_q=(maximum_q/resolution_in_q)+1
			!thus, the bin centres are:
			IF (ALLOCATED(q_modulus)) CALL report_error(0)
			ALLOCATE(q_modulus(0:number_of_bins_in_q),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			DO n=0,number_of_bins_in_q
				q_modulus(n)=DFLOAT(n)*resolution_in_q
			ENDDO
			PRINT *,"number_of_bins_in_q=",number_of_bins_in_q
			PRINT *,"Using number_of_wavevectors=",number_of_wavevectors
			! !Now, finally, we need to calculate the q averages again - because we could only allocate the memory now.
			! N_averages_perqbin(:)=0
			! q_modulus(:)=0
			! DO current_index=1,number_of_wavevectors
				! current_qbin=wavevectors(current_index)%q_entry
				! N_averages_perqbin(current_qbin)=N_averages_perqbin(current_qbin)+1
				! q_modulus(current_qbin)=q_modulus(current_qbin)+wavevectors(current_index)%q_modulus
			! ENDDO
			! q_modulus(:)=q_modulus(:)/FLOAT(N_averages_perqbin(:))

		CONTAINS

			!The following subroutine is a quicksort algorithm to sort the donor_vs_species_list.
			RECURSIVE SUBROUTINE sort_wavevectors(left,right,wavevectors_inout)
			IMPLICIT NONE
			INTEGER,INTENT(IN) :: left,right
			INTEGER :: a,b
			TYPE(wavevector),DIMENSION(*),INTENT(INOUT) :: wavevectors_inout
			REAL(KIND=WORKING_PRECISION) :: pivot
			TYPE(wavevector) :: wavevectors_clip	
				IF (left<right) THEN
					pivot=wavevectors_inout(left)%q_modulus
					a=left
					b=right
					DO
						DO WHILE (wavevectors_inout(a)%q_modulus<pivot)
							a=a+1
						ENDDO
						DO WHILE (wavevectors_inout(b)%q_modulus>pivot)
							b=b-1
						ENDDO
						IF (a>=b) EXIT
						!swap elements, unless they're the same
						IF (wavevectors_inout(a)%q_modulus==wavevectors_inout(b)%q_modulus) THEN
							b=b-1
							IF (a==b) EXIT
						ELSE
							wavevectors_clip=wavevectors_inout(a)
							wavevectors_inout(a)=wavevectors_inout(b)
							wavevectors_inout(b)=wavevectors_clip
						ENDIF
					ENDDO
					CALL sort_wavevectors(left,b,wavevectors_inout)
					CALL sort_wavevectors(b+1,right,wavevectors_inout)
				ENDIF
			END SUBROUTINE sort_wavevectors

		END SUBROUTINE choose_wavevectors

		SUBROUTINE initialise_total_xray_sfactor()
		IMPLICIT NONE
		INTEGER :: allocstatus,site_counter,site_member_counter,current_qbin
		CHARACTER(LEN=2),DIMENSION(64) :: element_list
			IF (ALLOCATED(total_structure_factor)) CALL report_error(0)
			ALLOCATE(total_structure_factor(0:number_of_bins_in_q),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			!How many sites are there?
			number_of_sites=give_elements_in_simulation_box(element_list)
			WRITE (*,ADVANCE="NO",FMT='("There are ",I0," sites (")') &
			&number_of_sites
			DO site_counter=1,number_of_sites
				WRITE (*,ADVANCE="NO",FMT='(A)') TRIM(element_list(site_counter))
				IF (site_counter<number_of_sites) WRITE (*,ADVANCE="NO",FMT='(",")')
			ENDDO
			WRITE (*,'(")")')
			IF (ALLOCATED(sites)) CALL report_error(0)
			ALLOCATE(sites(number_of_sites),STAT=allocstatus)
			IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
			DO site_counter=1,number_of_sites
				sites(site_counter)%N_site_members=give_number_of_specific_atoms(TRIM(element_list(site_counter)))
				IF (ALLOCATED(sites(site_counter)%atoms)) CALL report_error(0)
				ALLOCATE(sites(site_counter)%atoms(sites(site_counter)%N_site_members,2),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				CALL give_indices_of_specific_atoms(TRIM(element_list(site_counter)),sites(site_counter)%atoms)
				sites(site_counter)%N_atoms=0
				DO site_member_counter=1,sites(site_counter)%N_site_members,1
					sites(site_counter)%N_atoms=sites(site_counter)%N_atoms+&
					&give_number_of_molecules_per_step(sites(site_counter)%atoms(site_member_counter,1))
				ENDDO
				PRINT *,"Site "//TRIM(element_list(site_counter))//" with N=",sites(site_counter)%N_atoms
				IF (ALLOCATED(sites(site_counter)%form_factor)) CALL report_error(0)
				ALLOCATE(sites(site_counter)%form_factor(0:number_of_bins_in_q),STAT=allocstatus)
				IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
				DO current_qbin=0,number_of_bins_in_q
					sites(site_counter)%form_factor(current_qbin)=&
					&formfactor(TRIM(element_list(site_counter)),q_modulus(current_qbin))
				ENDDO
			ENDDO
		END SUBROUTINE initialise_total_xray_sfactor

		!finalises the structure factor module.
		SUBROUTINE finalise_sfactor()
		IMPLICIT NONE
		INTEGER :: deallocstatus
			IF (ALLOCATED(wavevectors)) THEN
				DEALLOCATE(wavevectors,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(total_structure_factor)) THEN
				DEALLOCATE(total_structure_factor,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(sites)) THEN
				DEALLOCATE(sites,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
			IF (ALLOCATED(q_modulus)) THEN
				DEALLOCATE(q_modulus,STAT=deallocstatus)
				IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
			ELSE
				CALL report_error(0)
			ENDIF
		END SUBROUTINE finalise_sfactor

		SUBROUTINE trajectory_total_sfactor_analysis_parallel()
		IMPLICIT NONE
	 !$ INTERFACE
	 !$ 	FUNCTION OMP_get_num_threads()
	 !$ 	INTEGER :: OMP_get_num_threads
	 !$ 	END FUNCTION OMP_get_num_threads
	 !$ END INTERFACE
		INTEGER :: stepcounter,allocstatus,deallocstatus,site_counter,site_member_counter,molecule_index,current_index,qbin1,qbin2
		REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: fourier_component_re,fourier_component_im
		REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: fourier_component_SITE_re,fourier_component_SITE_im
		REAL(KIND=WORKING_PRECISION),DIMENSION(:),ALLOCATABLE :: N_averages_perqbin
		REAL(KIND=WORKING_PRECISION) :: normalisation(0:number_of_bins_in_q),reclip,imclip
		REAL(KIND=WORKING_PRECISION) :: r(3),kr,share1,share2
		total_structure_factor(:)=0.0d0
		!$OMP PARALLEL IF((PARALLEL_OPERATION).AND.(.NOT.(READ_SEQUENTIAL)))&
		!$OMP PRIVATE(site_counter,site_member_counter,molecule_index,current_index,r,kr,qbin1,qbin2,N_averages_perqbin)&
		!$OMP PRIVATE(fourier_component_re,fourier_component_im,fourier_component_SITE_re,fourier_component_SITE_im)&
		!$OMP PRIVATE(share1,share2,reclip,imclip)
		!$OMP SINGLE
	 !$ IF ((VERBOSE_OUTPUT).AND.(PARALLEL_OPERATION)) THEN
	 !$ 	WRITE(*,'(A,I0,A)') " ### Parallel execution on ",OMP_get_num_threads()," threads (structure factor)"
	 !$ 	CALL timing_parallel_sections(.TRUE.)
	 !$ ENDIF
		IF (VERBOSE_OUTPUT) CALL print_progress(MAX((nsteps-1+sampling_interval)/sampling_interval,0))
		!$OMP END SINGLE
		!allocate memory for microscopic density fourier components
		IF (ALLOCATED(fourier_component_re)) CALL report_error(0)
		ALLOCATE(fourier_component_re(0:number_of_bins_in_q),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(179,exit_status=allocstatus)
		IF (ALLOCATED(fourier_component_im)) CALL report_error(0)
		ALLOCATE(fourier_component_im(0:number_of_bins_in_q),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(179,exit_status=allocstatus)
		IF (ALLOCATED(fourier_component_SITE_re)) CALL report_error(0)
		ALLOCATE(fourier_component_SITE_re(0:number_of_bins_in_q),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(179,exit_status=allocstatus)
		IF (ALLOCATED(fourier_component_SITE_im)) CALL report_error(0)
		ALLOCATE(fourier_component_SITE_im(0:number_of_bins_in_q),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(179,exit_status=allocstatus)
		IF (ALLOCATED(N_averages_perqbin)) CALL report_error(0)
		ALLOCATE(N_averages_perqbin(0:number_of_bins_in_q),STAT=allocstatus)
		IF (allocstatus/=0) CALL report_error(22,exit_status=allocstatus)
		!average over timesteps, a few should be enough.
		!$OMP DO SCHEDULE(STATIC,1)
		DO stepcounter=1,nsteps,sampling_interval
			fourier_component_re(:)=0.0d0
			fourier_component_im(:)=0.0d0
			N_averages_perqbin(:)=0.0d0
			DO site_counter=1,number_of_sites
				fourier_component_SITE_re(:)=0.0d0
				fourier_component_SITE_im(:)=0.0d0
				DO current_index=1,number_of_wavevectors
					!need to calculate the averages and the bin position
					qbin1=wavevectors(current_index)%q_modulus/resolution_in_q
					qbin2=qbin1+1
					share2=(wavevectors(current_index)%q_modulus-q_modulus(qbin1))/resolution_in_q
					share1=(q_modulus(qbin2)-wavevectors(current_index)%q_modulus)/resolution_in_q
					N_averages_perqbin(qbin1)=N_averages_perqbin(qbin1)+share1
					N_averages_perqbin(qbin2)=N_averages_perqbin(qbin2)+share2
					DO site_member_counter=1,sites(site_counter)%N_site_members,1
						DO molecule_index=1,give_number_of_molecules_per_step(sites(site_counter)%atoms(site_member_counter,1))
							!Here we have a specific atom. We now need to calculate exp(-iqr) for all the qs...
							!Maybe I should think about using FFT
							!get the position of this atom
							r=give_atom_position(stepcounter,sites(site_counter)%atoms(site_member_counter,1),&
							&molecule_index,sites(site_counter)%atoms(site_member_counter,2))
							!iterate over the wavevectors, sort into q bins.
							kr=DOT_PRODUCT(wavevectors(current_index)%q(:),r(:))
							!They should have fixed complex numbers in the 90s but I won't use them anyways
							!real part
							reclip=DCOS(kr)
							fourier_component_SITE_re(qbin1)=fourier_component_SITE_re(qbin1)+share1*reclip
							fourier_component_SITE_re(qbin2)=fourier_component_SITE_re(qbin2)+share2*reclip
							!imaginary part
							imclip=DSIN(kr)
							fourier_component_SITE_im(qbin1)=fourier_component_SITE_im(qbin1)+share1*imclip
							fourier_component_SITE_im(qbin2)=fourier_component_SITE_im(qbin2)+share2*imclip
						ENDDO
					ENDDO
				ENDDO
				!We are now done with this site. Thus, we need to multiply with the form factors and collect the result
				fourier_component_re(:)=fourier_component_re(:)+fourier_component_SITE_re(:)*sites(site_counter)%form_factor(:)
				fourier_component_im(:)=fourier_component_im(:)+fourier_component_SITE_im(:)*sites(site_counter)%form_factor(:)
			ENDDO
			!now we have the full fourier component for this step!
			!I abuse fourier_component_re to store rho(k)*rho(-k)=Re(rho(k))^2+Im(rho(k))^2...
			WHERE (N_averages_perqbin>0.001d0)
				fourier_component_re=(fourier_component_re**2+fourier_component_im**2)/N_averages_perqbin
			ELSEWHERE
				fourier_component_re=0.0d0
			END WHERE
			!... the reason for this is to keep the OMP CRITICAL as limited as possible.
			!Could even be ATOMIC I think
			!$OMP CRITICAL(sfactor_updates)
			total_structure_factor(:)=total_structure_factor(:)+fourier_component_re(:)
			!$OMP END CRITICAL(sfactor_updates)
			!$OMP CRITICAL(progress)
			IF (VERBOSE_OUTPUT) CALL print_progress()
			!$OMP END CRITICAL(progress)
		ENDDO
		!$OMP END DO
		IF (ALLOCATED(fourier_component_SITE_re)) THEN
			DEALLOCATE(fourier_component_SITE_re,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(N_averages_perqbin)) THEN
			DEALLOCATE(N_averages_perqbin,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(fourier_component_SITE_im)) THEN
			DEALLOCATE(fourier_component_SITE_im,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(fourier_component_re)) THEN
			DEALLOCATE(fourier_component_re,STAT=deallocstatus)
			IF (deallocstatus/=0) CALL report_error(23,exit_status=deallocstatus)
		ELSE
			CALL report_error(0)
		ENDIF
		IF (ALLOCATED(fourier_component_im)) THEN
			DEALLOCATE(fourier_component_im,STAT=deallocstatus)
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
		!we need to normalise the whole structure factor: divide by volume, and divide by number of timesteps averaged
		total_structure_factor(:)=total_structure_factor(:)/&
		&(give_box_volume()*FLOAT(MAX((nsteps-1+sampling_interval)/sampling_interval,0)))
		!also, provide the normalisation, which is the square of the sum of the form factors
		DO qbin1=0,number_of_bins_in_q
			normalisation(qbin1)=0.0d0
			DO site_counter=1,number_of_sites
				normalisation(qbin1)=normalisation(qbin1)+&
				&sites(site_counter)%form_factor(qbin1)*DFLOAT(sites(site_counter)%N_atoms)
			ENDDO
			normalisation(qbin1)=normalisation(qbin1)**2
			WRITE(*,*) q_modulus(qbin1),total_structure_factor(qbin1),&
			&normalisation(qbin1)
		ENDDO
		END SUBROUTINE trajectory_total_sfactor_analysis_parallel

		SUBROUTINE perform_sfactor_analysis()
		IMPLICIT NONE
		INTEGER :: n_acceptor
			PRINT *,"call wavevectors"
			CALL choose_wavevectors()
			PRINT *,"call initialisation"
			CALL initialise_total_xray_sfactor()
			CALL trajectory_total_sfactor_analysis_parallel()
			CALL finalise_sfactor()
		END SUBROUTINE perform_sfactor_analysis

		REAL FUNCTION formfactor(element_name_input,q_in)
		IMPLICIT NONE
		CHARACTER(LEN=*),INTENT(IN) :: element_name_input
		REAL(KIND=WORKING_PRECISION),INTENT(IN) :: q_in
		REAL(KIND=WORKING_PRECISION) :: redk
			redk=(q_in/(pi*4.0))**2
			SELECT CASE (TRIM(ADJUSTL(element_name_input)))
			CASE ("H")
				formfactor&
				&=0.489918*EXP(-20.6593*redk)+0.262003*EXP(-7.74039*redk)+0.196767*EXP(-49.5519*redk)+0.049879*EXP(-2.20159*redk)+0.001305
			CASE ("He")
				formfactor=0.8734*EXP(-9.1037*redk)+0.6309*EXP(-3.3568*redk)+0.3112*EXP(-22.9276*redk)+0.178*EXP(-0.9821*redk)+0.0064
			CASE ("Li")
				formfactor=1.1282*EXP(-3.9546*redk)+0.7508*EXP(-1.0524*redk)+0.6175*EXP(-85.3905*redk)+0.4653*EXP(-168.261*redk)+0.0377
			CASE ("Be")
				formfactor=1.5919*EXP(-43.6427*redk)+1.1278*EXP(-1.8623*redk)+0.5391*EXP(-103.483*redk)+0.7029*EXP(-0.542*redk)+0.0385
			CASE ("B")
				formfactor=2.0545*EXP(-23.2185*redk)+1.3326*EXP(-1.021*redk)+1.0979*EXP(-60.3498*redk)+0.7068*EXP(-0.1403*redk)-0.1932
			CASE ("C")
				formfactor=2.31*EXP(-20.8439*redk)+1.02*EXP(-10.2075*redk)+1.5886*EXP(-0.5687*redk)+0.865*EXP(-51.6512*redk)+0.2156
			CASE ("N")
				formfactor=12.2126*EXP(-0.0057*redk)+3.1322*EXP(-9.8933*redk)+2.0125*EXP(-28.9975*redk)+1.1663*EXP(-0.5826*redk)-11.529
			CASE ("O")
				formfactor=3.0485*EXP(-13.2771*redk)+2.2868*EXP(-5.7011*redk)+1.5463*EXP(-0.3239*redk)+0.867*EXP(-32.9089*redk)+0.2508
			CASE ("F")
				formfactor=3.5392*EXP(-10.2825*redk)+2.6412*EXP(-4.2944*redk)+1.517*EXP(-0.2615*redk)+1.0243*EXP(-26.1476*redk)+0.2776
			CASE ("Ne")
				formfactor=3.9553*EXP(-8.4042*redk)+3.1125*EXP(-3.4262*redk)+1.4546*EXP(-0.2306*redk)+1.1251*EXP(-21.7184*redk)+0.3515
			CASE ("Na")
				formfactor=4.7626*EXP(-3.285*redk)+3.1736*EXP(-8.8422*redk)+1.2674*EXP(-0.3136*redk)+1.1128*EXP(-129.424*redk)+0.676
			CASE ("Mg")
				formfactor=5.4204*EXP(-2.8275*redk)+2.1735*EXP(-79.2611*redk)+1.2269*EXP(-0.3808*redk)+2.3073*EXP(-7.1937*redk)+0.8584
			CASE ("Al")
				formfactor=6.4202*EXP(-3.0387*redk)+1.9002*EXP(-0.7426*redk)+1.5936*EXP(-31.5472*redk)+1.9646*EXP(-85.0886*redk)+1.1151
			CASE ("Si")
				formfactor=6.2915*EXP(-2.4386*redk)+3.0353*EXP(-32.3337*redk)+1.9891*EXP(-0.6785*redk)+1.541*EXP(-81.6937*redk)+1.1407
			CASE ("P")
				formfactor=6.4345*EXP(-1.9067*redk)+4.1791*EXP(-27.157*redk)+1.78*EXP(-0.526*redk)+1.4908*EXP(-68.1645*redk)+1.1149
			CASE ("S")
				formfactor=6.9053*EXP(-1.4679*redk)+5.2034*EXP(-22.2151*redk)+1.4379*EXP(-0.2536*redk)+1.5863*EXP(-56.172*redk)+0.8669
			CASE ("Cl")
				formfactor=11.4604*EXP(-0.0104*redk)+7.1964*EXP(-1.1662*redk)+6.2556*EXP(-18.5194*redk)+1.6455*EXP(-47.7784*redk)-9.5574
			CASE ("Ar")
				formfactor=7.4845*EXP(-0.9072*redk)+6.7723*EXP(-14.8407*redk)+0.6539*EXP(-43.8983*redk)+1.6442*EXP(-33.3929*redk)+1.4445
			CASE ("K")
				formfactor=8.2186*EXP(-12.7949*redk)+7.4398*EXP(-0.7748*redk)+1.0519*EXP(-213.187*redk)+0.8659*EXP(-41.6841*redk)+1.4228
			CASE ("Ca")
				formfactor=8.6266*EXP(-10.4421*redk)+7.3873*EXP(-0.6599*redk)+1.5899*EXP(-85.7484*redk)+1.0211*EXP(-178.437*redk)+1.3751
			CASE ("Sc")
				formfactor=9.189*EXP(-9.0213*redk)+7.3679*EXP(-0.5729*redk)+1.6409*EXP(-136.108*redk)+1.468*EXP(-51.3531*redk)+1.3329
			CASE ("Ti")
				formfactor=9.7595*EXP(-7.8508*redk)+7.3558*EXP(-0.5*redk)+1.6991*EXP(-35.6338*redk)+1.9021*EXP(-116.105*redk)+1.2807
			CASE ("V")
				formfactor=10.2971*EXP(-6.8657*redk)+7.3511*EXP(-0.4385*redk)+2.0703*EXP(-26.8938*redk)+2.0571*EXP(-102.478*redk)+1.2199
			CASE ("Cr")
				formfactor=10.6406*EXP(-6.1038*redk)+7.3537*EXP(-0.392*redk)+3.324*EXP(-20.2626*redk)+1.4922*EXP(-98.7399*redk)+1.1832
			CASE ("Mn")
				formfactor=11.2819*EXP(-5.3409*redk)+7.3573*EXP(-0.3432*redk)+3.0193*EXP(-17.8674*redk)+2.2441*EXP(-83.7543*redk)+1.0896
			CASE ("Fe")
				formfactor=11.7695*EXP(-4.7611*redk)+7.3573*EXP(-0.3072*redk)+3.5222*EXP(-15.3535*redk)+2.3045*EXP(-76.8805*redk)+1.0369
			CASE ("Co")
				formfactor=12.2841*EXP(-4.2791*redk)+7.3409*EXP(-0.2784*redk)+4.0034*EXP(-13.5359*redk)+2.3488*EXP(-71.1692*redk)+1.0118
			CASE ("Ni")
				formfactor=12.8376*EXP(-3.8785*redk)+7.292*EXP(-0.2565*redk)+4.4438*EXP(-12.1763*redk)+2.38*EXP(-66.3421*redk)+1.0341
			CASE ("Cu")
				formfactor=13.338*EXP(-3.5828*redk)+7.1676*EXP(-0.247*redk)+5.6158*EXP(-11.3966*redk)+1.6735*EXP(-64.8126*redk)+1.191
			CASE ("Zn")
				formfactor=14.0743*EXP(-3.2655*redk)+7.0318*EXP(-0.2333*redk)+5.1652*EXP(-10.3163*redk)+2.41*EXP(-58.7097*redk)+1.3041
			CASE ("Ga")
				formfactor=15.2354*EXP(-3.0669*redk)+6.7006*EXP(-0.2412*redk)+4.3591*EXP(-10.7805*redk)+2.9623*EXP(-61.4135*redk)+1.7189
			CASE ("Ge")
				formfactor=16.0816*EXP(-2.8509*redk)+6.3747*EXP(-0.2516*redk)+3.7068*EXP(-11.4468*redk)+3.683*EXP(-54.7625*redk)+2.1313
			CASE ("As")
				formfactor=16.6723*EXP(-2.6345*redk)+6.0701*EXP(-0.2647*redk)+3.4313*EXP(-12.9479*redk)+4.2779*EXP(-47.7972*redk)+2.531
			CASE ("Se")
				formfactor=17.0006*EXP(-2.4098*redk)+5.8196*EXP(-0.2726*redk)+3.9731*EXP(-15.2372*redk)+4.3543*EXP(-43.8163*redk)+2.8409
			CASE ("Br")
				formfactor=17.1789*EXP(-2.1723*redk)+5.2358*EXP(-16.5796*redk)+5.6377*EXP(-0.2609*redk)+3.9851*EXP(-41.4328*redk)+2.9557
			CASE ("Kr")
				formfactor=17.3555*EXP(-1.9384*redk)+6.7286*EXP(-16.5623*redk)+5.5493*EXP(-0.2261*redk)+3.5375*EXP(-39.3972*redk)+2.825
			CASE ("Rb")
				formfactor=17.1784*EXP(-1.7888*redk)+9.6435*EXP(-17.3151*redk)+5.1399*EXP(-0.2748*redk)+1.5292*EXP(-164.934*redk)+3.4873
			CASE ("Sr")
				formfactor=17.5663*EXP(-1.5564*redk)+9.8184*EXP(-14.0988*redk)+5.422*EXP(-0.1664*redk)+2.6694*EXP(-132.376*redk)+2.5064
			CASE ("Y")
				formfactor=17.776*EXP(-1.4029*redk)+10.2946*EXP(-12.8006*redk)+5.72629*EXP(-0.125599*redk)+3.26588*EXP(-104.354*redk)+1.91213
			CASE ("Zr")
				formfactor=17.8765*EXP(-1.27618*redk)+10.948*EXP(-11.916*redk)+5.41732*EXP(-0.117622*redk)+3.65721*EXP(-87.6627*redk)+2.06929
			CASE ("Nb")
				formfactor=17.6142*EXP(-1.18865*redk)+12.0144*EXP(-11.766*redk)+4.04183*EXP(-0.204785*redk)+3.53346*EXP(-69.7957*redk)+3.75591
			CASE ("Mo")
				formfactor=3.7025*EXP(-0.2772*redk)+17.2356*EXP(-1.0958*redk)+12.8876*EXP(-11.004*redk)+3.7429*EXP(-61.6584*redk)+4.3875
			CASE ("Ru")
				formfactor=19.2674*EXP(-0.80852*redk)+12.9182*EXP(-8.43467*redk)+4.86337*EXP(-24.7997*redk)+1.56756*EXP(-94.2928*redk)+5.37874
			CASE ("Rh")
				formfactor=19.2957*EXP(-0.751536*redk)+14.3501*EXP(-8.21758*redk)+4.73425*EXP(-25.8749*redk)+1.28918*EXP(-98.6062*redk)+5.328
			CASE ("Pd")
				formfactor=19.3319*EXP(-0.698655*redk)+15.5017*EXP(-7.98929*redk)+5.29537*EXP(-25.2052*redk)+0.605844*EXP(-76.8986*redk)+5.26593
			CASE ("Ag")
				formfactor=19.2808*EXP(-0.6446*redk)+16.6885*EXP(-7.4726*redk)+4.8045*EXP(-24.6605*redk)+1.0463*EXP(-99.8156*redk)+5.179
			CASE ("Cd")
				formfactor=19.2214*EXP(-0.5946*redk)+17.6444*EXP(-6.9089*redk)+4.461*EXP(-24.7008*redk)+1.6029*EXP(-87.4825*redk)+5.0694
			CASE ("In")
				formfactor=19.1624*EXP(-0.5476*redk)+18.5596*EXP(-6.3776*redk)+4.2948*EXP(-25.8499*redk)+2.0396*EXP(-92.8029*redk)+4.9391
			CASE ("Sn")
				formfactor=19.1889*EXP(-5.8303*redk)+19.1005*EXP(-0.5031*redk)+4.4585*EXP(-26.8909*redk)+2.4663*EXP(-83.9571*redk)+4.7821
			CASE ("Sb")
				formfactor=19.6418*EXP(-5.3034*redk)+19.0455*EXP(-0.4607*redk)+5.0371*EXP(-27.9074*redk)+2.6827*EXP(-75.2825*redk)+4.5909
			CASE ("Te")
				formfactor=19.9644*EXP(-4.81742*redk)+19.0138*EXP(-0.420885*redk)+6.14487*EXP(-28.5284*redk)+2.5239*EXP(-70.8403*redk)+4.352
			CASE ("I")
				formfactor=20.1472*EXP(-4.347*redk)+18.9949*EXP(-0.3814*redk)+7.5138*EXP(-27.766*redk)+2.2735*EXP(-66.8776*redk)+4.0712
			CASE ("Xe")
				formfactor=20.2933*EXP(-3.9282*redk)+19.0298*EXP(-0.344*redk)+8.9767*EXP(-26.4659*redk)+1.99*EXP(-64.2658*redk)+3.7118
			CASE ("Cs")
				formfactor=20.3892*EXP(-3.569*redk)+19.1062*EXP(-0.3107*redk)+10.662*EXP(-24.3879*redk)+1.4953*EXP(-213.904*redk)+3.3352
			CASE ("Ba")
				formfactor=20.3361*EXP(-3.216*redk)+19.297*EXP(-0.2756*redk)+10.888*EXP(-20.2073*redk)+2.6959*EXP(-167.202*redk)+2.7731
			CASE ("La")
				formfactor=20.578*EXP(-2.94817*redk)+19.599*EXP(-0.244475*redk)+11.3727*EXP(-18.7726*redk)+3.28719*EXP(-133.124*redk)+2.14678
			CASE ("Ce")
				formfactor=21.1671*EXP(-2.81219*redk)+19.7695*EXP(-0.226836*redk)+11.8513*EXP(-17.6083*redk)+3.33049*EXP(-127.113*redk)+1.86264
			CASE ("Pr")
				formfactor=22.044*EXP(-2.77393*redk)+19.6697*EXP(-0.222087*redk)+12.3856*EXP(-16.7669*redk)+2.82428*EXP(-143.644*redk)+2.0583
			CASE ("Nd")
				formfactor=22.6845*EXP(-2.66248*redk)+19.6847*EXP(-0.210628*redk)+12.774*EXP(-15.885*redk)+2.85137*EXP(-137.903*redk)+1.98486
			CASE ("Sm")
				formfactor=24.0042*EXP(-2.47274*redk)+19.4258*EXP(-0.196451*redk)+13.4396*EXP(-14.3996*redk)+2.89604*EXP(-128.007*redk)+2.20963
			CASE ("Eu")
				formfactor=24.6274*EXP(-2.3879*redk)+19.0886*EXP(-0.1942*redk)+13.7603*EXP(-13.7546*redk)+2.9227*EXP(-123.174*redk)+2.5745
			CASE ("Gd")
				formfactor=25.0709*EXP(-2.25341*redk)+19.0798*EXP(-0.181951*redk)+13.8518*EXP(-12.9331*redk)+3.54545*EXP(-101.398*redk)+2.4196
			CASE ("Tb")
				formfactor=25.8976*EXP(-2.24256*redk)+18.2185*EXP(-0.196143*redk)+14.3167*EXP(-12.6648*redk)+2.95354*EXP(-115.362*redk)+3.58324
			CASE ("Dy")
				formfactor=26.507*EXP(-2.1802*redk)+17.6383*EXP(-0.202172*redk)+14.5596*EXP(-12.1899*redk)+2.96577*EXP(-111.874*redk)+4.29728
			CASE ("Ho")
				formfactor=26.9049*EXP(-2.07051*redk)+17.294*EXP(-0.19794*redk)+14.5583*EXP(-11.4407*redk)+3.63837*EXP(-92.6566*redk)+4.56796
			CASE ("Er")
				formfactor=27.6563*EXP(-2.07356*redk)+16.4285*EXP(-0.223545*redk)+14.9779*EXP(-11.3604*redk)+2.98233*EXP(-105.703*redk)+5.92046
			CASE ("Tm")
				formfactor=28.1819*EXP(-2.02859*redk)+15.8851*EXP(-0.238849*redk)+15.1542*EXP(-10.9975*redk)+2.98706*EXP(-102.961*redk)+6.75621
			CASE ("Yb")
				formfactor=28.6641*EXP(-1.9889*redk)+15.4345*EXP(-0.257119*redk)+15.3087*EXP(-10.6647*redk)+2.98963*EXP(-100.417*redk)+7.56672
			CASE ("Lu")
				formfactor=28.9476*EXP(-1.90182*redk)+15.2208*EXP(-9.98519*redk)+15.1*EXP(-0.261033*redk)+3.71601*EXP(-84.3298*redk)+7.97628
			CASE ("Hf")
				formfactor=29.144*EXP(-1.83262*redk)+15.1726*EXP(-9.5999*redk)+14.7586*EXP(-0.275116*redk)+4.30013*EXP(-72.029*redk)+8.58154
			CASE ("Ta")
				formfactor=29.2024*EXP(-1.77333*redk)+15.2293*EXP(-9.37046*redk)+14.5135*EXP(-0.295977*redk)+4.76492*EXP(-63.3644*redk)+9.24354
			CASE ("W")
				formfactor=29.0818*EXP(-1.72029*redk)+15.43*EXP(-9.2259*redk)+14.4327*EXP(-0.321703*redk)+5.11982*EXP(-57.056*redk)+9.8875
			CASE ("Re")
				formfactor=28.7621*EXP(-1.67191*redk)+15.7189*EXP(-9.09227*redk)+14.5564*EXP(-0.3505*redk)+5.44174*EXP(-52.0861*redk)+10.472
			CASE ("Os")
				formfactor=28.1894*EXP(-1.62903*redk)+16.155*EXP(-8.97948*redk)+14.9305*EXP(-0.382661*redk)+5.67589*EXP(-48.1647*redk)+11.0005
			CASE ("Ir")
				formfactor=27.3049*EXP(-1.59279*redk)+16.7296*EXP(-8.86553*redk)+15.6115*EXP(-0.417916*redk)+5.83377*EXP(-45.0011*redk)+11.4722
			CASE ("Pt")
				formfactor=27.0059*EXP(-1.51293*redk)+17.7639*EXP(-8.81174*redk)+15.7131*EXP(-0.424593*redk)+5.7837*EXP(-38.6103*redk)+11.6883
			CASE ("Au")
				formfactor=16.8819*EXP(-0.4611*redk)+18.5913*EXP(-8.6216*redk)+25.5582*EXP(-1.4826*redk)+5.86*EXP(-36.3956*redk)+12.0658
			CASE ("Hg")
				formfactor=20.6809*EXP(-0.545*redk)+19.0417*EXP(-8.4484*redk)+21.6575*EXP(-1.5729*redk)+5.9676*EXP(-38.3246*redk)+12.6089
			CASE ("Tl")
				formfactor=27.5446*EXP(-0.65515*redk)+19.1584*EXP(-8.70751*redk)+15.538*EXP(-1.96347*redk)+5.52593*EXP(-45.8149*redk)+13.1746
			CASE ("Pb")
				formfactor=31.0617*EXP(-0.6902*redk)+13.0637*EXP(-2.3576*redk)+18.442*EXP(-8.618*redk)+5.9696*EXP(-47.2579*redk)+13.4118
			CASE ("Bi")
				formfactor=33.3689*EXP(-0.704*redk)+12.951*EXP(-2.9238*redk)+16.5877*EXP(-8.7937*redk)+6.4692*EXP(-48.0093*redk)+13.5782
			CASE ("Th")
				formfactor=35.5645*EXP(-0.563359*redk)+23.4219*EXP(-3.46204*redk)+12.7473*EXP(-17.8309*redk)+4.80703*EXP(-99.1722*redk)+13.4314
			CASE ("Pa")
				formfactor=35.8847*EXP(-0.547751*redk)+23.2948*EXP(-3.41519*redk)+14.1891*EXP(-16.9235*redk)+4.17287*EXP(-105.251*redk)+13.4287
			CASE ("U")
				formfactor=36.0228*EXP(-0.5293*redk)+23.4128*EXP(-3.3253*redk)+14.9491*EXP(-16.0927*redk)+4.188*EXP(-100.613*redk)+13.3966
			CASE DEFAULT
				formfactor = 0.0  ! Unknown element
			END SELECT
		END FUNCTION formfactor

END MODULE SFACTOR
!--------------------------------------------------------------------------------------------------------------------------------!