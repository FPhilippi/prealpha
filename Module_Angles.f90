
!This Module can be used to perform rotations and also to turn a set of coordinates into a dihedral angle.
MODULE ANGLES ! Copyright (C) !RELEASEYEAR! Frederik Philippi
    USE SETTINGS
	IMPLICIT NONE
	PRIVATE
	!default values
	REAL(KIND=WORKING_PRECISION),PARAMETER :: tolerance=0.0001d0!tolerance in degrees used to define identity of angles. Also used in some other case as distance criterion.
	!variables
    REAL(KIND=WORKING_PRECISION) :: uvector(3) !uvector is the unit vector around which the rotation will occur
	REAL(KIND=WORKING_PRECISION) :: sina,cosa !sina and cosa are actually sinus(alpha) and cosinus(alpha), with alpha being the angle of the rotation.
    REAL(KIND=WORKING_PRECISION) :: rotation(3,3) !the rotation matrix
	!PRIVATE/PUBLIC declarations
	PUBLIC :: prepare_rotation,rotate,dihedral_angle,rotate_random!dihedral_angle changes the rotation matrix!
	PRIVATE :: make_rotation_matrix,uvector,sina,cosa,tolerance,rotation!Thou shalt not temper with the matrix.
	CONTAINS

		!brute-force generation of a uniformly distributed random rotation matrix. needs srand to be called before
		SUBROUTINE rotate_random(vector)
		IMPLICIT NONE
		REAL(KIND=GENERAL_PRECISION),INTENT(INOUT) ::  vector(3)
		REAL(KIND=WORKING_PRECISION) :: x,y,z
		REAL(KIND=WORKING_PRECISION) :: cosg,sing,cosb,sinb,alpha,beta,gamma,sina2,cosa2
		REAL(KIND=WORKING_PRECISION) :: rotation_local(3,3) !the rotation matrix - needs to be a local variable for parallelisation
			alpha=RAND()*twopi
			beta=RAND()*twopi
			gamma=RAND()*twopi
			sina2=SIN(alpha)
			cosa2=COS(alpha)
			sinb=SIN(beta)
			cosb=COS(beta)
			sing=SIN(gamma)
			cosg=COS(gamma)
			rotation_local(1,1)=cosb*cosg
			rotation_local(1,2)=-sing*cosb
			rotation_local(1,3)=sinb
			rotation_local(2,1)=sina2*sinb*cosg+sing*cosa2
			rotation_local(2,2)=-sina2*sinb*sing+cosa2*cosg
			rotation_local(2,3)=-sina2*cosb
			rotation_local(3,1)=-sinb*cosa2*cosg+sina2*sing
			rotation_local(3,2)=sinb*sing*cosa2+sina2*cosg
			rotation_local(3,3)=cosa2*cosb
			x=vector(1)
			y=vector(2)
			z=vector(3)
			vector(1)=x*rotation_local(1,1)+y*rotation_local(1,2)+z*rotation_local(1,3)
			vector(2)=x*rotation_local(2,1)+y*rotation_local(2,2)+z*rotation_local(2,3)
			vector(3)=x*rotation_local(3,1)+y*rotation_local(3,2)+z*rotation_local(3,3)
		END SUBROUTINE rotate_random

		SUBROUTINE make_rotation_matrix()! This subroutine builds the rotation matrix. Not accessible globally.
		IMPLICIT NONE
			rotation(1,1)=cosa+(uvector(1)*uvector(1))*(1.0d0-cosa)
			rotation(1,2)=uvector(1)*uvector(2)*(1.0d0-cosa)-uvector(3)*sina
			rotation(1,3)=uvector(1)*uvector(3)*(1.0d0-cosa)+uvector(2)*sina
			rotation(2,1)=uvector(1)*uvector(2)*(1.0d0-cosa)+uvector(3)*sina
			rotation(2,2)=cosa+(uvector(2)*uvector(2))*(1.0d0-cosa)
			rotation(2,3)=uvector(2)*uvector(3)*(1.0d0-cosa)-uvector(1)*sina
			rotation(3,1)=uvector(1)*uvector(3)*(1.0d0-cosa)-uvector(2)*sina
			rotation(3,2)=uvector(3)*uvector(2)*(1.0d0-cosa)+uvector(1)*sina
			rotation(3,3)=cosa+(uvector(3)*uvector(3))*(1.0d0-cosa)
		END SUBROUTINE make_rotation_matrix

		SUBROUTINE prepare_rotation(startvector,targetvector,aligned)!prepares the rotation matrix that maps 'startvector' onto 'targetvector' by rotation around an axis perpendicular to both these vectors.
		IMPLICIT NONE
		LOGICAL,INTENT(OUT),OPTIONAL :: aligned
		REAL(KIND=GENERAL_PRECISION),INTENT(IN) :: startvector(3),targetvector(3)
		REAL(KIND=WORKING_PRECISION) :: a(3),b(3),angle,dummy_axis(3)
			a=startvector
			b=targetvector
			CALL normalize3D(a)
			CALL normalize3D(b)
			cosa=DOT_PRODUCT(a,b)!get the angle between a and b
			sina=SQRT(1-cosa*cosa)!also sinus(angle), needed for the rotation matrix later.
			angle=(ACOS(cosa)*degrees)!calculate angle in degrees, mainly for debugging purposes
			IF (PRESENT(aligned)) aligned=.FALSE. !Some external procedures need to know whether the vectors are aligned, hence the optional variable.
			IF (angle<=tolerance) THEN !angle is zero - vectors are aligned!
				IF (DEVELOPERS_VERSION) WRITE(*,'(A,E7.1)') "  ! vectors are aligned, angle = ",angle
				IF (PRESENT(aligned)) aligned=.TRUE.
				!The rotation matrix is the identity matrix, because nothing should change.
				rotation(:,:)=0.0d0
				rotation(1,1)=1.0d0
				rotation(2,2)=1.0d0
				rotation(3,3)=1.0d0
			ELSE
				IF ((180.0-angle)<=tolerance) THEN
					IF (DEVELOPERS_VERSION) WRITE(*,'(A,F7.3)') "  ! vectors are antiparallel, angle = ",angle
					IF ((a(1)<=(1.0d0-tolerance)).AND.(a(1)>=(-1.0d0+tolerance))) THEN !This part is just to handle the rare case of antiparallel a and b along the x axis...
						dummy_axis(:)=0.0d0
						dummy_axis(1)=1.0d0
					ELSE
						dummy_axis(:)=0.0d0
						dummy_axis(3)=1.0d0
					ENDIF
					uvector=crossproduct(a,dummy_axis)
				ELSE !The two vectors are neither aligned nor antiparallel. This should be the normal case.
					uvector=crossproduct(a,b)
				ENDIF
				CALL normalize3D(uvector)
				CALL make_rotation_matrix()
			ENDIF
		END SUBROUTINE prepare_rotation

		SUBROUTINE rotate(vector)!subroutine that rotates the given vector, using the rotation matrix ("Rotation")
		IMPLICIT NONE
		REAL(KIND=GENERAL_PRECISION),INTENT(INOUT) ::  vector(3)
		REAL(KIND=WORKING_PRECISION) :: x,y,z
				x=vector(1)
				y=vector(2)
				z=vector(3)
				!I don't like MATMUL.
				vector(1)=x*Rotation(1,1)+y*Rotation(1,2)+z*Rotation(1,3)
				vector(2)=x*Rotation(2,1)+y*Rotation(2,2)+z*Rotation(2,3)
				vector(3)=x*Rotation(3,1)+y*Rotation(3,2)+z*Rotation(3,3)
		END SUBROUTINE rotate

		!dihedral_angle returns the dihedral from 0 to 360° if minusplus=FALSE (default), or from -180 to +180 if minusplus=TRUE. minusplus is an optional variable.
		REAL(KIND=GENERAL_PRECISION) FUNCTION dihedral_angle(dihedral_members,minusplus) ! dihedral_members is an array containing the four positions of the atoms.
		IMPLICIT NONE
		LOGICAL,INTENT(IN),OPTIONAL :: minusplus
		REAL(KIND=WORKING_PRECISION) :: dihedral_members(4,3),connection_vector(3),uvector1(3),uvector2(3)
		REAL(KIND=WORKING_PRECISION) :: projection1(2),projection2(2),ccw(2),clip,length
			connection_vector=(dihedral_members(2,:)-dihedral_members(3,:)) !vector connecting the second and third atom
			uvector1=(dihedral_members(2,:)-dihedral_members(1,:)) !vector from first to second atom
			uvector2=(dihedral_members(3,:)-dihedral_members(4,:)) !vector from fourth to third atom
			CALL normalize3D(connection_vector)
			CALL normalize3D(uvector1)
			CALL normalize3D(uvector2)
			!INITIALIZE ROTATION, RX=0
			uvector(1)=0.0d0 !X component is zero
			sina=SQRT(connection_vector(2)**2+connection_vector(3)**2)! 'sina' is here the length of the projection of the connection_vector along the x axis.
			IF (sina>=tolerance) THEN
					uvector(2)=connection_vector(3)/sina
					uvector(3)=-connection_vector(2)/sina!now, uvector is normalized.
					cosa=connection_vector(1)
					CALL make_rotation_matrix()
					CALL rotate(uvector1)
					CALL rotate(uvector2)
			ELSE!there are no components in y or z direction in the connection vector.
				IF (DEVELOPERS_VERSION) WRITE(*,*) "  ! dihedral_angle: already preoriented"
			ENDIF
			!Project to plane
			projection1=uvector1(2:3)
			projection2=uvector2(2:3)
			CALL normalize2D(projection1)
			CALL normalize2D(projection2)
			projection1=(projection1*projection2)
			length=(projection1(1)+projection1(2))
			IF (length<(-1.0d0)) THEN
				IF (length<(-1.0d0-tolerance)) CALL report_error(2)
				length=(-1.0d0)
			ELSEIF (length>(1.0d0)) THEN
				IF (length>(1.0d0+tolerance)) CALL report_error(2)
				length=(1.0d0)
			ENDIF
			clip=ACOS(length)
			ccw(1)=+COS(clip)*uvector1(2)-SIN(clip)*uvector1(3)
			ccw(2)=+SIN(clip)*uvector1(2)+COS(clip)*uvector1(3)
			CALL normalize2D(ccw)
			dihedral_angle=(clip*degrees)
			ccw=(ccw*projection2)
			clip=ccw(1)+ccw(2)
			IF (clip<(-1.0d0)) THEN
				IF (clip<(-1.0d0-tolerance)) CALL report_error(2)
				clip=(-1.0d0)
			ELSEIF (clip>(1.0d0)) THEN
				IF (clip>(1.0d0+tolerance)) CALL report_error(2)
				clip=(1.0d0)
			ENDIF
			IF ((ACOS(clip)*degrees)<=tolerance) THEN
				IF (PRESENT(minusplus)) THEN
					IF (minusplus) THEN 
						dihedral_angle=(-dihedral_angle)
					ELSE
						dihedral_angle=(360.0d0-dihedral_angle)
					ENDIF
				ELSE
					dihedral_angle=(360.0d0-dihedral_angle)
				ENDIF
			ENDIF
		END FUNCTION dihedral_angle

END MODULE ANGLES
!--------------------------------------------------------------------------------------------------------------------------------!