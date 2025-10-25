Include 'CompilerVersion.f90'
Subroutine Inertia(N,C,AtomMass,CM,PMOI,PAxes,kRotTyp)
!USE Numerical_Libraries
!USE UNITS

Implicit Real(8) (A-H,O-Z)

Real(8)	C(3,N),			&	! Atom coordinates (in Angstroems)
		AtomMass(N),	&	! Atomic masses (in amu)
		CM(3),			&	! Center of masses (in the initial coordinate system)
		PMOI(3),		&	! Principal moments of inertia (in a.u.)
		PAxes(3,3),		&	! Principal axes (in the initial coordinate system)
		T(3,3)				! Tensor of inertia (in a.u.)

Real(8) CIPA(3,N)			! Temporary array (coordinates reduced to the princpal axes)

! Subroutine INERTIA calulates center of masses (CM, in Angstroems), 
! principal moments of inertia (PMOI, in a.u.),
! principal axes (PAxes), and the type of rotor (kRotTyp). 
! Input parameters: number of atoms N, atomic coordinates C(in Angstroems),
! and atomic masses AtomMass (in AMU).
! At outup, C is replaced by the coordinates centered to CM and reduced to principal axes.
! kRotTyp	= 1 - Spheric Top
!			= 2 - Prolate Symmetric Top
!			= 3 - Oblate  Symmetric Top
!			= 4 - Asymmetric Top
!			= 5 - Linear Top
!			= 6 - Single atom molecule

au=1.d0/0.529177249d0

PMOI=0.d0
PAxes=0.d0

! Check for single atom molecule
If (N==1) Then
	CM(1:3)=C(1:3,1)
	C(1:3,1)=0.d0
	TotMass=AtomMass(1)
	kRotTyp=6
	Return
Endif

!
! Center of masses
!

TotMass=0.d0
CM=0.d0
Do i=1,N
	TotMass=TotMass+AtomMass(i)
	Do k=1,3
		CM(k)=CM(k)+AtomMass(i)*C(k,i)
	Enddo
Enddo
CM=CM/TotMass

Do i=1,N
	Do k=1,3
		C(k,i)=C(k,i)-CM(k)
	Enddo
Enddo
	
!
! Tensor of Inertia
!
T=0.d0
Do i=1,N
	T(1,1)=T(1,1)+AtomMass(i)*(C(2,i)*C(2,i)+C(3,i)*C(3,i))
	T(2,2)=T(2,2)+AtomMass(i)*(C(1,i)*C(1,i)+C(3,i)*C(3,i))
	T(3,3)=T(3,3)+AtomMass(i)*(C(1,i)*C(1,i)+C(2,i)*C(2,i))
	T(1,2)=T(1,2)-AtomMass(i)*C(1,i)*C(2,i)
	T(1,3)=T(1,3)-AtomMass(i)*C(1,i)*C(3,i)
	T(2,3)=T(2,3)-AtomMass(i)*C(2,i)*C(3,i)
Enddo
T(2,1)=T(1,2)
T(3,1)=T(1,3)
T(3,2)=T(2,3)

T=T*au*au

!
! Calculate PMOI and PAxes
!

CALL Jacobi(T,3,3,PMOI,PAxes,NRotation)
!Call DEVCSF(3,T,3,PMOI,PAxes,3)

!
! Reduce to principal axes
!
do i=1,N
    do j=1,3
	  	CIPA(j,i)=(C(1,i)*PAxes(1,j)+ &
		           C(2,i)*PAxes(2,j)+ &
				   C(3,i)*PAxes(3,j))
	enddo
enddo

C=CIPA

!
! Rotor type determination
!
kRotTyp=0
ep=1.d-3
epAbs=1.d-3
dm12=DABS(DMIN1(PMOI(1),PMOI(2)))
dm23=DABS(DMIN1(PMOI(2),PMOI(3)))
IF (DABS(PMOI(1)-PMOI(2))/dm12<ep .and. DABS(PMOI(2)-PMOI(3))/dm23<ep) Then
	kRotTyp=1    ! Spherical top
!	If (IPrint==1) Write(6,'('' Spherical top'')')
ELSE IF (DABS(PMOI(1)-PMOI(2))/dm12<ep .and. DABS(PMOI(3))<epAbs) Then
    kRotTyp=5    ! Linear top
!	If (IPrint==1) Write(6,'('' Linear top'')')
ELSE IF (DABS(PMOI(1)-PMOI(2))/dm12<ep .and. PMOI(3)<PMOI(2)) Then
    kRotTyp=2    ! Prolate symmetric top
!	If (IPrint==1) Write(6,'('' Prolate symmetric top'')')
ELSE IF (DABS(PMOI(2)-PMOI(3))/dm23<ep .and. PMOI(1)>PMOI(2)) Then
    kRotTyp=3    ! Oblate symmetric top
!	If (IPrint==1) Write(6,'('' Oblate symmetric top'')')
ELSE IF (DABS(PMOI(1)-PMOI(2))/dm12>ep .and. DABS(PMOI(2)-PMOI(3))/dm23>ep) Then
    kRotTyp=4    ! Asymmetric top
!	If (IPrint==1) Write(6,'('' Asymmetric top'')')
ENDIF

End
!******************************************************************************************
Subroutine InertiaNew(mode,iSort,iRHS,N,C,AtomMass,CM,PMOI,PAxes,kRotTyp)
!DEC$ IF (_COMPILER_==1)
USE Numerical_Libraries
!DEC$ ENDIF

Implicit Real(8) (A-H,O-Z)

Real(8)	C(3,N),			&	! Atom coordinates (in Angstroems)
		AtomMass(N),	&	! Atomic masses (in amu)
		CM(3),			&	! Center of masses (in the initial coordinate system)
		PMOI(3),		&	! Principal moments of inertia (in a.u.)
		PAxes(3,3),		&	! Principal axes (in the initial coordinate system)
		T(3,3),T0(3,3)		! Tensor of inertia (in a.u.)

Real(8) CIPA(3,N)			! Temporary array (coordinates reduced to the princpal axes)
Real(8) X(3),Y(3),Z(3),P(3)	! Temporary arrays

! Subroutine INERTIA calulates center of masses (CM), principal moments of inertia (PMOI),
! principal axes (PAxes), type of rotor (kRotTyp) and reduces molecule to CM and PA if needed.
!
! At INPUT: number of atoms N, atomic coordinates C(in Angstroems or Bohrs),
! and atomic masses AtomMass (in AMU).
!
! At OUTPUT:
! CM	- center of masses (in the same units as C)
! PMOI  - principal moments of inertia in AMU*(units of C)**2
! PAxes - principal axes (normalized to 1)
! kRotTyp	= 1 - Spheric Top
!			= 2 - Prolate Symmetric Top
!			= 3 - Oblate  Symmetric Top
!			= 4 - Asymmetric Top
!			= 5 - Linear Top
!			= 6 - Single atom molecule
!
! Additional control:
!    Mode =	0	- C is not changed
!			1	- C is reduced to CM
!			2	- C is reduced to CM and PAxes
!			3	- C is reduced to PAxes without shift of CM (rotated at the initial position)
!
!    iSort= 0	- PMOI and PAxes are not sorted
!			1	- PMOI and PAxes are sorted in increasing order of PMOI
!		   -1	- PMOI and PAxes are sorted in decreasing order of PMOI
!
!	 iRHS = 1	- PAxes should form right-hand system
!		   -1	- PAxes should form left-hand system
!			0	- remain PAxes as is


!au=1.d0/0.529177249d0

PMOI=0.d0
PAxes=0.d0
Do k=1,3
	PAxes(k,k)=1.d0
Enddo

! Check for single atom molecule
If (N==1) Then
	CM(1:3)=C(1:3,1)
	If (Mode==1.or.Mode==2) C(1:3,1)=0.d0
	TotMass=AtomMass(1)
	kRotTyp=6
	Return
Endif

!
! Center of masses
!

TotMass=Sum(AtomMass(1:n))
CM=0.d0
Do i=1,N
    CM(1:3)=CM(1:3)+AtomMass(i)*C(1:3,i)
Enddo
CM=CM/TotMass

Do i=1,N
	Do k=1,3
		C(k,i)=C(k,i)-CM(k)
	Enddo
Enddo


!
! Tensor of Inertia
!
T=0.d0
Do i=1,N
	T(1,1)=T(1,1)+AtomMass(i)*(C(2,i)*C(2,i)+C(3,i)*C(3,i))
	T(2,2)=T(2,2)+AtomMass(i)*(C(1,i)*C(1,i)+C(3,i)*C(3,i))
	T(3,3)=T(3,3)+AtomMass(i)*(C(1,i)*C(1,i)+C(2,i)*C(2,i))
	T(1,2)=T(1,2)-AtomMass(i)*C(1,i)*C(2,i)
	T(1,3)=T(1,3)-AtomMass(i)*C(1,i)*C(3,i)
	T(2,3)=T(2,3)-AtomMass(i)*C(2,i)*C(3,i)
Enddo
T(2,1)=T(1,2)
T(3,1)=T(1,3)
T(3,2)=T(2,3)

!T=T*au*au

!
! Calculate PMOI and PAxes
!
T0=T
!DEC$ IF (_COMPILER_==1)
Call DEVCSF(3,T,3,PMOI,PAxes,3)				! Sorted PMOI (decreasing order)
!DEC$ ELSE
Call JacobiSorted(T,3,3,PMOI,Paxes,itmp,-1)	
!DEC$ ENDIF 


!
! Rotor type determination
!
kRotTyp=0
ep=1.d-3
epAbs=1.d-3
dm12=DABS(DMIN1(PMOI(1),PMOI(2)))
dm23=DABS(DMIN1(PMOI(2),PMOI(3)))
IF (DABS(PMOI(1)-PMOI(2))/dm12<ep .and. DABS(PMOI(2)-PMOI(3))/dm23<ep) Then
	kRotTyp=1    ! Spherical top
!	If (IPrint==1) Write(6,'('' Spherical top'')')
ELSE IF (DABS(PMOI(1)-PMOI(2))/dm12<ep .and. DABS(PMOI(3))<epAbs) Then
    kRotTyp=5    ! Linear top
!	If (IPrint==1) Write(6,'('' Linear top'')')
ELSE IF (DABS(PMOI(1)-PMOI(2))/dm12<ep .and. PMOI(3)<PMOI(2)) Then
    kRotTyp=2    ! Prolate symmetric top
!	If (IPrint==1) Write(6,'('' Prolate symmetric top'')')
ELSE IF (DABS(PMOI(2)-PMOI(3))/dm23<ep .and. PMOI(1)>PMOI(2)) Then
    kRotTyp=3    ! Oblate symmetric top
!	If (IPrint==1) Write(6,'('' Oblate symmetric top'')')
ELSE IF (DABS(PMOI(1)-PMOI(2))/dm12>ep .and. DABS(PMOI(2)-PMOI(3))/dm23>ep) Then
    kRotTyp=4    ! Asymmetric top
!	If (IPrint==1) Write(6,'('' Asymmetric top'')')
ENDIF

!
! Unsorted or increasing PMOI 
!
If (iSort>0) Then
	tmp=PMOI(1)
	PMOI(1)=PMOI(3)
	PMOI(3)=tmp
	Do k=1,3
		tmp=PAxes(k,1)
		PAxes(k,1)=PAxes(k,3)
		PAxes(k,3)=tmp
	Enddo
ElseIf (iSort==0) Then
	CALL Jacobi(T0,3,3,PMOI,PAxes,NRotation)
	Do j=1,3
		tmp=1.d0/DSQRT(Paxes(1,j)**2+PAxes(2,j)**2+PAxes(3,j)**2)
		PAxes(1:3,j)=PAxes(1:3,j)*tmp
	Enddo
Endif


!
! Check PAxes to be right-hand system
!
If (irhs/=0) Then
	X(1:3)=PAxes(1:3,1)
	Y(1:3)=PAxes(1:3,2)
	Z(1:3)=PAxes(1:3,3)
	Call CrossProduct(X,Y,P,PZ)
	PZ=P(1)*Z(1)+P(2)*Z(2)+P(3)*Z(3)
	If ((PZ<0.d0.and.irhs>0).or.(PZ>0.d0.and.irhs<0)) PAxes(1:3,3)=-Z(1:3)
Endif


!
! Reduce to principal axes
!
If (mode==2.or.mode==3) Then
	do i=1,N
		do j=1,3
	  		CIPA(j,i)=(C(1,i)*PAxes(1,j)+ &
				       C(2,i)*PAxes(2,j)+ &
					   C(3,i)*PAxes(3,j))
		enddo
	enddo
	C=CIPA
Endif

!
! Shift back from CM
!
If (mode==0.or.mode==3) ForAll(k=1:3,i=1:N) C(k,i)=C(k,i)+CM(k)


End
