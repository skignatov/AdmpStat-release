Subroutine ReadNXYZ(iUnit,Mode,iSkip,LeftDel,iForm,MaxAt,NA,C,NumAt,SearchStr)
Use StringMod, Only:ToUpperCase,SubString,SubStr,LenStrDef
Implicit Real(8) (A-H,O-Z)

Character(*) SearchStr
Real(8) C(3,MaxAt)
Integer(4) NA(MaxAt)
Character(10) AName(MaxAt)
Character(LenStrDef) String   !,ToUpperCase
Character(LenStrDef) Str132   !,SubStr(10)

! Subroutine ReadNXYZ reads from file iUnit the atomic numbers NA and cartesian
! coordinates of atoms C(3,MaxAt). Number of atoms read in is NumAt.
! The read format is determined by iForm:
!	Mode==1		search first occurence of SearchStr in the file iUnit 
!               starting from the beginning of file
!	Mode==<+n>	search n-th  occurence of SearchStr in the file
!	Mode==<-n>	search n-th  occurence of SearchStr from the end of file 
!               (i.e. Mode==-1 corresponds to the last occurrence)
!	Mode==0		start reading coordinates from the (iSkip+1)th line of the file
!               without search of SearchStr in the file
!	if iUnit<0 -reading from the current position
!	iSkip		skip iSkip lines before reading the XYZ coordinates
!	iForm==1	the coordinates are in NXYZ format
!	iForm==2	the coordinates are in NDXYZ format
!	iForm==3	the coordinates are in DNXYZ
!	LeftDel		ignore first LeftDel symbols in the string with coordinates

!ierr=-1
!
! Find the starting position
!
If (iUnit>=0) Then
	Rewind(iUnit)
	iu=iUnit
Else
	iu=-iUnit
Endif

ifound=0
line=0
If (Mode>0) Then
	Do While (.not.EOF(iu))
		Read(iu,'(a255)') String
		line=line+1
		String=ToUpperCase(String)
		If (INDEX(String,Trim(SearchStr))>0) Then
			ifound=ifound+1
			linefound=line
			If (ifound==Mode) Goto 10
		Endif
	Enddo
ElseIf (Mode<0) Then
	Do While (.not.EOF(iu))
		Read(iu,'(a255)') String
		line=line+1
		String=ToUpperCase(String)
		If (INDEX(String,Trim(SearchStr))>0) Then
			ifound=ifound+1
			linefound=line
		Endif
	Enddo
	Rewind(iu)
	If (Mode==-1) Then
		Do i=1,linefound
			Read(iu,*)
		Enddo
		Goto 10
	Endif
	ifound1=ifound+Mode+1
	line=0
	ifound=0
	Do While (.not.EOF(iu))
		Read(iu,'(a255)') String
		line=line+1
		String=ToUpperCase(String)
		If (INDEX(String,Trim(SearchStr))>0) Then
			ifound=ifound+1
			linefound=line
			If (ifound==ifound1) Goto 10
		Endif
	Enddo
Endif

10 Continue
Do i=1,iSkip
	If (.not.EOF(iu)) Read(iu,*)
Enddo

!
! Read in the coordinates
!
!ierr=0
NumAt=0
Do While (.not.EOF(iu))
	Read(iu,'(a255)')String
	String=ToUpperCase(String)
	If (Len_Trim(String)==0) Exit
	If (INDEX(String,'------')>0) Exit
	If (INDEX(String,'......')>0) Exit
	If (INDEX(String,'$END')>0) Exit
	NumAt=NumAt+1
	Str132=String(LeftDel+1:132)
	Call ParseXYZ(iForm,NNuc,X,Y,Z,ierr,Str132)
	!Call ReadWordLeft(Str132,AName(NumAt))
!    Call SubString(Str132,10,nsubstr,SubStr)
!    Read(SubStr(1),FMT=*,ERR=20)Nuc
!    Read(SubStr(3),FMT=*,ERR=20)x
!    Read(SubStr(4),FMT=*,ERR=20)y
!    Read(SubStr(5),FMT=*,ERR=20)z
	NA(NumAt)=NNuc
	C(1,NumAt)=X
	C(2,NumAt)=Y
	C(3,NumAt)=Z
!    Cycle
!20  Continue
!	NA(NumAt)=99
!	C(1,NumAt)=Dble(Numat)
!	C(2,NumAt)=Dble(Numat)
!	C(3,NumAt)=Dble(Numat)
!    ierr=1
Enddo

Return
End
!**************************************************************************
Subroutine ReadWordLeft(String,Str2)

Character(len=*) String,Str2
Character(len=Len(String)) Str1

! Subroutine ReadWordLeft reads a word (sequence of non-blank characters)
! from the beginning of Str1 and puts it to Str2

Str1=AdjustL(String)
l1=Len(Str1)
l2=Len(Str2)

If (l1==0.or.Len_Trim(Str1)==0) Then
	Str2(1:l2)=' '
	Return
Endif	

k=0
Do i=1,l1
	If (Str1(i:i)==' ') Exit
	k=k+1
	If (k>l2) Exit
	Str2(k:k)=Str1(i:i)
Enddo

If (k<l2) Str2(k+1:l2)=' '

Return
End
!***********************************************************************************
Subroutine ReadNXYZvib(iu,Mode,iSkip,LeftDel,iForm,MaxAt,NA,C,NumAt,nvibr,f1,f2,ierr,SearchStr)
Use StringMod, Only:ToUpperCase,SubString,SubStr,LenStrDef
Implicit Real(8) (A-H,O-Z)

Character(*) SearchStr
Real(8) C(3,MaxAt),Freq(3,MaxAt),Vibr1(3,MaxAt)
Integer(4) NA(MaxAt)
Character(10) AName(MaxAt)
Character(LenStrDef) Str,Str1   !,ToUpperCase
Character(LenStrDef) line           !Str132   !,SubStr(10)
Character(30) buf30

ierr=-1
Numat=0
nfreq=0
nvibr=0
nneg=0

Do While (.not.EOF(iu))
    Read(iu,'(a255)')Str
    Str1=AdjustL(ToUpperCase(Str))
    If (INDEX(Str1,'STOP')==1) Then
        ierr=-99
        Return
    Endif
    If (INDEX(Trim(Str1),Trim(ToUpperCase(SearchStr)))>0) Goto 10
Enddo
Return

10 Continue
line=repeat(' ',LenStrDef)
Do While (.not.EOF(iu))
    Read(iu,'(a255)')Str
    Str1=ToUpperCase(Str)
    If (Len_Trim(Str)==0) Exit
    Call SubString(Str,nSubStr,SubStr)
    line=Trim(SubStr(1))//' '//Trim(SubStr(2))//' '//Trim(SubStr(3))//' '//Trim(SubStr(4))
    Call ParseXYZ(1,NNuc,xx,yy,zz,ierr,line)
    Numat=Numat+1
    NA(Numat)=NNuc
    C(1,Numat)=xx
    C(2,Numat)=yy
    C(3,Numat)=zz
    If (nSubStr>=7) Then
        Do k=1,3
            nfreq=nfreq+1
            buf30=SubStr(k+4)
            Read(buf30,'(f30.10)')ff
            Freq(k,Numat)=ff
            If (ff<0.d0) nneg=nneg+1
        Enddo
    Endif
    If (nSubStr>=10) Then
        Do k=1,3
            nvibr=nvibr+1
            buf30=SubStr(k+7)
            Read(buf30,'(f30.10)')ff
            Vibr1(k,Numat)=ff
        Enddo
    Endif
Enddo

! Get min nonnegative and max frequencies
If (nfreq>0) Then
A: Do i=1,Numat
    Do k=1,3
        If (Freq(k,i)>0.d0) Then
			f1=Freq(k,i)
            Exit A
        Endif
    Enddo
Enddo A
Do i=Numat,1,-1
	Do k=3,1,-1
		If (Freq(k,i)>0.d0) Then
			f2=Freq(k,i)
            Return
        Endif
    Enddo
Enddo
Endif

End
!*********************************************************************************************
Subroutine ParseXYZ(Mode,NNuc,X,Y,Z,ierr,line)
Use Elements, Only: ElName
Use StringMod, Only: LenStrDef,ToUpperCase
Implicit Real(8) (A-H,O-Z)
Character(LenStrDef) Line,lin(4),lin1   !,ToUpperCase
Character(5) Anam,Anam1

! Subroutine ParseXYZ parses the input string Line extracting the
! atom number NNuc, and cartesian coords X, Y, Z.
! The line format is determined by Mode
! Mode==1 : AN, X, Y, Z
! Mode==2 : AN, dummy, X, Y, Z
! Mode==3 : dummy, AN, X, Y, Z
! where dummy is a variable to be skipped
! AN - element symbol, element symbol followed by number, or atomic number.

lin=repeat(' ',LenStrDef)
lin1=repeat(' ',LenStrDef)
k=0
Do i=1,5
	If (Mode==1.and.i==5) Cycle
	line=AdjustL(line)
	ll=Len_Trim(line)
	If (ll==0) Cycle
	iend=INDEX(line,' ')-1
	If (iend<=0.or.iend>ll) iend=ll
	lin1=line(1:iend)
	line(1:iend)=' '
	If (Mode==2.and.i==2) Cycle
	If (Mode==3.and.i==1) Cycle
	k=k+1
	lin(k)=lin1
	lin1=repeat(' ',LenStrDef)
Enddo

lin1=lin(2)
Read(lin1,FMT=*,ERR=10)x
lin1=lin(3)
Read(lin1,FMT=*,ERR=10)y
lin1=lin(4)
Read(lin1,FMT=*,ERR=10)z

lin1=AdjustL(Lin(1))
ll=Len_Trim(lin1)
ic=ICHAR(lin1(1:1))
If (ic>=48.and.ic<=57) Then
	iend1=INDEX(lin1,' ')-1
	iend2=INDEX(lin1,'.')-1
	If (iend2<=0) iend2=999999
	iend=MIN0(iend1,iend2)
	lin(1)=lin1(1:iend)
	Read(lin(1),FMT=*,ERR=20)NNuc
Else
	Do i=1,ll
		ic=ICHAR(lin1(i:i))
		If (ic>=48.and.ic<=57) Exit
	Enddo
	iend=MIN0(i-1,ll)
	anam=lin1(1:iend)
	Anam=ToUpperCase(Anam)
	Do NNuc=1,99
		Anam1=ToUpperCase(Elname(Nnuc))
		If (Anam==Anam1) Exit
	Enddo
Endif
ierr=0
Return

10 Continue
ierr=1
x=0.d0
y=0.d0
z=0.d0
Return

20 Continue
ierr=1
Nuc=99
Anam='        '
   
End