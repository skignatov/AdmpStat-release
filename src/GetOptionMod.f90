Module GetOptionMod
Use StringMod, Only: MaxSubStrLength=>LenStrDef,MaxSubStr,ToUpperCase
Implicit Real(8) (A-H,O-Z)

Integer(4),parameter::MaxRes=50,MaxSubStr1=MaxSubStr

Integer(4) lRes(MaxRes),iRes(MaxRes)
Real(8) xRes(MaxRes)
Character(255) sRes(MaxRes)

Integer(4) lRes1(MaxRes),iRes1(MaxRes)
Real(8) xRes1(MaxRes)
Character(255) sRes1(MaxRes)

End Module
!******************************************************************************
Function iGetOption(OptionName,String,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
Use StringMod, Only: ToUpperCase,SubString,SubStr,MaxSubStr
Use GetOptionMod, Only:MaxRes,MaxSubStr,MaxSubStr1,ToUpperCase,MaxSubStrLength
Use Vars, Only: iPrint6
Implicit Real(8) (A-H,O-Z)

!Integer(4),parameter::MaxSubStrLength=255,MaxSubStr=100,MaxSubStr1=50

Character(*) String
Character(*) OptionName
Character(len=Len(String)) Str,Opt,Buf,Buf1,SubStr1(MaxSubStr1)  !,SubStr(MaxSubStr),
!Character(255) ToUpperCase
Character(1) ch

Integer(4) lRes(MaxRes),iRes(MaxRes)
Real(8) xRes(MaxRes)
Character(255) sRes(MaxRes)

! Function iGetOption reads string of keywords determining if there are keyword OptionName
! in forms [Option], [Option=Val], [Option=(Val1,Val2,Val3...)], [Option(Val1,Val2,Val3)]
! If option is found (iGetOption=1, otherwise=0), it determines:
!	nRes - number of values
!   lRes - type of values (1 - integer, 2 - float,3 - string, 4-string w substr, -1 -error)
!   iRes,xRes,sRes - arrays of results
! If iPrintUnit>0 the program print out the keywords found to Unit=iPrintUnit
!
! WARNING! Make sure that strings and arrays in formal and actual arguments coincide! (Otherwise, it will work incorrect)
!


!Call PrintHelp(-1,OptionName)	! Save all keywords (uncomment this if PrintHelp system is present)

iGetOption=0
nRes=0
lRes=0
lErr=0
iRes=0
xRes=0.d0
!ForAll(i=1:50) sRes(i)=Repeat(' ',255)
Do i=1,MaxRes
	sRes(i)=Repeat(' ',255)
Enddo

Str=ToUpperCase(AdjustL(String))        ! String to be converted to upper case
Opt=ToUpperCase(AdjustL(OptionName))    ! OptionName converted to upper case
lOpt=Len(OptionName)
lStr=Len(Str)
ltStr=Len_Trim(Str)

Call SubString(Str,NSubStr,SubStr)
Do i=1,NSubStr
	If (INDEX(SubStr(i),Trim(Opt))==1) Then
		iGetOption=1
		If (iPrintUnit>0.and.iPrint6>0) Write(iPrintUnit,'('' Keyword found:  '',255a1)')(SubStr(i)(j:j),j=1,Len_Trim(SubStr(i)))
		Buf=SubStr(i)

		! Find external parentheses (if any)
		jbeg=INDEX(Buf,'(')
		jbeg1=INDEX(Buf,'=')
		If (jBeg>0) Then
			jbeg=jbeg+1
			jend=INDEX(Buf,')',Back=.TRUE.)-1
			If (jend<=1) jend=lOpt
			nRes=2
		ElseIf (jbeg1>0) Then
		! Find '='
			jbeg=jbeg1+1
			jend=ltStr
			nRes=1
		Else
			nRes=0
			Return
		Endif
		
		Buf1=AdjustL(Buf(jbeg:jend))
		lBuf1=Len_Trim(Buf1)
		Do j=1,lBuf1
			ch=Buf1(j:j)
			If (ch==','.or.ch==';') Buf1(j:j)=' '
		Enddo
		Call SubString(Buf1,nRes,SubStr1)
		Do j=1,nRes
			Buf=SubStr1(j)
            sRes(j)=AdjustL(Buf)
			inb=IsNumber(Buf,ir,xr)
			If (inb==1.or.inb==2) Then
				iRes(j)=ir
				xRes(j)=xr
			ElseIf (inb==0) Then
				inb=3
			Endif
			lRes(j)=inb
			Cycle
101			lRes(j)=-1
		Enddo
		Return
	Endif

Enddo

End
!******************************************************************************
Integer(4) Function IsNumber(String,iiRes,xRes)
Implicit Real(8) (A-H,O-Z)

Character(*) String
Character(len=Len(String)) Str,IntPart,FrcPart,OrdPart
Character(1) ch
Integer(8) i1,i2	! Just comment this line if the compiler does not support I8 (in this case the maximal number of significant digits will be less)

! Function IsNumber determines whether String contains integer number (IsNumber=1)
! float number (IsNumber=2), or just text string (IsNumber=0)

IsNumber=0
Str=Trim(AdjustL(String))
ls=Len_Trim(Str)

iSign=1
ch=Str(1:1)
If (ch=='-') iSign=-1
ibeg=1
If (ch=='-'.or.ch=='+') ibeg=2

iend=ls
If (iend<ibeg) Return
ipt=0
idg=0
isg=0
Do i=ibeg,iend
	ch=Str(i:i)
	If (ch=='.') Then
		If (ipt>0) Return	! several points found
		If (ipt>idg) Return ! Improper position of point
		ipt=i
		Cycle
	Endif
	If (ch=='d'.or.ch=='D'.or.ch=='e'.or.ch=='E') Then
		If (idg>0) Return !several degrees found
		If (idg<ipt) Return	! Improper position of degree
		idg=i
		Cycle
	Endif
	If (ch=='-'.or.ch=='+') Then
		If (isg>0) Return		! Several signs found (except main sign)
		If (i==ibeg+1) Return	! Repeated sign after main one
		If (i<=ipt) Return		! Improper position
		If (i/=idg+1) Return	! Improper position of degree sign
		isg=i
		Cycle
	Endif
	If (IsDigit(ch)==0) Return	! ch is Not a digit
Enddo

If (ipt>0.and.idg>0.and.ipt>idg) Return	! Additional check
If (idg>0.and.isg>0.and.idg>isg) Return

If (ipt>0) Then
    ie=idg
	If (idg==0) ie=iend
	IntPart=AdjustR(Str(ibeg:ipt-1))
	FrcPart=AdjustR(Str(ipt+1:ie))
	Read(IntPart,'(i100)')i1
	Read(FrcPart,'(i100)')i2
	n=Len_Trim(AdjustL(FrcPart))
	xRes=Dble(i1)+Dble(i2)*10.d0**(-n)
	xRes=Dble(iSign)*xRes
	If (idg>0) Then
		OrdPart=AdjustR(Str(idg+1:iend))
		Read(OrdPart,'(i100)')i3
		xRes=xRes*10.d0**i3
	Endif
	iiRes=INT4(xRes)
	IsNumber=2
ElseIf (ipt==0.and.idg==0.and.isg==0) Then
	IntPart=AdjustR(Str(ibeg:iend))
	Read(IntPart,'(i100)')iiRes
	If (iSign==-1) iiRes=-iiRes
	xRes=Dble(iiRes)
	IsNumber=1
ElseIf (ipt==0.and.idg/=0) Then
	IntPart=AdjustR(Str(ibeg:idg-1))
	OrdPart=AdjustR(Str(idg+1:iend))
	Read(IntPart,'(i100)')i1
	Read(OrdPart,'(i100)')i3
	xRes=Dble(iSign)*Dble(i1)*10.d0**i3
	iiRes=INT4(xRes)
	IsNumber=2
Endif

End
!******************************************************************************
Integer(4) Function IsDigit(ch)
Character(1) ch

IsDigit=0
ic=ICHAR(ch)
If (ic>=48.and.ic<=57) IsDigit=1

End
!******************************************************************************
Function iGetOptionNoUCase(OptionName,String,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
Use StringMod, Only: ToUpperCase,SubString,SubStr,MaxSubStr
Use GetOptionMod, Only:MaxRes,MaxSubStr,MaxSubStr1,ToUpperCase,MaxSubStrLength
Use Vars, Only: iPrint6
Implicit Real(8) (A-H,O-Z)

!Integer(4),parameter::MaxSubStrLength=255,MaxSubStr=100,MaxSubStr1=50

Character(*) String
Character(*) OptionName
Character(len=Len(String)) Str,Opt,SubStr1(MaxSubStr1),Buf,Buf1 !,SubStr(MaxSubStr)
!Character(255) ToUpperCase
Character(1) ch

Integer(4) lRes(MaxRes),iRes(MaxRes)
Real(8) xRes(MaxRes)
Character(255) sRes(MaxRes)

! This is the version of iGetOption whci does not convert substring to upper case (in order to keep file names given in a command intact)

! Function iGetOption reads string of keywords determining if there are keyword OptionName
! in forms [Option], [Option=Val], [Option=(Val1,Val2,Val3...)], [Option(Val1,Val2,Val3)]
! If option is found (iGetOption=1, otherwise=0), it determines:
!	nRes - number of values
!   lRes - type of values (1 - integer, 2 - float,3 - string, 4-string w substr, -1 -error)
!   iRes,xRes,sRes - arrays of results
! If iPrintUnit>0 the program print out the keywords found to Unit=iPrintUnit
!
! WARNING! Make sure that strings and arrays in formal and actual arguments coincide! (Otherwise, it will work incorrect)
!

!Call PrintHelp(-1,OptionName)	! Save all keywords (uncomment this if PrintHelp system is present)

iGetOptionNoUcase=0
nRes=0
lRes=0
lErr=0
iRes=0
xRes=0.d0
!ForAll(i=1:50) sRes(i)=Repeat(' ',255)
Do i=1,MaxRes
	sRes(i)=Repeat(' ',255)
Enddo

Str=String
lStr=Len(Str)
ltStr=Len_Trim(Str)
lOpt=Len(OptionName)
Opt=ToUpperCase(AdjustL(OptionName))    ! OptionName converted to upper case
!Str=ToUpperCase(AdjustL(String))        ! String to be converted to upper case


Call SubString(Str,NSubStr,SubStr)
Do i=1,NSubStr
	If (INDEX(ToUpperCase(SubStr(i)),Trim(Opt))==1) Then
		iGetOptionNoUcase=1
		If (iPrintUnit>0.and.iPrint6>0) Write(iPrintUnit,'('' Keyword found:  '',255a1)')(SubStr(i)(j:j),j=1,Len_Trim(SubStr(i)))
		Buf=SubStr(i)

		! Find external parentheses (if any)
		jbeg=INDEX(Buf,'(')
		jbeg1=INDEX(Buf,'=')
		If (jBeg>0) Then
			jbeg=jbeg+1
			jend=INDEX(Buf,')',Back=.TRUE.)-1
			If (jend<=1) jend=lOpt
			nRes=2
		ElseIf (jbeg1>0) Then
		! Find '='
			jbeg=jbeg1+1
			jend=ltStr
			nRes=1
		Else
			nRes=0
			Return
		Endif
		
		Buf1=AdjustL(Buf(jbeg:jend))
		lBuf1=Len_Trim(Buf1)
		Do j=1,lBuf1
			ch=Buf1(j:j)
			If (ch==','.or.ch==';') Buf1(j:j)=' '
		Enddo
		Call SubString(Buf1,nRes,SubStr1)
		Do j=1,nRes
			Buf=SubStr1(j)
            sRes(j)=AdjustL(Buf)
			inb=IsNumber(Buf,ir,xr)
			If (inb==1.or.inb==2) Then
				iRes(j)=ir
				xRes(j)=xr
			ElseIf (inb==0) Then
				inb=3
			Endif
			lRes(j)=inb
			Cycle
101			lRes(j)=-1
		Enddo
		Return
	Endif

Enddo

End


