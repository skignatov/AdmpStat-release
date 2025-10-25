Module StringMod
Implicit Real(8) (A-H,O-Z)

Integer(4),parameter::MaxSubStr=50,LenStrDef=255

Character(LenStrDef) SubStr(MaxSubStr)

! Module:
! Function   ToUpperCase(String)                - converts all symbols of String into upper case
! Subroutine SubString(String,NSubstr,SubStr)   - divide String into substrings SubStr
! Subroutine CleanStr(Str)                      - remove all non-visible symbols from strings and replace them with spaces
! Subroutine RemoveStrComment(mode,symb,Str)    - remove comments (marked with symb) from the beginning (mode<0) or to the end (mode>0)
! Subroutine ReadWordLeftShift(Str1,Str2)       - read the word from Str1, remove it from Str1, place it to Str2, and adjust Str1 to the left

CONTAINS
!*****************************************************
Character(len=LenStrDef) Function ToUpperCase(string)
Character(len=*) string

ll=Len_Trim(string)
!If (ll>132) ll=132
ToUpperCase=''
Do i=1,ll
   ic=ICHAR(string(i:i))
   If (ic>=97 .and. ic<=122) Then
		ToUpperCase(i:i)=CHAR(ic-32)
	Else
	    ToUpperCase(i:i)=CHAR(ic)
	Endif
Enddo
Return
End function
!*****************************************************
Subroutine SubString(String,NSubstr,SubStr)

Character(*) string
Character(*) substr(MaxSubStr)
		
ls=Len(string)
lss=Len(substr)
ll=Len_Trim(string)

do i=1,MaxSubStr
	substr(i)=repeat(' ',lss)
enddo

nsubstr=0

     i=0
1	 i=i+1
        If (i>ls) Return
	    If (string(i:i)==' ') goto 1
		k=0
		nsubstr=nsubstr+1
Cyc2:	Do j=i,ll
			If (string(j:j)==' ') exit Cyc2
			k=k+1
		SubStr(nsubstr)(k:k)=string(j:j)
		Enddo cyc2
	    i=i+k
		If (i>=ll) Return
	  goto 1

End subroutine
!*****************************************************
Subroutine CleanStr(Str)
Implicit Real(8) (A-H,O-Z)
Character(len=*) Str

! Subroutine CleanStr replaces all non-printable ASCII characters in Str by SPACE

lStr=Len(Str)
Do i=1,lStr
    ich=ICHAR(Str(i:i))
    If (ich>=32.and.ich<=126) Cycle
    If (ich>=128.and.ich<=254) Cycle
    Str(i:i)=CHAR(32)
Enddo

End subroutine
!******************************************************************
Subroutine ReadWordLeftShift(Str1,Str2)

Character(len=*) Str1,Str2

! Subroutine ReadWordLeftShift reads a word (sequence of non-blank characters)
! from the beginning of Str1 and puts it to Str2. Str1 is shifted to the left

l1=Len(Str1)
l2=Len(Str2)

Str1=AdjustL(Str1)
Call CleanStr(Str1)
If (Len_Trim(Str1)==0) Then
    Str2(1:l2)=' '
    Return
Endif

k2=INDEX(Str1,' ')-1
Str2(1:k2)=Str1(1:k2)
Str2(k2+1:l2)=' '
Str1(1:k2)=' '
Str1=AdjustL(Str1)

End subroutine
!******************************************************************
Subroutine RemoveStrComment(mode,symb,Str)
    
Implicit Real(8) (A-H,O-Z)

Character(*) Str
Character(1) Symb

!Subroutine RemoveStrComment looks through Str and removes comments marked with Symb.
! If mode>0 it removes comments from symb position to the end of Str
! If mode<0 it removes comments from the beginning of Str to the Symb position (inclusively with symb)

is=INDEX(Str,symb)
If (is==0) Return

If (mode>=0) Then
        ls=Len(Str)
        Str(is:ls)=repeat(' ',ls-is+1)
Else
        Str(1:is)=repeat(' ',is)
Endif

End subroutine
!******************************************************************
Function Fread(string)
Implicit Real(8) (A-H,O-Z)

Character(len=*),intent(in):: string
Character(Len=Len(String)) str

str=AdjustL(string)
str=AdjustR(str)

Read(str,'(g<LenStrDef>.0)')Fread

End function Fread
!******************************************************************
Function Iread(string)
Implicit Real(8) (A-H,O-Z)

Character(len=*),intent(in):: string
Character(Len=Len(String)) str

str=AdjustL(string)
str=AdjustR(str)

Read(str,'(i<LenStrDef>)')Iread

End function Iread

End module
