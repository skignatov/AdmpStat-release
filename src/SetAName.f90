Subroutine SetAName(na,n,AName)
USE Elements, Only: ElName
Character(*) AName
Character(6) Buf

ll=Len(AName)
AName=repeat(' ',ll)

Write(Buf,'(a2,i4)')ElName(na),n
l=0
Do i=1,Len_Trim(Buf)
	If (buf(i:i)==' ') Cycle
	l=l+1
	AName(l:l)=Buf(i:i)
Enddo

End

