Module GraphConnMod
Use Vars, Only: MaxAt,Numat,C,RbondMax,NA

Integer(4) Marked(MaxAt)
Integer(4) iConnectivity,nmarked
Logical(1) lConn(MaxAt,MaxAt)

End module
!******************************************************
Subroutine AdjMatrix(nbndx)
Use Vars, Only: iAMx,C,Numat,rBondMax,rBondScale
Use GraphConnMod
Use Elements, Only: RA
Implicit Real(8) (A-H,O-Z)

If (rBondMax>0.d0) rbond=rBondMax
iAMx=0
nbndx=0

Do i=2,Numat
    Do j=1,i-1
        rij=Distance(i,j,Numat,C)
        If (rBondMax<0.d0) rbond=(RA(NA(i))+RA(NA(j)))*rBondScale
        If (rij<rbond) Then
            iAMx(i,j)=1
            iAMx(j,i)=1
            nbndx=nbndx+1
        Endif
    Enddo
Enddo

End
!******************************************************
Subroutine GraphConnectivity(iConnectivity0)
Use Vars, Only: iAMx,nbndx,StrDiameter,rBondScale
Use GraphConnMod
Use Elements, Only: RA
Implicit Real(8) (A-H,O-Z)

If (rBondMax>0.d0) rbond=rBondMax
lConn=0
iAMx=0
nbndx=0
StrDiameter=0.d0

Do i=2,Numat
    Do j=1,i-1
        rij=Distance(i,j,Numat,C)
        If (rij>StrDiameter) StrDiameter=rij
        If (rBondMax<0.d0) rbond=(RA(NA(i))+RA(NA(j)))*rBondScale
        If (rij<rbond) Then
            lConn(i,j)=.true.
            lConn(j,i)=.true.
            iAMx(i,j)=1
            iAMx(j,i)=1
            nbndx=nbndx+1
        Endif
    Enddo
Enddo

iConnectivity=0
Marked=0
nMarked=0
Do i=1,Numat
    If (Marked(i)==0) Then
        iConnectivity=iConnectivity+1
        Call MarkAtoms(i)
    Endif
    If (nMarked==Numat) Exit
Enddo

iConnectivity0=iConnectivity

End
!***************************************************************
Recursive Subroutine MarkAtoms(iHere)
Use GraphConnMod
Implicit Real(8) (A-H,O-Z)

If (Marked(iHere)>0) Return

Marked(iHere)=iConnectivity
nmarked=nmarked+1
Do iNext=1,Numat
    If (iNext==iHere) Cycle
    If (lConn(iHere,iNext).and.Marked(iNext)==0) Then
        Call MarkAtoms(iNext) 
        If (nmarked==Numat) Return
    Endif
Enddo

End
!************************************************************************
Subroutine AM2g6(n,iAM,g6str)
Implicit Real(8) (A-H,O-Z)

Integer(4) iAM(n,n)
Character(*) g6str
Integer(4) iX(6,1000)
Character(6) buf

! Generation of g6 string from adjacency matrix
! WARNING! Current code can be applied to n<=63 only!

k=0
kg=1
Do j=2,n
    Do i=1,j-1
        k=k+1
        If (k>6) Then
            k=1
            kg=kg+1
        Endif
        iX(k,kg)=iAM(i,j)
    Enddo
Enddo
If (k<6) iX(k+1:6,kg)=0

g6str(1:1)=CHAR(n+63)
Do i=1,kg
    Write(buf,'(6i1)')iX(1:6,i)
    Read(buf,'(b6)')itmp
    itmp=itmp+63
    g6str(i+1:i+1)=CHAR(itmp)
Enddo

End
!************************************************************************
Subroutine g62am(g6str,n,iAM)
Implicit Real(8) (A-H,O-Z)

Integer(4) iAM(n,n)
Character(1000) g6str
Integer(4) iX(1000)
Character buf*6,ch*1

! Generation of adjacency matrix from g6 string
! WARNING! Current code can be applied to n<=63 only!


!n=ICHAR(g6str(1:1))-63
ll=Len_Trim(g6str)
k=0
Do i=2,ll
    ch=g6str(i:i)
    ich=ICHAR(ch)-63
    Write(buf,'(b6)')ich
    Do j=1,6
        ch=buf(j:j)
        Read(ch,'(i1)')itmp
        k=k+1
        iX(k)=itmp
    Enddo
Enddo

iAM=0
k=0
Do j=2,n
    Do i=1,j-1
        k=k+1
        iAM(i,j)=iX(k)
        iAM(j,i)=iX(k)
    Enddo
Enddo

End