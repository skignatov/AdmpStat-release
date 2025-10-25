Subroutine FPassign(mode,nref,Numat,C,ires,epsmin)
Use Vars, Only: nstr,FPref,NAref,Cref,NA,istr,MaxFP

Implicit Real(8) (A-H,O-Z)

!Real(8) FP(Numat*(Numat-1)/2),D(Numat*(Numat-1)/2),Dmin(Numat*(Numat-1)/2)
Real(8) FP(MaxFP),D(MaxFP),Dmin(MaxFP)
Real(8) C(3,Numat),C0(3,Numat)

Character(255) Str

If (mode==0) Then
    Open(NewUnit=iu,File='refstr-FP.dat')
    Do iref=1,nref
        C=Cref(1:3,1:Numat,iref)
        NA(1:Numat)=NAref(1:Numat)
        Call GeoMetrics(Numat,n4,NA,C,FP)
        FPref(1:n4,iref)=FP(1:n4)
        Write(iu,'(i5,<n4>f10.5)')iref,(FP(j),j=n4,1,-1)
    Enddo        
    Close(iu)
Else
    Write(41,'(''****'')')
    Call GeoMetrics(Numat,n4,NA,C,FP)
    epsmin=999999.d99
    Write(41,'(i5,5x,<n4>f10.5)')istr,(FP(j),j=n4,1,-1)
    Do i=1,nref
        D(1:n4)=FP(1:n4)-FPref(1:n4,i)
        eps=DSQRT(dot_product(D(1:n4),D(1:n4)))
        If (eps<epsmin) Then
            epsmin=eps
            ires=i
            Dmin(1:n4)=D(1:n4)
        Endif

    Enddo
    Write(41,'(i5,5x,<n4>f10.5)')ires,(FPref(j,ires),j=n4,1,-1)
    Write(41,'(f10.5,<n4>f10.5)')epsmin,(Dmin(j),j=n4,1,-1)
    Write(41,*)

    !!! Debug
    iDebug=1
    If (iDebug==0) Return
    Do i=1,nref
        D(1:n4)=FP(1:n4)-FPref(1:n4,i)
        eps=DSQRT(dot_product(D(1:n4),D(1:n4)))
        Write(41,'(i5,5x,<n4>f10.5)')i,(FPref(j,i),j=n4,1,-1)
        Write(41,'(f10.5,<n4>f10.5)')eps,(D(j),j=n4,1,-1)
        Write(41,*)
    Enddo
    Write(Str,'(''str'',i0.4)')istr
    Call PrintNXYZnew(41,4,Numat,NA,C,15,6,1,Str)
    C0(1:3,1:Numat)=Cref(1:3,1:Numat,ires)
    Write(Str,'(''str'',i0.4,''ref'')')istr
    Call PrintNXYZnew(41,4,Numat,NA,C0,15,6,1,Str)
    Write(41,*)
    Write(41,*)
    !!!
    
Endif
    
End
!************************************************************************************
Subroutine GeoMetrics(Numat,n4,NA,C0,D)

Use Vars, Only: MaxAt,MaxFP,iFPSubmethod   !,Amass
Use Elements, Only: RA

Implicit Real(8) (A-H,O-Z)

Integer(4) NA(Numat),IA(Numat*(Numat-1)/2)
Real(8) C(3,Numat),CM(3),C0(3,Numat)
Real(8) S(MaxFP,MaxFP),D(MaxFP),V(MaxFP,MaxFP)

Real(8) Alp(Numat)

C(1:3,1:Numat)=C0(1:3,1:Numat)

! Set Gaussian exponents
!    xn=Dble(Numat)
!    bs=((Product(BoxSize))**(1.d0/3.d0))    ! Geometric mean of a,b,c
!    !Alp(1:120)=1.d0/RA(1:120)**2                       ! Set up as in Goedecker2013
!    !Alp=(Product(BoxSize)/Dble(Numat))**(-1.d0/3.d0)   !good choice
!    !Alp=2.d0/(bs*bs)
!    rb=0.5d0*bs/xn**(1.d0/3.d0)      ! Average bond radius of N atoms in Vol (average closest distance ra=(Vol/n)^(-1/3), rb=0.5*ra)
!    Alp=1/rb**2

    !Alp(1:Numat)=1.d0/RA(NA(1:Numat))**2
Do i=1,Numat
    Alp(i)=1.d0/RA(NA(i))**2
Enddo

Method=iFPSubMethod

!DEC$ DEFINE _MKL=1
If (Method==1) Then         ! Fingerprints
    n4=Numat*4
    Call MassCenter(Numat,C)
    Call Sint(Numat,Alp,C,S)
!DEC$ IF DEFINED (_MKL)
    Call DiagMKL(n4,S,MaxFP,D,V,nrot,-1)
!DEC$ ELSE
    Call JacobiSorted(S,n4,MaxFP,D,V,nrot,-1)
!DEC$ ENDIF

Else If (Method==2) Then        ! Diagonalized distances matrix
    n4=Numat
    Do i=2,Numat
        S(i,i)=0.d0
        Do j=1,i-1
            S(i,j)=Distance(i,j,Numat,C)
            S(j,i)=S(i,j)
        Enddo
    Enddo
!DEC$ IF DEFINED (_MKL)
    Call DiagMKL(n4,S,MaxFP,D,V,nrot,-1)
!DEC$ ELSE
    Call JacobiSorted(S,n4,MaxFP,D,V,nrot,-1)
!DEC$ ENDIF
    
ElseIf (Method==3) Then         ! Sorted Distances
    n4=Numat*(Numat-1)/2
    k=0
    Do i=2,Numat
        Do j=1,i-1
            k=k+1
            D(k)=-Distance(i,j,Numat,C)
        Enddo
    Enddo
    Call HeapSort1(n4,D,IA)
    D(1:n4)=-D(1:n4)
    
Endif


End
!************************************************************************************
Subroutine Sint(Numat,Alp,C,S)
!Use Elements, Only: RA

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,Numat),S(4*Numat,4*Numat),Alp(Numat)
Real(8), parameter:: A2au=0.529177d0

n4=Numat*4
S=0.d0

Do i=2,Numat
!    nai=NA(i)
!    ai=Alp(nai)
    ai=Alp(i)
    i0=(i-1)*4
    Do j=1,i-1
!        naj=NA(j)
!        aj=Alp(naj)
        aj=Alp(j)
        a1=2.d0*DSQRT(ai*aj)/(ai+aj)
        a2=(ai*aj)/(ai+aj)
        rij=Distance(i,j,Numat,C)
        ss=DSQRT(a1*a1*a1)*DEXP(-a2*rij*rij)
        a1ss=a1*ss
        spx=-a1ss*DSQRT(aj)*(C(1,i)-C(1,j))
        spy=-a1ss*DSQRT(aj)*(C(2,i)-C(2,j))
        spz=-a1ss*DSQRT(aj)*(C(3,i)-C(3,j))
        pxpx=a1ss*(1.d0-2.d0*a2*(C(1,i)-C(1,j))**2)
        pypy=a1ss*(1.d0-2.d0*a2*(C(2,i)-C(2,j))**2)
        pzpz=a1ss*(1.d0-2.d0*a2*(C(3,i)-C(3,j))**2)
        pxpy=-a1ss*2.d0*a2*(C(1,i)-C(1,j))*(C(2,i)-C(2,j))
        pxpz=-a1ss*2.d0*a2*(C(1,i)-C(1,j))*(C(3,i)-C(3,j))
        pypz=-a1ss*2.d0*a2*(C(2,i)-C(2,j))*(C(3,i)-C(3,j))
        j0=(j-1)*4
!
        S(i0+1,j0+1)=ss
        S(i0+1,j0+2)=spx
        S(i0+1,j0+3)=spy
        S(i0+1,j0+4)=spz
! - minus?
        S(i0+2,j0+1)=-spx
        S(i0+3,j0+1)=-spy
        S(i0+4,j0+1)=-spz
!
        S(i0+2,j0+2)=pxpx
        S(i0+3,j0+3)=pypy
        S(i0+4,j0+4)=pzpz
!
        S(i0+2,j0+3)=pxpy
        S(i0+3,j0+2)=pxpy
!
        S(i0+2,j0+4)=pxpz
        S(i0+4,j0+2)=pxpz
!
        S(i0+3,j0+4)=pypz
        S(i0+4,j0+3)=pypz
    Enddo
Enddo

Do i=1,n4
    S(i,i)=1.d0
    Do j=i,n4
        S(i,j)=S(j,i)
    Enddo
Enddo
        
End
