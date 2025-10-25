Subroutine Histo(iu,mode,iref)
Use Vars, Only: Cstr,Numat,MaxAt,NA,iAMref,nstr,istrbeg,iType

Implicit Real(8) (A-H,O-Z)

Character(1) symb
Real(8) C(3,Numat)
Real(8), allocatable:: xBin(:)
Integer(4), allocatable:: iBin(:)

! Histogram initialization
    rmin=999999.d0
    rmax=0.d0
        Do istr=istrbeg,nstr
            If (iType(6,istr)/=iref) Cycle
            C(1:3,1:Numat)=Cstr(1:3,1:Numat,istr)
            Do i=2,Numat
                Do j=1,i-1
                    If (mode==2.and.iAMref(i,j,iref)==0) Cycle
                    rr=Distance(i,j,MaxAt,C)
                    If (rr<rmin) rmin=rr
                    If (rr>rmax) rmax=rr
                Enddo
            Enddo
        Enddo

    dx=0.05d0
    dx2=dx*0.5d0
    xmin=rmin-10.d0*dx
    xmax=rmax+10.d0*dx
    nbin=INT((xmax-xmin)/dx)+1
    Allocate(xBin(nbin),iBin(nbin))
    xx=xmin
    Do i=1,nbin
        xBin(i)=xx
        xx=xx+dx
        iBin(i)=0
    Enddo

! Accumulation of histo data
    nval=0
    kstr=0
    Do istr=istrbeg,nstr
        If (iType(6,istr)/=iref) Cycle
        C(1:3,1:Numat)=Cstr(1:3,1:Numat,istr)
        kstr=kstr+1
        Do i=2,Numat
            Do j=1,i-1
                If (mode==2.and.iAMref(i,j,iref)==0) Cycle
                rr=Distance(i,j,MaxAt,C)
                ib=INT((rr-xmin)/dx)+1
                iBin(ib)=iBin(ib)+1
                nval=nval+1
            Enddo
        Enddo
    Enddo
    
! Printing 
Write(iu,'(//'' Bond lengths histogram for the structures assigned to REF'',i4\)')iref
If (mode==1) Write(iu,'(2x,''(all interatomic distances are included)'')')
If (mode==2) Write(iu,'(2x,''(only the distances corresponding to the bonds in the REF structure are included)'')')
Write(iu,'(  '' NumStr:'',i10  )')kstr
Write(iu,'(  '' Nvals :'',i10  )')nval
Write(iu,'(  '' Nbins :'',i10  )')nbin
Write(iu,'(  '' Xmin  :'',f10.4)')xmin
Write(iu,'(  '' Xmax  :'',f10.4)')xmax
Write(iu,'(  '' dX    :'',f10.4)')dx
nmax=0
Do i=1,nbin
    If (iBin(i)>nmax) nmax=iBin(i)
Enddo
!Write(iu,'(160x,''1'')')
!Write(iu,'(51x,30i10)')(j,j=0,9),100
Write(iu,'(/''  Bin      Xmin      Xmax   Xcenter   Count Normalized Scaled 0'',10i10)')(j,j=1,9),100
!Write(iu,'(60x,300i1)')((j,j=0,9),k=1,10),0
Write(iu,'(1x,162(''-''),''>'')')
Do i=1,nbin
    binbeg=xBin(i)
    binend=binbeg+dx
    binmid=binbeg+dx2
    ii=iBin(i)
    symb='*'
!    If (ii>250) Then
!        ii=250
!        symb='>'
!    Endif
!    iff=Dble(ii) ! unscaled
    iff=0
    If (nmax>0) iff=IDNINT(Dble(ii*100)/Dble(nmax))  ! scaled histo
    xii=Dble(ii)/Dble(nval)
    Write(iu,'(i5,3f10.5,i8,f10.4,i6,3x,''|'',<iff+1>a1)')i,binbeg,binend,binmid,ii,xii,iff,(symb,j=2,iff+1)
Enddo
Write(iu,*)

End
!********************************************************************************
Subroutine HistoRef(iu,mode,iref)
Use Vars, Only: Cref,Numat,MaxAt,NA,iAMref,nstr,istrbeg,iType

Implicit Real(8) (A-H,O-Z)

Character(1) symb
Real(8) C(3,Numat)
Real(8), allocatable:: xBin(:)
Integer(4), allocatable:: iBin(:)

! Histogram for Reference structure initialization
    rmin=999999.d0
    rmax=0.d0
            C(1:3,1:Numat)=Cref(1:3,1:Numat,iref)
            Do i=2,Numat
                Do j=1,i-1
                    If (mode==2.and.iAMref(i,j,iref)==0) Cycle
                    rr=Distance(i,j,MaxAt,C)
                    If (rr<rmin) rmin=rr
                    If (rr>rmax) rmax=rr
                Enddo
            Enddo

    dx=0.05d0
    dx2=dx*0.5d0
    xmin=rmin-10.d0*dx
    xmax=rmax+10.d0*dx
    nbin=INT((xmax-xmin)/dx)+1
    Allocate(xBin(nbin),iBin(nbin))
    xx=xmin
    Do i=1,nbin
        xBin(i)=xx
        xx=xx+dx
        iBin(i)=0
    Enddo

! Accumulation of histo data
    nval=0
    kstr=0
        C(1:3,1:Numat)=Cref(1:3,1:Numat,iref)
        kstr=kstr+1
        Do i=2,Numat
            Do j=1,i-1
                If (mode==2.and.iAMref(i,j,iref)==0) Cycle
                rr=Distance(i,j,MaxAt,C)
                ib=INT((rr-xmin)/dx)+1
                iBin(ib)=iBin(ib)+1
                nval=nval+1
            Enddo
        Enddo
    
! Printing 
Write(iu,'(//'' Bond lengths histogram for the reference structure'',i4\)')iref
If (mode==1) Write(iu,'(2x,''(all interatomic distances are included)'')')
If (mode==2) Write(iu,'(2x,''(only the distances corresponding to the bonds in the REF structure are included)'')')
Write(iu,'(  '' NumStr:'',i10  )')kstr
Write(iu,'(  '' Nvals :'',i10  )')nval
Write(iu,'(  '' Nbins :'',i10  )')nbin
Write(iu,'(  '' Xmin  :'',f10.4)')xmin
Write(iu,'(  '' Xmax  :'',f10.4)')xmax
Write(iu,'(  '' dX    :'',f10.4)')dx
nmax=0
Do i=1,nbin
    If (iBin(i)>nmax) nmax=iBin(i)
Enddo
!Write(iu,'(160x,''1'')')
Write(iu,'(/''  Bin      Xmin      Xmax   Xcenter   Count Normalized Scaled 0'',10i10)')(j,j=1,9),100
!Write(iu,'(60x,300i1)')((j,j=0,9),k=1,10),0
Write(iu,'(1x,160(''-''),''>'')')
Do i=1,nbin
    binbeg=xBin(i)
    binend=binbeg+dx
    binmid=binbeg+dx2
    ii=iBin(i)
    symb='*'
!    If (ii>250) Then
!        ii=250
!        symb='>'
!    Endif
!    iff=Dble(ii) ! unscaled
    iff=0
    If (nmax>0) iff=IDNINT(Dble(ii*100)/Dble(nmax))  ! scaled histo
    xii=Dble(ii)/Dble(nval)
    Write(iu,'(i5,3f10.5,i8,f10.4,i6,3x,''|'',<iff+1>a1)')i,binbeg,binend,binmid,ii,xii,iff,(symb,j=2,iff+1)
Enddo
Write(iu,*)

End
!********************************************************************************
Subroutine HistoE(iu,iref)
Use Vars, Only: xStep,Numat,MaxAt,NA,iAMref,nstr,istrbeg,iType,ebase,au2kcalmol

Implicit Real(8) (A-H,O-Z)

Character(1) symb
Real(8), allocatable:: xBin(:)
Integer(4), allocatable:: iBin(:)

! Histogram initialization
    emin=999999.d0
    emax=0.d0
        Do istr=istrbeg,nstr
            If (iref/=0.and.iType(6,istr)/=iref) Cycle
            ei=xStep(2,istr)
            erel=(ei-ebase)*au2kcalmol
!           If (mode==2.and.iAMref(i,j,iref)==0) Cycle
            If (erel<emin) emin=erel
            If (erel>emax) emax=erel
        Enddo

    dx=0.5d0
    dx2=dx*0.5d0
    xmin=emin-2.d0*dx
    xmax=emax+2.d0*dx
    If (xmin<0.d0) xmin=0.d0
    nbin=INT((xmax-xmin)/dx)+1
    Allocate(xBin(nbin),iBin(nbin))
    xx=xmin
    Do i=1,nbin
        xBin(i)=xx
        xx=xx+dx
        iBin(i)=0
    Enddo

! Accumulation of histo data
    nval=0
    kstr=0
    Do istr=istrbeg,nstr
        If (iref/=0.and.iType(6,istr)/=iref) Cycle
        ei=xStep(2,istr)
        erel=(ei-ebase)*au2kcalmol
        kstr=kstr+1
        ib=INT((erel-xmin)/dx)+1
        iBin(ib)=iBin(ib)+1
        nval=nval+1
    Enddo
    
! Printing 
If (iref>0) Then
    Write(iu,'(//'' Erel histogram for the structures assigned to REF'',i4)')iref
Else
    Write(iu,'(//'' Erel histogram for ALL structures'')')
Endif
Write(iu,'(  '' NumStr:'',i10  )')kstr
Write(iu,'(  '' Nvals :'',i10  )')nval
Write(iu,'(  '' Nbins :'',i10  )')nbin
Write(iu,'(  '' Xmin  :'',f10.4)')xmin
Write(iu,'(  '' Xmax  :'',f10.4)')xmax
Write(iu,'(  '' dX    :'',f10.4)')dx
nmax=0
Do i=1,nbin
    If (iBin(i)>nmax) nmax=iBin(i)
Enddo
!Write(iu,'(160x,''1'')')
!Write(iu,'(51x,30i10)')(j,j=0,9),100
Write(iu,'(/''  Bin      Xmin      Xmax   Xcenter   Count Normalized Scaled 0'',10i10)')(j,j=1,9),100
!Write(iu,'(60x,300i1)')((j,j=0,9),k=1,10),0
Write(iu,'(1x,162(''-''),''>'')')
Do i=1,nbin
    binbeg=xBin(i)
    binend=binbeg+dx
    binmid=binbeg+dx2
    ii=iBin(i)
    symb='*'
!    If (ii>250) Then
!        ii=250
!        symb='>'
!    Endif
!    iff=Dble(ii) ! unscaled
    iff=0
    If (nmax>0) iff=IDNINT(Dble(ii*100)/Dble(nmax))  ! scaled histo
    xii=Dble(ii)/Dble(nval)
    Write(iu,'(i5,3f10.5,i8,f10.4,i6,3x,''|'',<iff+1>a1)')i,binbeg,binend,binmid,ii,xii,iff,(symb,j=2,iff+1)
Enddo
Write(iu,*)

End

