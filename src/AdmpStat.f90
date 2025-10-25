Module Vars
Implicit Real(8) (A-H,O-Z)

Integer(4), parameter::MaxAt=18,MaxRef=300,MaxStr=100100     ! WARNING! MaxAt>63 can result in wrong g6 string generation! (see Am2g6 and g62am subroutines)

Character(255) FilDat,FilOut,FilRef(10),FilCan,FilAssign(10),FilSmt
Character(255), allocatable::FilInp(:)
Character(50) StrLabel0
Character(20) EtotLabel,EtotLabel0/'SCF Done:           '/
Character(20) RefLabel/'@GEO                '/,StrLabel

Integer(4) iRefLabelSkip/0/,iRefLabelFormat/1/,nFilRef/1/,iReadAssign/0/,iScaling/0/,MethodAssign/1/,iFPsubmethod/0/
Integer(4) iPrint6/1/,Numat,NA(MaxAt),nFilInp,NumatRef,NAref(MaxRef),nfact/1/,IAbuf(0:MaxAt,MaxRef)
Integer(4) nv,ne,iConn,nbndx,nref,nstr
Integer(4),allocatable::iAM(:,:),iAMx(:,:),iAMref(:,:,:)
Integer(4),allocatable::iStep(:,:),iType(:,:),Jump(:,:),iAssigned(:,:),IAbest(:,:),JumpII(:)

Real(8) C(3,Maxat),Cref(3,MaxAt,MaxRef),Cstr(3,MaxAt,MaxStr),EtotStr(MaxStr),EpsBuf(MaxRef),Freq(2,MaxRef)
Real(8) RbondMax/-1.d0/,RbondMin/-1.d0/,StrDiameter,rBondScale/1.15d0/,cscale,CoordScaleFactor
Real(8),allocatable::xStep(:,:)
Real(8),allocatable::StatUni(:,:),StatRef(:,:),RbondRef(:),DiamRef(:)
Real(8) au2kcalmol,Ebase/1.d0/,epsTresh/999999.d0/,tstep,xLTminTresh/40.d0/
Real(8),allocatable::FPref(:,:)

Character(63), allocatable:: g6ref(:),g6str(:),g6uni(:)

Logical(4) lWin,lexist

Integer(4) iref,istr,istrbeg
Real(8) Cmin0(3,MaxAt),Cref0(3,MaxAt)

Integer(4) iRecalc(3)/1000,100,10/,iSkipADMP/0/,iAddLT/0/,iRefRMSDcalc/0/
Real(8) epsrecalc/1.d-3/

Integer(4),allocatable::iRefRot(:,:),irr(:,:,:)

Integer(4) MaxFP    ! FP array size during FP assignments in GeoMetrics
Real(8) epst(maxref)
    
End module
!****************************************************************
Program AdmpStat
Use Vars

Implicit Real(8) (A-H,O-Z)

Real(8) v(8)
!
! Program AdmpStat analyzes the ADMP trajectory (converted by admp2xyz) and calculates statistics for the structures and their lifetimes
!

FilDat='admpstat.inp'
FilOut='admpstat.out'
RefLabel='@GEO'
StrLabel='!GEOSTEP'
au2kcalmol=27.21d0*23.06d0

Open(6,File=FilOut)
!Write(6,'(//''  *****  ADMP statistics v.1.0 *****''/)')
Call BuildNumber(nBuild)
Call PrintTitle(6,nBuild)
Call Timer(0,0,6,TT)
Call Timer(1,2,6,TT)

ioutopen=1
iPrintUnit=6

Inquire(DIRECTORY='c:\windows',EXIST=lWin)
If (lWin) Then
    Write(6,'(/'' Detected OS type: Windows'')')
Else
    Write(6,'(/'' Detected OS type: Linux'')')
Endif

Write(6,'(/'' Program limitations:'')')
Write(6,'('' MaxAt       '',i10)')MaxAt
Write(6,'('' MaxStr      '',i10)')MaxStr
Write(6,*)

!
! Read main command file
!
Call ReadCommandFile


! Read reference structures
Call ReadRefStr

Call Timer(1,2,6,TT)

! Read reference structures
Call ReadADMPstr

Call Timer(1,2,6,TT)
Call Timer(1,3,6,TT)

End
!****************************************************************************
Subroutine ReadADMPstr
!Use Vars, Only: nFilInp,FilInp,Numat,NA,C,StrLabel,NumatRef,nv,ne,iConn,g6str,lWin,g6ref,iAMx,FilCan,nbndx,iStep,xStep,g6uni,iType,nref,nstr
Use Vars
Use StringMod
Use Elements, Only: AMS

Implicit Real(8) (A-H,O-Z)

Character(nv*(nv-1)) g6string,g6strX
Character(14) g6type
Character(1) ch,symb,symb1
Character(10) buf
Character(255) Str,Str1
Integer(4) iStepX(10),iEminUni(1),iEmaxUni(1),iEminStr(1),iEmaxStr(1),IAtmp(MaxRef)
Real(8) xStepX(10)

Real(8),allocatable::tuni(:),Ctmp(:,:)
Integer(4),allocatable::iuni(:)

Integer(4),allocatable::iTL(:,:)
Real(8),allocatable::TL(:,:)

Integer(4),allocatable::iTSdat(:,:)
Real(8),allocatable::xTSdat(:,:),Esmt(:)

Real(8) Amass(MaxAt),CM(3),PMI(3),Paxes(3,3)

Real(8),allocatable::X(:)
Integer(4),allocatable::IXA(:)


! Read ADMP structures
Open(7,File='admp-str.dat')
Open(8,File='admp-str.g6')
!StrLabel='@GEO'
nStr=0
Write(7,'( ''File istr  istep  idump  itime      t,fs       Etot,au    nv   ne conn   g6'')')
Do i=1,nFilInp
    Open(5,File=FilInp(i))
    Do While(.not.EOF(5))
        Read(5,'(a255)')Str
        If (Len_Trim(Str)==0) Cycle
        Str1=ToUpperCase(Str)
        If (Index(Str1,Trim(StrLabel))==0) Cycle
        nstr=nstr+1
        Call SubString(Str,nSubStr,SubStr)
        iStepX(1)=iRead(SubStr(2))
        iStepX(2)=iRead(SubStr(4))
        If (INDEX(Str,'TimeStep:')>0) Then
            iStepX(3)=iRead(SubStr(6))
            xStepX(1)=fRead(SubStr(9))
            xStepX(2)=fRead(SubStr(10))
        Else
            xStepX(1)=fRead(SubStr(7))
            xStepX(2)=fRead(SubStr(9))
            iStepX(3)=IDNINT(xStepX(1)*100.d0)
        Endif
        Backspace(5)
        Call ReadNXYZ(-5,1,0,0,1,MaxAt,NA,C,NumAt,StrLabel)

        If (Numat==0) Then
            nstr=nstr-1
            Exit
        Endif
        
        If (Numat/=NumatRef) Then
            Write(6,'('' ERROR! Numat is wrong in the admp structure:'',i2)')nStr
            Stop
        Endif
    
        nv=Numat
        Call GraphConnectivity(iConn)
        ne=nbndx
!        Call AdjMatrix(ne)
        g6string=repeat(' ',Len(g6string))
        Call AM2g6(nv,iAMx,g6string)
        
        ns=Len_Trim(g6string)
        !Write(6,'(i4,2i5,2i8,f10.2,f16.6,i4,i5,i4,4x,a<ns>)')i,nStr,iStepX(1:3),xStepX(1:2),nv,ne,iConn,g6string(1:ns)
        Write(7,'(i4,2i8,2i10,f12.2,f16.6,i4,i5,i4,4x,a<ns>)')i,nStr,iStepX(1:3),xStepX(1:2),nv,ne,iConn,g6string(1:ns)
        Write(8,'(a<ns>)')g6string(1:ns)
        
        Cstr(1:3,1:Numat,nstr)=C(1:3,1:Numat)
    Enddo
    Close(5)
Enddo
Close(8)

Allocate(iStep(5,nstr),xStep(3,nStr),g6str(nstr),g6uni(nstr),iType(13,nstr))
iType=0
Rewind(7)
Read(7,*)
Do i=1,nstr
    Read(7,'(i4,2i8,2i10,f12.2,f16.6,i4,i5,i4,4x,a<ns>)')iStep(4,i),istr,iStep(1:3,i),xStep(1:2,i),iType(1:3,i),g6str(i)
Enddo
Close(7)

If (EBase==1.d0) Then
    EBase=xStep(2,1)
    Write(6,'(//'' *** Base energy is '',f16.6,''  (energy of first structure)'')')Ebase
Else
    Write(6,'(//'' *** Base energy is '',f16.6,''  (user-supplied value)'')')Ebase
Endif
    
!
! Try to canonicalize ADMP structures invoking labelg program
!
If (.not.lWin) Then
    Write(6,'(/'' Trying to invoke labelg to canonicalize the g6 strings (search in $PWD and $PATH directories only)'')')
    Inquire(File='labelg',EXIST=lexist)
    If (lexist) Then
        Call System('./labelg admp-str.g6 admp-str-canonic.g6')
    Else
        Call System('labelg admp-str.g6 admp-str-canonic.g6')
    Endif
Endif

g6type='(non-canonic) '
Inquire(File='admp-str-canonic.g6',EXIST=lexist)
If (lexist) Then
    Open(9,File='admp-str-canonic.g6')
    g6type='(canonic)     '
    Write(6,'(/'' Reading the canonicalized structure info from admp-str-canonic.g6'')')
Else
    Open(9,File='admp-str.g6')
    Write(6,'(/'' Reading the structure info from admp-str.g6'')')
Endif
nuniq=0
If (iPrint6>1) Then
    Write(6,'(/'' Raw ADMP structures'')')
    Write(6,'( ''File istr  istep  idump  itime      t,fs       Etot,au    nv   ne conn   g6 '',a14)')g6type
Endif
Do i=1,nStr
    Read(9,'(a<ns>)')g6str(i)
    If (iPrint6>1) Write(6,'(i4,2i5,2i8,f10.2,f16.6,i4,i5,i4,4x,a<ns>)')iStep(4,i),istr,iStep(1:3,i),xStep(1:2,i),iType(1:3,i),g6str(i)
Enddo


!
! Determine unique ADMP structures
!
nuniq=1
g6uni(1)=g6str(1)
iType(4,1)=1
A1: Do i=2,nstr
        Do j=1,nuniq
            If (Trim(g6str(i))==Trim(g6uni(j))) Then
                iType(4,i)=j
                Cycle A1
            Endif
        Enddo
        nuniq=nuniq+1
        g6uni(nuniq)=g6str(i)
        iType(4,i)=nuniq
Enddo A1

!
! Assign ADMP to reference structures by graph correspondence
!
A2: Do i=1,nstr
        Do j=1,nref
            If (Trim(g6str(i))==Trim(g6ref(j))) Then
                iType(5,i)=j
                Cycle A2
            Endif
        Enddo
        iType(5,i)=0
Enddo A2

! Time intervals
Do i=1,nstr-1
    xStep(3,i)=xStep(1,i+1)-xStep(1,i)
Enddo
xStep(3,nstr)=xStep(3,nstr-1)
tstep=xStep(3,nstr-1)   ! typical time step

iEminStr=MinLoc(xStep(2,1:nstr))
iEmaxStr=MaxLoc(xStep(2,1:nstr))
EminStr=xStep(2,iEminStr(1))
EmaxStr=xStep(2,iEmaxStr(1))
Write(6,'(//'' Minimum Etot value among ADMP structures:'',f16.6,''  structure'',i5)')EminStr,iEminStr
Write(6,'(  '' Maximum Etot value among ADMP structures:'',f16.6,''  structure'',i5)')EmaxStr,iEmaxStr
ErelAverStr=Sum((xStep(2,1:nstr)-EminStr)*au2kcalmol/Dble(nstr))
Write(6,'(  '' Average Erel value among ADMP structures:'',f16.2,'' kcal/mol'')')ErelAverStr

! Prnting 
Write(6,'(/''Processed ADMP structures (assigned to unique and reference structures):'')')
Write(6,'( ''File istr  istep  idump  itime      t,fs       Etot,au       dt,fs  nv   ne conn iUniq iRef Erel,kcal/mol  g6 '',a14)')g6type
Do i=1,nstr
    erel=(xStep(2,i)-EminStr)*au2kcalmol
    buf='          '
    If (i==iEminStr(1)) buf=' <-- Min E'
    If (i==iEmaxStr(1)) buf=' <-- Max E'
    Write(6,'(i4,2i6,2i7,f10.2,f16.6,f10.2,i4,i5,i4,2i5,f12.2,6x,a<ns>,a10)')iStep(4,i),i,iStep(1:3,i),xStep(1:3,i),iType(1:5,i),erel,g6str(i),buf
Enddo




!
! Allocate statistics arrays for unique and reference data
!
Allocate(StatUni(10,nuniq),StatRef(16,0:nref))
StatUni=0.d0
StatRef=0.d0

!
! Statistics for unique structures
!
TotTime=Sum(xStep(3,1:nstr))
Do i=1,nuniq
    Do j=1,nstr
        If (iType(4,j)==i) Then
            StatUni(1,i)=StatUni(1,i)+1.d0          ! Number of structures
            StatUni(2,i)=StatUni(2,i)+xStep(3,j)    ! Total time for uniq structure
            StatUni(3,i)=StatUni(3,i)+xStep(2,j)    ! Average energy
        Endif
    Enddo
Enddo
StatUni(3,1:nuniq)=StatUni(3,1:nuniq)/StatUni(1,1:nuniq)    ! Energy averaging 
Do i=1,nuniq
    eaver=StatUni(3,i)
    tmp=0.d0
    Do j=1,nstr
        If (iType(4,j)==i) Then
            tmp=tmp+(xStep(2,j)-eaver)**2
        Endif
    Enddo
    su1m1=StatUni(1,i)-1.d0
    If (su1m1>1.d-8) Then
        StatUni(4,i)=DSQRT(tmp/su1m1) ! Standard deviation of energy 
    Else
        StatUni(4,i)=0.d0
    Endif
Enddo

iEminUni=MinLoc(StatUni(3,1:nuniq))
iEmaxUni=MaxLoc(StatUni(3,1:nuniq))
EminUni=StatUni(3,iEminUni(1))
EmaxUni=StatUni(3,iEmaxUni(1))
Write(6,'(//'' Minimum Etot value among unique structures:'',f16.6,''  structure'',i5)')EminUni,iEminUni
Write(6,'('' Maximum Etot value among unique structures:'',f16.6,''  structure'',i5)')EmaxUni,iEmaxUni

Allocate(tuni(nuniq),iuni(nuniq))
tuni(1:nuniq)=StatUni(2,1:nuniq)
Call HeapSort1(nuniq,tuni,iuni)

Open(12,File='uniq-str-sorted.xyz')
Write(6,'(/'' Unique structures statistics (SORTED by dt). Total time, fs ='',f10.2)')TotTime
Write(6,'( ''    i Uniq ADMP Nstr    dt,fs      dt,%    Sum(dt,%)    EtotAver    +-   SD         ErelUni,kcal/mol'')')
Write(12,'(11x,''i Uniq ADMP Nstr    dt,fs      dt,%    Sum(dt,%)    EtotAver    +-   SD         ErelUni,kcal/mol'')')
tpdt=0.d0
ii=0
Do j=nuniq,1,-1
    ii=ii+1
    i=iuni(j)
    erel=(StatUni(3,i)-EminUni)*au2kcalmol
    dt=StatUni(2,i)                             ! Life time for unique str
    pdt=dt/TotTime*100.d0                       ! Life time in percents
    tpdt=tpdt+pdt                               ! Integrated life time for best structure
    Do iadmp=1,nstr
        If (iType(4,iadmp)==i) Exit
    Enddo
    Write(6,'(4i5,3f10.2,2f16.6,f12.2)')ii,i,iadmp,INT(StatUni(1,i)),dt,pdt,tpdt,StatUni(3:4,i),erel
    C(1:3,1:Numat)=Cstr(1:3,1:Numat,iadmp)
    Write(12,'(''@GeoUni'',4i5,3f10.2,2f16.6,f12.2)')ii,i,iadmp,INT(StatUni(1,i)),dt,pdt,tpdt,StatUni(3:4,i),erel
    Call PrintNXYZnew(12,4,Numat,NA,C,16,6,0,Str)
Enddo


!
! Statistics for reference structures
!
StatRef=0.d0
Do i=0,nref
    Do j=1,nstr
        If (iType(5,j)==i) Then
            StatRef(1,i)=StatRef(1,i)+1.d0          ! Number of structures
            StatRef(2,i)=StatRef(2,i)+xStep(3,j)    ! Total time for uniq structure
            StatRef(3,i)=StatRef(3,i)+xStep(2,j)    ! Average energy
        Endif
    Enddo
    If (StatRef(1,i)>0.d0) StatRef(3,i)=StatRef(3,i)/StatRef(1,i)    ! Energy averaging 
Enddo

Do i=0,nref
    If (StatRef(1,i)<1.d0) Cycle
    eaver=StatRef(3,i)
    tmp=0.d0
    Do j=1,nstr
        If (iType(4,j)==i) Then
            tmp=tmp+(xStep(2,j)-eaver)**2
        Endif
    Enddo
    sr1m1=StatRef(1,i)-1.d0
    If (sr1m1>1.d-8) Then
        StatRef(4,i)=DSQRT(tmp/sr1m1) ! Standard deviation of energy 
    Else
        StatRef(4,i)=0.d0
    Endif
Enddo

Write(6,'(//'' Statistics for the reference structures identified. Total time, fs ='',f10.2)')TotTime
Write(6,'(  ''       i    Nstr     dt,fs         dt,%       EtotAver      +-    SD'')')
Do i=0,nref
    inum=INT(StatRef(1,i))
    If (inum==0) Cycle
    Write(6,'(2i8,2f12.2,2f16.6)')i,inum,StatRef(2,i),StatRef(2,i)/TotTime*100.d0,StatRef(3:4,i)
Enddo

!
! Assign ADMP to REF structures using the Iterative Closest Point algorithm (ICP) 
!
nfact=1
Do i=2,Numat
    nfact=nfact*i
Enddo

! First, reduce reference structures to CM and principal axes
Open(12,File='ref-str.xyz')
Write(12,'(''! REF structures - reduced to CM and principal axes''/)')
Amass(1:Numat)=AMS(NA(1:Numat))
Do iref=1,nref
    C(1:3,1:Numat)=Cref(1:3,1:Numat,iref)
    Call InertiaNew(2,1,1,Numat,C,Amass,CM,PMI,PAxes,kRotTyp)
    Cref(1:3,1:Numat,iref)=C(1:3,1:Numat)
    Write(Str,'(''@georef'',i0.4)')iref
    Call PrintStr(12,Numat,C,Str)
Enddo    
Close(12)

! Jumps and permutations arrays
Allocate(Jump(0:nref,0:nref),iAssigned(2,nstr),IAbest(Numat,nstr),JumpII(nref))
iAssigned=0

! Permutations of REF structures which are rotations
Allocate(iRefRot(-1:0,nref))  
Allocate(irr(Numat,100,nref))   ! Max first index can be reduced to 60 (max number of rotations for Ih symmetry group). 100 is just in case if some atoms coincide due to numerical errors
iRefRot=0
irr=0
!!! temp debug - permutations
Open(76,File='ref-rot.tmp')
!!!


! Read old assignment if any
If (iReadAssign>0) Then
    Open(27,File='str-admp2ref4-tmp.dat')
    Do i=1,iReadAssign
        Open(25,File=FilAssign(i))
        Call ReadAssignment
        Close(25)
        Rewind(27)
    Enddo
Endif

Open(20,File='str-admp2ref2.dat')
Open(12,File='str-admp2ref.xyz')
Open(24,File='str-admp2ref4.dat')
Open(25,File='str-admp2ref5.dat')
Open(17,File='str-admp2-ts.xyz')

Open(41,File='str-admp2-FP.dat')

Allocate(Ctmp(3,Numat))

lim1=iRecalc(1)
lim2=iRecalc(2)
mref=Min(iRecalc(3),nref)
If (lim1<=0) lim1=500
If (lim2<=0) lim2=50
If (mref<=0) mref=10

Write(6,'(//'' [ADMPassign] Assignment of ADMP structures to RefStr using the ICP algorithm. eps is the RMSD estimate between ADMPstr and REFstr after ICP (eps=DSQRT(SSE/Numat))'')')
Write(6,'(  '' See file str-admp2ref.xyz for details and assigned geometries'')')
Write(6,'(  '' Full  ICP calculation periods (RECALC=(LIM1,LIM2,MREF) command):'',3i5)')lim1,lim2,mref
Write(6,'(  '' (i.e. ICP for all REFs - every'',i5,'' steps, ICP for the best structures only - every'',i5,'' steps. The number of best structures is'',i5)')lim1,lim2,mref
Write(12,'( '' Assignment of ADMP structures to RefStr using the ICP algorithm. eps is the RMSD estimate between ADMPstr and REFstr after ICP (eps=DSQRT(SSE/Numat))'')')
Write(20,'(''  ADMP iresold  REF   LTS        eps     Time,fs       LT,fs'')')
Write(6,'( ''  ADMP      t,fs   REF   Iter      Eps      REF(final)      Etot          Erel       RMSD   CoordScaling  Atom correspondence ADMPstr->REFstr(1234...)     RminMin  RminMax   Diam  AverRefBond  dRefBond'')')
Write(24,'(''  ADMP      t,fs   REF   Iter      Eps      REF(final)      Etot          Erel'')')
Write(24,'(''  iref       nf        nftot    Atom mapping ADMP->REF(1234...)             eps'')')
Jump=0
iter=0
iresold=-1
nts=0
istrbeg=iSkipADMP+1
istrold=istrbeg-1
Call Timer(2,0,6,CalcTime)
Do istr=istrbeg,nstr        
    C(1:3,1:Numat)=Cstr(1:3,1:Numat,istr)
!    Call ICPassign(Numat,C,Amass,ires,iter,eps)
    If (iAssigned(1,istr)==0) Then
        If (MethodAssign==1) Then
            If (istr==istrbeg.or.istr/lim1*lim1==istr) Then
                Call ICPassignFull(0,nref,Numat,C,Amass,ires,eps)
                !Write(25,'(i6,i4,<nref>f10.4)')istr,ires,epst(1:nref)
            ElseIf (istr/lim2*lim2==istr) Then
                Call ICPassignFull(1,mref,Numat,C,Amass,ires,eps)
            Else
                Call ICPassignShort(mref,Numat,C,Amass,ires,eps)
            Endif
            Write(25,'(i6,i4,<nref>f10.4)')istr,ires,epst(1:nref)
            IAbest(1:Numat,istr)=IAbuf(1:Numat,1)
        Else
            Call FPassign(1,nref,Numat,C,ires,eps)
            Cref0(1:3,1:Numat)=Cref(1:3,1:Numat,ires)
            Do i=1,Numat
                Cmin0(1:3,i)=C(1:3,i)
            Enddo
            !Call LSQR(Numat,Cref0,Cmin0,C,eps1)
        Endif
    Else
        ires=iAssigned(1,istr)
        itmp=iAssigned(2,istr)
        eps=Dble(itmp)*1.d-4
        Read(27,'(i6,i4,<Numat>i3)')is,ir,IAbest(1:Numat,istr)
        
        If (is/=istr.or.ir/=ires) IAbest(1:Numat,istr)=-IAbest(1:Numat,istr)
        
        Cref0(1:3,1:Numat)=Cref(1:3,1:Numat,ires)
        Do i=1,Numat
            ii=IAbest(i,istr)
            Cmin0(1:3,i)=C(1:3,ii)
        Enddo
        Call LSQR(Numat,Cref0,Cmin0,C,eps1)
        
!        If (istr==4869) Then
!            Call PrintStr(17,Numat,Cref0,'@geo-cref0-4869')
!            Call PrintStr(17,Numat,C,'@geo-cmin0-4869')
!        Endif
!        If (istr==4870) Then
!            Call PrintStr(17,Numat,Cref0,'@geo-cref0-4870')
!            Call PrintStr(17,Numat,C,'@geo-cmin0-4870')
!        Endif
        
        Cmin0(1:3,1:Numat)=C(1:3,1:Numat)
        C(1:3,1:Numat)=Cstr(1:3,1:Numat,istr)
        If (DABS(eps1-eps)/eps>epsrecalc) Then
            Write(6,'(//'' ERROR! During recalculation, can not properly restore the structure'',i8,'' REF ='',i8)')istr,ires
            Write(6,'(  '' Current value of recognition eps'',f12.6,''  is different  from  the source value:'',f12.6,'' rel.diff:'',f12.6)')eps1,eps,DABS(eps1-eps)/eps
            Write(6,'(  '' You can weak the accuracy requirements increasing EPSRECALC option, current value:'',f12.6/)')epsrecalc
            Stop
        Endif
        CoordScaleFactor=cscale
        
    Endif
    Call Timer(2,1,6,CalcTime)
    tt=xStep(1,istr)    !(istr-1)*tstep
    ei=xStep(2,istr)
    erel=(ei-ebase)*au2kcalmol
    Call CalcSSE(Numat,Cmin0,Cref0,sse)
    Call CalcDist(Numat,Cmin0,rmin,rmax,diam)
    Call CalcBonds(Numat,Cmin0,ires,nb,arbx)   ! Number of bonds and bond lengths
    drbx=arbx-RbondRef(ires)
    
    symb1=' '
    If (ires/=iresold) Then
        symb1='*'
    Else
        If (iRefRot(-1,ires)==0) Call RefPerm(ires)
        Call CheckRot(istr,ires,irot)
        If (irot==0) symb1='X'
    Endif
    
    lx=2
    If (Numat<=15) lx=3*(15-Numat)+2
    If (eps<=epsTresh) Then
!        Write(6, '('' Str'',i5,'' assigned to Ref'',i5,''  iter:'',i5,''  eps:'',f12.4,i6)')istr,ires,iter,eps,ires
        Write(6, '(i6,f10.2,i6,i6,f12.4,i6,f20.6,f12.2,2f12.4,1x,a1,<Numat>i3,<lx>x,5f10.4)') &
            istr,tt,ires,iter,eps,ires,xStep(2,istr),erel,sse,CoordScaleFactor,symb1,IAbest(1:Numat,istr),rmin,rmax,diam,arbx,drbx
        Write(24,'(i6,f10.2,i6,i6,f12.4,i6,f20.6,f12.2,6x,''time,s:'',f10.2)')istr,tt,ires,iter,eps,ires,xStep(2,istr),erel,CalcTime
        Write(12,'(''@str'',i0.6,'' assigned to Ref'',i0.4,''  iter:'',i5,''  eps:'',f12.4,i6\)')istr,ires,iter,eps,ires
        If (iresold/=ires) Then
            If (istr>1) Write(20,'(/4i6,3f12.2/)')0,0,0,(istr-istrold-1),0.d0,0.d0,Dble(istr-istrold-1)*tstep
            it=0
            istrold=istr
        Endif
        it=it+1
        Write(20,'(4i6,3f12.2)')istr,iresold,ires,it,eps,Dble(istr-1)*tstep,Dble(it)*tstep
    Else
        Write(6, '(i6,f10.2,i6,i6,f12.4,i6,f20.6,f12.2,3f12.4,2x,<Numat>i3,2x,2f10.4,''  *** eps>'',f6.3,'' => REF set to 0 (unrecognized)'')')istr,tt,ires,iter,eps,0,xStep(2,istr),erel,epsTresh,sse,CoordScaleFactor,IAbest(1:Numat,istr),rmin,rmax
        Write(24,'(i6,f10.2,i6,i6,f12.4,i6,f20.6,f12.2,''  *** eps>'',f6.3,'' => REF set to 0 (unrecognized)'',6x,''time,s:'',f10.2)')istr,tt,ires,iter,eps,0,xStep(2,istr),erel,epsTresh,CalcTime
        !Write(6, '('' Str'',i5,'' assigned to Ref'',i5,''  iter:'',i5,''  eps:'',f12.4,i6,''  *** eps>'',f6.3,'' => REF set to 0 (unrecognized)'')')istr,ires,iter,eps,0,epsTresh
        Write(12,'(''@str'',i0.6,'' assigned to Ref'',i0.4,''  iter:'',i5,''  eps:'',f12.4,i6,''  *** eps>'',f6.3,'' => REF set to 0 (unrecognized)''\)')istr,ires,iter,eps,0,epsTresh
        ires=0        
    Endif
!    Write(6,'(3i5,f12.4)')istr,ires,iter,eps
!    C(1:3,1:Numat)=Cmin0(1:3,1:Numat)
    Call PrintStr(12,Numat,Cmin0,' ')
!    C(1:3,1:Numat)=Cref0(1:3,1:Numat)
!    Call PrintStr(12,Numat,C,' BestStrREF:')
    iType(6,istr)=ires
    iType(7,istr)=iter
    iType(8,istr)=IDNINT(eps*1.d6)
    
    iType(9,istr)=IDNINT(rmin*1.d6)
    iType(10,istr)=IDNINT(rmax*1.d6)
    iType(11,istr)=IDNINT(diam*1.d6)
    iType(12,istr)=IDNINT(arbx*1.d6)
    iType(13,istr)=IDNINT(drbx*1.d6)
    
    
    If (istr==istrbeg) iresold=ires
    If (istr==1) Then
        idiff=0
    Else
        idiff=Sum(IABS(IAbest(1:Numat,istr)-IAbest(1:Numat,istr-1)))
    Endif
    If (ires/=iresold.or.idiff/=0) Jump(iresold,ires)=Jump(iresold,ires)+1      ! Jumps
    If (ires/=iresold .and. MethodAssign==1) Then                                                     ! Transition states
        nts=nts+1

        Do i=1,Numat
            If (i/2*2==i) NA(i)=13  ! Make odd atoms Al to make bonds more recognizable
        Enddo
        
        Ctmp(1:3,1:Numat)=Cmin0(1:3,1:Numat)
        Call ICPassignSingle(iresold,Numat,Ctmp,Amass,iresold,eps,IAtmp)
        Do i=1,Numat
            ii=IAtmp(i)
            Ctmp(1:3,ii)=Cref(1:3,i,iresold)
        Enddo
        Call LSQR(Numat,Cmin0,Ctmp,C,epstmpold)
        Write(17,'(''@str'',i0.3,''-ref'',i0.6,'' structure  REF1 '',f10.4\)')nts,iresold,epstmpold
        Call PrintStr(17,Numat,C,' ')
        
        Do ix=1,19
            xx=Dble(ix)*0.05d0
            Ctmp(1:3,1:Numat)=(1.d0-xx)*C(1:3,1:Numat)+xx*Cmin0(1:3,1:Numat)
            Write(17,'(''@str'',i0.3,''-intrpl'',i0.2\)')nts,ix
            Call PrintStr(17,Numat,Ctmp,' ')
        Enddo
        
        !Do i=1,Numat
        !    ii=IAbest(i,istr-1)
        !    C(1:3,i)=Cstr(1:3,ii,istr-1)
        !Enddo
        !Do i=1,Numat
        !    ii=IAtmp(i)
        !    Ctmp(1:3,ii)=C(1:3,i)
        !Enddo
        !Call MassCenter(Numat,Ctmp)
        !Call LSQR(Numat,Cmin0,Ctmp,C,epstmp)
        !Write(17,'(''@str-TS-1'',i0.4,'' istr'',i0.8,'' TS: REF'',i0.4,'' -> REF'',i0.4,2x,''Etot:'',f15.6,2x,''Erel:'',f8.2\)')nts,istr,iresold,ires,xStep(2,istr),erel
        !Call PrintStr(17,Numat,C,' ')
        
        Write(17,'(''@str'',i0.3,''-TS istr'',i0.8,'' TS: REF'',i0.4,'' -> REF'',i0.4,2x,''Etot:'',f15.6,2x,''Erel:'',f10.4\)')nts,istr,iresold,ires,xStep(2,istr),erel
        Call CalcSSE(Numat,Cref0,Cmin0,epstmp)
        Call PrintStr(17,Numat,Cmin0,' ')
        
        Do ix=1,19
            xx=Dble(ix)*0.05d0
            Ctmp(1:3,1:Numat)=(1.d0-xx)*Cmin0(1:3,1:Numat)+xx*Cref0(1:3,1:Numat)
            Write(17,'(''@str'',i0.3,''-intrpl'',i0.2\)')nts,ix+20
            Call PrintStr(17,Numat,Ctmp,' ')
        Enddo
        
        Write(17,'(''@str'',i0.3,''-ref'',i0.6,'' structure  REF2 '',f10.4\)')nts,ires,epstmp
        Call PrintStr(17,Numat,Cref0,' ')
        
        Do i=1,Numat
            If (i/2*2==i) NA(i)=12  ! Make odd atoms Al to make bonds more recognizable
        Enddo

    
    Endif
    iresold=ires
Enddo
Close(12)
Close(17)
Close(24)
If (iReadAssign>0) Close(27)
Rewind(20)

Open(22,File='str-admp2ref3.dat')
Read(20,'(a255)')Str
Write(22,'(a255)')Str
istr1=0
istr2=0
irefold=0
it=0
Do While (.not.EOF(20))
    Read(20,'(a255)')Str
    If (Len_Trim(Str)==0) Cycle
    Read(Str,'(4i6,3f12.2)')istr,itmp,iref,itmp2,xtmp1,xtmp2
    If (istr==0) Then
        !Write(22,'(4i6,3f12.2)')istr,itmp,iref,itmp2,xtmp1,xtmp2
        Write(22,'(4i6,3f12.2)')istr1,istr2,irefold,it,0.d0,xtmp2
    Endif
    If (itmp==0) Then
        it=0
        istr1=istr
    Else
        it=it+1
        irefold=iref
        istr2=istr
        xtmp2old=xtmp2
    Endif
Enddo
Write(22,'(4i6,3f12.2)')istr1,istr2,irefold,it,0.d0,xtmp2
Close(22)
Close(20)

! iType(1,:) - 
! iType(2,:) - 
! iType(3,:) - 
! iType(4,:) - 
! iType(5,:) - 
! iType(6,:) - REF str identified by ICP
! iType(7,:) - iterations in ICP 
! iType(8,:) - eps of ICP (IDNINT(eps*1.d6)
!    iType(9,istr)=IDNINT(rmin*1.d6)
!    iType(10,istr)=IDNINT(rmax*1.d6)
!    iType(11,istr)=IDNINT(diam*1.d6)
!    iType(12,istr)=IDNINT(arbx*1.d6)
!    iType(13,istr)=IDNINT(drbx*1.d6)


! Statistics for ADMP to REF after ICP assignment
Do i=0,nref
    Do j=istrbeg,nstr
        If (iType(6,j)==i) Then
            StatRef(6,i)=StatRef(6,i)+1.d0          ! Number of structures
            StatRef(7,i)=StatRef(7,i)+xStep(3,j)    ! Total time for uniq structure
            StatRef(8,i)=StatRef(8,i)+xStep(2,j)    ! Average energy
            StatRef(9,i)=StatRef(9,i)+Dble(iType(8,j))*1.d-6    ! Average eps (SSE for given REF after ICP assignment)
            
            StatRef(12,i)=StatRef(12,i)+Dble(iType(9,j))*1.d-6    ! RminMin
            StatRef(13,i)=StatRef(13,i)+Dble(iType(10,j))*1.d-6    ! RminMax
            StatRef(14,i)=StatRef(14,i)+Dble(iType(11,j))*1.d-6    ! Diam
            StatRef(15,i)=StatRef(15,i)+Dble(iType(12,j))*1.d-6    ! AverBondRef
            StatRef(16,i)=StatRef(16,i)+Dble(iType(13,j))*1.d-6    ! dAverBondRef
            
        Endif
    Enddo
    sr6=StatRef(6,i)
    If (sr6>0.d0) Then
        StatRef(8,i)=StatRef(8,i)/sr6   ! Energy averaging 
        StatRef(9,i)=StatRef(9,i)/sr6   ! eps averaging 
        StatRef(12:16,i)=StatRef(12:16,i)/sr6  ! 
    Endif
Enddo

! Standard deviations for Etot and eps
Do i=0,nref
    Do j=istrbeg,nstr
        If (iType(6,j)==i) Then
            StatRef(10,i)=StatRef(10,i)+((xStep(2,j)-StatRef(8,i))**2)/(StatRef(6,i))                ! RMS for energy
            StatRef(11,i)=StatRef(11,i)+((Dble(iType(8,j))*1.d-6-StatRef(9,i))**2)/(StatRef(6,i))    ! RMS for eps
        Endif
    Enddo
Enddo
StatRef(10:11,1:nref)=DSQRT(StatRef(10:11,1:nref))

If (iPrint6>1) Then
    Write(6,'(//'' Statistics for the reference structures identified by ICP algorithm. Total time, fs ='',f10.2)')TotTime
    Write(6,'( /''    i  Nstr    dt,fs      dt,%       EtotAver    +-     RMS             eps      +-     RMS'')')
    Do i=0,nref
        inum=INT(StatRef(6,i))
        Write(6,'(2i5,2f10.2,4f16.6)')i,inum,StatRef(7,i),StatRef(7,i)/TotTime*100.d0,StatRef(8,i),StatRef(10,i),StatRef(9,i),StatRef(11,i)
    Enddo
Endif

! Sorted and condensed printing of ADMP/REF after ICP
Deallocate(tuni,iuni)
Allocate(tuni(nref),iuni(nref))
Do i=1,nref
    tuni(i)=StatRef(6,i)
Enddo
Call HeapSort1(nref,tuni,iuni)
Write(6,'(//'' Statistics for the SORTED reference structures identified by ICP algorithm. Total time, fs ='',f10.2,''  Base energy for Erel ='',f16.6)')TotTime,Ebase
Write(6,'(  '' Nstr - number of trajectory points assigned to REF structure; dt - their total duration; EtotAver - their average energy +- RMS;'')')
Write(6,'(  '' eps - their average deviation (SSE) from REF structure +- RMS; Erel - their energy relatively to the base energy Ebase'')')
Write(6,'( /''    i  Ref  Nstr     dt,fs      dt,%   Sum(dt,%)     EtotAver    +-     RMS             eps      +-     RMS      Erel,kcal/mol  +- RMS'',20x,''RminMin   RminMax      Diam AverRefBond dRefBond  DiamRef'')')
ttim=0.d0
k=0
Open(19,File='str-admp-histo1.dat')
Open(29,File='str-admp-histo2.dat')
Open(39,File='str-admp-histoE.dat')
Do j=nref,1,-1
    i=iuni(j)
    inum=INT(StatRef(6,i))
    If (inum==0) Cycle
    tim=StatRef(7,i)/TotTime*100.d0
    ttim=ttim+tim
    k=k+1
    EtotAver=StatRef(8,i)
    EtotSD=StatRef(10,i)
    Erel=(EtotAver-Ebase)*au2kcalmol
    Write(6,'(2i5,i6,3f10.2,4f16.6,2f12.2,17x,6f10.4)')k,i,inum,StatRef(7,i),tim,ttim,EtotAver,EtotSD,StatRef(9,i),StatRef(11,i),Erel,EtotSD*au2kcalmol,StatRef(12:16,i),DiamRef(i)
    Do mode=1,2
        iu=10*mode+9
        If (j/=nref) Write(iu,'(200(''#''))')
        Call HistoRef(iu,mode,i)
        Call Histo(iu,mode,i)
    Enddo
    If (j/=nref) Write(iu,'(200(''#''))')
    Call HistoE(39,i)
Enddo
Call HistoE(39,0)
Close(19)
Close(29)
Close(39)

! Statistics for jumps between REF structures
jjump=Sum(Jump)
Write(6,'(//'' ***** Jumps statistics (including jumps between the same-type REF structures)'')')
Write(6,'(  '' Statistics for jumps between REF structures identified by ICP algorithm. Total number of jumps ='',i5)')jjump
!Write(6,'(6x,<nref+1>i3)')(j,j=0,nref)
!Do i=0,nref
!    Write(6,'(i3,3x,<nref+1>i3)')i,Jump(i,0:nref)
!Enddo
!stop
nref2=(nref+1)*(nref+1)
Deallocate(tuni,iuni)
Allocate(tuni(nref2),iuni(nref2))
ij=0
Do i=0,nref
    Do j=0,nref
        ij=ij+1
        tuni(ij)=Dble(Jump(i,j))
    Enddo
Enddo
Call HeapSort1(nref2,tuni,iuni)
ij=0
k=0
txj=0.d0
nr1=nref+1
Do ii=nref2,1,-1
    jmp=IDNINT(tuni(ii))
    If (jmp==0) Cycle
    ij=iuni(ii)

    i=ij/nr1+1
    j=ij-(i-1)*nr1
    If (j==0) Then
        i=i-1
        j=nr1
    Endif
    
    k=k+1
    xjump=Dble(jmp)/Dble(jjump)
    txj=txj+xjump
    Write(6,'(i5''  Ref'',i3,''  -->  Ref'',i3,''   jumps:'',i5,''   Pij:'',f10.4,''   SumPij:'',f10.4)')k,i-1,j-1,jmp,xjump,txj
Enddo

njmp=Sum(Jump(1:nref,1:nref))
Write(6,'(//'' Jumps sorted by REF structures (among identified structures only). Number of jumps among identified structures N ='',i6)')njmp
Write(6,'(  ''   No.   k  REFi->REFj Nij   Nij/Ni      Sum    Prob=Nij/N     Graph for drawing at https://programforyou.ru/graph-redactor'')')
ij=0
ttjmp=Dble(njmp)
Do i=1,nref
    tuni(1:nref)=Dble(Jump(i,1:nref))
    ttuni=Sum(tuni(1:nref))
    Call HeapSort1(nref,tuni,iuni)
    xtot=0.d0
    jj=0
    Do j=nref,1,-1
        jmp=IDNINT(tuni(j))
        If (jmp==0) Cycle
        ij=ij+1
        jj=jj+1
        symb=' '
        If (j==nref) symb='*'
        xjmp=tuni(j)/ttuni
        xtot=xtot+xjmp
        ttij=tuni(j)/ttjmp
        Write(6,'(5i5,3f10.4,2x,a1,5x,i3,'' -> '',i3,f6.3)')ij,jj,i,iuni(j),jmp,xjmp,xtot,ttij,symb,i,iuni(j),ttij
    Enddo
Enddo

! Jumps statistics between different refs
njmpold=njmp
Do i=1,nref
    JumpII(i)=Jump(i,i)
    Jump(i,i)=0
Enddo
njmp=Sum(Jump(1:nref,1:nref))
Write(6,'(//'' Jumps sorted by REF structures (between DIFFERENT REFs only). Number of jumps among identified structures N ='',i6)')njmp
Write(6,'(  ''   No.   k  REFi->REFj Nij   Nij/Ni      Sum    Prob=Nij/N     Graph for drawing at https://programforyou.ru/graph-redactor'')')
ij=0
ttjmp=Dble(njmp)
Do i=1,nref
    tuni(1:nref)=Dble(Jump(i,1:nref))
    ttuni=Sum(tuni(1:nref))
    Call HeapSort1(nref,tuni,iuni)
    xtot=0.d0
    jj=0
    Do j=nref,1,-1
        jmp=IDNINT(tuni(j))
        If (jmp==0) Cycle
        ij=ij+1
        jj=jj+1
        symb=' '
        If (j==nref) symb='*'
        xjmp=tuni(j)/ttuni
        xtot=xtot+xjmp
        ttij=tuni(j)/ttjmp
        Write(6,'(5i5,3f10.4,2x,a1,5x,i3,'' -> '',i3,f6.3)')ij,jj,i,iuni(j),jmp,xjmp,xtot,ttij,symb,i,iuni(j),ttij
    Enddo
Enddo
Do i=1,nref
    Jump(i,i)=JumpII(i)
Enddo
njmp=njmpold


Call LifeTimes1(0,0.d0)

Call LifeTimes2(0,0.d0)

If (iAddLT>0) Then
    Do i=1,10
        xLTminTresh=Dble(i)*10.d0    
        Call LifeTimes1(1,xLTminTresh)
    Enddo
Endif


!
! Additional analysis of TSs
!
Allocate(iTSdat(7,nts),xTSdat(2,nts))
its=0
Open(17,File='str-admp2-ts.xyz',Action='READ')
Do While(.not.EOF(17))
    Read(17,'(a255)')Str
    If (Index(Str,'TS:')==0) Cycle
    Read(Str,'(4x,i3,8x,i8,8x,i4,7x,i4,7x,f15.6,7x,f10.4)')jts,istrts,iref1,iref2,etotts,erelts
    its=its+1
    iTSdat(1,its)=iref1
    iTSdat(2,its)=iref2
    iTSdat(3,its)=istrts
    iTSdat(4,its)=its
    xTSdat(1,its)=etotts
    xTSdat(2,its)=erelts
Enddo
Close(17)

! Read smoothed energies if present
Allocate(Esmt(nstr))
FilSmt=FilInp(1)
Call FileExtension(FilSmt,'.smt')
INQUIRE(File=FilSmt,EXIST=lexist)
If (lexist) Then
    Open(17,File=FilSmt)
    Read(17,*)
    Do i=1,nstr
        Read(17,'(35x,f15.5)')Esmt(i)
    Enddo
    Close(17)
Else 
    Esmt=0.d0
Endif
!
Open(17,File='str-admp2-ts.dat')
Write(17,'('' Energies of all'',i4,'' transition states located at the ADMP trajectory (between different REF structures only), see file str-admp2-ts.xyz'')')nts
Write(17,'('' Erel - energy of TS relative to BaseEnergy. No.TS - ordering number of TS in the file str-admp2-ts.xyz'')')
Do i=2,nref
    Do j=1,i-1
        etottsmin=999999.d0
        ereltsmin=999999.d0
        etottsaver=0.d0
        ereltsaver=0.d0
        etm12=999999.d0
        erm12=999999.d0
        etm21=999999.d0
        erm21=999999.d0
        eta12=0.d0
        eta21=0.d0
        era12=0.d0
        era21=0.d0
        esa=0.d0
        esa12=0.d0
        esa21=0.d0
        mts=0
        m12=0
        m21=0
        Do it=1,nts
            ir1=iTSdat(1,it)
            ir2=iTSdat(2,it)
            If(ir1==i.and.ir2==j.or.ir1==j.and.ir2==i) Then
                mts=mts+1
                et=xTSdat(1,it)
                er=xTSdat(2,it)
                If (et<etottsmin) Then
                    etottsmin=et
                    mmin=mts
                Endif
                If (er<ereltsmin) ereltsmin=er
                etottsaver=etottsaver+et
                ereltsaver=ereltsaver+er
                ii=ii+1
                If (mts==1) Write(17,'(/''   i    Ref1    Ref2  No.TS Ref1   Ref2  iPoint     Etot, a.u.      Erel, kcal/mol    Erel-Esmooth, kcal/mol'')')
                itsd3=iTSdat(3,it)
                es=xTSdat(2,it)-Esmt(itsd3)
                esa=esa+es
                Write(17,'(i4,i6,''  <=>'',i4,i7,i4,''  ->'',i4,i8,2f16.4,4x,f16.4)')mts,i,j,it,iTSdat(1:3,it),xTSdat(1:2,it),es
                If (ir1<ir2) Then
                    If (et<etm12) Then
                        etm12=et
                        mmin12=mts
                    Endif
                    If (er<erm12) erm12=er
                    eta12=eta12+et
                    era12=era12+er
                    esa12=esa12+es
                    m12=m12+1
                Else
                    If (et<etm21) Then
                        etm21=et
                        mmin21=mts
                    Endif
                    If (er<erm21) erm21=er
                    eta21=eta21+et
                    era21=era21+er
                    esa21=esa21+es
                    m21=m21+1
                Endif
            Endif
        Enddo
        If (mts==0) Cycle
        etottsaver=etottsaver/Dble(mts)
        ereltsaver=ereltsaver/Dble(mts)
        esa=esa/Dble(mts)
        Write(17,'(/3x,''TS:'',i4,''  <=>'',i4,2x,''numTS ='',i5,2x,'' EtotMin,EtotAver ='',2f16.8,4x,''ErelMin,ErelAver ='',2f12.2,5x,''imin ='',i4,6x,''(Er-Esm)aver ='',f10.4)')i,j,mts,etottsmin,etottsaver,ereltsmin,ereltsaver,mmin,esa
        If (m12>0) Then
            eta12=eta12/Dble(m12)
            era12=era12/Dble(m12)
            esa12=esa12/Dble(m12)
            Write(17,'(3x,''TS:'',i4,''  -->'',i4,2x,''numTS ='',i5,2x,'' EtotMin,EtotAver ='',2f16.8,4x,''ErelMin,ErelAver ='',2f12.2,5x,''imin ='',i4,6x,''(Er-Esm)aver =''f10.4)')j,i,m12,etm12,eta12,erm12,era12,mmin12,esa12
        Endif
        If (m21>0) Then
            eta21=eta21/Dble(m21)
            era21=era21/Dble(m21)
            esa21=esa21/Dble(m21)
            Write(17,'(3x,''TS:'',i4,''  -->'',i4,2x,''numTS ='',i5,2x,'' EtotMin,EtotAver ='',2f16.8,4x,''ErelMin,ErelAver ='',2f12.2,5x,''imin ='',i4,6x,''(Er-Esm)aver ='',f10.4)')i,j,m21,etm21,eta21,erm21,era21,mmin21,esa21
        Endif
        Write(17,*)        
    Enddo
Enddo
Close(17)

!
! Dump interatomic distances for further cluster analysis
!
iu=17
nn=Numat*(Numat-1)/2
Allocate(X(nn),IXA(nn))
Open(iu,File='str-admp-clust1.txt')
Do ia=1,Numat-1
    Do ja=ia+1,Numat
        Write(iu,'(3x,''r'',2i0.2\)')ia,ja
    Enddo
Enddo
Write(iu,*)
Do istr=istrbeg,nstr
    C(1:3,1:Numat)=Cstr(1:3,1:Numat,istr)
    k=0
    Do ia=1,Numat-1
        Do ja=ia+1,Numat
            k=k+1
            rij=Distance(ia,ja,Numat,C)
            X(k)=rij
        Enddo
    Enddo
    Call HeapSort1(nn,X,IXA)
    Write(iu,'(<nn>f8.4)')X(1:nn)
Enddo
Close(iu)            
        
End
!******************************************************************************
Subroutine LifeTimes1(irun,xLTminTresh)
Use Vars, xLTMinTresh1=>xLTminTresh

Implicit Real(8) (A-H,O-Z)
   
Real(8),allocatable::TL(:,:),tuni(:)
Integer(4),allocatable::iuni(:)
Integer(8),allocatable::iTL(:,:)
Integer(8) lti

Character(255) File14,File16,File18,File15

If (irun==0) Then
    File14='str-lifetimes1a.dat'
    File15='str-lifetimes1d.dat'
    File16='str-lifetimes1b.dat'
    File18='str-lifetimes1c.dat'
ElseIf (irun>0) Then
    File14='str-lifetimes1a_x.dat'
    File15='str-lifetimes1d_x.dat'
    File16='str-lifetimes1b_x.dat'
    File18='str-lifetimes1c_x.dat'
    Write(6,'(//'' *** Repeat the lifetime statistics using the minimum LifeTime treshold'',f10.2,'' fs'')')xLTMinTresh
Endif

! Lifetime statistics
Allocate(TL(4,nref),iTL(10,nref))
!tstep=xStep(3,nstr-1)   ! typical time step
iTL=0
TL=0.d0
tti=0.d0
Open(14,File=File14)
Write(6,'(//'' Lifetime distributions are dumped  to the file str-ref-lifetimes.dat.'')')
Write(6,'(  '' Lifetime statistics    are written to the file str-ref-lifetimes3.dat (for additional analysis on tstep influence).'')')
Write(14,'(''  Ref  Int   Steps Lifetime,fs  ibeg    iend'')')
Open(18,File=File18)
!If (iPrint6>1) Then
    Write(18,'(//'' Lifetime statistics for REF structures identified by ICP. TimeStep,fs = '',f12.4)')tstep
    Write(18,'( /''  REF    n    N  Nmin Nmax       N/n    LifeTime,fs     SD,fs   Sum(LifeTime*n),fs'')')
!Endif
irec14=0
Do i=1,nref
    ltimin=999999
    ltimax=0
    Do j=istrbeg,nstr
        it=iType(6,j)
        If (j==istrbeg) Then
            jbeg=j
            If (it==i) Then
                lti=1
            Else
                lti=0
            Endif
            Cycle
        Endif
        itold=iType(6,j-1)
        
        ! counting of jumps between different REFs only        
        If (it==i.and.itold/=i) Then
            lti=1
            jbeg=j
        ElseIf (it==i.and.itold==i) Then
            lti=lti+1
        ElseIf (it/=i.and.itold==i) Then
            iTL(1,i)=iTL(1,i)+1
            iTL(2,i)=iTL(2,i)+lti
            iTL(5,i)=iTL(5,i)+lti*lti   ! for dispersion
            If (lti<ltimin) ltimin=lti
            If (lti>ltimax) ltimax=lti
            Write(14,'(2i5,i8,f10.2,2i8)')i,iTL(1,i),lti,Dble(lti)*tstep,jbeg,j
        Endif
        If (j==nstr.and.it==i) Then
            iTL(1,i)=iTL(1,i)+1
            iTL(2,i)=iTL(2,i)+lti
            iTL(5,i)=iTL(5,i)+lti*lti   ! for dispersion
            If (lti<ltimin) ltimin=lti
            If (lti>ltimax) ltimax=lti
            Write(14,'(2i5,i8,f10.2,2i8)')i,iTL(1,i),lti,Dble(lti)*tstep,jbeg,j
            irec14=irec14+1
        Endif            
    Enddo
    If (iTL(1,i)==0) Then
        Write(18,'(5i5,f12.4,4f12.2)')i,0,0,0,0,0.d0,0.d0,0.d0,tti
        Cycle
    Endif
    iTL(3,i)=ltimin
    iTL(4,i)=ltimax
    xtl1=Dble(iTL(1,i))
    xtl2=Dble(iTL(2,i))
    xtl5=Dble(iTL(5,i))
    tli=xtl2/xtl1
    TL(1,i)=tli*tstep
    tti=tti+TL(1,i)*xtl1
    tdisp=xtl5/xtl1-tli**2
    If (iTL(1,i)>1) Then
        tdisp=tdisp*xtl1/(xtl1-1.d0)
    Else
        tdisp=0.d0
    Endif
    TL(2,i)=DSQRT(tdisp)*tstep
    !If (iPrint6>1) 
    Write(18,'(5i5,f12.4,4f12.2)')i,iTL(1:4,i),tli,TL(1,i),TL(2,i),tti
Enddo
Close(18)

! iType(6,:) - REF str identified by ICP
! iType(7,:) - iterations in ICP 
! iType(8,:) - eps of ICP (IDNINT(eps*1.d6)

! iTL(1,i) - number of traj intervals for REF i (n)
! iTL(2,i) - sum of interval durations in steps (N)
! iTL(3:4,i) - Nmin and Nmax
! iTL(5,i) - lifetime**2, then LT disperison
! iTL(6:10,i) - copy of old data iTL(1:5,i) after exclusion of short intervals from iTL(1:5,i) 

! TL(1,i) - average lifetime of trj intervals for the given REF i (N/n*tstep) 
! TL(2,i) - SD of lifetimes for traj intervals for the given REF i


! Exclude short intervals
iTL(6:10,1:nref)=iTL(1:5,1:nref)
TL(3:4,1:nref)=TL(1:2,1:nref)
If (irun>0) Then
    iTL(1:5,1:nref)=0
    TL(1:2,1:nref)=0.d0
    iTL(6:10,1:nref)=iTL(1:5,1:nref)
    TL(3:4,1:nref)=TL(1:2,1:nref)
    Rewind(14)
    Read(14,*)
    Open(15,File=File15)
    Write(15,'(''  Ref    i  Steps Lifetime,fs'')')
    irefold=0
    ij=0
    nr=0
    Do While(.not.EOF(14))
        Read(14,'(2i5,i8,f10.2,2i8)')iref,itmp,ilt,xlt
        If (iref==0) Cycle
        If (xlt<xltminTresh) Cycle
        If (iref/=irefold) Then
            nr=0
            iltmin=999999
            iltmax=0
            ilttot=0
        Endif
        nr=nr+1
        ij=ij+1
        Write(15,'(2i5,i8,f10.2,2i8)')iref,nr,ilt,xlt
        If (ilt<iltmin) iltmin=ilt
        If (ilt>iltmax) iltmax=ilt
        ilttot=ilttot+ilt
        iTL(1,iref)=nr
        iTL(2,iref)=ilttot
        iTL(3,iref)=iltmin
        iTL(4,iref)=iltmax
        xtl=Dble(ilt)*tstep
        TL(1,iref)=TL(1,iref)+xtl
        TL(2,iref)=TL(2,iref)+xtl*xtl
        irefold=iref
    Enddo
    TL(1,1:nref)=TL(1,1:nref)/Dble(iTL(1,1:nref))
    TL(2,1:nref)=DSQRT(TL(2,1:nref)/Dble(iTL(1,1:nref))-TL(1,1:nref)**2)   ! dispersion (RMSD), not SD 
    Close(14)
    Close(15)
    Open(14,File=File15)
Endif       
    

! Sort and condence data from file14 to file16
Rewind(14)
Open(16,File=File16)
Write(16,'(''  Ref    i     LTS  Nint    LT,fs     Histo  <--Condensed and sorted data of lifetime intervals: Ref - REFstructure ID; i - ord.number; LTS - life time in steps; Nint - number of trajectory intervals with the given LTS; LT=LTS*TimeStep'')')
Read(14,*)
Do i=1,nref
    nr=iTL(1,i)
    If (nr==0) Cycle
    If (Allocated(tuni)) Deallocate(tuni,iuni)
    Allocate(tuni(nr),iuni(nr))
    Do j=1,nr
        Read(14,'(18x,f10.2)')tuni(j)
    Enddo
    Call HeapSort1(nr,tuni,iuni)
    il=0
    ts=0.d0
    Do lts=1,iTL(4,i)
        nint=0
        Do j=1,nr
            lt=IDNINT(tuni(j)/tstep)
            If (lt==lts) nint=nint+1
        Enddo
        If (nint==0) Cycle
        il=il+1
        symb=' '
        If (il==1) symb='*'
        tx=Dble(lts)*tstep
        tt=Dble(iTL(2,i))*tstep
        ts=ts+Dble(nint*lts)*tx/tt
        Write(16,'(2i5,i8,i5,2f10.2,1x,a1)')i,il,lts,nint,tx,ts,symb
    Enddo
Enddo
Close(16)
Close(14)
        
! Final statistics on lifetimes
If (Allocated(tuni)) Deallocate(tuni,iuni)
Allocate(tuni(nref),iuni(nref))
tuni(1:nref)=TL(1,1:nref)
Call HeapSort1(nref,tuni,iuni)
Write(6,'(//'' **** Lifetimes statistics for the jumps between the different REFs only ***'')')
Write(6,'(  '' Lifetime statistics for REF structures identified by ICP. SORTED by NON-ZERO lifetimes. TimeStep,fs = '',f8.4)')tstep
Write(6,'(  '' (n - number of trj intervals assigned to REF; Nsum - their total duration in steps; Naver=Nsum/n; Lifetime=Naver*TimeStep; SD - standard deviation of LifeTime)'')')
Write(6,'(  '' Detailed statistics see in files str-ref-lifetimes.dat and str-ref-lifetimes2.dat'')')
Write(6,'( /''    i  REF    n  Nsum  Nmin  Nmax      Naver    LifeTime,fs +-  SD,fs    LTmin,fs    LTmax,fs   Sum(LifeTime*n),fs      Nvib1     Nvib2      Freq1,Freq2,cm-1      Tvib1,Tvib2,fs'')')
ii=0
tti=0.d0
ltmin=IDNINT(xLTminTresh/tstep)
Do j=nref,1,-1
    i=iuni(j)
    If (iTL(1,i)==0) Cycle
    ii=ii+1
    tli=TL(1,i)/tstep
    tti=tti+TL(1,i)*Dble(iTL(1,i))
    ! Compare LT with vibration periods
    Call SpectralConversion(i,f1,f2,t1,t2)
    xnv1=TL(1,i)/t1
    xnv2=TL(1,i)/t2
!    Write(6,'(3i5,3i6,f12.4,5f12.2)')ii,i,iTL(1:4,i),tli,TL(1,i),TL(2,i),Dble(iTL(3:4,i))*tstep,tti
    Write(6,'(3i5,3i6,f12.4,5f12.2,10x,3(2f10.2,2x))')ii,i,iTL(1:4,i),tli,TL(1,i),TL(2,i),Dble(iTL(3:4,i))*tstep,tti,xnv1,xnv2,f1,f2,t1,t2
Enddo

End
!****************************************************************
Subroutine ICPassign(n,C,Amass,irefmin,epsmin)
Use Vars, Only: nref,Cref,iref,istr,Cmin0,Cref0

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n),C0(3,n),Amass(n),CM(3),PMI(3),PAxes(3,3),Crot(3,n)

! Try all the REF structures and determine the best RMS deviation between ADMP and REF 

Call InertiaNew(2,1,1,n,C,Amass,CM,PMI,PAxes,kRotTyp)
epsmin=999999.d99
Do iref=1,nref
    C0(1:3,1:n)=Cref(1:3,1:n,iref)
    Do k=0,3
        Crot(1:3,1:n)=C(1:3,1:n)
        If (k==1) Then
            Crot(2:3,1:n)=-Crot(2:3,1:n)
        ElseIf (k==2) Then
            Crot(1,1:n)=-Crot(1,1:n)
            Crot(3,1:n)=-Crot(3,1:n)
        ElseIf (k==3) Then
            Crot(1:2,1:n)=-Crot(1:2,1:n)
        Endif
        Call DoICP(n,Crot,C0,ires,iter,eps)    ! Find the best coincidence between current REF and the given ADMP
        If (eps<epsmin) Then
            irefmin=iref
            epsmin=eps
            Cmin0(1:3,1:n)=Crot(1:3,1:n)       ! Remember the best ADMP/REF pair
            Cref0(1:3,1:n)=C0(1:3,1:n)
        Endif
    Enddo
Enddo

!Write(6,'(//2i5,f12.4)')irefmin,itermin,epsmin

End
!****************************************************************
Subroutine ICPassignFull(mode,mref,n,C,Amass,irefmin,epsmin)
Use Vars, Only: nref,Cref,iref,istr,Cmin0,Cref0,nfact,IAbuf,EpsBuf,cscale,CoordScaleFactor,epst

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n),C0(3,n),Amass(n),CM(3),PMI(3),PAxes(3,3),Crot(3,n),C1(3,n),Cmin(3,n),Eps0(mref)
Integer(4) p(0:n),s,IA(n),IA0(0:n,mref),IX(mref)


! Try all the REF structures and determine the best RMS deviation between ADMP and REF 

Call InertiaNew(2,1,1,n,C,Amass,CM,PMI,PAxes,kRotTyp)
epsmin=999999.d99
epst(1:nref)=999.d0

Do ir=1,mref
    
    iref=ir
    If (mode>0) iref=IAbuf(0,ir)
    
    C0(1:3,1:n)=Cref(1:3,1:n,iref)
    
    ! Do all permutation of atom indices
    epsrefmin=999999.d0
    ipbest=0
    Call InitializePermutation(n,IA,nf,p,s)
    Do ip=1,nfact
        
        Do i=1,n
            C1(1:3,i)=C(1:3,IA(i))
        Enddo
        
        Call LSQR(n,C0,C1,Crot,eps)

        If (eps<epsrefmin) Then
            epsrefmin=eps
            Cmin(1:3,1:n)=Crot(1:3,1:n)       ! Remember the best permutation
            IA0(1:n,ir)=IA(1:n)
            IA0(0,ir)=iref
            ipbest=ip
            cbest=cscale
        Endif

        Call Permutation(n,IA,nf,p,s)
        
    Enddo
    
    If (epsrefmin<epsmin) Then
        irefmin=iref
        epsmin=epsrefmin
        Cmin0(1:3,1:n)=Cmin(1:3,1:n)
        Cref0(1:3,1:n)=C0(1:3,1:n)
        CoordScaleFactor=cbest
    Endif
     Write(24,'(i5,2i12,2x,<n>i3,''  epsmin='',f10.4)')iref,ipbest,nf,IA0(1:n,ir),epsrefmin
     epst(iref)=epsrefmin
    Eps0(ir)=epsrefmin
    
Enddo

Call HeapSort1(mref,Eps0,IX)
Do i=1,mref
    EpsBuf(i)=Eps0(IX(i))
    IAbuf(1:n,i)=IA0(1:n,IX(i))
    IAbuf(0,i)=IA0(0,IX(i))
Enddo

End
!****************************************************************
Subroutine ICPassignSingle(mref,n,C,Amass,irefmin,epsmin,IAtmp)
Use Vars, Only: nref,Cref,iref,istr,Cmin0,Cref0,nfact,IAbuf,EpsBuf,cscale,CoordScaleFactor

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n),C0(3,n),Amass(n),CM(3),PMI(3),PAxes(3,3),Crot(3,n),C1(3,n),Cmin(3,n),Eps0(mref)
Integer(4) p(0:n),s,IA(n),IA0(0:n,mref),IX(mref),IAtmp(n)


! Try all the REF structures and determine the best RMS deviation between ADMP and REF 

Call InertiaNew(2,1,1,n,C,Amass,CM,PMI,PAxes,kRotTyp)
epsmin=999999.d99

!Do ir=1,mref
ir=mref
    
    iref=ir
    !If (mode>0) iref=IAbuf(0,ir)
    
    C0(1:3,1:n)=Cref(1:3,1:n,iref)
    
    ! Do all permutation of atom indices
    epsrefmin=999999.d0
    ipbest=0
    Call InitializePermutation(n,IA,nf,p,s)
    Do ip=1,nfact
        
        Do i=1,n
            C1(1:3,i)=C(1:3,IA(i))
        Enddo
        
        Call LSQR(n,C0,C1,Crot,eps)

        If (eps<epsrefmin) Then
            epsrefmin=eps
            Cmin(1:3,1:n)=Crot(1:3,1:n)       ! Remember the best permutation
!            IA0(1:n,ir)=IA(1:n)
!            IA0(0,ir)=iref
            IAtmp(1:n)=IA(1:n)
            ipbest=ip
            cbest=cscale
        Endif

        Call Permutation(n,IA,nf,p,s)
        
    Enddo
    
    
    !If (epsrefmin<epsmin) Then
        irefmin=iref
        epsmin=epsrefmin
        !Cmin0(1:3,1:n)=Cmin(1:3,1:n)
        !Cref0(1:3,1:n)=C0(1:3,1:n)
        !CoordScaleFactor=cbest

        return
    
        
    !Endif
     !Write(24,'(i5,2i12,2x,<n>i3,''  epsmin='',f10.4)')iref,ipbest,nf,IA0(1:n,ir),epsrefmin
    Eps0(ir)=epsrefmin
    
!Enddo

!    
Call HeapSort1(mref,Eps0,IX)
Do i=1,mref
    EpsBuf(i)=Eps0(IX(i))
    IAbuf(1:n,i)=IA0(1:n,IX(i))
    IAbuf(0,i)=IA0(0,IX(i))
Enddo

End
!****************************************************************
Subroutine CalcSSE(n,C1,C2,sse)
Implicit Real(8) (A-H,O-Z)

Real(8) C1(3,n),C2(3,n)

sse=0.d0
Do i=1,n
    sse=sse+Sum((C1(1:3,i)-C2(1:3,i))**2)        ! Direct RMSD calculation between two structures
Enddo
sse=DSQRT(sse/Dble(n))

End
!****************************************************************
Subroutine ICPassignShort(mref,n,C,Amass,irefmin,epsmin)
Use Vars, Only: nref,Cref,iref,istr,Cmin0,Cref0,nfact,IAbuf,EpsBuf,cscale,CoordScaleFactor,epst

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n),C0(3,n),Amass(n),CM(3),PMI(3),PAxes(3,3),Crot(3,n),C1(3,n),Cmin(3,n),Eps0(mref)
Integer(4) p(0:n),s,IA(n),IA0(0:n,mref),IX(mref)


! Try all the REF structures and determine the best RMS deviation between ADMP and REF 

Call InertiaNew(2,1,1,n,C,Amass,CM,PMI,PAxes,kRotTyp)
epsmin=999999.d99
epst(1:nref)=999.d0
ipbest=0


Do ir=1,mref

    iref=IAbuf(0,ir)
    IA(1:n)=IAbuf(1:n,ir)
    
    C0(1:3,1:n)=Cref(1:3,1:n,iref)
    Do i=1,n
        C1(1:3,i)=C(1:3,IA(i))
    Enddo
        
    Call LSQR(n,C0,C1,Crot,eps)

    Write(24,'(i5,2i12,2x,<n>i3,''  epsmin='',f10.4)')iref,ipbest,nfact,IA(1:n),eps
    Eps0(ir)=eps
    epst(iref)=eps
    
    If (eps<epsmin) Then
        epsmin=eps
        irefmin=iref
        Cmin0(1:3,1:n)=Crot(1:3,1:n)
        Cref0(1:3,1:n)=C0(1:3,1:n)
        CoordScaleFactor=cscale
    Endif

Enddo

Call HeapSort1(mref,Eps0,IX)
Do i=1,mref
    EpsBuf(i)=Eps0(IX(i))
    IA0(0:n,i)=IAbuf(0:n,IX(i))
Enddo
IAbuf(0:n,1:mref)=IA0(0:n,1:mref)


End
!**************************************************************************************        
Subroutine InitializePermutation(n,IA,nf,p,s)
Implicit Real(8) (A-H,O-Z)

Integer(4) IA(n),p(0:n),s

    Do i=1,n
        IA(i)=i
        p(i)=-i
    Enddo
    
    p(0)=0
    s = 1
    nf=0

End
!**************************************************************************************        
Subroutine Permutation(n,IA,nf,p,s)
Implicit Real(8) (A-H,O-Z)

Integer(4) IA(n),p(0:n),s

        
    nf=nf+1
    !    Write(7,'(i8,<n>i2)')nf,IABS(p(1:n))
    k = 0
    Do i=2,n
        If (p(i)<0 .and. (iabs(p(i)) > iabs(p(i-1))) .and. (iabs(p(i)) > iabs(p(k)))) k = i
    Enddo
    Do i=1,n-1
        If (p(i)>0 .and. (iabs(p(i)) > iabs(p(i+1))) .and. (iabs(p(i)) > iabs(p(k)))) k = i
    Enddo
    If (k/=0) then
        Do i=1,n    !reverse elements > k
            If (iabs(p(i)) > iabs(p(k))) p(i)=-p(i)
        Enddo
        If (p(k)<0) then 
            i = k-1 
        else 
            i = k+1
        Endif
        itmp = p(k)
        p(k) = p(i)
        p(i) = itmp
        s = -s
            
        ! Swap atom numbers
        itmp=iA(i)
        iA(i)=iA(k)
        iA(k)=itmp

    Endif
!    If (k == 0) Exit

End   
!****************************************************************
Subroutine PrintStr(iu,n,C,Str)
Use Vars, Only: iref,istr,NA

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n)
Character(*) Str
Character(6) aname

Write(iu,'(a<Len_Trim(Str)>)')Trim(Str)
Do i=1,n
    Call SetAName(NA(i),i,aname)
    Write(iu,'(a5,3f12.6)')aname,C(1:3,i)
Enddo
Write(iu,*)
    
End    
!****************************************************************
Subroutine DoICP(n,C,C0,ires,iter,eps)
Use Vars, Only: iref,istr

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n),C0(3,n),C1(3,n),D(n,n),Cmin(3,n)
Integer(4) IA(n)

! First step of ICP algorithm - find iteratively the best atom matching between ADMP and REF structures

MaxIter=100
eps=999999.d0
epsmin=999999.d0
method=2    ! matching method

Do iter=1,MaxIter
    
    ! Matching atoms
    If (method==1) Then
        
        IA=0
        Do i=1,n
            rmin=999999.d0
            Do j=1,n
                If (IA(j)>0) Cycle
                rij=DSQRT(Sum((C(1:3,j)-C0(1:3,i))**2))
                If (rij<rmin) Then
                    rmin=rij
                    jmin=j
                Endif
            Enddo
    !        Write(6,'(2i5,g20.4,3x,3f10.4,3x,3f10.4,3x,3f10.4,f15.4)')i,jmin,rmin,C0(1:3,i),C(1:3,i),C(1:3,jmin),rmin
            C1(1:3,i)=C(1:3,jmin)
            IA(jmin)=i
        Enddo
    
    ElseIf (method==2) Then
    
        D=0.d0
        Do i=1,n
            Do j=1,n
                rij=DSQRT(Sum((C(1:3,j)-C0(1:3,i))**2))
                D(i,j)=rij
            Enddo
        Enddo
        Do i1=1,n
            rmin=99999.d0
            Do i=1,n
                Do j=1,n
                    If (D(i,j)<rmin) then
                        rmin=D(i,j)
                        imin=i
                        jmin=j
                    Endif
                Enddo
            Enddo
            D(imin,1:n)=999999.d0
            D(1:n,jmin)=999999.d0
            C1(1:3,imin)=C(1:3,jmin)
    !        Write(6,'(2i5,g20.4,3x,3f10.4,3x,3f10.4,3x,3f10.4,f15.4)')i1,jmin,rmin,C0(1:3,imin),C(1:3,i1),C(1:3,jmin),rmin
        Enddo
    Endif

    ! Make the best linear transformation (in a least-squares sense) between REF and ADMP (with atoms swapped for best coincidence)
    epsold=eps
    Call LSQR(n,C0,C1,C,eps)
    
!    Write(6,'(''istr:'',i5,''  iref:'',i5,''  iter:'',i5,2f15.6)')istr,iref,iter,eps,epsold
    !Write(6,'(/''@coord-C0'')')
    !Do i=1,n
    !    Write(6,'(i3,3f10.4)')12,C0(1:3,i)
    !Enddo
    !Write(6,'(/''@coord-Ctransformed'')')
    !Do i=1,n
    !    Write(6,'(i3,3f10.4)')12,C(1:3,i)
    !Enddo
    !Write(6,*)
    
    ! Remember the best transformation among all iterations
    If (eps<epsmin) Then
        epsmin=eps
        Cmin=C
    Endif
    
    If (DABS((eps-epsold)/epsold)<1.d-3) Exit
    
Enddo

C=Cmin
eps=epsmin
sse=0.d0
Do i=1,n
    sse=sse+Sum((C(1:3,i)-C0(1:3,i))**2)        ! Additional check for RMSD (direct RMSD calculation)
Enddo
sse=DSQRT(sse/Dble(n))

!Write(6,'(/''*** istr:'',i5,''  iref:'',i5,''  iter:'',i5,''  eps:'',f15.6,''  sse:''f15.6)')istr,iref,iter,eps,sse
!Write(6,'(''@coord-C0'')')
!Do i=1,n
!    Write(6,'(i3,3f10.4)')12,C0(1:3,i)
!Enddo
!Write(6,'(/''@coord-Ctransformed'')')
!Do i=1,n
!    Write(6,'(i3,3f10.4)')12,C(1:3,i)
!Enddo
!Write(6,*)

End
!***************************************************************************
Subroutine LSQR(n,Y,X,X1,eps)    !,rmsd,c,T,R)
Use Vars, Only: iScaling,c=>cscale

Implicit Real(8) (A-H,O-Z)

Integer(4),parameter::m=3   ! dimension of vectors
Real(8) X(m,n),Y(m,n),X1(m,n),R(m,m),D(m),U(m,m),V(m,m),DD(m,m),XK(n,n),Sxy(m,m),S(m,m),Rx(m),Ry(m),T(m),X1min(m,n),Rmin(m,n),dtmp(m),utmp(3,3),vtmp(3,3)
Integer(4) kkk(1)

! Optimal rotation (in a least-squares sence) of point set B to achive best coincidense with poin set A (algorithm by Umeyama1991)
! This realizes the exact formulas of Umeyama exept setting c=1 (no scaling)
! Y - reference vectors (not moved) , X - moved vectors; X1 - transformed X
! R - rotation matrix

!iScaling=1

!rmsdmin=999999.d0
!Do k=1,3

!Y=0.d0
!Y(1,2)=-1.d0
!Y(2,3)=2.d0
!X=0.d0
!X(1,2)=1.d0
!X(2,3)=2.d0

m1=m-1
xn=1.d0/Dble(n)

XK=0.d0
Do i=1,n
    XK(i,i)=1.d0
Enddo
XK=XK-xn


! Mean values of X, Y
Rx=0.d0
Ry=0.d0
Do ia=1,n
    Rx(1:m)=Rx(1:m)+X(1:m,ia)
    Ry(1:m)=Ry(1:m)+Y(1:m,ia)
Enddo

!! Variances
!sx2=0.d0
!sy2=0.d0
!Do i=1,n
!    sx2=sx2+Sum(X(1:3,i)-Rx(1:3))**2
!    sy2=sy2+Sum(Y(1:3,i)-Ry(1:3))**2
!Enddo
!sx2=sx2*xn
!sy2=sy2*xn

X1=matmul(X,XK)
sx2=Sum(X1**2)*xn   ! Frobenius norm (squared) of x.XK     - this is a varianse of X vectors
X1=matmul(Y,XK)
sy2=Sum(X1**2)*xn   ! Frobenius norm (squared) of y.XK     - this is a variance of Y vectors


! Covariance, its SVD decomposition, rank and determinants
Sxy=matmul(Y,matmul(XK,Transpose(X)))*xn   
U=Sxy
! SVD decomposition Sxy=U.d.V^T w/o sorting 
Call svdcmp(U,m,m,m,m,d,V)  
kkk=MinLoc(d)
k3=kkk(1)

!Sorting of d,U,V in descending order of d. (Strange, but it seems that it works without it...)
kkk=MaxLoc(d)
k1=kkk(1)
Do k=1,3
    If (k==k1.or.k==k3) Cycle
    k2=k
Enddo
dtmp(1)=d(k1)
dtmp(2)=d(k2)
dtmp(3)=d(k3)
d=dtmp
utmp(1:3,1)=u(1:3,k1)
utmp(1:3,2)=u(1:3,k2)
utmp(1:3,3)=u(1:3,k3)
u=utmp
utmp(1:3,1)=v(1:3,k1)
utmp(1:3,2)=v(1:3,k2)
utmp(1:3,3)=v(1:3,k3)
v=utmp

! check for sorted decomposition
!utmp=0.d0
!utmp(1,1)=d(1)
!utmp(2,2)=d(2)
!utmp(3,3)=d(3)
!vtmp=matmul(u,matmul(utmp,Transpose(v)))
!tmp=Sum((Sxy-vtmp)**2)
!write(35,*)tmp


irankS=0
Do i=1,m
    If(DABS(D(i))>1.d-14) irankS=irankS+1
Enddo
ds=Det23(m,Sxy)


! Matrix S
S=0.d0
Do i=1,m
    S(i,i)=1.d0
Enddo
If (irankS==m) Then
    If (ds<0.d0) S(3,3)=-1.d0
ElseIf (irankS==m1) Then
    du=Det23(m,U)
    dv=Det23(m,V)
    If (du*dv==-1) S(3,3)=-1.d0
Else
    Write(6,'('' WARNING! Rank of Sxy matrix is less than m-1:'',2i5,'' Transformation is probably erroneous'')')irankS,m-1
Endif
    
! Transformation parameters R, t, c

! Optimal rotation matrix
R=matmul(U,matmul(S,Transpose(V)))


! Optimal scaling factor
trDS=0.d0
Do i=1,m
    trDS=trDS+D(i)*S(i,i)
Enddo

c=trDS/sx2      ! c is optimal scaling factor for structure X
If (iScaling==0.or.iScaling==2) c=1.d0     ! No scaling

! Best SSE
eps2=sy2+c*c*sx2-2*c*trDS
If (eps2<0.d0.and.eps2<1.d-8) eps2=DABS(eps2)
eps=DSQRT(eps2)

ctmp=c
c=1.d0  ! now, after best scaling and eps are estimated to make the assignment, set c=1 to obtain unscaled structure (to avoid ADMP structure distortion)

! Optimal translation vector
T=(Ry-c*matmul(R,Rx))*xn

! Transformed vector X
Do i=1,n
    X1(1:m,i)=c*matmul(R,X(1:m,i))+T
Enddo

! Check for assignment validity 
!Write(3,'(//''@XY1'')')
!Do i=1,n
!    Write(3,'(i3,3f10.6,5x,3f10.6,5x,f10.4)')4,Y(1:3,i)
!Enddo
!Write(3,*)
!Write(3,'(//''@XY2'')')
!Do i=1,n
!    Write(3,'(i3,3f10.6,5x,3f10.6,5x,f10.4)')4,X1(1:3,i)    ,X1(1:3,i),Sum((Y(1:3,i)-X1(1:3,i))**2)
!Enddo
!Write(3,*)
!iu=3
!Call CalcRMSD(iu,0,n,Y,X1,rmsd)
!Write(3,*)
!
!R(1,1:3)=(/1.d0,        0.d0 ,         0.d0/)
!R(2,1:3)=(/0.d0,DCOSD(180.d0),-DSIND(180.d0)/)
!R(3,1:3)=(/0.d0,DSIND(180.d0), DCOSD(180.d0)/)
!Do i=1,n
!    X(1:m,i)=matmul(R,X1(1:m,i))
!Enddo
!X1=X
!Write(3,'(//''@XY3'')')
!Do i=1,n
!    Write(3,'(i3,3f10.6,5x,3f10.6,5x,f10.4)')4,X1(1:3,i)    ,X1(1:3,i),Sum((Y(1:3,i)-X1(1:3,i))**2)
!Enddo
!Write(3,*)


c=ctmp  ! restore proper c to print it out

!Call CalcRMSD(0,n,Y,X1,rmsd)
!If (rmsd<rmsdmin) Then
!    kmin=k
!    rmsdmin=rmsd
!    Rmin=R
!    X1min=X1
!Endif
!enddo
!X1=X1min
!R=Rmin

rmsd=eps
! Check whether rmsd=eps, just in case. Actually, RMSD should be equal to eps
!iu=3   !33
!Call CalcRMSD(iu,0,n,Y,X1,rmsd)
!Write(33,'(f12.6)')eps-rmsd
!Write(3,*)

End
!***************************************************************************
Subroutine LSQRold(n,Y,X,X1,eps)
Use Vars, Only: iScaling,c=>cscale

Implicit Real(8) (A-H,O-Z)

Integer(4),parameter::m=3   ! dimension of vectors
Real(8) X(m,n),Y(m,n),X1(m,n),R(m,m),D(m),U(m,m),V(m,m),DD(m,m),XK(n,n),Sxy(m,m),S(m,m),Rx(m),Ry(m),T(m)

! Optimal rotation (in a least-squares sence) of point set B to achive best coincidense with poin set A (algorithm by Umeyama1991)
! This realizes the exact formulas of Umeyama exept setting c=1 (no scaling)
! Y - reference vectors (not moved) , X - moved vectors; X1 - transformed X
! R - rotation matrix

!Y=0.d0
!Y(1,2)=-1.d0
!Y(2,3)=2.d0
!X=0.d0
!X(1,2)=1.d0
!X(2,3)=2.d0

m1=m-1
xn=1.d0/Dble(n)

XK=0.d0
Do i=1,n
    XK(i,i)=1.d0
Enddo
XK=XK-xn


! Mean values of X, Y
Rx=0.d0
Ry=0.d0
Do ia=1,n
    Rx(1:m)=Rx(1:m)+X(1:m,ia)
    Ry(1:m)=Ry(1:m)+Y(1:m,ia)
Enddo

!! Variances
!sx2=0.d0
!sy2=0.d0
!Do i=1,n
!    sx2=sx2+Sum(X(1:3,i)-Rx(1:3))**2
!    sy2=sy2+Sum(Y(1:3,i)-Ry(1:3))**2
!Enddo
!sx2=sx2*xn
!sy2=sy2*xn

X1=matmul(X,XK)
sx2=Sum(X1**2)*xn   ! Frobenius norm (squared) of x.XK     - this is a varianse of X vectors
X1=matmul(Y,XK)
sy2=Sum(X1**2)*xn   ! Frobenius norm (squared) of y.XK     - this is a variance of Y vectors


! Covariance, its SVD decomposition, rank and determinants
Sxy=matmul(Y,matmul(XK,Transpose(X)))*xn   
U=Sxy
Call svdcmp(U,m,m,m,m,d,V)

irankS=0
Do i=1,m
    If(DABS(D(i))>1.d-6) irankS=irankS+1
Enddo
ds=Det23(m,Sxy)

! Matrix S
S=0.d0
Do i=1,m
    S(i,i)=1.d0
Enddo
If (irankS==m) Then
    If (ds<0.d0) S(1,1)=-1.d0
ElseIf (irankS==m1) Then
    du=Det23(m,U)
    dv=Det23(m,V)
    If (du*dv==-1) S(1,1)=-1.d0
Else
    Write(6,'('' WARNING! Rank of Sxy matrix is less than m-1:'',2i5,'' Transformation is probably erroneous'')')irankS,m-1
Endif
    
! Transformation parameters R, t, c

! Optimal rotation matrix
R=matmul(U,matmul(S,Transpose(V)))


! Optimal scaling factor
trDS=0.d0
Do i=1,m
    trDS=trDS+D(i)*S(i,i)
Enddo
If (iScaling==0) Then
    c=1.d0          ! No scaling
Else
    c=trDS/sx2      ! c is optimal scaling factor for structure X
Endif

! Best SSE
eps2=sy2+c*c*sx2-2*c*trDS
If (eps2<0.d0.and.eps2<1.d-8) eps2=DABS(eps2)
eps=DSQRT(eps2)

ctmp=c
c=1.d0  ! now, after best scaling and eps are estimated to make the assignment, set c=1 to obtain unscaled structure (to avoid ADMP structure distortion)

! Optimal translation vector
T=(Ry-c*matmul(R,Rx))*xn

! Transformed vector X
Do i=1,n
    X1(1:m,i)=c*matmul(R,X(1:m,i))+T
Enddo

c=ctmp  ! restore proper c to print it out

End
!****************************************************************
Function Det23(n,A)
Implicit Real(8) (A-H,O-Z)

Real(8) A(n,n)
   
If (n==2) Then
    Det23=A(1,1)*A(2,2)-A(2,1)*A(1,2)
ElseIf (n==3) Then
    Det23=      a(1,1)*a(2,2)*a(3,3)  &
              + a(1,2)*a(2,3)*a(3,1)  &
              + a(1,3)*a(2,1)*a(3,2)  &
              - a(1,3)*a(2,2)*a(3,1)  &
              - a(1,1)*a(2,3)*a(3,2)  &
              - a(1,2)*a(2,1)*a(3,3) 
Endif

End
!****************************************************************
Subroutine ReadCommandFile
Use Vars
Use StringMod
Use GetOptionMod

Implicit Real(8) (A-H,O-Z)

Character(255) Str,Str1

iPrintUnit=6

Open(3,File=FilDat)

Do While (.not.EOF(3))
    Read(3,'(a255)')Str
 
    If (Len_Trim(Str)==0) Cycle
    Call RemoveStrComment(0,'!',Str)
    If (Len_Trim(Str)==0) Cycle
    Str1=ToUpperCase(AdjustL(Str))
!    If (INDEX(Str1,'!')==1) Cycle
    If (INDEX(Str1,'*STRUCT')==1) Goto 10
    
    ii=iGetOptionNoUCase('FILREF',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
    If (ii>0) Then
        Do i=1,nRes
            FilRef(i)=sRes(i)
        Enddo
        nFilRef=nRes
    Endif

    ii=iGetOption('REFLABEL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)   ! Structure label in the input files
    If (ii>0) Then
        RefLabel=sRes(1)
        If (nRes>1) iRefLabelSkip=iRes(2)
        If (nRes>2.and.INDEX(sRes(3),'NXYZ')>0) iRefLabelFormat=1
        If (nRes>2.and.INDEX(sRes(3),'NDXYZ')>0) iRefLabelFormat=2
        If (nRes>2.and.INDEX(sRes(3),'DNXYZ')>0) iRefLabelFormat=3
    Endif

    ii=iGetOption('PRINTLEVEL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Print level for additional printing
    If (ii>0) iPrint6=iRes(1)

    ii=iGetOption('RECALC',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Recalc=(lim1,lim2,n) - steps to be assigned by full permutation search
    If (ii>0) Then                                                        !  lim1 - make full search each lim1 steps for all REF structures
        Do i=1,nres                                                       !  lim2 - make search for n structures only each lim2 steps
            iRecalc(i)=iRes(i)
        Enddo
    Endif

    ii=iGetOptionNoUcase('READASSIGNMENT',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Read assignments made in a previous calculation
    If (ii>0) Then                                                          ! To recalc the structure reading the old assignment you need:
        iReadAssign=nRes                                                    ! (1) old file admpstat.out - rename it to admpstat.out_source 
        Do i=1,nRes                                                         ! (2) Indicate it in the command ReadAssignment=admpstat.out_source
            FilAssign(i)=sRes(i)                                            ! (3) old file str-admp2ref4.dat - rename it to str-admp2ref4.dat_source
        Enddo                                                               ! (4) old *.g6 files (Windows)
    Endif

    ii=iGetOption('ADDL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Additional LifeTimes calculation
    If (ii>0) iAddLT=1

    ii=iGetOption('SCALING',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Additional LifeTimes calculation
    If (ii>0) iScaling=1
    
    ii=iGetOption('SKIPADMP',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Skip some starting ADMP structure (equilibration steps)
    If (ii>0) iSkipADMP=iRes(1)
    
    ii=iGetOption('EPSMAX',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Accuracy threshold for REF analysis (eps after ICP)
    If (ii>0) epsTresh=xRes(1)

    ii=iGetOption('LTMIN',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Accuracy threshold for REF analysis (eps after ICP)
    If (ii>0) xLTminTresh=xRes(1)
    
    ii=iGetOption('RBONDMIN',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)      ! Maximum bond length between two atoms to be considered as bonded atoms
    If (ii>0) RbondMin=xRes(1)
    ii=iGetOption('RBONDMAX',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)      ! Maximum bond length between two atoms to be considered as bonded atoms
    If (ii>0) RbondMax=xRes(1)

    !ii=iGetOption('EPSE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)       ! Accuracy threshold for energy analysis
    !If (ii>0) epsE=xRes(1)
    !ii=iGetOption('EPSFP',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)      ! Accuracy threshold for FP analysis
    !If (ii>0) epsFP=xRes(1)
    !ii=iGetOption('EPSSTR',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)     ! Accuracy threshold for the analysis on the basis of sorted interatomic distances 
    !If (ii>0) epsStr=xRes(1)
    !ii=iGetOption('EPSPMOI',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)    ! Accuracy threshold for PMOI analysis
    !If (ii>0) epsPMOI=xRes(1)
    !ii=iGetOptionNoUcase('ETOTLABEL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)   ! Total Energy label in the input files of QC program. Default: 'SCF Done:'
    !If (ii>0) EtotLabel0=sRes(1)

    ii=iGetOption('RBONDSCALE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)   ! Structure label in the input files
    If (ii>0) rBondScale=xRes(1)

    ii=iGetOption('STRLABEL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)   ! Structure label in the input files
    If (ii>0) StrLabel0=sRes(1)

    ii=iGetOption('EWORD',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)      ! Number of word/value on the input line where energy is given
    If (ii>0) iEword0=iRes(1)

    ii=iGetOption('E1VALUE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)    ! Energy of a single atom
    If (ii>0) Then
        E1value0=xRes(1)
        iE1value=1
    Endif

    ii=iGetOption('METHOD',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)    ! Energy of a single atom
    If (ii>0) Then
!        If (INDEX(sRes(1),'RMSD')==1) Then  ! Default assignment method
!            MethodAssign=1
!            iFPsubmethod=0
!        Endif
        If (INDEX(sRes(1),'FP')==1) Then
            MethodAssign=2
            iFPsubmethod=1
        Endif
        If (INDEX(sRes(1),'DD')==1) Then
            MethodAssign=2
            iFPsubmethod=2
        Endif
        If (INDEX(sRes(1),'SD')==1) Then
            MethodAssign=2
            iFPsubmethod=3
        Endif
    Endif

    ii=iGetOption('EPSRECALC',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)    ! Energy of a single atom
    If (ii>0) Then
        epsrecalc=xRes(1)
    Endif

    ii=iGetOption('REFRMSD',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)    ! Calculate RMSD of REF structures relative to other REFs
    If (ii>0) Then
        iRefRMSDcalc=1
    Endif

    ii=iGetOption('EBASE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)      ! Base energy to calculate Erel (if not provided, the first Etot value will be used)
    If (ii>0) Ebase=xRes(1)

    
    !ii=iGetOption('EINPUT',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)     ! Which kind of energies are supplied? (1 - Etot, 2 - Eb, 3 - Eb per atom)
    !If (ii>0) iEinput=iRes(1)
    !ii=iGetOption('EUNITS=AU',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)  ! Units of energy in the input file(s)
    !If (ii>0) iEunits=1
    !ii=iGetOption('EUNITS=EV',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
    !If (ii>0) iEunits=2
    !ii=iGetOption('EUNITS=KCAL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
    !If (ii>0) iEunits=3
    !ii=iGetOption('EUNITS=KJ',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
    !If (ii>0) iEunits=4
    !
    !ii=iGetOption('SHOWALL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)    ! Draw all clusters on a single page: ShowAll(nRow,iUp2Down,HorizontalShift,VerticalShift,iShowPerPage)
    !If (ii>0) Then
    !    iShowAll=1
    !!    If (nRes>=1.and.iRes(1)>0) nRow=iRes(1)
    !!    If (nRes>=2.and.iRes(2)>0) iUp2Down=iRes(2)
    !!    If (nRes>=3.and.iRes(3)>0) HorizontalShift=xRes(3)
    !!    If (nRes>=4.and.iRes(4)>0) VerticalShift=xRes(4)
    !!    If (nRes>=5.and.iRes(5)>0) iShowPerPage=iRes(5)
    !
    !    If (nRes>=1) Then
    !        Do j=1,nRes
    !            Str1=sRes(j)
    !            ii=iGetOption('DOWN2UP',Str1,iPrintUnit,nRes1,lRes1,iRes1,xRes1,sRes1)  ! Down2Up - Draw structures upwards(from bottom to top, row by row)
    !            If (ii>0) iUp2Down=-1
    !            ii=iGetOption('COL',Str1,iPrintUnit,nRes1,lRes1,iRes1,xRes1,sRes1)      ! COL= - Number of columns on a page
    !            If (ii>0.and.iRes1(1)>0) nCol=iRes1(1)
    !            ii=iGetOption('NCOL',Str1,iPrintUnit,nRes1,lRes1,iRes1,xRes1,sRes1)     ! The same as COL/COLS
    !            If (ii>0.and.iRes1(1)>0) nCol=iRes1(1)
    !            ii=iGetOption('HSHIF',Str1,iPrintUnit,nRes1,lRes1,iRes1,xRes1,sRes1)    ! Horizontal shift, in Angstroems
    !            If (ii>0.and.iRes1(1)>0) HorizontalShift=xRes1(1)
    !            ii=iGetOption('VSHIF',Str1,iPrintUnit,nRes1,lRes1,iRes1,xRes1,sRes1)    ! Vertical shift, in Angstroems
    !            If (ii>0.and.iRes1(1)>0) VerticalShift=xRes1(1)
    !            ii=iGetOption('PAGE',Str1,iPrintUnit,nRes1,lRes1,iRes1,xRes1,sRes1)     ! Number of structures per page
    !            If (ii>0.and.iRes1(1)>0) iShowPerPage=iRes1(1)
    !        Enddo
    !    Endif
    !
    !Endif
    
Enddo

Write(6,'(//'' ERROR! *STRUCT keyword followed by the structures file list must be present in the input file!'')')
Stop ' Can not find *STRUCT keyword!'

10 Continue

!
! Read input file name(s)
!
nFilInp=0
Open(1,File='tmp.tmp')
Do While(.not.EOF(3))
    Read(3,'(a255)')Str
    If (Len_Trim(Str)==0) Exit
    Str=AdjustL(Str)
    If (Str(1:1)=='!') Cycle
    nFilInp=nFilInp+1
    Write(1,'(a<Len_Trim(Str)>)')Trim(Str)
Enddo
Close(3)

Allocate(FilInp(nFilInp))
!EtotLabel=repeat(' ',Len(EtotLabel))
!ll=Len_Trim(EtotLabel0)
!EtotLabel(1:ll)=EtotLabel0(1:ll)
!StrLabel=StrLabel0
!iEword=iEword0
!E1value=E1value0
iFilGroup=1
nGroups=1
Rewind(1)
Do i=1,nFilInp
    Read(1,'(a255)')Str
    Call SubString(Str,nSubStr,SubStr)
    FilInp(i)=Trim(SubStr(1))
!    Str=ToUpperCase(Str)
!    ii=iGetOption('STRLABEL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
!    If (ii>0) StrLabel(i)=sRes(1)
!    ii=iGetOption('EWORD',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
!    If (ii>0) iEword(i)=iRes(1)
!    ii=iGetOption('E1VALUE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
!    If (ii>0) Then
!        E1value(i)=xRes(1)
!        iE1value=1
!    Endif
!    ii=iGetOption('GROUPB',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)      ! These files will be considered as Group B during comparison between two groups of structures
!    If (ii>0) Then
!        iFilGroup(i)=-1
!        nGroups=2
!    Endif
Enddo
Close(1,Status='DELETE')
!If (nGroups==2) Allocate(iGroupStr(MaxStr))

Write(6,*)
If (rBondMin>0.d0) Then
    Write(6,'('' RbondMin    '',f10.5)')rBondMin
Else
    Write(6,'('' RbondMin not set -- will be used by default (0)'')')
Endif
If (rBondMax>0.d0) Then
    Write(6,'('' RbondMax    '',f10.5)')rBondMax
Else
    Write(6,'('' RbondMax not set -- will be used by default as rBondScale*(Ra+Rb)'')')
Endif

Write(6,'('' RbondScale = '',f10.4)')rBondScale


Write(6,'(//'' Structures will be read in from the files using the structure labels (case-insensitive):'',a<Len_Trim(StrLabel)>)')StrLabel
Do i=1,nFilInp
    Write(6,'(i4,2x,a<Len_Trim(FilInp(i))>)')i,Trim(FilInp(i))   !,Trim(StrLabel(i))
Enddo
!Write(6,'(/'' Structures will be selected from file using the label (case-insensitive): '',a<Len_Trim(StrLabel)>)')Trim(StrLabel)
Write(6,*)

End
!*********************************************************************
Subroutine ReadRefStr
Use Vars
Use Elements, Only: AMS

Implicit Real(8) (A-H,O-Z)

Character(63) g6string
Character(14) g6type
Integer(4),allocatable::iequiv(:)
Real(8),allocatable::C1(:,:),C2(:,:),C3(:,:),Amass(:)
Integer(4),allocatable::IAtmp(:)

Open(8,File='ref-str.g6')
Open(10,File='ref-str.dat')

Freq=0.d0

nRef=0
Do ifr=1,nFilRef
    Open(7,File=FilRef(ifr))
    nfstr=0
    Do While (.not.EOF(7))
    
        !Call ReadNXYZ(-7,1,iRefLabelSkip,0,iRefLabelFormat,MaxAt,NA,C,NumAt,RefLabel)
        Call ReadNXYZvib(7,1,iRefLabelSkip,0,iRefLabelFormat,MaxAt,NA,C,NumAt,nvibr,f1,f2,ierr,RefLabel)
        
        If (ierr==-99) Exit ! STOP command found in FilRef 
        
        If (Numat==0) Exit
        nfstr=nfstr+1
        nRef=nRef+1
        If (nRef>MaxRef) Then
            Write(6,'(//'' Too many reference structures! Increase MaxRef. Stop.'')')
            Stop 'Too many Ref structs!'
        Endif
    
        If (nRef==1) Then
            NumatRef=Numat
            NAref(1:Numat)=NA(1:Numat)
            Allocate(iAM(Numat,Numat),iAMx(Numat,Numat))
        Else
            If (Numat/=NumatRef) Then
                Write(6,'('' ERROR! Numat is wrong in the reference structure:'',i2)')nRef
                Stop
            Endif
        Endif
        Cref(1:3,1:Numat,nRef)=C(1:3,1:Numat)
        
        If (nvibr>0) Freq(1:2,nref)=(/f1,f2/)
    
        nv=Numat
        Call GraphConnectivity(iConn)
        ne=nbndx
        !Call AdjMatrix(ne)
        g6string=repeat(' ',Len(g6string))
        Call AM2g6(nv,iAMx,g6string)
        
        arbr=0.d0
        diam=0.d0
        Do i=2,Numat
            Do j=1,i-1
                rb=Distance(i,j,Numat,C)
                If (rb>diam) diam=rb
                If (iAMx(i,j)>0) arbr=arbr+rb
            Enddo
        Enddo
        arbr=arbr/Dble(ne)
    
        ns=Len_Trim(g6string)
    !    Write(6,'(''Ref:'',i4,2x,''nv:'',i3,''  ne:'',i5,''  iConn:'',i2,''  g6: '',a<ns>)')nRef,nv,ne,iConn,g6string(1:ns)
        ll=Numat*(Numat-1)/2
        Write(10,'(i3,2i4,i3,i5,i2,a<ns>,5x,2f12.4,2x,<ll>i1)')ifr,nfstr,nRef,nv,ne,iConn,g6string(1:ns),arbr,diam,(iAMx(i,1:i-1),i=2,Numat)
        Write(8,'(a<ns>)')g6string(1:ns)
    
    Enddo
    Close(7)
Enddo
Numat=NumatRef

Allocate(g6ref(nRef))
lg6=Len(g6ref(1))
Rewind(8)
Do i=1,nRef
    g6ref(i)=repeat(' ',lg6)
    Read(8,'(a<ns>)')g6ref(i)
Enddo
Close(8)    
    
    
! Try to invoke labelg program
If (.not.lWin) Then
    Write(6,'(/'' Trying to invoke labelg to canonicalize the g6 strings (search in $PWD and $PATH directories only)'')')
    Inquire(File='labelg',EXIST=lexist)
    If (lexist) Then
        Call System('./labelg ref-str.g6 ref-str-canonic.g6')
    Else
        Call System('labelg ref-str.g6 ref-str-canonic.g6')
    Endif
Endif

! Read canonic structures if the file exists
g6type='(non-canonic) '
FilCan='ref-str-canonic.g6'
Inquire(File=FilCan,EXIST=lexist)
If (lexist) Then
    Open(8,File=FilCan)
    Do i=1,nRef
        Read(8,'(a<ns>)')g6string
        g6ref(i)=repeat(' ',lg6)
        g6ref(i)=g6string(1:ns)
    Enddo
    Close(8)    
    Write(6,'(i10,'' g6 strings have been read in.'',a<Len_Trim(FilCan)>)')nRef,FilCan
    g6type='(canonic)     '
Endif

! Check for equivalent structure graphs
Allocate(iequiv(nref))
iequiv=0
Do i=2,nref
    Do j=1,i-1
        If (Trim(g6ref(i))==Trim(g6ref(j))) iequiv(i)=j
    Enddo
Enddo

Allocate(iAMref(Numat,Numat,nref),RbondRef(nref),DiamRef(nref))
iAMref=0

 
! Printing ref structures info
Rewind(10)
Write(6, '(/'' Reference structures:'')')
Write(6, '( '' File Fstr iStr   nv   ne conn    g6 '',a14,t58,''AverRbond    Diameter'')')g6type
Do i=1,nRef
    Read(10,'(i3,2i4,i3,i5,i2,<ns+5>x,2f12.4,2x,<Numat*(Numat-1)/2>i1)')ifr,nfstr,itmp,nv,ne,iConn,RbondRef(i),DiamRef(i),(iAMref(k,1:k-1,i),k=2,Numat)
!   Write(10,'(i3,i4,i3,i5,i2,a<ns>)')ifr,nRef,nv,ne,iConn,g6string(1:ns)
    If (iequiv(i)==0) Then
        Write(6,'(6i5,4x,a<ns>,t55,2f12.4,5x)')ifr,nfstr,i,nv,ne,iConn,g6ref(i),RbondRef(i),DiamRef(i)
    Else
        Write(6,'(6i5,4x,a<ns>,2x''--> Ref'',i3,''  (Graphs are equivalent but the geometry can be different. Check it manually)'')')ifr,nfstr,i,nv,ne,iConn,g6ref(i),iequiv(i)
    Endif
Enddo
Close(10)

! Adiitional informaton for REFS -- RMSD of the given REF relative to all other REFs
!Write(10, '(/'' File Fstr iStr   nv   ne conn    g6 '',a14,t58,''AverRbond    Diameter'')')g6type
!Write(10,'(/'' -------------'')')
If (iRefRMSDcalc>0) Then
    
    Open(10,File='ref-str-rmsd.dat')
    Write(10,'(''  REF   Ref RMSD relative to other REFs'')')
    Write(10,'(5x,<nref>i10)')(j,j=1,nref)
    Allocate(C1(3,Numat),C2(3,Numat),C3(3,Numat),Amass(Numat),IAtmp(Numat))

    nfact=1
    Do ia=1,Numat
        Amass(ia)=AMS(NAref(ia))
        nfact=nfact*ia
    Enddo

    Do i=1,nref
        C1(1:3,1:Numat)=Cref(1:3,1:Numat,i)
        Call MassCenter(Numat,C1)
        Cref(1:3,1:Numat,i)=C1(1:3,1:Numat)
    Enddo

    Do i=1,nRef
        C1(1:3,1:Numat)=Cref(1:3,1:Numat,i)
        Do j=1,nref
            Call ICPassignSingle(j,Numat,C1,Amass,iref,eps,IAtmp)
            epst(j)=eps
        Enddo
        Write(10,'(i5,<nref>f10.4)')i,epst(1:nref)
    Enddo

    Close(10)
    Write(6,'(//'' Calculation of RMSD between REF structures (RefRMSD keyword) has been completed. See results in file ref-str-rmsd.dat -- Calculation finished'')')
    
    Call Timer(1,2,6,TT)
    Call Timer(1,3,6,TT)
    Stop
    
Endif


! Calculate fingerprints for REF structs
If (MethodAsssign>1) Then
    MaxFP=Max(4*Numat,Numat*(Numat-1)/2)
    Allocate(FPref(MaxFP,nref))
    Call FPassign(0,nref,Numat,C,ires,epsmin)
Endif
    
End
!**************************************************************************
Subroutine PrintNXYZnew(iU,iForm,N,NA,C,nf,nd,iLabel,Str)
Implicit Real(8) (A-H,O-Z)

Real(8) C(3,N)
Integer(4) NA(N)
Character(10) AName
Character(*) Str


!
! Subroutine PrintNXYZ print out the coordinates C(3,N) from array C(3,MaxAt) to the file iU 
! by format iForm = 0(DNXYZ), 1(NXYZ), 2(NDXYZ)
! Title is Str if iLabel>0
!


If (iLabel/=0) ls=Len_Trim(Str)
If (iLabel>0) Then
	Write(iU,'(<ls>a1)')(Str(i:i),i=1,ls)
ElseIf (iLabel<0) Then
	Write(iU,'(<ls>a1,i4)')(Str(i:i),i=1,ls),-iLabel
Endif

If (nf<=0.or.nd<=0.or.nf-2<=nd.or.nf>30) Then
    nf=15
    nd=8
Endif

Do i=1,N
	nai=na(i)
	Call SetAName(NAi,i,AName)
	If		(iForm==1) Then					! NXYZ
		Write(iU,'(i4,2x,3f<nf>.<nd>)')NA(i),(C(k,i),k=1,3)
	ElseIf	(iForm==2) Then					! NDXYZ
		Write(iU,'(i4,2x,a5,2x,3f<nf>.<nd>)')NA(i),AName,(C(k,i),k=1,3)
	ElseIf	(iForm==3.or.iForm==0) Then		! DNXYZ (Default format)
		Write(iU,'(a5,i4,2x,3f<nf>.<nd>)')AName,NA(i),(C(k,i),k=1,3)
	ElseIf	(iForm==4) Then		! NXYZ with ANAME
		Write(iU,'(a5,2x,3f<nf>.<nd>)')AName,(C(k,i),k=1,3)
	Endif
Enddo
Write(iU,*)

    End
!***********************************************************************
Function Quantile(p,n,X,iSorted)
Implicit Real(8) (A-H,O-Z)

Real(8) X(n),X0(0:n-1)
Integer(4) IA(n)

! Function Quantile gives the interpolated value for quantile position p (0<=p<=1) of discrete sample X. E.g. 75-th percentile is p=0.75, median is p=0.5
! X(n) is the data. If it is sorted, set iSorted=0, otherwise iSorted=0
! For p=0, Quantile gives the Xmin, for p=1 it gives Xmax.
! This code is based on GSL library https://www.gnu.org/software/gsl/doc/html/statistics.html#median-and-percentiles

X0(0:n-1)=X(1:n)
If (iSorted==0) Call HeapSort1(n,X0,IA)

xn1=Dble(n-1)
i=INT(xn1*p)
d=xn1*p-Dble(i)

Quantile=(1.d0-d)*X0(i)+d*X0(i+1)

End
!**************************************************************************
Subroutine ReadAssignment
Use Vars, Only:iAssigned,nstr,nref,Numat
Use StringMod
Implicit Real(8) (A-H,O-Z)

Character(255) Str,Str1
Character(10) buf,buf1
Integer(4) IAtmp(Numat)
Logical lex

INQUIRE(File='str-admp2ref4.dat_source',EXIST=lex)
If (lex) Open(28,File='str-admp2ref4.dat_source') 

Do While (.not.EOF(25))
    Read(25,'(a255)')Str
    If (INDEX(Str,'ADMP structures to RefStr using the ICP')==0) Cycle
    Do While(.not.EOF(25))
        Read(25,'(a255)')Str
        If (INDEX(Str,'t,fs   REF   Iter')>0) Exit
    Enddo
    If (EOF(25)) Stop 'Can not find OLD ASSIGNMENTS in the file admpstat.out_source'
    Do While(.not.EOF(25))
        Read(25,'(a255)')Str
        If (Len_Trim(Str)>80) Then
!            Read(Str,'(i6,10x,i6,6x,f12.4,52x,<Numat>i3)')ii,iref,eps,IAtmp(1:Numat)
            Read(Str,'(i6,10x,i6,6x,f12.4,64x,<Numat>i3)')ii,iref,eps,IAtmp(1:Numat)
        Else
            Read(Str,'(i6,10x,i6,6x,f12.4)')ii,iref,eps
            Do While (.not.EOF(28))
                If (ii==1) Then
                    Read(28,*)
                    Read(28,*)
                    Call ReadIREF(iref,Numat,IAtmp)
                    Exit
                Endif    
                Read(28,'(a255)')Str1
                If (Len_Trim(Str1)>80) Then
                    Call SubString(Str1,nsubstr,SubStr)
                    buf1=AdjustL(SubStr(1))
                    buf=AdjustR(buf1)
                    Read(buf,'(i10)')ii1
                    If (ii1==ii-1) Then
                        Call ReadIREF(iref,Numat,IAtmp)
                        Exit
                    Endif
                Endif
            Enddo
        Endif
        If (ii==0.or.ii>nstr) Exit
        If (iref>nref) Cycle
        iAssigned(1,ii)=iref
        iAssigned(2,ii)=INT(eps*1.d4)
        Write(27,'(i6,i4,<Numat>i3)')ii,iref,IAtmp(1:Numat)
    Enddo
Enddo

If (lex) Close(28)

End
!**************************************************************
Subroutine ReadIREF(iref,Numat,IAtmp)
Implicit Real(8) (A-H,O-Z)

Integer(4) IAtmp(Numat)
Character(255) Str1

Do While (.not.EOF(28))
    Read(28,'(a255)')Str1
    Read(Str1,'(i5)')iref1
    If (iref1==0) Then
        BackSpace(28)
        Exit
    Endif
    If (iref1==iref) Then
        Read(Str1,'(31x,<Numat>i3)')IAtmp(1:Numat)
        Exit
    Endif
Enddo
    
End
!**************************************************************
Subroutine LifeTimes2(irun,xLTminTresh)
Use Vars, xLTMinTresh1=>xLTminTresh

Implicit Real(8) (A-H,O-Z)
   
Real(8),allocatable::TL(:,:),tuni(:)
Integer(4),allocatable::iuni(:)
Integer(8),allocatable::iTL(:,:)
Integer(8) lti

Character(255) File14,File16,File18,File15

If (irun==0) Then
    File14='str-lifetimes2a.dat'
    File15='str-lifetimes2d.dat'
    File16='str-lifetimes2b.dat'
    File18='str-lifetimes2c.dat'
ElseIf (irun>0) Then
    File14='str-lifetimes2a_x.dat'
    File15='str-lifetimes2d_x.dat'
    File16='str-lifetimes2b_x.dat'
    File18='str-lifetimes2c_x.dat'
    Write(6,'(//'' *** Repeat the lifetime statistics using the minimum LifeTime treshold'',f10.2,'' fs'')')xLTMinTresh
Endif


! Lifetime statistics
Allocate(TL(4,nref),iTL(10,nref))
!tstep=xStep(3,nstr-1)   ! typical time step
iTL=0
iTL(3,1:nref)=999999
TL=0.d0
tti=0.d0
Open(14,File=File14)
Write(6,'(//'' Lifetime distributions are dumped  to the file str-ref-lifetimes.dat.'')')
Write(6,'(  '' Lifetime statistics    are written to the file str-ref-lifetimes3.dat (for additional analysis on tstep influence).'')')
Write(14,'(''  Ref  Int   Steps Lifetime,fs  ibeg    iend'')')
Open(18,File=File18)
!If (iPrint6>1) Then
    Write(18,'(//'' Lifetime statistics for REF structures identified by ICP. TimeStep,fs = '',f12.4)')tstep
    Write(18,'( /''  REF    n    N  Nmin Nmax       N/n    LifeTime,fs     SD,fs   Sum(LifeTime*n),fs'')')
!Endif
irec14=0
    Do j=istrbeg,nstr
        If (j==istrbeg) Then
            lti=1
            jbeg=istrbeg
            Cycle
        Endif
        it=iType(6,j)
        itold=iType(6,j-1)

        idiff1=Sum(IABS(IAbest(1:Numat,j)-IAbest(1:Numat,j-1)))
        idiff2=IABS(it-itold)
        idiff=idiff1+idiff2

        If (idiff1>0.and.idiff2==0) Then    ! jumps within the same REF
            Call CheckRot(j,it,irot)
            If (irot>0) idiff=0             ! if irot>0, the difference in IAbest is due to rotation, not due to the swap of atoms
        Endif
        
        ! counting of jumps 
        If (idiff>0) Then
            iTL(1,itold)=iTL(1,itold)+1
            iTL(2,itold)=iTL(2,itold)+lti
            iTL(5,itold)=iTL(5,itold)+lti*lti   ! for dispersion
            If (lti<iTL(3,itold)) iTL(3,itold)=lti
            If (lti>iTL(4,itold)) iTL(4,itold)=lti
            Write(14,'(2i5,i8,f10.2,2i8)')itold,iTL(1,itold),lti,Dble(lti)*tstep,jbeg,j
            irec14=irec14+1
            lti=1
            jbeg=j
        ElseIf (idiff==0) Then
            lti=lti+1
        Endif
        If (j==nstr) Then
            iTL(1,it)=iTL(1,it)+1
            iTL(2,it)=iTL(2,it)+lti
            iTL(5,it)=iTL(5,it)+lti*lti   ! for dispersion
            If (lti<iTL(3,it)) iTL(3,it)=lti
            If (lti>iTL(4,it)) iTL(4,it)=lti
            Write(14,'(2i5,i8,f10.2,2i8)')it,iTL(1,it),lti,Dble(lti)*tstep,jbeg,j
            irec14=irec14+1
        Endif            

    Enddo
    
Do i=1,nref    
    If (iTL(1,i)==0) Then
        Write(18,'(5i5,f12.4,4f12.2)')i,0,0,0,0,0.d0,0.d0,0.d0,tti
        Cycle
    Endif
    xtl1=Dble(iTL(1,i))
    xtl2=Dble(iTL(2,i))
    xtl5=Dble(iTL(5,i))
    tli=xtl2/xtl1
    TL(1,i)=tli*tstep
    tti=tti+TL(1,i)*xtl1
    tdisp=xtl5/xtl1-tli**2
    If (iTL(1,i)>1) Then
        tdisp=tdisp*xtl1/(xtl1-1.d0)
    Else
        tdisp=0.d0
    Endif
    TL(2,i)=DSQRT(tdisp)*tstep
    Write(18,'(5i5,f12.4,4f12.2)')i,iTL(1:4,i),tli,TL(1,i),TL(2,i),tti
Enddo
Close(18)

! iType(6,:) - REF str identified by ICP
! iType(7,:) - iterations in ICP 
! iType(8,:) - eps of ICP (IDNINT(eps*1.d6)

! iTL(1,i) - number of traj intervals for REF i (n)
! iTL(2,i) - sum of interval durations in steps (N)
! iTL(3:4,i) - Nmin and Nmax
! iTL(5,i) - lifetime**2, then LT disperison
! iTL(6:10,i) - copy of old data iTL(1:5,i) after exclusion of short intervals from iTL(1:5,i) 

! TL(1,i) - average lifetime of trj intervals for the given REF i (N/n*tstep) 
! TL(2,i) - SD of lifetimes for traj intervals for the given REF i


! Exclude short intervals
iTL(6:10,1:nref)=iTL(1:5,1:nref)
TL(3:4,1:nref)=TL(1:2,1:nref)
If (irun>0) Then
    iTL(1:5,1:nref)=0
    TL(1:2,1:nref)=0.d0
    iTL(6:10,1:nref)=iTL(1:5,1:nref)
    TL(3:4,1:nref)=TL(1:2,1:nref)
    Rewind(14)
    Read(14,*)
    Open(15,File=File15)
    Write(15,'(''  Ref    i  Steps Lifetime,fs'')')
    irefold=0
    ij=0
    nr=0
    Do While(.not.EOF(14))
        Read(14,'(2i5,i8,f10.2,2i8)')iref,itmp,ilt,xlt
        If (iref==0) Cycle
        If (xlt<xltminTresh) Cycle
        If (iref/=irefold) Then
            nr=0
            iltmin=999999
            iltmax=0
            ilttot=0
        Endif
        nr=nr+1
        ij=ij+1
        Write(15,'(2i5,i8,f10.2,2i8)')iref,nr,ilt,xlt
        If (ilt<iltmin) iltmin=ilt
        If (ilt>iltmax) iltmax=ilt
        ilttot=ilttot+ilt
        iTL(1,iref)=nr
        iTL(2,iref)=ilttot
        iTL(3,iref)=iltmin
        iTL(4,iref)=iltmax
        xtl=Dble(ilt)*tstep
        TL(1,iref)=TL(1,iref)+xtl
        TL(2,iref)=TL(2,iref)+xtl*xtl
        irefold=iref
    Enddo
    TL(1,1:nref)=TL(1,1:nref)/Dble(iTL(1,1:nref))
    TL(2,1:nref)=DSQRT(TL(2,1:nref)/Dble(iTL(1,1:nref))-TL(1,1:nref)**2)   ! dispersion (RMSD), not SD 
    Close(14)
    Close(15)
    Open(14,File=File15)
Endif       
    

! Sort and condence data from file14 to file16
Rewind(14)
Open(16,File=File16)
Write(16,'(''  Ref    i     LTS  Nint    LT,fs     Histo  <--Condensed and sorted data of lifetime intervals: Ref - REFstructure ID; i - ord.number; LTS - life time in steps; Nint - number of trajectory intervals with the given LTS; LT=LTS*TimeStep'')')
Read(14,*)
Do i=1,nref
    nr=iTL(1,i)
    If (nr==0) Cycle
    If (Allocated(tuni)) Deallocate(tuni,iuni)
    Allocate(tuni(nr),iuni(nr))
    Do j=1,nr
        Read(14,'(18x,f10.2)')tuni(j)
    Enddo
    Call HeapSort1(nr,tuni,iuni)
    il=0
    ts=0.d0
    Do lts=1,iTL(4,i)
        nint=0
        Do j=1,nr
            lt=IDNINT(tuni(j)/tstep)
            If (lt==lts) nint=nint+1
        Enddo
        If (nint==0) Cycle
        il=il+1
        symb=' '
        If (il==1) symb='*'
        tx=Dble(lts)*tstep
        tt=Dble(iTL(2,i))*tstep
        ts=ts+Dble(nint*lts)*tx/tt
        Write(16,'(2i5,i8,i5,2f10.2,1x,a1)')i,il,lts,nint,tx,ts,symb
    Enddo
Enddo
Close(16)
Close(14)
        
! Final statistics on lifetimes
If (Allocated(tuni)) Deallocate(tuni,iuni)
Allocate(tuni(nref),iuni(nref))
tuni(1:nref)=TL(1,1:nref)
Call HeapSort1(nref,tuni,iuni)
Write(6,'(//'' **** Lifetimes statistics including the jumps between the same REFs ***'')')
Write(6,'(  '' Lifetime statistics for REF structures identified by ICP. SORTED by NON-ZERO lifetimes. TimeStep,fs = '',f8.4)')tstep
Write(6,'(  '' (n - number of trj intervals assigned to REF; Nsum - their total duration in steps; Naver=Nsum/n; Lifetime=Naver*TimeStep; SD - standard deviation of LifeTime)'')')
Write(6,'(  '' Detailed statistics see in files str-ref-lifetimes.dat and str-ref-lifetimes2.dat'')')
Write(6,'( /''    i  REF    n  Nsum  Nmin  Nmax      Naver    LifeTime,fs +-  SD,fs    LTmin,fs    LTmax,fs   Sum(LifeTime*n),fs      Nvib1     Nvib2      Freq1,Freq2,cm-1      Tvib1,Tvib2,fs'')')
ii=0
tti=0.d0
ltmin=IDNINT(xLTminTresh/tstep)
Do j=nref,1,-1
    i=iuni(j)
    If (iTL(1,i)==0) Cycle
    ii=ii+1
    tli=TL(1,i)/tstep
    tti=tti+TL(1,i)*Dble(iTL(1,i))
    ! Compare LT with vibration periods
    Call SpectralConversion(i,f1,f2,t1,t2)
    xnv1=TL(1,i)/t1
    xnv2=TL(1,i)/t2
    Write(6,'(3i5,3i6,f12.4,5f12.2,10x,3(2f10.2,2x))')ii,i,iTL(1:4,i),tli,TL(1,i),TL(2,i),Dble(iTL(3:4,i))*tstep,tti,xnv1,xnv2,f1,f2,t1,t2
Enddo


End
!*******************************************************************************************
Subroutine RefPerm(it)
Use Vars

Implicit Real(8) (A-H,O-Z)

Integer(4) IA(Numat),P(0:Numat),s
Real(8) C0(3,Numat),C1(3,Numat),Crot(3,Numat)

! Determination of the permutations in REF atom ordering which can be reduced to the simple rotation


iref=it
n=Numat
epsrefmin=0.1d0
irot=0
eps=0.d0

    C0(1:3,1:n)=Cref(1:3,1:n,iref)
    Crot=0.d0
    
    ! Do all permutation of atom indices
    Call InitializePermutation(n,IA,nf,p,s)
    Do ip=1,nfact
        
        Do i=1,n
            C1(1:3,i)=Cref(1:3,IA(i),iref)
        Enddo
        
        Call LSQR(n,C0,C1,Crot,eps)

        If (eps<epsrefmin) Then
            irot=irot+1
            irr(1:n,irot,it)=IA(1:n)
            Write(76,'(2i8,f10.4,i4,2x,<n>i2)')it,ip,eps,irot,irr(1:n,irot,iref)
        Endif

        Call Permutation(n,IA,nf,p,s)
        
    Enddo

iRefRot(0,it)=irot
iRefRot(-1,it)=1


End
!*******************************************************************************************
Subroutine SpectralConversion(iref,f1,f2,t1,t2)    
Use Vars, Only: Freq
Implicit Real(8) (A-H,O-Z)

cc=299792458.d0

f1=Freq(1,iref)
f2=Freq(2,iref)

cm12nm1=1.d7/f1
cm12nm2=1.d7/f2

fHz1=cc/cm12nm1*1.d9
fHz2=cc/cm12nm2*1.d9

t1=1.d15/fHz1
t2=1.d15/fHz2

End
!*******************************************************************************************
Subroutine CalcDist(n,C,rmin,rmax,Diameter)
Implicit Real(8) (A-H,O-Z)

Real(8) C(3,n)

Diameter=0.d0
rmin=999999.d0
rmax=0.d0
Do i=2,n
    rjmin=999999.d0
    Do j=1,i-1
        rij=Distance(i,j,n,C)
        If (rij<rmin) rmin=rij
        If (rij>Diameter) Diameter=rij
        If (rij<rjmin) Then
            jmin=j
            rjmin=rij
        Endif
    Enddo
    If (rjmin>rmax) rmax=rjmin
Enddo

End
!*******************************************************************************************
Subroutine CalcBonds(Numat,C,iref,nb,arbx)   ! Number of bonds and bond lengths
Use Vars, Only: iAMref

Implicit Real(8) (A-H,O-Z)
Real(8) C(3,Numat)

nb=0
arbx=0.d0
Do i=2,Numat
    Do j=1,i-1
        ir=iAMref(i,j,iref)
        If (ir==0) Cycle
        rb=Distance(i,j,Numat,C)
        arbx=arbx+rb
        nb=nb+1
    Enddo
Enddo
arbx=arbx/Dble(nb)

End
!*********************************************************************************************
Subroutine CheckRot(istr,it,irot)
Use Vars, Only: n=>Numat,IAbest,iRefRot,irr

Implicit Real(8) (A-H,O-Z)

Integer(4) Inow(n),Iold(n),inow1(n)

irot=1
If (istr==1) Return
idiff1=Sum(IABS(IAbest(1:n,istr)-IAbest(1:n,istr-1)))
If (idiff1==0) Return

Iold(1:n)=IAbest(1:n,istr-1)    
Inow(1:n)=IAbest(1:n,istr)
Do i=1,n
    ia1=iold(i)
    Do j=1,n
        If (inow(j)==ia1) inow1(j)=i
    Enddo
Enddo

irot=0
nrefrot=iRefRot(0,it)
Do i=1,nrefrot
    If (Sum(IABS(inow1(1:n)-irr(1:n,i,it)))==0) Then
        irot=i
        Exit
    Endif
Enddo

    End
!!**********************************************************    
!Subroutine imul(n,imat,iv1,iv2)
!Implicit Real(8) (A-H,O-Z)
!
!Integer(4) iv1(n),iv2(n),imat(n,n)
!
!Do i=1,n
!    iv2(i)=0
!    Do k=1,n
!        iv2(i)=iv2(i)+imat(i,k)*iv1(k)
!    Enddo
!Enddo
!
!    End

    
!iWordX=0
!iwordtmp=0

!imat=0
!Do i=1,n
!    j=iold(i)
!    imat(j,i)=1
!Enddo

!Call imul(n,imat,iold,iv1)
!Call imul(n,Transpose(imat),inow,iv2)

    
!Do k=1,n
!    ia1=IAbest(k,istr-1)
!    iwordtmp=iwordtmp+10**(n-k)*ia1
!    ia2=IAbest(ia1,istr)
!    ia2=IAbest(k,istr)
!    iWordX=iWordX+10**(n-k)*ia2
!Enddo
!

!If (irot==0) Then
!    iatp(1:n,it)=inow(1:n)
!Endif

    
            
!            iWordX=0
!            Do k=1,Numat
!                ia1=IAbest(k,j-1)
!                ia2=IAbest(ia1,j)
!                iWordX=iWordX+10**(Numat-k)*ia2
!            Enddo
!            
!!            If (iRefRot(-1,it)==0) Call RefPerm(it)
!            nrefrot=iRefRot(0,it)
!            Do i=1,nrefrot
!                If (iWordX==iRefRot(i,it)) Then
!                    irot=i
!                    idiff=0
!                    Exit
!                Endif
!            Enddo
    