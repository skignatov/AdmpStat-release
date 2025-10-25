Subroutine Timer(it,mode,iu,Time)
Use IFPORT, Only: RTC
Implicit Real(8) (A-H,O-Z)

Integer(4), parameter:: MaxTT=10    ! Number of timers
Real(8), save:: TT0(MaxTT)
Integer(4) idt(MaxTT)

! mode=0  - initalize timer it or all timers (it=0)
! mode=1  - get value of timer it (seconds since last initialization) to TT
! mode=2  - print current time on unit IU
! mode=3  - print calculation time on unit IU

If (it<0.or.it>MaxTT) it=1

If (mode==0) Then ! Initialize timers
    TT=RTC()
    If (it==0) Then
        TT0=TT
    Else
        TT0(it)=TT
    Endif
    Return
ElseIf (mode==1.or.mode==3) Then
    TT=RTC()
    If (it==0) it=1
    Time=TT-TT0(it)
    If (mode==3) Then
        Write(iu,901)Time,Time/60.d0,Time/3600.d0,Time/86400.d0    
    Endif
ElseIf (mode==2) Then
    Call DATE_AND_TIME(VALUES=idt)
    iY=idt(1)
    iM=idt(2)
    iD=idt(3)
    iHrs=idt(5)
    iMin=idt(6)
    iSec=idt(7)
    iHS=idt(8)/10
    WRITE(iu,900)ID,IM,IY,IHrs,IMin,ISec,IHS
Endif

900 FORMAT(//' Date:',2(I0.2,'/'),I0.4,T63,'Time:',I0.2,2(':',I0.2),'.',I0.2/)
901 Format(' Calculation time:',f10.2,' sec   = ',f8.2,' min   = ',f6.2,' hrs  =',f6.2,' days')

End
