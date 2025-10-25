Subroutine BuildNumber(nbuild)

Logical res,res1
Character(255) FilObj,FilInc
!Character(36) str

Include 'BuildNumber.fi'

!!DEC$ IF DEFINED (_DEBUG)    ! Determine if the operation in the Debug configuration (this line is not obligatory)

! BuildNumber saves the build number of the current program version
! using the BuilNumber.fi file
!
! Usage in the program (no additional code is needed):
!
!		.....
!		Call Build Number(nBuild)
!		Write(*,'('' Build number: '',i5')nBuild
!		....
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Before use, adjust the RELATIVE paths where two your signal files are located! (relatively to the project directory)
! WARNING! Program changes BuildNumber if your executable is compiled in the Debug configuration only
! WARNING! It does not change BuildNumber if your working directory is different from the initial project directory (even if you run it under VisualStudio API!)
! WARNING! If absolute paths are indicated, the Debug-executable will change BuildNumber being executed from any directory!
FilObj='.\x64\Debug\BuildNumber.obj'
FilInc='.\BuildNumber.fi'		
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

INQUIRE(file=FilObj,EXIST=res)
If (res) Then                   ! Program started under VisualStudio in the project directory (or compiled earlier but was not running and transfreed to this dir)
    INQUIRE(file=FilInc,EXIST=res1)
    If (res1) Then              ! Include file *.fi exists in the current durectory
        Open(98,file=FilInc)
        Read(98,'(10x,i6)')nBuild0
        If (nBuild0==nbuild) Then
            nbuild=nbuild+1
            Rewind(98)
            Write(98,'(''nBuild   ='',i6)')nbuild
            Close(98)
        Endif
    Endif
    Open(98,file=FilObj)
    Close(98,status='delete')
Else                ! Program started from any other dir not under VS (in a normal working mode)
    nBuild=nBuild+1 ! Correct included code which is less by 1 than the just compiled nBuild 
    
Endif

!!DEC$ ENDIF

End



