Subroutine FileExtension(Fil,Ext)
Character(*) Fil,Ext

! Subroutine FileExtension changes or adds the file extension Ext
! to the file name Fil 
! Fil should have a format:  filename.old or filename
! Ext should have a format:  .ext
 
ind=INDEX(Fil,'.',Back=.true.)
If (ind==0) Then
	Fil=Trim(Fil)//Trim(Ext)
Else
	le=Len_Trim(Ext)
	Fil(ind:ind+le-1)=Trim(Ext)
Endif

    End
!*******************************************************************
Subroutine GetFileName(mode,Fil,FilName)

Character(*) Fil,FilName

! Subroutine GetFileName extracts FileName from Fil (dot-separated FileName.Extension)
! If mode >= 0 dot is sought from the beginning of Fil 
!    mode <  0 dot is sought from the end of Fil

If (mode>=0) Then
    ind=INDEX(Fil,'.',Back=.true.)
Else
    ind=INDEX(Fil,'.',Back=.true.)
Endif
FilName=Fil(1:ind-1)

End

