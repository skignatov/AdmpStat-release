!
! Compiler version:
!  _WIN_=1    : OS type: 1-Windows, 0-LINUX
!  _OMP_=1 : OpenMP for Intel Compiler (both Win and Linux): 0 -off, 1 - on  (requires /Qopenmp or --qopenmp options Professional)
!  _MPI_=0 : Intel Fortran 12.1 or higher for Windows or Linux
!  _MPI_=1 : Linux MPI compiler (MPIF90 -fc=ifort -f90=ifort)
!  _IFORT_ : Compiler type: Intel Fortran - 1, Gfortran - 0
!  No IMSL, no QWIN
!
!DEC$ DEFINE _WIN_=1
!DEC$ DEFINE _OMP_=0
!DEC$ DEFINE _MPI_=0
!DEC$ DEFINE _IFORT_=1
!DEC$ DEFINE _IMSL_=0
!DEC$ DEFINE _QWIN_=0
