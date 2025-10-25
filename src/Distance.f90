!***********************************************************************
Function Distance(i,j,N,C)

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,N)

! Subroutine Distance calculates the interatomic distance between two atoms i and j in molecule C(3,MaxAt)

dx=C(1,i)-C(1,j)
dy=C(2,i)-C(2,j)
dz=C(3,i)-C(3,j)

Distance=DSQRT(dx*dx+dy*dy+dz*dz)

End
