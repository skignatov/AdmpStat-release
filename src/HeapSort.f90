!*************************************************************************
Subroutine HeapSort(n,ra,rb,ia)
Implicit Real(8) (A-H,O-Z)

Real(8) ra(n),rb(n,n),tmp(n,n)
Integer(4) ia(n)

! HeapSort algorithm from "Numerical Recipes in F90"
! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is
! input; ra is replaced on output by its sorted rearrangement.

! rb the second nxn array with columns sorted as ra. ia - indexes of permutations

! CITED REFERENCES AND FURTHER READING:
! Knuth, D.E. 1973, Sorting and Searching, vol. 3 of The Art of Computer Programming (Reading,
! MA: Addison-Wesley), x5.2.3. [1]
! Sedgewick, R. 1988, Algorithms, 2nd ed. (Reading, MA: Addison-Wesley), Chapter 11. [2]
! 8.4 Indexing and Ranking

if (n.lt.2) return

Do i=1,n
	ia(i)=i
Enddo

!The index l will be decremented from its initial value down to 1 during the \hiring" (heap
!creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
!down to 1 during the \retirement-and-promotion" (heap selection) phase.

l=n/2+1
ir=n
10 continue
if(l.gt.1)then		!Still in hiring phase.
	l=l-1
	rra=ra(l)
	iia=ia(l)			!*
else				!In retirement-and-promotion phase.
	rra=ra(ir)		!Clear a space at end of array.
	iia=ia(ir)			!*
	ra(ir)=ra(1)	!Retire the top of the heap into it.
	ia(ir)=ia(1)		!*
	ir=ir-1			!Decrease the size of the corporation.
	if(ir.eq.1)then	!Done with the last promotion.
		ra(1)=rra	!The least competent worker of all!
		ia(1)=iia		!*
		Goto 30
	endif
endif

i=l					!Whether in the hiring phase or promotion phase, we here
					!set up to sift down element 
j=l+l				!rra to its proper level.
20 if(j.le.ir)then	!Do while j.le.ir:"
if(j.lt.ir)then
	if(ra(j).lt.ra(j+1))j=j+1		!Compare to the better underling.
endif
if(rra.lt.ra(j))then				!Demote rra.
	ra(i)=ra(j)
	ia(i)=ia(j)		!*
	i=j
	j=j+j
else								!This is rra's level. Set j to terminate the sift-down.
	j=ir+1
endif
goto 20
endif
	
ra(i)=rra			!Put rra into its slot.
ia(i)=iia			!*
goto 10


! Permutations of the second aray
30 Continue
Do i=1,n
	j=ia(i)
	tmp(1:n,i)=rb(1:n,j)
Enddo
rb=tmp

END
!*************************************************************************
Subroutine HeapSort0(n,ra)
Implicit Real(8) (A-H,O-Z)

Real(8) ra(n)
Integer(4) ia(n)

Call HeapSort1(n,ra,ia)

End
!*************************************************************************
Subroutine HeapSortI(n,ira)
Implicit Real(8) (A-H,O-Z)

Real(8) ra(n)
Integer(4) ia(n),ira(n)

ra=Dble(ira)

Call HeapSort1(n,ra,ia)

ira=IDNINT(ra)

End
!*************************************************************************
Subroutine HeapSort1(n,ra,ia)
Implicit Real(8) (A-H,O-Z)

Real(8) ra(n)
Integer(4) ia(n)

! HeapSort algorithm from "Numerical Recipes in F90"
! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is
! input; ra is replaced on output by its sorted rearrangement.

! rb the second nxn array with columns sorted as ra. ia - indexes of permutations

! CITED REFERENCES AND FURTHER READING:
! Knuth, D.E. 1973, Sorting and Searching, vol. 3 of The Art of Computer Programming (Reading,
! MA: Addison-Wesley), x5.2.3. [1]
! Sedgewick, R. 1988, Algorithms, 2nd ed. (Reading, MA: Addison-Wesley), Chapter 11. [2]
! 8.4 Indexing and Ranking

if (n.lt.2) return

Do i=1,n
	ia(i)=i
Enddo

!The index l will be decremented from its initial value down to 1 during the \hiring" (heap
!creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
!down to 1 during the \retirement-and-promotion" (heap selection) phase.

l=n/2+1
ir=n
10 continue
if(l.gt.1)then		!Still in hiring phase.
	l=l-1
	rra=ra(l)
	iia=ia(l)			!*
else				!In retirement-and-promotion phase.
	rra=ra(ir)		!Clear a space at end of array.
	iia=ia(ir)			!*
	ra(ir)=ra(1)	!Retire the top of the heap into it.
	ia(ir)=ia(1)		!*
	ir=ir-1			!Decrease the size of the corporation.
	if(ir.eq.1)then	!Done with the last promotion.
		ra(1)=rra	!The least competent worker of all!
		ia(1)=iia		!*
		Goto 30
	endif
endif

i=l					!Whether in the hiring phase or promotion phase, we here
					!set up to sift down element 
j=l+l				!rra to its proper level.
20 if(j.le.ir)then	!Do while j.le.ir:"
if(j.lt.ir)then
	if(ra(j).lt.ra(j+1))j=j+1		!Compare to the better underling.
endif
if(rra.lt.ra(j))then				!Demote rra.
	ra(i)=ra(j)
	ia(i)=ia(j)		!*
	i=j
	j=j+j
else								!This is rra's level. Set j to terminate the sift-down.
	j=ir+1
endif
goto 20
endif
	
ra(i)=rra			!Put rra into its slot.
ia(i)=iia			!*
goto 10

30 Continue

END










