Subroutine PrintTitle(iUnit,nBuild)

Write(iUnit,101)nbuild

101 FORMAT(13x,'*******************************************************'/ &
		   13x,'*                   ADMPstat v.0.1                    *'/ &
		   13x,'*               (build number',i6,')                  *'/ &
		   13x,'*                                                     *'/ &
		   13x,'*             Statistics of ADMP trajectory           *'/ &
		   13x,'*                                                     *'/ &
		   13x,'*      Originally written by Stanislav Ignatov        *'/ &
		   13x,'*          Theoretical Chemistry Group                *'/ &
		   13x,'* N.I.Lobachevsky State University of Nizhny Novgorod *'/ &
		   13x,'*                 ignatov@ichem.unn.ru                *'/ &
		   13x,'*               http://ichem.unn.ru/tcg               *'/ &
		   13x,'*             Nizhny Novgorod, Russia, 2023           *'/ &
		   13x,'*******************************************************'/)


End