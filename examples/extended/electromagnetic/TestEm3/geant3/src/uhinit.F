
      SUBROUTINE UHINIT
*
*     To book the user's histograms
*
#include "calor.inc"
#include "histo.inc"
*
      character*25 title
*
* *** complete the booking of histograms
      do k=1,NbAbsor
        if (iexist(k).eq.0) then
	   title = 'Edep/Ebeam in absorber '//char(ichar('0')+(k-1))
	   call hbook1(k,title,100, 0.,1.,0.)
	   iexist(k) = 1
	endif
      enddo	   	   
*                                                                               
      END                                                                       
