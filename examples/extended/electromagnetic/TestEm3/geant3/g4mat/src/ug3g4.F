
      SUBROUTINE UG3G4                                                         
*
*                                                                               
#include "geant321/gcmulo.inc"
#include "geant321/gckine.inc"
#include "ugmate.inc"
#include "ggvalues.inc"
*
*
* *** fetch Geant3 tables      
      call gftmat(imat,ikine,'LOSS',nekbin,elow,g3value,dum,idum)
*
* *** histograms to plot the difference
      elow(nek1) = ekmax
      call HBOOKB (1,'(g4-g3)/g3 in percent',nekbin,elow,0.)
      call HBOOKB (2,'abs((g4-g3)/g3) in percent',nekbin,elow,0.)           
*
* *** compute difference
      do i=1,nekbin
        dif = 100*(g4value(i)-g3value(i))/g3value(i)
	call hfill (1, elow(i), 0., dif)
	call hfill (2, elow(i), 0., abs(dif))	
        print *,i,elow(i),g4value(i),g3value(i), 'dif(%)=',dif
      enddo      
*                                                                               
      END                                                                       
