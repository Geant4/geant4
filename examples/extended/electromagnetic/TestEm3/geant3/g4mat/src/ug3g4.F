
      SUBROUTINE UG3G4                                                         
*
*                                                                               
#include "geant321/gcmulo.inc"
#include "geant321/gckine.inc"
#include "ugmate.inc"
#include "ggvalues.inc"
*
      dimension dum(5)
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
      print 750
      do i=1,nekbin
        dif = 100*(g4value(i)-g3value(i))/g3value(i)
	call hfill (1, elow(i), 0., dif)
	call hfill (2, elow(i), 0., abs(dif))	
        print 751,i,elow(i),g4value(i),g3value(i),dif
      enddo
*
  750 format(/,6X,'Tkin (GeV)',3X,'G4 (MeV/cm)',3X,'G3 (MeV/cm)',
     +         4X,'G4-G3 (%)',/)
  751 format(1X,I3, E12.4, 2X, E12.4, 2X, E12.4, 2X, E12.4) 
*                                                                               
      END                                                                       
