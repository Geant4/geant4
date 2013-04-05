
      SUBROUTINE UHINIT
*
*
#include "histo.inc"
*
      CHARACTER*50 title(7)
*
      data title /
     1    'Edep (MeV/mm) along absorber',	      
     2    'total energy deposit in absorber',					 
     3    'true track length of the primary particle',
     4    'true step size of the primary particle',
     5    'projected range of the primary particle',
     6    'true track length of charged secondaries',
     7    'true step size of charged secondaries'       /
*
      if (idhist.gt.MaxHist) return
*      
      if (histo(idhist)) call hdelet(idhist)
*      	
      call hbook1(idhist,title(idhist),nbBins,valmin,valmax,0.)
*      
      histo   (idhist) = .true.
      binWidth(idhist) = (valmax-valmin)/nbBins
      if (valunit.le.0.) valunit = 1.                
      histUnit(idhist) = valunit        
*
      END
