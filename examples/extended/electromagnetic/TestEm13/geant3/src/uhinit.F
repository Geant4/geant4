
      SUBROUTINE UHINIT
*
*
#include "histo.inc"
*
      CHARACTER*50 title(6)
*
      data title /
     1    'scattered primary particle: energy spectrum',					 
     2    'scattered primary particle: costheta distribution',
     3    'charged secondaries: energy spectrum',					 
     4    'charged secondaries: costheta distribution',
     5    'neutral secondaries: energy spectrum',					 
     6    'neutral secondaries: costheta distribution'     /          			 	 
*
      if (histo(idhist)) call hdelet(idhist)
*      
      vmin = valmin
      vmax = valmax      	
      call hbook1(idhist,title(idhist),nbBins,vmin,vmax,0.)
*      
      histo   (idhist) = .true.
      binWidth(idhist) = (valmax-valmin)/nbBins
      if (valunit.le.0.) valunit = 1.                
      histUnit(idhist) = valunit        
*
      END
