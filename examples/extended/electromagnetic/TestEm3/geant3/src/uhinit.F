
      SUBROUTINE UHINIT
*
*     To book the user's histograms
*
#include "histo.inc"
*
      character*25 title
*
      if (histo(idhist)) call hdelet(idhist)      
*
      title = 'Edep in absorber '//char(ichar('0')+(idhist-1))
*      
      call hbook1(idhist,title,nbBins,valmin,valmax,0.)
*      
      histo   (idhist) = .true.
      binWidth(idhist) = (valmax-valmin)/nbBins
      if (valunit.le.0.) valunit = 1.                
      histUnit(idhist) = valunit              	   
*                                                                               
      END                                                                       
