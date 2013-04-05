
      SUBROUTINE UHINIT
*
*     To book the user's histograms
*
#include "histo.inc"
#include "calor.inc"
*
      character*50 title
*
      if (idhist.gt.MaxHist) return
*      
      if (histo(idhist)) call hdelet(idhist)      
*
      if (idhist.lt.MaxAbs)
     +  title = 'Edep in absorber '//char(ichar('0')+(idhist))
      if ((idhist.gt.MaxAbs).and.(idhist.lt.2*MaxAbs))
     +  title = 'Edep long. profile in absorber '
     +          //char(ichar('0')+(idhist-MaxAbs))//' (MeV/event)'
      if (idhist.eq.(2*MaxAbs+1))
     +  title = 'energy flow (MeV/event)'
      if (idhist.eq.(2*MaxAbs+2))
     +  title = 'lateral energy leak (MeV/event)'               
*      
      call hbook1(idhist,title,nbBins,valmin,valmax,0.)
*      
      histo   (idhist) = .true.
      binWidth(idhist) = (valmax-valmin)/nbBins
      if (valunit.le.0.) valunit = 1.                
      histUnit(idhist) = valunit              	   
*                                                                               
      END                                                                       
