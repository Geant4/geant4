
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*
#include "geant321/gcflag.inc" 
#include "geant321/gctrak.inc"
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
*
      CHARACTER*20 matnam
      CHARACTER*4 unitdm,unitdr,unitlm,unitlr
*
* *** mean energy deposit and track length
      PRINT 749
      PRINT 750,IEVENT
      PRINT 751 
      fnorm    = 1./IEVENT
      do k=1,NbAbsor
        call GFMATE (materAbs(k),matnam,dua,duz,dud,dur,dui,udu,idu)      
        aveEdep = fnorm*sumEdep(k)
        rmsEdep = fnorm*sqrt(abs(IEVENT*su2Edep(k)-sumEdep(k)**2))
        aveTrck = fnorm*sumTrck(k)
        rmsTrck = fnorm*sqrt(abs(IEVENT*su2Trck(k)-sumTrck(k)**2))
        aveElea = fnorm*sumElea(k)
        rmsElea = fnorm*sqrt(abs(IEVENT*su2Elea(k)-sumElea(k)**2))	
*     
        CALL GEVKEV (aveEdep,Edmean,unitdm)
        CALL GEVKEV (rmsEdep,Edrms ,unitdr)
        CALL GEVKEV (aveElea,Elmean,unitlm)
        CALL GEVKEV (rmsElea,Elrms ,unitlr)	
        PRINT 752,k-1,matnam,Edmean,unitdm,Edrms,unitdr,aveTrck,rmsTrck
     +                      ,Elmean,unitlm,Elrms,unitlr
      enddo
      PRINT 749
*                                             
* *** geant termination
      CALL GLAST
*
* *** close HIGZ
      CALL HPLEND
*                                                                               
* *** Save histo                                                  
      CALL HRPUT(0,fileName,'N')

*      
*
* *** formats
  749 FORMAT(/ ,60(1H-),/)            
  750 FORMAT(1X,'Nb of events:',I6)
  751 FORMAT(1X,'AbsNo  Material',5X,'Energy deposit',14X,'Trck length',
     +      12X,'E leakage by secondaries')
  752 FORMAT(I5,3X,A10,F8.3,A4,' +- ',F6.2,A4,F9.3,' cm +- ',F6.2,' cm',
     +                 F9.3,A4,' +- ',F6.2,A4)   
*              
      END
