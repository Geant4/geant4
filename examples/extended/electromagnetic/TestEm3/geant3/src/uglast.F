
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*                                                                               
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
#include "calor.inc"
#include "celoss.inc"
*
      CHARACTER*20 matnam
      CHARACTER*4 unitm,unitr
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
*     
        CALL GEVKEV (aveEdep,Emean,unitm)
        CALL GEVKEV (rmsEdep,Erms ,unitr)
        PRINT 752,k-1,matnam,Emean,unitm,Erms,unitr,aveTrck,rmsTrck
      enddo
      PRINT 749
*                                             
* *** geant termination
      CALL GLAST
*
* *** close HIGZ
      CALL HPLEND
*                                                                               
* *** Save ntuple                                                  
      IF (ISWIT(2).EQ.1) THEN                                                   
         CALL HRPUT(0,'testem3.paw','N')
      ENDIF                                                                     
*      
*
* *** formats
  749 FORMAT(/ ,60(1H-),/)            
  750 FORMAT(1X,'Nb of events:',I6)
  751 FORMAT(1X,'AbsNo  Material',5X,'Energy deposit',14X,'Trck length')
  752 FORMAT(I5,3X,A10,F8.3,A4,' +- ',F6.2,A4,F9.3,' cm +- ',F6.2,' cm')     
*              
      END
