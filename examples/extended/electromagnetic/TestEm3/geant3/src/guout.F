
      SUBROUTINE GUOUT 
*
*     User routine called at the end of each event
*     
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "calor.inc"
#include "celoss.inc"
*
      CHARACTER*4 unit
*                                                                               
* *** drawing
*                                                          
#ifndef batch
      IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'TestEm3',0.)
         CALL GDSHOW (3)
         CALL GDXYZ  (0)
      END IF
#endif      
*
*     histo: energy deposit per tracking medium
      do k=1,NbAbsor
        ratio = EdepAbs(k)/pkine(1)
	call hfill (k,ratio,0.,1.)
      enddo
*
*     cumul energy deposited and charged track length
      do k=1,NbAbsor
        sumEdep(k) = sumEdep(k) + EdepAbs(k)
        sumTrck(k) = sumTrck(k) + TrckAbs(k)
*
        su2Edep(k) = su2Edep(k) + EdepAbs(k)**2
        su2Trck(k) = su2Trck(k) + TrckAbs(k)**2
      enddo      
*
* *** print
      IF (ISWIT(1).NE.0) THEN
        do k=1,NbAbsor
          CALL GEVKEV (EdepAbs(k),Eunit,unit)
          PRINT 751,k,Eunit,unit,TrckAbs(k)
	enddo
	print 749
      ENDIF
          
  749 FORMAT(/) 
  751 FORMAT(1X,'Absorber',I2,' : Energy deposit: ',F8.3,A4,
     &                        '   Track length: '  ,F8.3,' cm')   
*
      END
