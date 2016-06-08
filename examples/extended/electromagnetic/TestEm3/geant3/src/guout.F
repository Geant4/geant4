
      SUBROUTINE GUOUT 
*
*     User routine called at the end of each event
*     
#include "geant321/gcflag.inc"
#include "calor.inc"
#include "celoss.inc"
*
      DIMENSION sumevent(2*MaxAbS)
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
*     ntuple: energies and track length per event
      do k=1,NbAbsor
        sumevent(2*k-1) = EdepAbs(k)
        sumevent(2*k  ) = TrckAbs(k)
      enddo
      CALL HFN(1,sumevent)     
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
