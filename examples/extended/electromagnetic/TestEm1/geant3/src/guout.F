                                                     
      SUBROUTINE GUOUT                                                          
*                                                                               
*     User routine called at the end of each event                              
*                                                                               
#include "geant321/gcflag.inc"
#include "celoss.inc"
*
      CHARACTER*4 unit                                                              
*                                                                               
* *** drawing
*
#ifndef batch                                                         
        IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'TestEm1',0.)
         CALL GDSHOW (3)
         CALL GDXYZ  (0)
      END IF
#endif      
*
*     energy deposited
*                                                                                
      CALL GEVKEV (ETOT,ETunit,unit)                                                 
      IF (ISWIT(1).NE.0) PRINT 750,ETunit,unit
*      
  750 FORMAT(1X,'Energy deposit: ',F8.3,A4)    
*                                                                                
      END                                                                       
