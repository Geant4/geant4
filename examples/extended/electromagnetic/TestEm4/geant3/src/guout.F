
      SUBROUTINE GUOUT                                                          
*                                                                               
*     User routine called at the end of each event                              
*                                                                               
#include "geant321/gcflag.inc"    
#include "celoss.inc"                                                               
*                                                                               
* *** drawing
*
#ifndef batch                                                          
        IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'TestEm4',0.)
         CALL GDSHOW (3)
         CALL GDXYZ  (0)
      END IF
#endif     
*
*     energy deposited
*                                                                                
      EMEV = 1000.*ETOT                                                 
      CALL HFILL(1,EMEV,0.,1.)
      IF (ISWIT(1).NE.0) PRINT 750,EMEV
*      
  750 FORMAT(1X,'Energy deposit: ',F6.3,' MeV')    
*                                                                                
      END                                                                       
