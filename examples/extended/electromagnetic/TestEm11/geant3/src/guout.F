                                                     
      SUBROUTINE GUOUT                                                          
*                                                                               
*     User routine called at the end of each event                              
*                                                                               
#include "geant321/gcflag.inc"
*                                                         
*                                                                               
* *** drawing
*
#ifndef batch                                                         
        IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'TestEm11',0.)
         CALL GDSHOW (3)
         CALL GDXYZ  (0)
      END IF
#endif      
*
*                                                                                
      END                                                                       
