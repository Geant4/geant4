
      SUBROUTINE GUTREV 
*
*     User routine to control tracking of one event
*     Called by GRUN
*
#include "calor.inc"
#include "celoss.inc"
*
      do k=1,NbAbsor
        EdepAbs(k) = 0.
        TrckAbs(k) = 0.
      enddo                                                               
*
      CALL GTREVE                                                               
* 
      END                                                                       
