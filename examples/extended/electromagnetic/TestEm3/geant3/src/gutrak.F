      SUBROUTINE GUTRAK
* 
*     User routine to control tracking of one track
*     Called by GTREVE 
#include "geant321/gctmed.inc"
#include "geant321/gcvolu.inc"

#include "calor.inc" 
#include "celoss.inc" 
*
* *** initialisation track per track
*     unique Id for absorber+layer
      Idold = 0
      if (numed.le.NbAbsor) then
        layer = number(nlevel-1)
	Idold = NbAbsor*(layer-1) + numed 
      endif 
*                                                                                     
      CALL GTRACK
*
      END 
