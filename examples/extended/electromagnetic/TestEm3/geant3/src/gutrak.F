      SUBROUTINE GUTRAK
* 
*     User routine to control tracking of one track
*     Called by GTREVE
*
#include "geant321/gckine.inc"
#include "geant321/gconst.inc"
#include "geant321/gctrak.inc" 
#include "geant321/gctmed.inc"
#include "geant321/gcvolu.inc"
#include "calor.inc" 
#include "celoss.inc" 
*
* *** initialisation of Eflow for primary track
*
      if ((itra.eq.1).and.(istak.eq.0)) then
*       unique Id for absorber+layer
        Idnow = 1
        if (numed.le.NbAbsor) then
          layer = number(nlevel-1)
	  Idnow = NbAbsor*(layer-1) + numed 
        endif 
*
        Eflow = gekin
        if (ipart.eq.2) Eflow = gekin + 2*emass
*
        do k=1,Idnow
          EnerFlow(Idnow) = EnerFlow(Idnow) + Eflow
        enddo
      endif	      
	                                                                                     
      CALL GTRACK
*
      END 
