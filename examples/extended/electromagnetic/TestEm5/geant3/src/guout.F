
      SUBROUTINE GUOUT
C

#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconst.inc"
#include "edepo.inc"
*                                                                               
      dimension runerr(2)
* *** drawing
*                                                          
#ifndef batch
      IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'TestEm5',0.)
         CALL GDSHOW (3)
         CALL GDXYZ  (0)
      END IF
#endif      
*
      call hfill(1,anbstep,0.,1.)
      sstep=sstep+anbstep
      s2step=s2step+anbstep*anbstep
      sanch = sanch+anch
      sanch2=sanch2+anch*anch
      sanne = sanne+anne
      sanne2=sanne2+anne*anne
      sanele=sanele+anele
      sanpos=sanpos+anpos
        SEDEP=SEDEP+EDEP
        SEDEP2=SEDEP2+EDEP**2
      stra=stra+tra
      sref=sref+ref
      eein=1000.*EDEP
      CALL HFILL(2,eein,0.,1.0)
C
      END
