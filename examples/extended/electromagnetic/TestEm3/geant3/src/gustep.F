
      SUBROUTINE GUSTEP                                                         
*                                                                               
*     User routine called at the end of each tracking step           
*                                                               
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconst.inc"
#include "geant321/gctrak.inc"
#include "geant321/gctmed.inc"
#include "calor.inc"
#include "celoss.inc"
*
*
* *** Debug event and store tracks for drawing
      IF (IDEBUG.NE.0)   CALL GPCXYZ      
      IF ((ISWIT(1).EQ.1).AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF  (ISWIT(1).EQ.2)                     CALL GSXYZ
*
* *** Something generated ?
      IF(NGKINE.GT.0) CALL GSKING(0)
*
* *** Energy deposited and charged track length
      IF (NUMED.LE.NbAbsor) THEN
        if (DESTEP.NE.0.)      EdepAbs(NUMED) = EdepAbs(NUMED) + DESTEP
        if (CHARGE*STEP.NE.0.) TrckAbs(NUMED) = TrckAbs(NUMED) + STEP
      ENDIF
*
* *** Energy leakage from secondaries
      IF ((NUMED.LE.NbAbsor).and.(INWVOL.EQ.2).and.(ISTAK.GT.0)) THEN
        Eleav = gekin
	if (ipart.eq.2) Eleav = gekin + 2*emass
	EleaAbs(NUMED) = EleaAbs(NUMED) + Eleav 
      ENDIF	     
*
      END
