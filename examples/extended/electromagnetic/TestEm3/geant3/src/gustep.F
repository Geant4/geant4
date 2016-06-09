
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
#include "geant321/gcvolu.inc"
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
*
      data aMeV /1.e-3/
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
*       longitudinal profile of Edep
        layer = number(nlevel-1)
	ih = MaxAbs + numed
	if (histo(ih)) call hfill (ih,0.+layer,0.,destep/(aMeV*nevent)) 	
      ENDIF
*
* *** Energy flow
      IF ((NUMED.LE.NbAbsor).and.(INWVOL.EQ.2)) THEN
        plane = 1. + NbAbsor*(layer-1) + numed
        Eleav = gekin
	if (ipart.eq.2) Eleav = gekin + 2*emass
        ih = 2*MaxAbs + 1
	if (istak.gt.0) ih = 2*MaxAbs + 2
	if (histo(ih)) call hfill (ih, plane, 0., Eleav/(aMeV*nevent))
      ENDIF	     
*
      END
