
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
* *** Energy flow : leaving an absorber ?
      IF ((NUMED.LE.NbAbsor).and.(INWVOL.eq.2)) THEN
        Idnow = NbAbsor*(layer-1) + numed
        Eflow = gekin
        if (ipart.eq.2) Eflow = gekin + 2*emass
*		      
        if (max(abs(vect(2)),abs(vect(3))).ge.(0.5*calorYZ)) then
	  EleakLat(Idnow)   = EleakLat(Idnow)   + Eflow
        elseif (vect(4).ge.0.) then
	  EnerFlow(Idnow+1) = EnerFlow(Idnow+1) + Eflow
        else
	  EnerFlow(Idnow)   = EnerFlow(Idnow)   - Eflow	  
        endif
      ENDIF	     
*
      END
