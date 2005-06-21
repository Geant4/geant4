
      SUBROUTINE GUSTEP
*                                                                               
*     User routine called at the end of each tracking step
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
*
      common/local1/ fnorm
*      
      data amm, aMeV /1.e-1, 1.e-3/
*
      if (ievent.eq.1) fnorm = amm/(nevent*aMeV*binWidth(1)*histUnit(1))                                                                         
*
* *** Debug event and store tracks for drawing
      IF (IDEBUG.NE.0)   CALL GPCXYZ
      IF (IDEBUG.NE.0)   CALL GPGKIN            
      IF ((ISWIT(1).EQ.1).AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF  (ISWIT(1).EQ.2)                     CALL GSXYZ      
*                                                                               
* *** Something generated ?                                                     
      IF(NGKINE.GT.0) CALL GSKING(0)
*
      if (destep.le.0.) return                                                   
*                                                                               
* *** Energy deposited per event                                                          
      Eevent = Eevent + destep
*
* *** histo edep longitudinal profile
      xlen = 0.5*sizeX + vect(1)
      ih = 1
      if (histo(ih))
     &   call hfill (ih,xlen/histUnit(ih),0.,fnorm*destep)
*
* *** histo step size
      ih = 0
      if (istak.eq.0) ih = 4
      if (istak*charge.ne.0.) ih = 7
      if ((ih.gt.0).and.(histo(ih)))
     &   call hfill (ih,step/histUnit(ih),0.,1.)
*               
      END                                                                       
