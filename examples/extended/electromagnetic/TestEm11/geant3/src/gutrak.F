
      SUBROUTINE GUTRAK                                                         
*                                                                               
*     User routine to control tracking of one track                 
*     Called by GTREVE                                                            
*
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "celoss.inc" 
#include "histo.inc"
*
* *** per track initialisations    
*
* *** process a track                                                                                     
      CALL GTRACK
*
* *** count track length of primary particle
      if (istak.eq.0) then 
        TrLen  = TrLen + sleng
	TrLen2 = TrLen2 + sleng*sleng
*	
	xLen    = vect(1) - vert(1)
	XTrLen  = XTrLen  + xLen
	XTrLen2 = XTrLen2 + xLen*xLen
*
        NbStep  = NbStep  + nstep
	NbStep2 = NbStep2 + nstep*nstep
*
        stepsize = sleng/nstep
	StepLen  = StepLen  + stepsize
	StepLen2 = StepLen2 + stepsize*stepsize			
      endif
*
* *** absorbed, transmited, reflected ?
      if (istak.eq.0) then 
        i = 1
	if (inwvol.eq.3) then
	  i = 2
	  if (vect(4).le.0.) i = 3
	endif
        kcoef(i) = kcoef(i) + 1
      endif
*
* *** histo track length
      ih = 0
      if (istak.eq.0) ih = 3
      if (istak*charge.ne.0.) ih = 6
      if ((ih.gt.0).and.(histo(ih)))
     &   call hfill (ih,sleng/histUnit(ih),0.,1.)
*
* *** histo projected range of primary
      if (istak.eq.0) then 
        ih = 5
        if (histo(ih)) call hfill (ih,xLen/histUnit(ih),0.,1.)
      endif                                                                                             
*                                                                               
      END                                                                            
