
      SUBROUTINE GUSTEP
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconst.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
*
#include "runinfo.inc"
#include "histo.inc"
*                                                                               
*                                                                               
* *** Debug event and store tracks for drawing
      IF (IDEBUG.NE.0)   CALL GPCXYZ      
      IF ((ISWIT(1).EQ.1).AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF  (ISWIT(1).EQ.2)                     CALL GSXYZ                        
*
* *** collect informations in the absorber only
*
      if (numed.ne.1) return
*
*     energy deposit
      edep = edep + destep
*
*     number of steps and total track length
      if (charge.ne.0.) then
        stpch = stpch + 1.
	trkch = trkch + step
      else
        stpne = stpne + 1.
	trkne = trkne + step
      endif	      	            
*
*     manage secondaries
      if (ngkine.gt.0) then
        do k = 1,ngkine
	  energy = gkin(4,k)
          id     = gkin(5,k)
          if (id.eq.1) then
	    ngamma = ngamma + 1
	    if (histo(3)) call hfill(3,log10(1000*energy),0.,1.)
	  else if (id.eq.2) then
	    nposit = nposit + 1
	    energy = energy - emass
	    if (histo(2)) call hfill(2,energy/histUnit(2),0.,1.)
	  else if (id.eq.3) then
	    nelect = nelect + 1
	    energy = energy - emass
	    if (histo(2)) call hfill(2,energy/histUnit(2),0.,1.)
	  endif  	    
	enddo
*       keep or kill all secondaries	
	if (.not.kill) call gsking(0)
      endif	
*
      END
