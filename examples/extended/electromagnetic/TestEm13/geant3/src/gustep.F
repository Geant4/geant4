
      SUBROUTINE GUSTEP                                                         
*                                                                               
*     User routine called at the end of each tracking step           
*                                                               
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
*
#include "process.inc"
#include "histo.inc"
*
      character*20 cdum                                                             
*                                                                                                                                                              
* *** Debug event and store tracks for drawing
      IF (IDEBUG.NE.0)   CALL GPCXYZ
      IF (IDEBUG.NE.0)   CALL GPGKIN            
      IF ((ISWIT(1).EQ.1).AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF  (ISWIT(1).EQ.2)                     CALL GSXYZ
*
* *** if no process: return
      IF (NMEC.EQ.0) return                                                                                    
*
* *** count nb of invoked processes
      DO IM = 1,NMEC
        IPROC = LMEC(IM)
        IF (IPROC.EQ.21) IPROC = 12
        IF (IPROC.LE.12) NBCALL(IPROC) = NBCALL(IPROC)+1
      ENDDO
*
* *** sum track length for discrete processes
      if ((iproc.gt.5)) then
        nbTot    = nbTot + 1
        sumTrak  = sumTrak  + sleng
        sumTrak2 = sumTrak2 + sleng*sleng
      endif
*
* *** plot final state
*
*     scattered primary particle (if still alive)
      if (istop.eq.0) then
        id = 1
        if (histo(id)) call hfill (id,gekin/histUnit(id),0.,1.)
	id = 2
	if (histo(id)) call hfill (id,vect(4),0.,1.)
      endif
*     
* *** secondaries  
      if (ngkine.gt.0) then
        do lp = 1,ngkine
	  ipar = gkin(5,lp) + 0.1
	  call gfpart(ipar,cdum,ndum,gmass,gcharg,dum,dum,ndum)
	  ekin = gkin(4,lp) - gmass
	  pc   = sqrt(gkin(1,lp)**2 + gkin(2,lp)**2 + gkin(3,lp)**2)
	  cost = gkin(1,lp)/pc
	  if (gcharg.ne.0.) id = 3
	  if (gcharg.eq.0.) id = 5
	  if (histo(id)) call hfill (id,ekin/histUnit(id),0.,1.)
	  id = id + 1
	  if (histo(id)) call hfill (id,cost,0.,1.)
        enddo
      endif            
*      
* *** stop the tracking
      istop = 1      		     
*
      END 
