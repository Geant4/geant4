
      SUBROUTINE GUOUT
*
*     called at end of event
*
#include "geant321/gcflag.inc"
*
#include "runinfo.inc"
#include "histo.inc"
*
* *** energy in absorber
      if (histo(1)) call hfill (1,edep/histUnit(1),0.,1.)
      sedep  = sedep  + edep
      sedep2 = sedep2 + edep**2      
*
* *** number of steps
      sstpch  = sstpch  + stpch
      sstpch2 = sstpch2 + stpch**2
*      
      sstpne  = sstpne  + stpne
      sstpne2 = sstpne2 + stpne**2            
*
* *** track length
      strkch  = strkch  + trkch
      strkch2 = strkch2 + trkch**2
*      
      strkne  = strkne  + trkne
      strkne2 = strkne2 + trkne**2
*
* *** transmission, reflection counters
      if (itransmi.eq.1) then
        nbtransm(1) = nbtransm(1) + 1
      else if (itransmi.eq.2) then
        nbtransm(1) = nbtransm(1) + 1
	nbtransm(2) = nbtransm(2) + 1
      endif
*      
      if (ireflect.eq.1) then
        nbreflec(1) = nbreflec(1) + 1
      else if (ireflect.eq.2) then
        nbreflec(1) = nbreflec(1) + 1
	nbreflec(2) = nbreflec(2) + 1
      endif      		      
*
#ifndef batch
*                                                                               
* *** drawing
*                                                          
      IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'TestEm5',0.)
         CALL GDSHOW (3)
         CALL GDXYZ  (0)
      END IF
#endif      
      END
