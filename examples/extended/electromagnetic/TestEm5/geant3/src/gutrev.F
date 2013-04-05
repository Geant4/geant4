
      SUBROUTINE GUTREV
*
*
*     User routine to control tracking of one event
*     Called by : GTRIG
*
#include "runinfo.inc"
*
* *** initialisation event by event
*
*     energy deposit in absorber
      edep = 0.
*
*     number of steps in absorber
      stpch = 0.
      stpne = 0.
*
*     track length in absorber
      trkch = 0.
      trkne = 0.
*
*     transmission, reflection flag
      itransmi = 0
      ireflect = 0                  
*
* *** perform the tracking of the event
*      
      CALL GTREVE
*
* *** the 'end of event' has been done in GUOUT
*
      END
