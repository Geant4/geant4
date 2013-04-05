
      SUBROUTINE GRUN
*.
*.    ******************************************************************
*.    *                                                                *
*.    *       Steering routine to process all the events               *
*.    *                                                                *
*.    *    ==>Called by : <USER>, main program                         *
*.    *       Author    R.Brun  *********                              *
*.    *                                                                *
*.    ******************************************************************
*.

#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctime.inc"

      SAVE IFIRST
      DATA IFIRST/0/
*.
*.    ------------------------------------------------------------------
*.
*             Keep starting time
*
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL TIMEL(TIMINT)
      ENDIF
*
   10 IF(IEVENT.LT.NEVENT) THEN
         IEVENT=IEVENT+1
*
*             Initialises event partition
*
         CALL GTRIGI
*	 
*             Process one event (trigger)
*
         CALL GTRIG
*
*             Clear event partition
*
         CALL GTRIGC
*
         IF(IEORUN.EQ.0) THEN
*
*             Check time left
*
            go to 10
            IF(ITIME.LE.0)GO TO 10
            IF(MOD(IEVENT,ITIME).NE.0)GO TO 10
            CALL TIMEL(TIMNOW)
            IF(TIMNOW.GT.TIMEND)GO TO 10
            WRITE(CHMAIL,10000)TIMEND
            CALL GMAIL(0,2)
            IEORUN = 1
         ENDIF
      ENDIF
*
10000 FORMAT(5X,'***** THE JOB STOPS NOW BECAUSE THE TIME LEFT IS LESS',
     +' THAN     ',F8.3,' SECONDS *****')
      END
