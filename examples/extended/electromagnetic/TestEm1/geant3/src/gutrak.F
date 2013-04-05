

      SUBROUTINE GUTRAK                                                         
*                                                                               
*     User routine to control tracking of one track                 
*     Called by GTREVE                                                            
*
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "celoss.inc"     
                                                                    
*                                                                                                                                         
*
* *** initialisation track per track
      NSTRK0 = 0
      NSTRK1 = 0
*                                                                                     
      CALL GTRACK
*
* *** count nb of tracks
      IF (CHARGE.NE.0.) THEN
         NTRAK1 = NTRAK1 + 1
      ELSE 
         NTRAK0 = NTRAK0 + 1
      ENDIF
*
* *** track length of primary particle
      IF ((ITRA.EQ.1).AND.(ISTAK.EQ.0)) THEN
	 CALL HFILL (1,SLENG*10,0.,1.)
	 nstrk = NSTRK0
	 if (CHARGE.NE.0.) nstrk = NSTRK1
         CALL HFILL (2,REAL(nstrk),0.,1.)
      ENDIF	                                                                                    
*                                                                               
      END                                                                            
