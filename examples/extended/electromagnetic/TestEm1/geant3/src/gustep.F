
      SUBROUTINE GUSTEP                                                         
*                                                                               
*     User routine called at the end of each tracking step           
*                                                               
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "celoss.inc"                                                             
*                                                                               
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
* *** Energy deposited                                                          
      ETOT = ETOT + DESTEP
*
* *** count number of steps
      IF (CHARGE.NE.0.) THEN
           NSTEP1 = NSTEP1 + 1
	   NSTRK1 = NSTRK1 + 1
	   CALL HFILL (3,STEP*10,0.,1.)
      ELSE 
           NSTEP0 = NSTEP0 + 1
	   NSTRK0 = NSTRK0 + 1
      ENDIF                                                          
*
* *** count nb of invoked processes
      IF (NMEC.GT.0) THEN
        DO IM = 1,NMEC
          IPROC = LMEC(IM)
          IF (IPROC.EQ.21) IPROC = 12
          IF (IPROC.LE.12) NBCALL(IPROC) = NBCALL(IPROC)+1
        ENDDO
      ENDIF
         
      END                                                                       
