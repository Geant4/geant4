
      SUBROUTINE GUSTEP                                                         
*                                                                               
*     User routine called at the end of each tracking step           
*                                                               
#include "geant321/gcflag.inc"
#include "geant321/gconst.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "pvolum.inc"
#include "celoss.inc"                                                             
*
*                                                                               
* *** Debug event and strore track for drawing
      IF (IDEBUG.NE.0) CALL GPCXYZ
      IF (ISWIT(1).EQ.1.AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF (ISWIT(1).EQ.2)                    CALL GSXYZ                                                         
*                                                                               
* *** Something generated ?                                                     
      IF(NGKINE.GT.0) CALL GSKING(0)                                            
*                                                                               
* *** Energy deposited                                                          
      IF (DESTEP.GT.0.)THEN                                                     
         radius = sqrt(vect(1)**2 + vect(2)**2)
	 NR = 1 + radius/(DRX0*X0)
	 NL = 1 + (Z1 + vect(3))/(DLX0*X0)
         DEDR(NR) = DEDR(NR) + DESTEP                                                
         DEDL(NL) = DEDL(NL) + DESTEP                                                                                    
      ENDIF
*
* *** track length and total energy deposit                                                                         
      IF (CHARGE.NE.0.) THEN
         STRCH  = STRCH  + STEP
	 EDEPCH = EDEPCH + DESTEP 
      ELSE
         STRNE  = STRNE  + STEP
	 EDEPNE = EDEPNE + DESTEP
      ENDIF	 	                                                                                                    
*
* *** energy of particles contributing to edep
      IF ((DESTEP).GT.0.)THEN
         elog = log10((gekin+destep)/pkine(3))
	 call hfill (21,elog,0.,destep)
      ENDIF 	                                                                                 
      END                                                                       
