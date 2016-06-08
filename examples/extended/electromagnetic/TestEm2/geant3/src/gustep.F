
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
      SAVE NLOLD                                                                               
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
         NR  = NUMBER(NLEVEL-1)                                                  
         NL  = NUMBER(NLEVEL)
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
* ***  Particle's flux                                                          
       NL  = NUMBER(NLEVEL)                                                     
       IF(SLENG.LE.0.) NLOLD = NL                                               
       IF(NL.NE.NLOLD) THEN                                                     
         NPL = (NL + NLOLD)/2                                               
	 IF (IPART.LE.3) FNPAT(NPL,IPART) = FNPAT(NPL,IPART) + 1.                                   
         NLOLD = NL                                                             
       ENDIF                                                                    
*
* *** energy of particles contributing to edep
      IF ((DESTEP).GT.0.)THEN
         elog = log10((gekin+destep)/pkine(3))
	 call hfill (21,elog,0.,destep)
      ENDIF 	                                                                                 
      END                                                                       
