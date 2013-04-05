
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
      IF ((ISWIT(1).EQ.1).AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF  (ISWIT(1).EQ.2)                     CALL GSXYZ                        
*                                                                               
* *** Something generated ?                                                     
      IF(NGKINE.GT.0) CALL GSKING(0)                                            
*                                                                               
* *** Energy deposited                                                          
      IF(DESTEP.GT.0.) ETOT = ETOT + DESTEP                                                     
*
      END                                                                       
