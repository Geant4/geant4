
      SUBROUTINE GUOUT                                                          
*                                                                               
*     User routine called at the end of each event                              
*                                                                               
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "pvolum.inc"
#include "celoss.inc"     
*
* *** drawing
*
#ifndef batch                                                          
      IF (ISWIT(1).NE.0) THEN
         CALL GDHEAD (110110,'testem2',0.)
         CALL GDSHOW (1)
         CALL GDXYZ  (0)
      END IF
#endif      
*
*
* *** statistic
*                                                                              
      DLC = 0.                                                                  
      DRC = 0.                                                                  

*     longitudinal profile
*                                                                                
      DO 2 I = 1,NLTOT
         SEL1 (I) = SEL1 (I) + DEDL(I)                                          
         SEL2 (I) = SEL2 (I) + DEDL(I)**2                                       
         DLC      = DLC      + DEDL(I)                                          
         SEL1C(I) = SEL1C(I) + DLC                                              
         SEL2C(I) = SEL2C(I) + DLC**2
         
	 BIN = (FLOAT(I)-0.5)*DLX0
         CALL HFILL(4,BIN,100*DEDL(I)/(DLX0*PKINE(3)),1.)
         BIN = FLOAT(I)*DLX0
         CALL HFILL(5,BIN,100*DLC    /PKINE(3),1.)                                                                                                                 
   2  CONTINUE                                                                  

*     radial profile
*                                                                                     
      DO 3 I = 1,NRTOT                                                                     
         SER1 (I) = SER1 (I) + DEDR(I)                                          
         SER2 (I) = SER2 (I) + DEDR(I)**2                                       
         DRC      = DRC      + DEDR(I)                                          
         SER1C(I) = SER1C(I) + DRC                                              
         SER2C(I) = SER2C(I) + DRC**2
	 
         BIN = (FLOAT(I)-0.5)*DRX0
         CALL HFILL(10,BIN,100*DEDR(I)/(DRX0*PKINE(3)),1.)
         BIN = FLOAT(I)*DRX0
         CALL HFILL(11,BIN,100*DRC    /PKINE(3),1.)                                                  
   3  CONTINUE                                                                  

*     particle flux
*                                                                                
      DO 14 IPAT = 1,3                                                          
      DO 14 NPL  = 1,NLTOT
         BIN = FLOAT(NPL)*DLX0
         CALL HFILL(6+IPAT,BIN,FNPAT(NPL,IPAT),1.)                                                   
   14 CONTINUE                                                                  

*     energy deposited and track length
*                                                                                
      ESEEN = 100.*DLC/PKINE(3)                                                 
      CALL HFILL(1, ESEEN,0.,1.)                                                
      CALL HFILL(2,STRCH/X0,0.,1.)                                            
      CALL HFILL(3,STRNE/X0,0.,1.)
*
      STRCH1 = STRCH1 + STRCH
      STRCH2 = STRCH2 + STRCH**2
      STRNE1 = STRNE1 + STRNE
      STRNE2 = STRNE2 + STRNE**2              
*                                                                              
      END                                                                       
