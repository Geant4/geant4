
      SUBROUTINE GUTREV                                                         
*                                                                               
*     User routine to control tracking of one event                 
*     Called by GRUN                                                            
*                                                                     
#include "celoss.inc"                                                                     
*                                                                          
      CALL VZERO(DEDL,5*NBIN)                                                                                                         
      STRCH = 0. 
      STRNE = 0.                                                              
*                                                                               
      CALL GTREVE                                                               
*                                                                               
      END                                                                       
