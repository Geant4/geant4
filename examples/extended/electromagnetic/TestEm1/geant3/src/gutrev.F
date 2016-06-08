
      SUBROUTINE GUTREV                                                         
*                                                                               
*     User routine to control tracking of one event                 
*     Called by GRUN                                                            
*                                                                     
#include "celoss.inc"                                                                       
*                                                                          
      ETOT = 0.                                                               
*                                                                               
      CALL GTREVE                                                               
*                                                                               
      END 
