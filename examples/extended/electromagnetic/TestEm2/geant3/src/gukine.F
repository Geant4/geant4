
      SUBROUTINE GUKINE                                                         
*                                                                               
*    Generates Kinematics for primary track
*
*    Data card Kine : Itype  Ekine
*                     (pkine(3) is used internaly to store Etot)                                       
*                                                                               
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "pvolum.inc"                                                                     
*                                                                               
      DIMENSION VERTEX(3),PLAB(3)                                               
      DATA VERTEX/3*0./                                                         
      DATA PLAB  /3*0./                                                         
*                                                            
      VERTEX(3) = - Z1     + 0.0001                                 
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)                          
*                                                         
      JPA = LQ(JPART-IKINE)                              
      XMASS = Q(JPA+7)                                     
      PKINE(3) = XMASS + PKINE(1)                          
      PLAB(3)  = SQRT(PKINE(1)*(PKINE(3)+XMASS))                                 
*                                                          
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)                           
*                                                                   
* *** Kinematics debug                      
      IF(IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)    
*                                                                               
      END                                                                       
