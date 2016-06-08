
      SUBROUTINE GUKINE                                                         
*                                                                               
*    Generates Kinematics for primary track
*
*    Data card Kine : Itype  Ekine                                       
*                                                                               
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"                                                                                                                          
*                                                                               
      DIMENSION VERTEX(3),PLAB(3)
*                                                     
      DATA VERTEX/3*0./
      DATA PLAB  /3*0./                                                                                                                  
*
*                                                                  
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
*
*
      JPA     = LQ(JPART-IKINE)
      XMASS   =  Q(JPA+7)
      PLAB(1) = SQRT(PKINE(1)*(PKINE(1)+2*XMASS))                                       
*                                                          
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)                           
*                                                                   
* *** Kinematics debug                     
      IF (IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)                                                                     
*                                                                               
      END                                                                       
