
      SUBROUTINE GUKINE
*                                                                               
*    Generates Kinematics for primary track
*
*    Data card Kine : Itype  Ekine                                       
*                                                                               
#include "geant321/gconst.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"                                                                                                                           
*                                                                               
      DIMENSION VERTEX(3),PLAB(3)
      DIMENSION RNDM(2)
*                                                     
      DATA VERTEX/3*0./                                                                                                                  
*
*                                                                  
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
*
      CALL GRNDM(RNDM,2)
      cosTheta = 2*RNDM(1) - 1.
      sinTheta = sqrt(1. - cosTheta**2)
      Phi   = TWOPI*RNDM(2)                                
*                                                         
      PLAB(1)  = PKINE(1)*sinTheta*COS(Phi)
      PLAB(2)  = PKINE(1)*sinTheta*SIN(Phi)
      PLAB(3)  = PKINE(1)*cosTheta                                       
*                                                          
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)                           
*                                                                   
* *** Kinematics debug
      IF (IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)
*                                                                               
      END
