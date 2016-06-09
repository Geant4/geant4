
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
#include "detector.inc"
*
      DIMENSION vertex(3),Plab(3)
      dimension rndm(2)
* 
      DATA vertex/3*0./
      DATA Plab  /3*0./                                                                                                                  
*
      vertex(1) = -0.5*BoxSize
*      
*     random in YZ 
      beam = 0.4*BoxSize
      call GRNDM (rndm,2)
*
      vertex(2) = (2*rndm(1)-1.)*beam     
      vertex(3) = (2*rndm(2)-1.)*beam      
*
      CALL GSVERT(vertex,0,0,0,0,NVERT)
*
      JPA     = LQ(JPART-IKINE)
      XMASS   =  Q(JPA+7)
      Plab(1) = SQRT(PKINE(1)*(PKINE(1)+2*XMASS))
*                                                          
      CALL GSKINE(Plab,IKINE,NVERT,0,0,NT)
*      
* *** Kinematics debug
      IF (IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)
*                                                                               
      END                            
