
      SUBROUTINE GUKINE                                                         
*                                                                               
*    Generates Kinematics for primary track
*
*    Data card Kine : Itype  Ekine  xbeam  rbeam                                      
*                                                                               
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "calor.inc"                                                                                                                          
*                                                                               
      DIMENSION VERTEX(3),PLAB(3)
      dimension rndm(3)
*                                                     
      DATA VERTEX/3*0./
      DATA PLAB  /3*0./                                                                                                                  
*
      IF (PKINE(2).NE.0.) vertex(1) = -0.5*SizeX
*      
*     random in YZ ?
      if ((pkine(3).lt.0.).or.(pkine(3).gt.1.)) pkine(3) = 0.
      rbeam = pkine(3)*0.5*SizeYZ
      call GRNDM (rndm,2)            
*
      VERTEX(2) = (2*rndm(1)-1.)*rbeam     
      VERTEX(3) = (2*rndm(2)-1.)*rbeam      
*                                                                  
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
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
