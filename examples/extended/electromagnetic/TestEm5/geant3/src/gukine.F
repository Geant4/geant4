
      SUBROUTINE GUKINE
*  
*    Generates Kinematics for primary track
*
*    Data card Kine : Itype  Ekine theta phi rbeam
*
*    rbeam : size of the random beam extension in fraction of 0.5*calorYZ      
*
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconst.inc"
#include "geant321/gckine.inc"
*
#include "geomate.inc"
#include "runinfo.inc"
*
      DIMENSION VERTEX(3),PLAB(3)
      dimension rndm(2)
*
*     random in YZ ?
      if ((pkine(4).lt.0.).or.(pkine(4).gt.1.)) pkine(4) = 0.
      rbeam = pkine(4)*0.5*YZABSO
      call GRNDM (rndm,2)            
*
      VERTEX(1)=XVERT     
      VERTEX(2)=YVERT + (2*rndm(1)-1.)*rbeam     
      VERTEX(3)=ZVERT + (2*rndm(2)-1.)*rbeam
*      
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)           
*
      JPA     = LQ(JPART-IKINE)
      XMASS   =  Q(JPA+7)
      PMOM    = SQRT(PKINE(1)*(PKINE(1)+2*XMASS))
      THETA   = PKINE(2)*DEGRAD
      PHI     = PKINE(3)*DEGRAD   
C
      PLAB(1)=PMOM*SIN(THETA)*COS(PHI)
      PLAB(2)=PMOM*SIN(THETA)*SIN(PHI)
      PLAB(3)=PMOM*COS(THETA)
*      
*     store the charge of primary particle in a free location
      pkine(10) = Q(JPA+8)         
C
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
*                                                                   
* *** Kinematics debug                     
      IF (IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)            
*
      END
