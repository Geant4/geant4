
      SUBROUTINE GUKINE
C
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconst.inc"
#include "geant321/gckine.inc"
#include "edepo.inc"
C
      DIMENSION VERTEX(3),PLAB(3)
C
      VERTEX(1)=XVERT     
      VERTEX(2)=YVERT     
      VERTEX(3)=ZVERT
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
C
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
*                                                                   
* *** Kinematics debug                     
      IF (IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)            
C
*
* *** initialisation per event
      anbstep=1.e-6
      anch=0.
      anne=0.
      anele=0.
      anpos=0.
      EDEP=0.
      tra=0.
      ref=0.
*      
      END
