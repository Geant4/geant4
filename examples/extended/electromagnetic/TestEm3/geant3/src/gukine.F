
      SUBROUTINE GUKINE
*  
*    Generates Kinematics for primary track
*
*    Data card Kine : Itype  Ekine x0  y0  z0  ux  uy  uz
*
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "calor.inc"                                                          
*                                                                               
      DIMENSION VERTEX(3),PLAB(3)
*
* *** set/reset default kinematic
      IF ((IKINE.LT.0).or.(abs(PKINE(2)).GT.0.5*worldX)) THEN
          IKINE = ABS(IKINE)
          PKINE(2) = -0.5*worldX
          PKINE(3) = 0.
          PKINE(4) = 0.
          PKINE(5) = 1.
          PKINE(6) = 0.
          PKINE(7) = 0.
      ENDIF
*
      VERTEX(1) = PKINE(2)
      VERTEX(2) = PKINE(3)
      VERTEX(3) = PKINE(4)
*
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
*
*
      JPA     = LQ(JPART-IKINE)
      XMASS   =  Q(JPA+7)
      PMOM    = SQRT(PKINE(1)*(PKINE(1)+2*XMASS))
      PLAB(1) = PMOM*PKINE(5)
      PLAB(2) = PMOM*PKINE(6)
      PLAB(3) = PMOM*PKINE(7)
*
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
*
* *** Kinematics debug
      IF (IEVENT.EQ.1.OR.IDEBUG.NE.0) CALL GPRINT('KINE',0)
*                                                                               
      END                                                                       
