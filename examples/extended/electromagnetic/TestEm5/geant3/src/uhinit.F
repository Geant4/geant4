
      SUBROUTINE UHINIT
C
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "edepo.inc"
C
      IF(NHI1.GT.0) THEN
        CALL HBOOK1(1,'number of steps/event'
     +               ,NHI1,HIS1(1),HIS1(2),0.)
      ENDIF
      IF(NHI2.GT.0) THEN
        CALL HBOOK1(2,'energy deposit in absorber (MeV)'
     +  ,NHI2,HIS2(1),HIS2(2),0.)
      ENDIF
      IF(NHI3.GT.0) THEN
        CALL HBOOK1(3,'angle distribution at exit(deg)'
     +  ,NHI3,HIS3(1),HIS3(2),0.)
      ENDIF
      IF(NHI4.GT.0) THEN
        CALL HBOOK1(4,'lateral distribution at exit(mm)'
     +  ,NHI4,HIS4(1),HIS4(2),0.)
      ENDIF
      IF(NHI5.GT.0) THEN
        CALL HBOOK1(5,'kinetic energy of the primary at exit(MeV)'
     +  ,NHI5,HIS5(1),HIS5(2),0.)
      ENDIF
      IF(NHI6.GT.0) THEN
        CALL HBOOK1(6,'angle distr. of backscattered primaries(deg)'
     +  ,NHI6,HIS6(1),HIS6(2),0.)
      ENDIF
      IF(NHI7.GT.0) THEN
        CALL HBOOK1(7,'kinetic energy of backscattered primaries(deg)'
     +  ,NHI7,HIS7(1),HIS7(2),0.)
      ENDIF
      IF(NHI8.GT.0) THEN
        CALL HBOOK1(8,'kinetic energy of the charged secondaries(deg)'
     +  ,NHI8,HIS8(1),HIS8(2),0.)
      ENDIF
      IF(NHI9.GT.0) THEN
        CALL HBOOK1(9,'z of secondary charged vertices (mm)'
     +  ,NHI9,HIS9(1),HIS9(2),0.)
      ENDIF
      IF(NHI10.GT.0) THEN
        CALL HBOOK1(10,'energy of gammas escaping the absorber(MeV)'
*    +  ,NHI10,HIS10(1),HIS10(2),0.)
     +  ,NHI10,log10(HIS10(1)),log10(HIS10(2)),0.)
      ENDIF
*
      CALL HIDOPT(0,'STAT')
*
      END
