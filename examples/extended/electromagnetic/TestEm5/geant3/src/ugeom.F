
      SUBROUTINE UGEOM
C
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "edepo.inc"
C
      DIMENSION AWAT(2),ZWAT(2),WWAT(2)
      DIMENSION AAIR(2),ZAIR(2),WAIR(2)
      DIMENSION PALL(3),PVOL(3)
C
C     COMPOUND/MIXTURE PARAMETERS
C
      DATA AWAT/1.01,16.00/
      DATA ZWAT/1.,8./
      DATA WWAT/2.,1./
      DATA AAIR/14.01,16.00/
      DATA ZAIR/7.,8./
      DATA WAIR/0.7,0.3/
C
C     DEFINE  MATERIALS
C
      CALL GSMATE(21,'BERYLLIUM',9.010,4.,1.848,35.30,0.,0,0)
      CALL GSMATE(22,'ALUMINIUM',26.98,13.,2.70,8.900,0.,0,0)
      CALL GSMATE(23,'SILICON',28.09,14.,2.33,9.36,45.49,0,0)
      CALL GSMATE(24,'LIQUID ARGON',39.95,18.,1.4,14.,83.71,0,0)
      CALL GSMATE(25,'IRON',55.85,26.,7.87,1.760,0.,0,0)
      CALL GSMATE(26,'COPPER',63.54,29.,8.96,1.430,0.,0,0)
      CALL GSMATE(27,'GOLD',196.967,79.,19.32,0.33,0.,0,0)
      CALL GSMATE(28,'LEAD',207.190,82.,11.35,0.560,0.,0,0)
      CALL GSMATE(31,'XenonGas',131.29,54.,5.858e-3,1447.8,0.,0,0)      
C
      CALL GSMIXT(29,'WATER',AWAT,ZWAT,1.00,-2,WWAT)
      CALL GSMIXT(30,'AIR',AAIR,ZAIR,1.205E-3,+2,WAIR)
C
C     DEFINE MEDIA
C
      IFIELD=0
      if (FieldW.ne.0.) IFIELD=3
      FIELDM=10*FieldW
*      
      TMAXFD=10.
      STEMAX=1.E+10
      DEEMAX=0.20
      EPSIL=0.0001
      STMIN=0.0010
C
      CALL GSTMED(99,'WORLD',MATWLD,0,IFIELD,
     *            FIELDM,TMAXFD,STEMAX,DEEMAX,EPSIL,STMIN,0,0)
C

      IFIELD=0
      if (FieldA.ne.0.) IFIELD=3
      FIELDM=10*FieldA
*      
      IF(IGAUTO.LE.0) THEN
        DEEMAX=DEEM
        STMIN =STMI
        STEMAX=STMA
      ENDIF
      CALL GSTMED(1,'ABSORBER',MATABS,0,IFIELD,
     *            FIELDM,TMAXFD,STEMAX,DEEMAX,EPSIL,STMIN,0,0)
C
C
C     DEFINE VOLUMES (WORLD+ABSORBER)
*
      if(XWORLD*YZWORL.le.0.) then
        XWORLD=1.5*XABSOR 
        YZWORL=1.2*YZABSO
      endif
*      	       
      PALL(1)=0.5*XWORLD   
      PALL(2)=0.5*YZWORL     
      PALL(3)=0.5*YZWORL
      CALL GSVOLU('worl','BOX ',99,PALL,3,IVOL)
C
      PVOL(1)=0.5*XABSOR
      PVOL(2)=0.5*YZABSO
      PVOL(3)=0.5*YZABSO 
      CALL GSVOLU('abso','BOX ', 1,PVOL,3,IVOL)
C
C     BUILD  GEOMETRY
C
      CALL GSPOS('abso',1,'worl',XPOSAB,0.,0.,0,'ONLY')
*                                                                               
* *** Close geometry banks. (mandatory system routine)                         
      CALL GGCLOS
*
*
* *** dessin
      CALL GSATT ('*','SEEN',1)
*                                                                     
      DO IX = 1,3
        CALL GDOPEN (IX)
        SCALE =   18./max(XWORLD,YZWORL)
        PAXIS =   0.
        SAXIS =   0.1*max(XWORLD,YZWORL)
        CALL GDRAWC ('worl',IX,0.,10.,9.3,SCALE,SCALE)
        CALL GDAXIS (PAXIS,PAXIS,PAXIS,SAXIS)
        CALL GDSCAL (10., 0.3)
        CALL GDCLOS
      END DO
*      
      END
