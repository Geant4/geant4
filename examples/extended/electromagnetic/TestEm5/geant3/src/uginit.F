
      SUBROUTINE UGINIT
C
#include "geant321/gcunit.inc"
#include "geant321/gclist.inc"
#include "geant321/gckine.inc"
#include "edepo.inc"

      CHARACTER*25 filnam
      CHARACTER*25 HF
      CHARACTER*4 KEY
      CHARACTER*2 SPACES 
      DATA HF/'testem5.paw'/
C
C     INITIALIZE GEANT STRUCTURES
C
      CALL GINIT
*
C     default values (materials,sizes)
*
      MATABS = 23
      XABSOR = 1.
      YZABSO = 2.
      FieldA = 0.
*
      MATWLD = 30      
      XWORLD = 0.
      YZWORL = 0.
      FieldW = 0.
*      
      XPOSAB = 0.
*      
      XVERT = -1000.
      YVERT = 0.
      ZVERT = 0.
*
      NAME=HF  
C
C     DEFINE DATA CARDS
C
      CALL FFKEY('WORLD',MATWLD,4,'mixed')
      CALL FFKEY('ABSOR',MATABS,4,'mixed')
      CALL FFKEY('XPOS' ,XPOSAB,1,'real')  
      CALL FFKEY('VERT' ,XVERT ,3,'real')  
*
      CALL FFKEY('HIS1',NHI1 ,3,'mixed')
      CALL FFKEY('HIS2',NHI2 ,3,'mixed')
      CALL FFKEY('HIS3',NHI3 ,3,'mixed')
      CALL FFKEY('HIS4',NHI4 ,3,'mixed')
      CALL FFKEY('HIS5',NHI5 ,3,'mixed')
      CALL FFKEY('HIS6',NHI6 ,3,'mixed')
      CALL FFKEY('HIS7',NHI7 ,3,'mixed')
      CALL FFKEY('HIS8',NHI8 ,3,'mixed')
      CALL FFKEY('HIS9',NHI9 ,3,'mixed')
      CALL FFKEY('HI10',NHI10,3,'mixed')
*
      CALL FFKEY('DEEM',DEEM, 1,'REAL')
      CALL FFKEY('STMA',STMA, 1,'REAL')
      CALL FFKEY('STMI',STMI, 1,'REAL')
*
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'run01.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*
*     filename should be 1st data card !
      READ(5,98)KEY,SPACES,NAME
98    FORMAT(A4,A2,A25)
*
      CALL GFFGO
*
        write(6,99) NAME
99      FORMAT(/,15x,'histogram file --> NAME:',A25)
C
C     INITIALIZE ZEBRA STRUCTURES
C
      CALL GZINIT
C
C     DEFINE STANDARD PARTICLES AND MATERIALS
C
      CALL GPART
      CALL GMATE
*      
      CALL GDINIT
C
C     DEFINE GEOMETRY
C
      CALL UGEOM
C
C     COMPUTE CROSS SECTIONS/ENERGY LOSSES
C
      CALL GPHYSI
C
C     PRINT MATERIALS,MEDIA,VOLUMES
C
      CALL GPRINT('MATE',0)
      CALL GPRINT('TMED',0)
      CALL GPRINT('VOLU',0)
C
C     DEFINE HISTOGRAMS
C
      CALL UHINIT
C
C     SOME MORE INITIALISATION
C
      SEDEP=0.
      SEDEP2=0.
      sanch=0.
      sanch2=0.
      sanne=0.
      sanne2=0.
      sanele=0.
      sanpos=0.
      sstep=0.
      s2step=0.
      stra=0.
      sref=0.
      nth0=0
*
      X1ABSO = XPOSAB-0.5*XABSOR
      X2ABSO = XPOSAB+0.5*XABSOR
*
      if(xvert.le.-1000.) XVERT=-0.5*XWORLD
*
      END
