
      SUBROUTINE UGINIT
*
#include "geomate.inc"
#include "runinfo.inc"
#include "histo.inc"

      CHARACTER*25 filnam  
      CHARACTER*4  key
      CHARACTER*2  spaces     
C
C     INITIALIZE GEANT STRUCTURES
C
      CALL GINIT
*
*     default values (materials,sizes,tracking parameters)
*
      matabs = 23
      Xabsor = 1.
      YZabso = 2.
      FieldA = 0.
*
      matwld = 30      
      Xworld = 0.
      YZworl = 0.
      FieldW = 0.
*      
      xposab = 0.
*      
      xvert = -1000.
      yvert = 0.
      zvert = 0.
*
      deem = 0.2
      stmi = 0.001
      stma = 1.e+10
      kill = .false.
*
*     histograms
      do ih = 1,MaxHist
        histo(ih) = .false.
      enddo	              
*
* *** Define data cards
*
*     geometry
      CALL FFKEY('WORLD',matwld,4,'mixed')
      CALL FFKEY('ABSOR',matabs,4,'mixed')
      CALL FFKEY('XPOS' ,xposab,1,'real')  
      CALL FFKEY('VERT' ,xvert ,3,'real')
*      
*     histograms
      CALL FFKEY('HISTO',idhist,5,'MIXED')        
*
*     tracking parameters
      CALL FFKEY('DEEM',deem, 1,'REAL')
      CALL FFKEY('STMA',stma, 1,'REAL')
      CALL FFKEY('STMI',stmi, 1,'REAL')
*      
      CALL FFKEY('KILL',kill, 1,'LOGICAL')
*
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'testem5.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*
*     filename should be 1st data card !
      fileName = 'testem5.paw'
      READ(5,98)key,spaces,fileName
98    FORMAT(A4,A2,A25)
*
      CALL GFFGO
*
      write(6,99) fileName
99    FORMAT(/,15x,'histogram file --> Name: ',A25)
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
*
*     global initialisation for the run
*
*     energy
      sedep  = 0.
      sedep2 = 0.
*
*     number of steps
      sstpch  = 0.
      sstpch2 = 0.
      sstpne  = 0.
      sstpne2 = 0.
*
*     track length
      strkch  = 0.
      strkch2 = 0.
      strkne  = 0.
      strkne2 = 0. 
*
*     number of secondaries
      ngamma = 0
      nelect = 0
      nposit = 0
*
*     transmission, reflection counters
      nbtransm(1) = 0
      nbtransm(2) = 0
      nbreflec(1) = 0
      nbreflec(2) = 0                        
*
      if (xvert.le.-1000.) xvert = -0.5*Xworld
*
      END
