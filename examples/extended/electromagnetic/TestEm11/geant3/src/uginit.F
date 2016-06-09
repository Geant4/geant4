
      SUBROUTINE UGINIT                                                         
*                                                                               
*     To initialise GEANT/USER  program and read data cards                     
*
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
*
      CHARACTER*20 filnam                                                                                                                                                                                                      
*                                                                               
* *** Define the GEANT parameters                                                          
      CALL GINIT 
*
* *** default values for histo
      do ih = 1,MaxHist
        histo(ih) = .false.
      enddo
*
* *** Calor definition                                               
      CALL FFKEY('CALOR', Imat,   4,'MIXED')
      CALL FFKEY('FIELD', Field,  1,'REAL')
      CALL FFKEY('STEPMX',stepmax,1,'REAL')
      CALL FFKEY('HISTO', idhist, 5,'MIXED')
*                                                                                      
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'run01.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*
*     fileName for histograms, must be 1st data card !
      fileName = 'testem11.paw'
      READ(5,98)key,spaces,fileName
98    FORMAT(A4,A2,A25)      
*
      CALL GFFGO
*                                                                               
      CALL GZINIT                                                               
      CALL GPART 
*
      CALL GDINIT                                                                     
*                                                                               
* *** Geometry and materials description                                        
      CALL UGEOM                                                                
*                                                                               
* *** Energy loss and cross-sections initialisations                            
      CALL GPHYSI                                                               
*                                                                               
      CALL GPRINT('MATE',0)                                                     
      CALL GPRINT('TMED',0)                                                     
      CALL GPRINT('VOLU',0)                                                     
*                                                                                                                                                              
* *** Define user histograms                                                    
      CALL UHINIT
*
* *** per run initialisations
      Etot     = 0.
      Etot2    = 0.
      TrLen    = 0.
      TrLen2   = 0.
      XTrLen   = 0.
      XTrLen2  = 0.
      StepLen  = 0.
      StepLen2 = 0.
      NbStep   = 0
      NbStep2  = 0
      kcoef(1) = 0
      kcoef(2) = 0
      kcoef(3) = 0
*                                                                               
      END
