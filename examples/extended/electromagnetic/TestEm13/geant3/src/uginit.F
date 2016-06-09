
      SUBROUTINE UGINIT                                                         
*                                                                               
*     To initialise GEANT/USER  program and read data cards                     
*
#include "detector.inc"
#include "process.inc"
#include "histo.inc"                                                            
*
      CHARACTER*20 filnam
      CHARACTER*4  key
      CHARACTER*2  spaces                                                                                                                                                                                                                 
*                                                                               
* *** Define the GEANT parameters                                                          
      CALL GINIT 
*
*     histograms
      do ih = 1,MaxHist
        histo(ih) = .false.
      enddo
*
* *** Detector definition                                               
      CALL FFKEY('DETECTOR',Imat,2,'MIXED')
*      
*     histograms
      CALL FFKEY('HISTO',idhist,5,'MIXED')
*      	                                                                          
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'allprocesses.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*
*     filename should be 1st data card !
      fileName = 'photonprocesses.paw'
      READ(5,98)key,spaces,fileName
98    FORMAT(A4,A2,A25)                 
*                                                                    
* *** read data cards                                                           
      CALL GFFGO
*
      write(6,99) fileName
99    FORMAT(/,15x,'histogram file --> Name: ',A25)       
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
* *** initialisation of /process/

      CALL VZERO (nbCall,12)
      nbTot    = 0
      sumTrak  = 0.
      sumTrak2 = 0.                                                               
*                                                                               
      END                                                                       
