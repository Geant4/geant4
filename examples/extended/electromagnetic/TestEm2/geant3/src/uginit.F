
      SUBROUTINE UGINIT                                                         
*                                                                               
*     To initialise GEANT/USER  program and read data cards                     
*                                                                                                                                                                                                                       
#include "pvolum.inc"
#include "celoss.inc"
*
      CHARACTER*20 filnam 
*                                                                                     
* *** Define the GEANT parameters                                                          
      CALL GINIT 
                                                                          
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'run01.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*                                              
* *** material definition                                                                                       
      CALL FFKEY('MATE',IMAT,1,'INTEGER') 
*
* *** volumes and bins definition                                               
      CALL FFKEY('BINS',NLTOT,4,'MIXED')
*                                                         
* *** read data cards                                                           
      CALL GFFGO 
                                                                     
* *** achieve initialization                                                    
      NLTOT = MIN(NLTOT,NBIN)                                                   
      NRTOT = MIN(NRTOT,NBIN)                                                   
      CALL VZERO(SEL1,13*NBIN+6)                                                  
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
      END                                                                       
