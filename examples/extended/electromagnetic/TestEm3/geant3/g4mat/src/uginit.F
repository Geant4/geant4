
      SUBROUTINE UGINIT                                                         
*                                                                               
*     To initialise GEANT/USER  program and read data cards                     
*
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
#include "ugmate.inc"
#include "ggvalues.inc"
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
* *** define data cards
      call FFKEY ('MATE' ,Imat,1,'INTEGER')
      call FFKEY ('G4VAL',g4value(1),200,'REAL')  
*                                                                    
* *** read data cards                                                           
      CALL GFFGO 
*
* *** overwrite some initialization
      iloss  = 4
      ppcutm = cutele
*                                                                                      
      CALL GZINIT                                                               
      CALL GPART
*                                                                               
* *** materials definition
      CALL UGIMAT                                                               
*                                                                               
* *** Energy loss and cross-sections initialisations
      CALL GPHYSI                                                               
*                                                                               
      CALL GPRINT('MATE',0)
*
* *** ready for geant3/geant4 comparison
      call UG3G4
*                                                                               
      END                                                                       
