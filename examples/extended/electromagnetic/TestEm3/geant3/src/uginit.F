
      SUBROUTINE UGINIT
*
*     To initialise GEANT/USER  program and read data cards
*
#include "calor.inc"
#include "celoss.inc"
#include "geant321/gckine.inc"
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
* *** Calor definition
      CALL FFKEY('CALOR',NbAbsor,4,'MIXED')
      CALL FFKEY('MATE' ,materAbs(1),MaxAbs,'INTEGER')
      CALL FFKEY('THICK',thickAbs(1),MaxAbs,'REAL')
*
* *** read data cards
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
      CALL UHINIT
*
* *** some initialisation
      do k=1,MaxAbs
        sumEdep(k) = 0.
        sumTrck(k) = 0.
        su2Edep(k) = 0.
        su2Trck(k) = 0.
      enddo                                                         
*                                                                               
      END                                                                       
