
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
* *** max allowed step size       
      CALL FFKEY('STEPMX',stepmax,1,'REAL')      
*                                                         
* *** read data cards                                                           
      CALL GFFGO 
*
* *** check size of arrays
      if (NLTOT.gt.NBIN) then
        write (6,51) NLTOT, NBIN
	NLTOT = NBIN
      endif
      if (NRTOT.gt.NBIN) then
        write (6,52) NRTOT, NBIN
	NRTOT = NBIN
      endif      	
51    FORMAT (/,5x,'warning (uginit): NLTOT= ',I3,' truncated to ',I3)
52    FORMAT (/,5x,'warning (uginit): NRTOT= ',I3,' truncated to ',I3)
*                                                                     
* *** achieve initialization                                                    
      CALL VZERO(SEL1,10*NBIN+8)                                                  
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
