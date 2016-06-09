
      SUBROUTINE UGINIT
*
*     To initialise GEANT/USER  program and read data cards
*
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
#include "geant321/gckine.inc"
#include "geant321/gcbank.inc"
*
      CHARACTER*20 filnam
*
* *** Define the GEANT parameters
      CALL GINIT
*
* *** default values
      do ih = 1,MaxHist
        histo(ih) = .false.
      enddo	                                                       
*
* *** Calor definition
      CALL FFKEY('CALOR',NbAbsor,4,'MIXED')
      CALL FFKEY('MATE' ,materAbs(1),MaxAbs,'INTEGER')
      CALL FFKEY('THICK',thickAbs(1),MaxAbs,'REAL')
* *** production cuts (bcute, dcute and ppcutm) for each absorber       
      CALL FFKEY('CUTPR',prodcut(1) ,4*MaxAbs,'REAL')      
* *** histograms
      CALL FFKEY('HISTO',idhist,5,'MIXED')
* *** max allowed step size       
      CALL FFKEY('STEPMX',stepmax,1,'REAL')
*                                                                          
* *** read data cards
      PRINT *, 'G3 > gives the filename of the data cards to be read:'
      READ (*,'(A)') filnam
      IF (filnam.EQ.' ') filnam = 'run01.dat'
      OPEN (unit=5,file=filnam,status='unknown',form='formatted')
*
      CALL GFFGO 
*
      CALL GZINIT
      CALL GPART
      CALL GPIONS
*
* *** overwrite ITRTYP for ion C12
      JPA = LQ(JPART-67)
      Q(JPA+6) = 8.0       
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
* *** some initialisation
      do k=1,MaxAbs
        sumEdep(k) = 0.
        sumTrck(k) = 0.
        su2Edep(k) = 0.
        su2Trck(k) = 0.
      enddo                                                         
*                                                                               
      END                                                                       
