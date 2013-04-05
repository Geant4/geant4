
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
* *** default values for histo
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
*     fileName for histograms, must be 1st data card !
      fileName = 'testem3.paw'
      READ(5,98)key,spaces,fileName
98    FORMAT(A4,A2,A25)

      CALL GFFGO
*
* *** check size of arrays
      if (NbAbsor.gt.MaxAbs) then
        write (6,51) NbAbsor, MaxAbs
	NbAbsor = MaxAbs
      endif
      if (NbLayer.gt.MaxLay) then
        write (6,52) NbLayer, MaxLay
	NbLayer = MaxLay
      endif      	
51    FORMAT (/,5x,'warning (uginit): NbAbsor= ',I2,' truncated to ',I2)
52    FORMAT (/,5x,'warning (uginit): NbLayer= ',I3,' truncated to ',I3)

      write(6,99) fileName
99    FORMAT(/,15x,'histogram file --> Name: ',A25)

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
      do k=1,MaxPlanes
        EnerFlow(k) = 0.
        EleakLat(k) = 0.
      enddo                                                         
*                                                                               
      END                                                                       
