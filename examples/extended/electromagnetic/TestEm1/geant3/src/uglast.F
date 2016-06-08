
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*                                                                                                                              
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
#include "celoss.inc"        
*
* *** nb of steps and tracks per event
      TOTEVT = IEVENT
      FNTRK0 = NTRAK0/TOTEVT
      FNTRK1 = NTRAK1/TOTEVT
      FNSTP0 = NSTEP0/TOTEVT
      FNSTP1 = NSTEP1/TOTEVT
      PRINT 750, FNTRK0,FNTRK1,FNSTP0,FNSTP1
*
* *** frequency of processes call
      CALL UCTOH('MUNU',NAMEC(12),4,4)
      PRINT 760,(NAMEC(I),I=1,12)
      PRINT 761,(NBCALL(I)/TOTEVT,I=1,12)      
*                                                                                                                        
* *** geant termination                                                                               
      CALL GLAST
*
* *** print  selected histo                                                      

*
* *** close HIGZ file
      CALL HPLEND
*                                                                               
* *** Save selected histograms                                                  
      IF (ISWIT(2).EQ.1) THEN                                                   
         CALL HRPUT(0,'testem1.paw','N')                                                                        
      ENDIF                                                                     
*
* *** formats            
  750 FORMAT(//,1X,'nb tracks/event neutral: ',F9.2,'  charged: ',F9.2,
     +        /,1X,'nb steps/event  neutral: ',F9.2,'  charged: ',F9.2)
     
  760 FORMAT(//,1X,'nb of process calls per event: ',
     +        /,3X,12A9)
  761 FORMAT(   3X,12F9.2,//)
*
      END
