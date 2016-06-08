
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*                                                                               
#include "geant321/gcflag.inc"                                                      
*                                                                                                                        
*                                                                                
      CALL GLAST
*
* *** print  selected histo                                                                               
*
* *** close HIGZ
      CALL HPLEND
*                                                                               
* *** Save selected histograms                                                  
      IF (ISWIT(2).EQ.1) THEN                                                   
         CALL HRPUT(0,'testem4.paw','N')                                                                        
      ENDIF                                                                     
*
      END
