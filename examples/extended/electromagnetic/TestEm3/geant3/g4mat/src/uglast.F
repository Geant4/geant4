
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*            
*                                                                                
* *** geant termination                                                                               
      CALL GLAST
*
* *** close HIGZ file
      CALL HPLEND 
*
* *** save histograms
      call HRPUT(0,'plmat.paw','N')                                                             
*
      END
