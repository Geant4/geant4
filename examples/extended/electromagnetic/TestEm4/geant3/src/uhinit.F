
      SUBROUTINE UHINIT                                                         
*                                                                               
*     To book the user's histograms                                             
*                                                                                                                                                              
      CALL HBOOK1(1,'total energy deposit in C6F6 (MeV)'             
     *,100, 0., 10., 0.)                                                     

      CALL HIDOPT(0,'STAT')                               
*                                                                               
      END                                                                       
