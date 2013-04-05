
      SUBROUTINE UHINIT                                                         
*                                                                               
*     To book the user's histograms                                             
*
      CALL HBOOK1(1,'track length (mm) of primary particle',
     +              100,0.,500.,0.)
      CALL HBOOK1(2,'Nb of steps of primary particle',
     +              100,0.,100.,0.)
      CALL HBOOK1(3,'step size (mm) of primary particle',
     +              100,0.,10. ,0.)                                                                                                                                                                                           
*                                                                               
      END                                                                       
