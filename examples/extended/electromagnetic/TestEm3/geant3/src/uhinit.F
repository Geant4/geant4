
      SUBROUTINE UHINIT
*
*     To book the user's histograms
*
#include "calor.inc"
*
      CHARACTER*5 chtags(20)
*
      DATA chtags /'EAbs0','LAbs0','EAbs1','LAbs1','EAbs2','LAbs2',
     +             'EAbs3','LAbs3','EAbs4','LAbs4','EAbs5','LAbs5',
     +             'EAbs6','LAbs6','EAbs7','LAbs7','EAbs8','LAbs8',
     +	           'EAbs9','LAbs9'/		   
*      
      CALL HBOOKN(1,'sum per event',2*NbAbsor,' ',4000,chtags)
*                                                                               
      END                                                                       
