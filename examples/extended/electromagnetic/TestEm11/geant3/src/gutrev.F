
      SUBROUTINE GUTREV                                                         
*                                                                               
*     User routine to control tracking of one event                 
*     Called by GRUN                                                            
*                                                                     
#include "celoss.inc"
#include "histo.inc"                                                                       
*
* *** par event initialisations                                                                          
      Eevent  = 0.
*
* *** process the event                                                                               
      CALL GTREVE
*
* *** returns if no energy deposition
      if (Eevent.le.0.) return      
*
* *** histo: total energy deposited
      ih = 2
      if (histo(ih)) call hfill (ih,Eevent/histUnit(ih),0.,1.)
*
* *** cumul energy deposited
      Etot  = Etot  + Eevent
      Etot2 = Etot2 + Eevent*Eevent
*                                                                                           
      END 
