
      SUBROUTINE UHINIT                                                         
*                                                                               
*     To book the user's histograms                                             
*                                                                               
#include "geant321/gckine.inc"
#include "pvolum.inc"                                                                     
*                                                                               
* *** Histograms for showers development                                                                                                                                                                                                                  
*                                                                                                                                         
      EDMIN=0.                                                              
      EDMAX=100.                                                             
      TRKMX=100.*PKINE(1)                                                                                                                         
*                                                                               
      CALL HBOOK1(1,'total energy deposition (in percent of E inc)'             
     *,100,EDMIN,EDMAX,0.0)                                                     
      CALL HBOOK1(2,'total charged tracklengh (in radl)'                       
     *,100, 0. , TRKMX, 0.0)
      CALL HBOOK1(3,'total neutral tracklengh (in radl)'                       
     *,100, 0. , 10*TRKMX, 0.0)                                                                                        
*                                                                               
      CALL HIDOPT(0,'STAT')                                                     
*
* *** Longitudinal profile
      ZMAX=NLTOT*DLX0
                                                                                           
      CALL HBPROF(4,'longit energy profile    (in percent of E inc)'          
     *, NLTOT, 0.,ZMAX, 0., 1000.,' ')
*      
      ZMIN = 0.5*DLX0
      ZMAX = ZMIN + NLTOT*DLX0
                                                                                                                          
      CALL HBPROF(5,'cumul longit energy dep. (in percent of E inc)'          
     *, NLTOT, ZMIN,ZMAX, 0., 100.,' ')
      CALL HBOOK1(6,'resolution: cumul L energy dep. (% of E inc)'          
     *, NLTOT, ZMIN,ZMAX, 0.0)
                                                                         
*                                                                               
* *** Particle flux                                                   
      CALL HBPROF(7,'nb of gamma per plane'                                    
     *, NLTOT, ZMIN,ZMAX, 0., 1.E+6,' ')                                                   
      CALL HBPROF(8,'nb of posit per plane'                                    
     *, NLTOT, ZMIN,ZMAX, 0., 1.E+6,' ')                                                   
      CALL HBPROF(9,'nb of elect per plane'                                    
     *, NLTOT, ZMIN,ZMAX, 0., 1.E+6,' ')
*
* *** Radial profile
      RMAX=NRTOT*DRX0
      
      CALL HBPROF(10,'radial energy profile    (in percent of E inc)'          
     *, NRTOT, 0.,RMAX, 0., 1000.,' ')
*
      RMIN = 0.5*DRX0
      RMAX = RMIN + NRTOT*DRX0
      
      CALL HBPROF(11,'cumul radial energy dep. (in percent of E inc)'          
     *, NRTOT, RMIN,RMAX, 0., 100.,' ')
      CALL HBOOK1(12,'resolution: cumul R energy dep. (% of E inc)'          
     *, NRTOT, RMIN,RMAX, 0.0)
     
*
* *** Origin of the deposited energy
      tminlog = -6.
      tmaxlog =  0.
      CALL HBOOK1(21,'edep (normalized 100 percent) versus log(ekin/e0)'          
     *, 100, tminlog, tmaxlog, 0.0)
      CALL HBOOK1(22,'cumul edep (100 percent) versus log(ekin/e0)'          
     *, 100, tminlog, tmaxlog, 0.0)
                                                                                                                                                     
*                                                                                                                                     
      END                                                                       
