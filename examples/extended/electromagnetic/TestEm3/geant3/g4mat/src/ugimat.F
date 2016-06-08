
      SUBROUTINE UGIMAT                                                          
*                                                                               
* *** Define user materials                                               
*
#include "geant321/gcbank.inc"                                                                                                                                                                                                                       
*
      DIMENSION Aair(2),Zair(2),Wair(2)
      DIMENSION ACO2(2),ZCO2(2),WCO2(2)
      DIMENSION AH2O(2),ZH2O(2),WH2O(2)
      DIMENSION AG10(4),ZG10(4),WG10(4)
      DIMENSION Asci(2),Zsci(2),Wsci(2)                                              
*      
* *** Air compound parameters          
      DATA Aair/14.01, 16.00/
      DATA Zair/ 7.  ,  8.  /
      DATA Wair/ 0.7 ,  0.3 /
*                                                                               
* *** CO2 compound parameters
      DATA ACO2/12.01, 16.00/
      DATA ZCO2/ 6.  ,  8.  /
      DATA WCO2/ 1.  ,  2.  /
*                                                                               
* *** Water compound parameters
      DATA AH2O/ 1.01, 16.00/
      DATA ZH2O/ 1.  ,  8.  /
      DATA WH2O/ 2.  ,  1.  /
*                                                                               
* *** G10 compound parameters
      DATA AG10/ 1.01, 12.00, 16.00, 28.00/
      DATA ZG10/ 1.  ,  6.  ,  8.  , 14.  /
      DATA WG10/ 3.  ,  3.  ,  2.  ,  1.  /      
*                                                                               
* *** Scintillator compound parameters
      DATA Asci/12.01,  1.01/
      DATA Zsci/ 6.  ,  1.  /
      DATA Wsci/ 9.  , 10.  /                                                                               

*
* *** Defines  materials
      CALL GSMIXT( 1,'Air'      , Aair ,Zair, 1.29E-3, 2   , Wair)
      CALL GSMIXT( 2,'CO2 gas'  , ACO2 ,ZCO2, 27.0E-3,-2   , WCO2)
      CALL GSMATE( 3,'H2 Liquid',  1.008,  1., 0.0708 , 865., 790., 0,0)
      CALL GSMIXT( 4,'Water'    , AH2O ,ZH2O, 1.0    ,-2   , WH2O)
      CALL GSMATE( 5,'Liquid Ar', 39.95, 18., 1.39   , 14.0, 84.0, 0,0)
      CALL GSMATE( 6,'Aluminium', 26.98, 13., 2.7    ,  8.9, 37.2, 0,0)
      CALL GSMATE( 7,'Iron     ', 55.85, 26., 7.87   , 1.76, 17.1, 0,0)
      CALL GSMATE( 8,'Lead     ',207.19, 82., 11.35  , 0.56, 18.5, 0,0)
      CALL GSMATE( 9,'Uranium  ',238.03, 92., 18.95  , 0.32, 12. , 0,0)
      CALL GSMATE(10,'Silicon  ', 28.09, 14.,  2.33  , 9.36, 45.5, 0,0)
      CALL GSMATE(11,'Tungsten ',183.85, 74., 19.30  , 0.35,  9.6, 0,0)
      CALL GSMIXT(12,'NemaG10'  , AG10 ,ZG10, 1.7    ,-4   , WG10)
      CALL GSMATE(13,'Copper   ', 63.55, 29., 8.96   , 1.43, 15.0, 0,0)
      CALL GSMIXT(14,'Scintilla', Asci ,Zsci, 1.032  ,-2   , Wsci)
      CALL GSMATE(15,'Gold     ',196.97, 79., 19.32  , 0.33,  9.6, 0,0)
*
*     total number of materals
      mtot = 15            
* 
* *** overwrite the computed radlength of some mixture
      JMA = LQ(JMATE-14)
      Q(JMA+9) = 42.549            
*                                               
*                                                                               
* *** Defines USER tracking media parameters
      IFIELD = 0                                                             
      FIELDM = 0.                                                           
      TMAXFD = 10.0                                                             
      STEMAX = 1000.                                                            
      DEEMAX = 0.20                                                            
      EPSIL  = 0.001                                                           
      STMIN  = 0.010                                                           
*
      do im =1,mtot                                                                               
        CALL GSTMED(im,'tracking medium',im, 0 ,IFIELD,FIELDM,TMAXFD,                
     *                 STEMAX,DEEMAX,EPSIL,STMIN, 0 , 0 )
      enddo        
*
      END                                                                       
