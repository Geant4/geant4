
      SUBROUTINE UGEOM                                                          
*                                                                               
* *** Define user geometry set up                                               
*                                                                                                                                                                                                                       
*
#include "calor.inc"  
*                                                                               
      DIMENSION PAR(3)
      
      DIMENSION Aair(2),Zair(2),Wair(2)
      DIMENSION AH2O(2),ZH2O(2),WH2O(2)                                        
*                                                                               
* *** Air compound parameters                                                   
      DATA Aair/14.01, 16.00/                                            
      DATA Zair/ 7.  ,  8.  /                                                     
      DATA Wair/ 0.7 ,  0.3 /
*                                                                               
* *** Air compound parameters                                                   
      DATA AH2O/ 1.01, 16.00/                                            
      DATA ZH2O/ 1.  ,  8.  /                                                     
      DATA WH2O/ 2.  ,  1.  /                                                           
*
* *** Defines USER perticular materials
      CALL GSMIXT( 1,'Air'      , Aair ,Zair, 1.29E-3, 2   , Wair)                                                
      CALL GSMATE( 2,'H2 Liquid',  1.01,  1., 0.0708 , 865., 790., 0,0)
      CALL GSMIXT( 3,'Water'    , AH2O ,ZH2O, 1.0    ,-2   , WH2O)      
      CALL GSMATE( 4,'Liquid Ar', 39.95, 18., 1.39   , 14.0, 84.0, 0,0)             
      CALL GSMATE( 5,'Aluminium', 26.98, 13., 2.7    ,  8.9, 37.2, 0,0)
      CALL GSMATE( 6,'Iron'     , 55.85, 26., 7.87   , 1.76, 17.1, 0,0)
      CALL GSMATE( 7,'Tungsten' ,183.85, 74., 19.30  , 0.35, 18.5, 0,0)                    
      CALL GSMATE( 8,'Lead'     ,207.19, 82., 11.35  , 0.56, 18.5, 0,0)                    
      CALL GSMATE( 9,'Uranium'  ,238.03, 92., 18.95  , 0.32, 12. , 0,0)                    
*                                               
*                                                                               
* *** Defines USER tracking media parameters
      IFIELD = 0                                                             
      IF (Field.GT.0.) IFIELD = 3
      FIELDM = 10*Field                                                               
      TMAXFD = 10.0                                                             
      STEMAX = 1000.                                                            
      DEEMAX = 0.20                                                            
      EPSIL  = 0.001                                                           
      STMIN  = 0.010                                                           
*                                                                               
      CALL GSTMED( 1,'Container',Imat, 0 ,IFIELD,FIELDM,TMAXFD,                
     *               STEMAX,DEEMAX,EPSIL,STMIN, 0 , 0 )        
*
*
* *** Geometry                                                                               
      PAR(1) = BoxSize/2.                                                               
      PAR(2) = BoxSize/2.                                                               
      PAR(3) = BoxSize/2.                                                               
      CALL GSVOLU('aBox','BOX ',1,PAR,3,IVOL)                     
*                                                                               
* *** Close geometry banks. (obligatory system routine)                         
      CALL GGCLOS                                                               
*
*
* *** dessin
      CALL GSATT ('*','SEEN',1)
*
      DO IX = 1,3
        CALL GDOPEN (IX)
        SCALE =   18./BoxSize
        PAXIS =   0.
        SAXIS =   0.1*BoxSize
        CALL GDRAWC ('aBox',IX,0.,10.,9.3,SCALE,SCALE)
        CALL GDAXIS (PAXIS,PAXIS,PAXIS,SAXIS)
        CALL GDSCAL (10. , 0.3)
        CALL GDCLOS
      END DO
*
                                                                               
      END                                                                       
