
      SUBROUTINE UGEOM                                                          
*                                                                               
* *** Define user geometry set up                                               
*                                                                               
#include "geant321/gcbank.inc"
#include "pvolum.inc"
*                                                                                                                                       
      DIMENSION ZAir (2),AAir (2),WAir (2)
      DIMENSION ZH2O (2),AH2O (2),WH2O (2)                                                  
      DIMENSION ZBGO (3),ABGO (3),WBGO (3)                                         
      DIMENSION ZPbWO(3),APbWO(3),WPbWO(3)
      
      DIMENSION PAR(3)      
                                                  
*                                                                               
* *** Air mixture parameters                                             
      DATA ZAir/  7.00,  8.00 /                  
      DATA AAir/ 14.01, 16.00 /              
      DATA WAir/  0.70,  0.30 /
*                                                                               
* *** H2O compound parameters                                             
      DATA ZH2O/  1.00,  8.00 /                  
      DATA AH2O/  1.01, 16.00 /              
      DATA WH2O/  2.  ,    1. /                
*                                                                               
* *** BGO compound parameters 
      DATA ZBGO/  8.00, 32.00,  83.00 /                                                         
      DATA ABGO/ 16.00, 72.59, 208.98 /                                                                                          
      DATA WBGO/ 12.  ,  3.  , 4.     /                                                      
*                                                                               
* *** PbWO4 compound parameters
      DATA ZPbWO/  8.00,  74.00,  82.00 /                                       
      DATA APbWO/ 16.00, 183.84, 207.19 /                                                                                                 
      DATA WPbWO/  4.  ,   1.  ,   1.   /                                          
*                                                                               
* *** Defines USER perticular materials
*
      CALL GSMIXT( 1,'Air'      , AAir , ZAir , 1.29E-3, 2,WAir)
      CALL GSMIXT( 2,'Water'    , AH2O , ZH2O , 1.0    ,-2,WH2O)
      CALL GSMATE( 3,'Ar Liquid', 40.00, 18.  , 1.39   ,14.0 ,84.0,0,0)      
      CALL GSMATE( 4,'Aluminium', 26.98, 13.  , 2.7    , 8.9 ,37.2,0,0)      
      CALL GSMATE( 5,'Iron'     , 55.85, 26.  , 7.87   , 1.76,17.1,0,0)
      CALL GSMIXT( 6,'BGO'      , ABGO , ZBGO , 7.1    ,-3,WBGO)
      CALL GSMIXT( 7,'PbWO'     , APbWO, ZPbWO, 8.28   ,-3,WPbWO)                 
      CALL GSMATE( 8,'Lead '    ,207.19, 82.  ,11.35   ,0.56,18.5,0,0)
      CALL GSMATE( 9,'Tungsten' ,183.85, 74.  ,19.30   ,0.35, 9.6,0,0)      
*                                                                               
* *** Defines USER tracking media parameters                                    
      FIELDM = 0.0                                                             
      IFIELD = 0                                                               
      TMAXFD = 10.0                                                             
      STEMAX = 1000.                                                            
      DEEMAX = 0.20                                                            
      EPSIL  = 0.001                                                           
      STMIN  = 0.010                                                           
*                                                                               
      CALL GSTMED( 1,'Absorber',IMAT, 0 ,IFIELD,FIELDM,TMAXFD,                
     *               STEMAX,DEEMAX,EPSIL,STMIN, 0 , 0 )        
*                                                                               
* ***  Defines USER'S VOLUMES                                                   
      JMA = LQ(JMATE-IMAT)                                                      
      X0 = Q(JMA + 9)                                                           
      R1     =  NRTOT*DRX0*X0                                                   
      Z1     =  NLTOT*DLX0*X0*0.5                                               
*                                                                               
      PAR(1) = 0.                                                               
      PAR(2) = R1                                                               
      PAR(3) = Z1                                                               
      CALL GSVOLU( 'ECAL' , 'TUBE' ,    1, PAR , 3 , IVOL )                     
*                                                                               
      CALL GSDVN( 'RING' , 'ECAL' ,   NRTOT , 1)                                
      CALL GSDVN( 'SLAB' , 'RING' ,   NLTOT , 3)                                
*                                                                               
* *** Close geometry banks. (obligatory system routine)                         
      CALL GGCLOS                                                               
**
* *** dessin
      CALL GSATT ('*'   ,'SEEN',1)
      CALL GSATT ('RING','SEEN',0)
*
      DO IX =1,3
         CALL GDOPEN (IX)
         SCALE =   9.5/Z1
         PAXIS =   0.
         SAXIS =   0.2*Z1
         CALL GDRAWC ('ECAL',IX,0.,10.,10.,SCALE,SCALE)
CCC         CALL GDAXIS (PAXIS,PAXIS,PAXIS,SAXIS)
         CALL GDSCAL ( 10., 0.3)
         CALL GDCLOS
      END DO
*
                                                                               
      END                                                                       
