
      SUBROUTINE UGEOM
*
* *** Define user geometry set up
*
#include "calor.inc"
#include "geant321/gcbank.inc"
*
      DIMENSION PAR(3)

      DIMENSION Aair(2),Zair(2),Wair(2)
      DIMENSION ACO2(2),ZCO2(2),WCO2(2)
      DIMENSION AH2O(2),ZH2O(2),WH2O(2)
      DIMENSION AG10(4),ZG10(4),WG10(4)
      DIMENSION Asci(2),Zsci(2),Wsci(2)
      DIMENSION ACsI(2),ZCsI(2),WCsI(2)                                                    
*
      CHARACTER*4  volnam
      CHARACTER*20 matnam
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
* *** CsI compound parameters
      DATA ACsI/ 126.90, 132.90/
      DATA ZCsI/ 53.   , 55.   /
      DATA WCsI/ 1.    , 1.    /     
*
* *** Defines USER perticular materials
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
      CALL GSMIXT(16,'CsI      ', ACsI ,ZCsI, 4.534  ,-2   , WCsI)            
* 
* *** overwrite the computed radlength of some mixture
      JMA = LQ(JMATE-14)
      Q(JMA+9) = 42.549            
*                                               
*                                                                               
* *** Defines USER tracking media parameters
      IFIELD = 0                                                             
      IF (Field.GT.0.) IFIELD = 3
      FIELDM = 10*Field
      TMAXFD = 10.0                                                             
      STEMAX = 1000.
      IF (stepmax.gt.0.) STEMAX = stepmax
      DEEMAX = 0.20                                                            
      EPSIL  = 0.001                                                           
      STMIN  = 0.010                                                           
*
      do k=1,NbAbsor
        CALL GSTMED( k,'absorber',materAbs(k), 0 ,IFIELD,FIELDM,TMAXFD,
     *                 STEMAX,DEEMAX,EPSIL,STMIN, 0 , 0 )
      enddo
* 
* *** set specific bcute/dcute (if any)     
      do k=1,4*NbAbsor,4
         itm = prodcut(k) + 0.01
	if(itm.ge.1) then
	   call GSTPAR(itm,'BCUTE' ,prodcut(k+1))
           call GSTPAR(itm,'BCUTM' ,prodcut(k+1))	   
	   call GSTPAR(itm,'DCUTE' ,prodcut(k+2))
           call GSTPAR(itm,'DCUTM' ,prodcut(k+2))
	   call GSTPAR(itm,'PPCUTM',prodcut(k+3))
	endif   
      enddo
*
      nudef = NbAbsor+1
      CALL GSTMED( nudef,'default' , 1    , 0 ,IFIELD,FIELDM,TMAXFD,
     *                 STEMAX,DEEMAX,EPSIL,STMIN, 0 , 0 )
*
*
* *** calor dimensions
      thLayer = 0.
      do k=1,NbAbsor
        thLayer = thLayer + thickAbs(k)
      enddo
      calorX  = NbLayer*thLayer
      worldX  = 1.2*calorX
      worldYZ = 1.2*calorYZ
*
* *** world
      PAR(1) = worldX /2.
      PAR(2) = worldYZ/2.
      PAR(3) = worldYZ/2.
      CALL GSVOLU('worl','BOX ',nudef,PAR,3,IVOL)
*
* *** calorimeter
      PAR(1) = calorX /2.
      PAR(2) = calorYZ/2.
      PAR(3) = calorYZ/2.
      CALL GSVOLU('calo','BOX ',nudef,PAR,3,IVOL)
      CALL GSPOS ('calo',1,'worl',0.,0.,0.,0,'ONLY')
*
* *** layers
      CALL GSDVN ('layr','calo',NbLayer,1)
*
* *** absorbers
      volnam = 'abs'
      xfront = -0.5*thLayer
      do k=1,NbAbsor
        PAR(1) = thickAbs(k)/2.
        PAR(2) = calorYZ/2.
        PAR(3) = calorYZ/2.
	volnam(4:4) = char(ichar('0')+k)
        CALL GSVOLU(volnam,'BOX ',k,PAR,3,IVOL)
        xcenter = xfront + 0.5*thickAbs(k)
        CALL GSPOS (volnam,1,'layr',xcenter,0.,0.,0,'ONLY')
	xfront = xfront + thickAbs(k)
      enddo                                      
*                                                                               
* *** Close geometry banks. (mandatory system routine)
      CALL GGCLOS
*
* *** print geometry
      PRINT 749
      PRINT 751,NbLayer 
      do k=1,NbAbsor
        call GFMATE (materAbs(k),matnam,dua,duz,dud,dur,dui,udu,idu)      
        PRINT 752,matnam,thickAbs(k)
      enddo
      PRINT 749                                                             
*      
  749 FORMAT(/ ,60(1H-),/)            
  751 FORMAT(1X,'The calorimeter is ',I2,' layers of:')
  752 FORMAT(5X,A10,': ',F8.4,' cm')     
*
* *** dessin
      CALL GSATT ('*'   ,'SEEN',1)
      CALL GSATT ('layr','SEEN',0)
*
      DO IX = 1,3
        CALL GDOPEN (IX)
        SCALE =   18./max(worldX,worldYZ)
        PAXIS =   0.
        SAXIS =   1.
        CALL GDRAWC ('worl',IX,0.,10.,9.3,SCALE,SCALE)
CCC        CALL GDAXIS (PAXIS,PAXIS,PAXIS,SAXIS)
        CALL GDSCAL (10. , 0.3)
        CALL GDCLOS
      END DO
*
      END                                                                       
