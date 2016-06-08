
      SUBROUTINE UGEOM                                                          
*                                                                               
* *** Define user geometry set up                                               
*                                                                                                                                                                                                                       
*                                                                               
      DIMENSION PAR(3)                                                         
      DIMENSION AC6F6(2),ZC6F6(2),WC6F6(2)
*                                                                               
* *** C6F6 compound parameters                                                   
      DATA AC6F6/12.02, 19.01/                                            
      DATA ZC6F6/ 6.  ,  9.  /                                                     
      DATA WC6F6/ 6.  ,  6.  /                                                                                                                
*                    
      CALL GSMIXT(1,'C6F6',AC6F6,ZC6F6,1.61,-2,WC6F6)
*                                               
*                                                                               
* *** Defines USER tracking media parameters
      IMAT   = 1                                    
      FIELDM = 0.0                                                             
      IFIELD = 0                                                               
      TMAXFD = 10.0                                                             
      STEMAX = 1000.                                                            
      DEEMAX = 0.20                                                            
      EPSIL  = 0.001                                                           
      STMIN  = 0.010                                                           
*                                                                               
      CALL GSTMED( 1,'C6F6',IMAT, 0 ,IFIELD,FIELDM,TMAXFD,                
     *               STEMAX,DEEMAX,EPSIL,STMIN, 0 , 0 )        
*
*
* *** Geometry                                                                               
      PAR(1) = 0.                                                               
      PAR(2) = 5.                                                               
      PAR(3) = 5.                                                               
      CALL GSVOLU('C6F6','TUBE',1,PAR,3,IVOL)                     
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
        SCALE =   1.75
        PAXIS =   0.
        SAXIS =   1.
        CALL GDRAWC ('C6F6',IX,0.,10.,9.5,SCALE,SCALE)
        CALL GDAXIS (PAXIS,PAXIS,PAXIS,SAXIS)
        CALL GDSCAL (10. , 0.3)
        CALL GDCLOS
      END DO
*
                                                                               
      END                                                                       
