 
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*                                                                               
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "pvolum.inc"
#include "celoss.inc"                                                               
*                                                                               
      DIMENSION XSEL1(NBIN),XSEL1C(NBIN),XSER1(NBIN),XSER1C(NBIN),          
     +          XSEL2(NBIN),XSEL2C(NBIN),XSER2(NBIN),XSER2C(NBIN)          

*                                          
*                                                                                
      CALL GLAST
*
* *** close HIGZ      
      CALL HPLEND                                                                      
*                                                                               
* *** Normalize and print energy distribution                                   
      XEVENT=IEVENT                                                             
      CNORM  = 100./(XEVENT*PKINE(3))                                           
*
* *** longitudinal profile                                                                               
      DO 2 I = 1,NLTOT                                                          
         XSEL1 (I) = CNORM * SEL1 (I)                                           
         XSEL2 (I) = CNORM*SQRT(ABS(XEVENT*SEL2 (I) - SEL1 (I)**2))             
         XSEL1C(I) = CNORM * SEL1C(I)                                           
         XSEL2C(I) = CNORM*SQRT(ABS(XEVENT*SEL2C(I) - SEL1C(I)**2))
    2 CONTINUE
         CALL HPAK (6,XSEL2C)                                                                  
*
* *** radial profile                                                                                
      DO 3 I = 1,NRTOT                                                          
         XSER1 (I) = CNORM * SER1 (I)                                           
         XSER2 (I) = CNORM*SQRT(ABS(XEVENT*SER2 (I) - SER1 (I)**2))             
         XSER1C(I) = CNORM * SER1C(I)                                           
         XSER2C(I) = CNORM*SQRT(ABS(XEVENT*SER2C(I) - SER1C(I)**2))             
    3 CONTINUE
         CALL HPAK (12,XSER2C)
*
* *** total track length
      CNORM  = 1./(XEVENT*X0)
      XTRCH1 = CNORM*STRCH1
      XTRCH2 = CNORM*SQRT(ABS(XEVENT*STRCH2 - STRCH1**2))
      XTRNE1 = CNORM*STRNE1
      XTRNE2 = CNORM*SQRT(ABS(XEVENT*STRNE2 - STRNE1**2))      	                                       
*
* *** Print profiles
*
      PRINT 749
      PRINT 750
      PRINT 751
      DO 15 I=1,NLTOT
         B0 = (I-1)*DLX0
         B1 = I*DLX0                                                                             
         PRINT 754,B0,B1,XSEL1(I),XSEL2(I),B1,XSEL1C(I),XSEL2C(I)                                         
   15 CONTINUE                                                

      PRINT 760
      PRINT 751
      DO 16 I=1,NRTOT
         B0 = (I-1)*DRX0
         B1 = I*DRX0                                                                             
         PRINT 754,B0,B1,XSER1(I),XSER2(I),B1,XSER1C(I),XSER2C(I)                                         
   16 CONTINUE
*
* *** print summary
      PRINT 770
      PRINT 771,XSEL1C(NLTOT),XSEL2C(NLTOT)
      PRINT 772,XTRCH1,XTRCH2
      PRINT 773,XTRNE1,XTRNE2   
*                                                                               
* *** Save selected histograms                                                  
      IF (ISWIT(2).EQ.1) THEN                                                   
         CALL HRPUT(0,'testem2.paw','N')                                                                        
      ENDIF
*
  749 FORMAT(///)                                                         
  750 FORMAT(15X,'LATERAL PROFILE',35X,'CUMULATIVE LATERAL PROFILE'/)                                               
  751 FORMAT( 8X,'Bin',12X,'  Mean ',5X,'  rms',
     *       19X,'Bin', 9X,'  Mean ',5X,'  rms',/)
  754 FORMAT( 3X,F5.2,'->',F5.2,' radl: ',F7.2,'%    ',F7.2,'%',                              
     *       13X,    '0->',F5.2,' radl: ',F7.2,'%    ',F7.2,'%')                                                 
  760 FORMAT(///,15X,'RADIAL PROFILE',35X,'CUMULATIVE RADIAL PROFILE'/)
  770 FORMAT(///,30X,'SUMMARY',/)
  771 FORMAT(    25X,'energy deposit : ',F7.2,' % E0 +- ',F7.2,' % E0')
  772 FORMAT(    25X,'charged traklen: ',F7.2,' radl +- ',F7.2,' radl')
  773 FORMAT(    25X,'neutral traklen: ',F7.2,' radl +- ',F7.2,' radl')                                                                           
*
      END

