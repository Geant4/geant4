 
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
     
      DIMENSION EDIST(100),EDISTC(100)                                                    
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
* *** mean total energy deposit by charged and by neutral
      EDEPCH = CNORM*EDEPCH
      EDEPNE = CNORM*EDEPNE                                                 
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
* *** normalize histo of energy ditribution of contributing particles
*     and compute cumulative distribution
      SUM = HSUM (21)
      CALL HUNPAK(21,EDIST,'HIST',1)
         EDIST( 1) = EDIST(1)*100/SUM
	 EDISTC(1) = EDIST(1)
      DO 17 I=2,100
         EDIST (I) = EDIST (I)*100/SUM
	 EDISTC(I) = EDISTC(I-1) + EDIST(I)
   17 CONTINUE
      CALL HPAK (21,EDIST)
      CALL HPAK (22,EDISTC)	 
*
* *** print summary
      PRINT 770
      PRINT 771,XSEL1C(NLTOT),XSEL2C(NLTOT)
      PRINT 774,EDEPCH
      PRINT 775,EDEPNE
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
  774 FORMAT(    25X,'edep by charged: ',F7.2,' % E0')
  775 FORMAT(    25X,'edep by neutral: ',F7.2,' % E0')                                                                            
*
      END

