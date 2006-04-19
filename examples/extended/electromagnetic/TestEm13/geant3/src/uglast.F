
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*                                                                                                                              
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
*
#include "detector.inc"  
#include "process.inc"
#include "histo.inc"        
*
      character*20 material, particle
      character*4  unit
      real massMFP, massAC
*
* *** run conditions
      call gfmate(Imat,material,dum,dum,density,dum,dum,dum,ndum)
      call gfpart(Ikine,particle,ndum,dum,dum,dum,dum,ndum)
      call gevkev(pkine(1),energy,unit)
      PRINT 750, ievent,particle,energy,unit,boxSize,material,density
*      
* *** frequency of processes call
      CALL UCTOH('MUNU',NAMEC(12),4,4)
      PRINT 760,(NAMEC(I),I=1,12)
      PRINT 761,(NBCALL(I),I=1,12)
*
* *** compute mean free path and related quantities
      avTrak  = sumTrak /nbTot
      avTrak2 = sumTrak2/nbTot
      rms = sqrt(abs(avTrak2 - avTrak*avTrak))
      AtteCoef = 1./avTrak
*
      massMFP = avTrak*density
      massAC  = 1./massMFP
      print 770, avTrak,rms,massMFP,AtteCoef,massAC                              
*                                                                                                                        
* *** geant termination                                                                               
      CALL GLAST
*
* *** close HIGZ file
      CALL HPLEND
*                                                                               
* *** histograms                                                                                                  
      CALL HRPUT(0,fileName,'N')
*
* *** formats            
*
  750 FORMAT(/,1X,'The run consists of ',I7,1X,A8,' of ',F7.2,A4,
     +      ' throught',F10.2,' cm of ',A8,' (density: ',F7.3,' g/cm3)')
          
  760 FORMAT(/,1X,'Frequency of process calls: ',
     +       /,1X,12A8)
  761 FORMAT(  1X,12I8,/)
  
  770 FORMAT(/,1X,'MeanFreePath:',F12.5,' cm +- ',F12.5,' cm',
     +          5X,'massic:',F12.5,' g/cm2', 
     +        /,1X,'CrossSection:',F12.5,' cm^-1 ',15X,
     +          5X,'massic:',F12.5,' cm2/g',/)  
*
      END
