
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc" 
#include "geant321/gctrak.inc"
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
*      
      CHARACTER*20 parnam,matnam
      CHARACTER*4  unit1,unit2
*
* *** run summary
      PRINT 749
      call GFPART (Ikine,parnam,ityp,dum,duc,dut,udu,idu)      
      call GFMATE (Imat ,matnam,dua,duz,dud,dur,dui,udu,idu)
      CALL GEVKEV (pkine(1),Ekine,unit1)
      PRINT 750, Ievent,parnam,Ekine,unit1,sizeX,matnam
      PRINT 749
*
      fnorm    = 1./Ievent      
*
* *** total energy deposit
      Etot  = fnorm*Etot
      Etot2 = fnorm*Etot2
      rms   = sqrt(abs(Etot2-Etot**2))
*
      CALL GEVKEV (Etot,aveE,unit1)
      CALL GEVKEV (rms ,rmsE ,unit2)
      PRINT 752, aveE,unit1,rmsE,unit2
*
* *** track length of primary track
      TrLen  = fnorm*TrLen
      TrLen2 = fnorm*TrLen2
      rms    = sqrt(abs(TrLen2-TrLen**2))            
      PRINT 753, TrLen,rms
*
* *** compare with csda range from table
      call GFTMAT(Imat,ikine,'RANG',1,pkine(1),csdaRang,dum,ixst)
      PRINT 763, csdaRang      
*     
* *** projected range of primary track
      XTrLen  = fnorm*XTrLen
      XTrLen2 = fnorm*XTrLen2
      rms     = sqrt(abs(XTrLen2-XTrLen**2))            
      PRINT 754, XTrLen,rms
*     
* *** nb of steps and step size of primary track
      fNbStep  = fnorm*NbStep
      fNbStep2 = fnorm*NbStep2
      rms1     = sqrt(abs(fNbStep2-fNbStep**2))
*      
      StepLen  = fnorm*StepLen
      StepLen2 = fnorm*StepLen2
      rms2     = sqrt(abs(StepLen2-StepLen**2))                  
      PRINT 755, fNbStep,rms1,StepLen,rms2
*
* *** transmission coefficients
      absorbed  = fnorm*100*kcoef(1)
      transmit  = fnorm*100*kcoef(2)
      reflected = fnorm*100*kcoef(3)
      PRINT 756, absorbed, transmit, reflected
*                  
      PRINT 749
*                                             
* *** geant termination
      CALL GLAST
*
* *** close HIGZ
      CALL HPLEND
*                                                                               
* *** Save histo
      if (iswit(2).eq.1) then
        write(6,99) fileName
99      format(/,15x,'histogram file --> Name: ',A25)                                                  
        CALL HRPUT(0,fileName,'N')
      endif	     
*
* *** formats
  749 FORMAT(/, 60(1H-),/)            
  750 FORMAT(  1X,'The run is',I8,1X,A10,' of',F8.2,A4,' through ',
     &             E12.4,' cm of ',A10)
  752 FORMAT(  1X,'Total Energy deposited        = ',
     &             F8.3, A4,' +- ', F8.3, A4)
  753 FORMAT(/,1X,'Track length of primary track = ',
     &             E12.5,' cm +- ', E12.5,' cm')
  763 FORMAT(  1X,'csda Range from table         = ',
     &             E12.5,' cm')     
  754 FORMAT(/,1X,'Projected range               = ',
     &             E12.5,' cm +- ', E12.5,' cm')      
  755 FORMAT(/,1X,'Nb of steps of primary track  = ',
     &             F6.2,' +- ', F6.2,
     &         3X,'Step size = ', E12.5,' cm +- ', E12.5,' cm')
  756 FORMAT(/,1X,'absorbed = ',F6.2,' %',3X,'transmit = ', F6.2,' %',
     &         3X,'reflected = ', F6.2,' %')              
*              
      END
