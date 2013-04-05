
      SUBROUTINE UGLAST
*
#include "geant321/gcflag.inc"
*
#include "runinfo.inc"
#include "histo.inc"
*
      CHARACTER*4 unitm,unitr
      DIMENSION transmit(2),reflect(2)
*
      CALL GLAST
*
      IF(ievent.le.0) return
*      
      PRINT 749
      PRINT 750,ievent
*             
      totevent = float(ievent)
*
*     energy deposit      
      fmean = sedep/totevent
      rms2 = sedep2/totevent-fmean**2
      rms = 0.
      if (rms2.gt.0.) rms = sqrt(rms2/totevent)
*
      call gevkev (fmean,Emean,unitm)
      call gevkev (rms  ,Erms ,unitr)      
      print 752, Emean,unitm,Erms,unitr           
*
*     track length (charged)      
      fmean = strkch/totevent
      rms2 = strkch2/totevent-fmean**2
      rms = 0.
      if (rms2.gt.0.) rms = sqrt(rms2/totevent)
      print 753, fmean,rms
*
*     track length (neutral)      
      fmean = strkne/totevent
      rms2 = strkne2/totevent-fmean**2
      rms = 0.
      if (rms2.gt.0.) rms = sqrt(rms2/totevent)
      print 754, fmean,rms
*
*     number of steps (charged)      
      fmean = sstpch/totevent
      rms2 = sstpch2/totevent-fmean**2
      rms = 0.
      if (rms2.gt.0.) rms = sqrt(rms2/totevent)
      print 755, fmean,rms
*
*     number of steps (neutral)      
      fmean = sstpne/totevent
      rms2 = sstpne2/totevent-fmean**2
      rms = 0.
      if (rms2.gt.0.) rms = sqrt(rms2/totevent)
      print 756, fmean,rms            
*
*     number of secondaries
      gamma = ngamma/totevent
      elect = nelect/totevent
      posit = nposit/totevent
      print 757, gamma,elect,posit
*
*     transmission & reflection coefficients
      transmit(1) = 100.*nbtransm(1)/totevent
      transmit(2) = 100.*nbtransm(2)/totevent
      reflect (1) = 100.*nbreflec(1)/totevent
      reflect (1) = 100.*nbreflec(2)/totevent                  
      print 761, transmit(2)
      print 762, transmit(1)
      print 763, reflect (2)
      print 764, reflect (1)
*
      print 749      
*
*
* *** histograms
*
      call HPLEND      
      call HRPUT(0,fileName,'N')
*      
* *** formats
*
  749 FORMAT(/ ,80(1H=),/)            
  750 FORMAT(' Number of events = ', I7)
*  
  752 FORMAT(/,' Total energy deposit in absorber per event = ',
     +       F8.3,A4,' +- ',F6.2,A4)
*     
  753 FORMAT(/,' Total track length (charged) in absorber per event = ',
     +       F8.4,' cm +- ',F8.4, ' cm')
  754 FORMAT(  ' Total track length (neutral) in absorber per event = ',
     +       F8.4,' cm +- ',F8.4, ' cm')
*     
  755 FORMAT(/,' Number of steps (charged) in absorber per event = ',
     +       F7.3,' +- ',F7.3)
  756 FORMAT(  ' Number of steps (neutral) in absorber per event = ',
     +       F7.3,' +- ',F7.3)
*     
  757 FORMAT(/,' Number of secondaries per event : ',
     +  'gammas = ',F8.3, '; electrons = ',F8.3,'; positrons =',F8.3)
*     
  761 FORMAT(/,' Nb of events with the primary particle transmitted = ',
     +       F6.2, ' %')
  762 FORMAT(  ' Nb of events with at least 1 particle transmitted ',
     +       '(same charge as primary) = ',F6.2, ' %')
  763 FORMAT(/,' Nb of events with the primary particle reflected = ',
     +       F6.2, ' %')
  764 FORMAT(  ' Nb of events with at least 1 particle reflected ',
     +       '(same charge as primary) = ',F6.2, ' %')                                                           
*
      END
