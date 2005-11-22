
      SUBROUTINE UGLAST                                                         
*                                                                               
*     Termination routine to print histograms and statistics                    
*
#include "geant321/gcflag.inc" 
#include "geant321/gctrak.inc"
#include "calor.inc"
#include "celoss.inc"
#include "histo.inc"
*
      DIMENSION EdepTot(MaxAbs)
*      
      CHARACTER*20 matnam
      CHARACTER*4 unit1, unit2
*      
      data aMeV /1.e-3/      
*
* *** mean energy deposit and track length
      PRINT 749
      PRINT 750,IEVENT
      PRINT 751 
      fnorm    = 1./IEVENT
      do k=1,NbAbsor
        call GFMATE (materAbs(k),matnam,dua,duz,dud,dur,dui,udu,idu)      
        aveEdep  = fnorm*sumEdep(k)
	aveEdep2 = fnorm*su2Edep(k)
        rmsEdep  = sqrt(abs(aveEdep2-aveEdep**2))
        aveTrck  = fnorm*sumTrck(k)
	aveTrck2 = fnorm*su2Trck(k)
        rmsTrck  = sqrt(abs(aveTrck2-aveTrck**2))
*     
        CALL GEVKEV (aveEdep,Edmean,unit1)
        CALL GEVKEV (rmsEdep,Edrms ,unit2)
        PRINT 752,k,matnam,Edmean,unit1,Edrms,unit2,aveTrck,rmsTrck
      enddo
      PRINT 749
*
* *** Energy flow
      Iplmax = NbAbsor*NbLayer + 1
      do k=1,Iplmax
        fk = k
	EnerFlow(k) = EnerFlow(k)*fnorm
	EleakLat(k) = EleakLat(k)*fnorm
        ih = 2*MaxAbs + 1
        if (histo(ih)) call hfill (ih, fk, 0., EnerFlow(k)/aMeV)
        ih = 2*MaxAbs + 2
        if (histo(ih)) call hfill (ih, fk, 0., EleakLat(k)/aMeV)	
      enddo 
*
* *** Energy deposit from energy flow balance
      do k=1,NbAbsor
        EdepTot(k)  = 0.
      enddo
*
      Idmax = NbAbsor*NbLayer
      do k=1,Idmax
        iAbsor = mod(k,NbAbsor)
        if (iAbsor.eq.0) iAbsor = NbAbsor
	Ebalance = EnerFlow(k) - EnerFlow(k+1) - EleakLat(k)
	EdepTot (iAbsor) = EdepTot (iAbsor) + Ebalance
      enddo
      
      PRINT 761
      do k=1,NbAbsor
        call GFMATE (materAbs(k),matnam,dua,duz,dud,dur,dui,udu,idu)
	CALL GEVKEV (abs(EdepTot (k)),Edtot,unit1)	      
	if (EdepTot (k).lt.0.) Edtot = -Edtot
        PRINT 762,k,matnam,Edtot,unit1
      enddo
      PRINT 749             
*                       
*                                             
* *** geant termination
      CALL GLAST
*
* *** close HIGZ
      CALL HPLEND
*                                                                               
* *** Save histo                                                  
      CALL HRPUT(0,fileName,'N')

*      
*
* *** formats
  749 FORMAT(/, 60(1H-),/)            
  750 FORMAT(1X,'Nb of events:',I6)
  751 FORMAT(1X,'AbsNo  Material',5X,'Energy deposit',14X,'Trck length')
  752 FORMAT(I5,3X,A10,F8.3,A4,' +- ',F6.2,A4,F9.3,' cm +- ',F6.2,' cm')
  761 FORMAT('  Energy deposition from energy flow balance:',/,/,
     &       1X,'AbsNo  Material',5X,'Total Edep')
  762 FORMAT(I5,3X,A10,F10.3,A4)     
*              
      END
