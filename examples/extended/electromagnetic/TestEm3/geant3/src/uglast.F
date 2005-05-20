
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
      DIMENSION EdepForw(MaxAbs), EdepBack(MaxAbs), EdepTot(MaxAbs)
*      
      CHARACTER*20 matnam
      CHARACTER*4 unit1,unit2,unit3
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
        aveEdep = fnorm*sumEdep(k)
        rmsEdep = fnorm*sqrt(abs(IEVENT*su2Edep(k)-sumEdep(k)**2))
        aveTrck = fnorm*sumTrck(k)
        rmsTrck = fnorm*sqrt(abs(IEVENT*su2Trck(k)-sumTrck(k)**2))
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
	Eforw(k) = Eforw(k)*fnorm
	Eback(k) = Eback(k)*fnorm
        ih = 2*MaxAbs + 1
        if (histo(ih)) call hfill (ih, fk, 0., Eforw(k)/aMeV)
        ih = 2*MaxAbs + 2
        if (histo(ih)) call hfill (ih, fk, 0., Eback(k)/aMeV)	
      enddo 
*
* *** Energy deposit from energy flow balance
      do k=1,NbAbsor
        EdepForw(k) = 0.
        EdepBack(k) = 0.
        EdepTot(k)  = 0.
      enddo
*
      Idmax = NbAbsor*NbLayer
      do k=1,Idmax
        forward  = Eforw(k)   - Eforw(k+1)
	backward = Eback(k+1) - Eback(k)
        iAbsor = mod(k,NbAbsor)
        if (iAbsor.eq.0) iAbsor = NbAbsor
	EdepForw(iAbsor) = EdepForw(iAbsor) + forward
        EdepBack(iAbsor) = EdepBack(iAbsor) + backward
	EdepTot (iAbsor) = EdepTot (iAbsor) + forward + backward
      enddo
      
      PRINT 761
      do k=1,NbAbsor
        call GFMATE (materAbs(k),matnam,dua,duz,dud,dur,dui,udu,idu)
	CALL GEVKEV (abs(EdepTot (k)),Edtot,unit1)	      
        CALL GEVKEV (abs(EdepForw(k)),Edfor,unit2)
	CALL GEVKEV (abs(EdepBack(k)),Edbak,unit3)
	if (EdepTot (k).lt.0.) Edtot = -Edtot
	if (EdepForw(k).lt.0.) Edfor = -Edfor
	if (EdepBack(k).lt.0.) Edbak = -Edbak
        PRINT 762,k,matnam,Edtot,unit1,Edfor,unit2,Edbak,unit3
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
     &       1X,'AbsNo  Material',5X,
     &       'Total Edep',5X,'Forward Edep',2X,'Backward Edep')
  762 FORMAT(I5,3X,A10,F10.3,A4,1X,F10.3,A4,1X,F10.3,A4)     
*              
      END
