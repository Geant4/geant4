
      SUBROUTINE UGLAST
C
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconst.inc"
#include "geant321/gckine.inc"
#include "edepo.inc"
C
      CALL GLAST
C
      IF(IEVENT.GT.0) THEN
        AEV=FLOAT(IEVENT)
        SEDEP=SEDEP/AEV
        SE2=SEDEP2/AEV-SEDEP**2
        SEDEP2=0.
        IF(SE2.GT.0.) SEDEP2=SQRT(SE2/AEV)
*
        WRITE(CHMAIL,1002)SEDEP,SEDEP2
        CALL GMAIL(0,0)
        sstep=sstep/aev
        s2step=s2step/aev-sstep*sstep
        sigst=0.
        if(s2step.gt.0.) sigst=sqrt(s2step/aev)
        write(6,122)sstep,sigst
122     format(/,'  mean number of steps=',f15.3,' +- ',f15.3) 
        sanch=sanch/aev
        sigch2=sanch2/aev-sanch*sanch
        sigch=0.
        if(sigch2.gt.0.) sigch=sqrt(sigch2/aev)
        sanne=sanne/aev
        signe2=sanne2/aev-sanne*sanne
        signe=0.
        if(signe2.gt.0.) signe=sqrt(signe2/aev)
        write(6,123) sanch,sigch,sanne,signe
123     format(/,'  mean nb of charged secondaries=',f15.4,' +- ',f15.4, 
     +         /,'  mean nb of neutral secondaries=',f15.4,' +- ',f15.4) 
        sanele=sanele/aev
        sanpos=sanpos/aev
        write(6,124)sanele,sanpos
124     format(/,'  mean nb of e-s = ',f15.4,'  and e+s=',f15.4)
        stra=stra/aev
        sref=sref/aev
        write(6,125)stra,sref
125     format(/,'  (number) transmission coeff=',e12.4,
     +           '  reflexion coeff=',e12.4)
        write(6,126) nth0
126     format(/,'  nb of (theta<0.001*bin) events=',i6)
*
CCC        CALL HISTDO
*
        CALL HPLEND
*
        CALL HRPUT(0,NAME,'N')
C
      ENDIF
C
1002  FORMAT(5X,'energy deposit=',E12.4,' +- ',E12.4)
C
      END
