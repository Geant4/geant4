
      SUBROUTINE GUSTEP
C

#include "geant321/gcbank.inc"
#include "geant321/gctmed.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gccuts.inc"
#include "geant321/gconst.inc"
#include "geant321/gcunit.inc"
#include "edepo.inc"
C
      parameter ( facdeg=180./3.1415927 )
      DATA IEVOLD/0/

*                                                                               
*                                                                               
* *** Debug event and store tracks for drawing
      IF (IDEBUG.NE.0)   CALL GPCXYZ      
      IF ((ISWIT(1).EQ.1).AND.(CHARGE.NE.0.)) CALL GSXYZ
      IF  (ISWIT(1).EQ.2)                     CALL GSXYZ                        
*
*   
      IF(IEVOLD.NE.IEVENT) THEN
        IEVOLD=IEVENT
        ITOLD=0
        NMOLD=0
        XO=VERT(1)
        RO=sqrt(vect(2)**2+vect(3)**2)
      ENDIF
C
C     TRACK FLAG
C
      ITFL=ISTAK+10000*ITRA+100000*IVERT
      IF(ITOLD.NE.IEVENT) THEN
        XO=VERT(1)
        RO=sqrt(vect(2)**2+vect(3)**2)
      ENDIF
C     
      X=VECT(1)
C
      IF(NGKINE.GT.0) THEN
        do 10 i=1,ngkine
        itypa=gkin(5,i)
*       all the charged secondaries
        if(itypa.ne.1)  call hfill(9,10.*GPOS(1,i),0.,1.)
*       in absorber only !
        if(NUMED.EQ.1) then
          if(itypa.eq.1) then
             anne=anne+1.
          else
             anch=anch+1.
             if(itypa.eq.2) anpos=anpos+1.
             if(itypa.eq.3) anele=anele+1.
             tsec=1000.*(GKIN(4,i)-510.999e-6)
             call hfill(8,tsec,0.,1.)
          endif 
        endif 
10      continue
      CALL GSKING(0)
      ENDIF
*
      R=SQRT(VECT(2)**2+VECT(3)**2)
      if((r.le.YZABSO).and.(xo.lt.X2ABSO).and.(x.ge.X2ABSO)) then
       if(vect(4).gt.0.) then
        if(ipart.eq.1) call hfill(10,log10(1000.*gekin),0.,1.)
        if(ipart.eq.ikine) then
          tra=1.
*
          call hfill(5,1000.*gekin,0.,1.)

          if((1-vect(4)).gt.1.e-6) then
            theta=acos(vect(4))
          else
            theta=sqrt(2.*(1.-vect(4)))
          endif
          if(nhi3.gt.0) then
            h1=his3(1)/facdeg
            h2=his3(2)/facdeg
            dth=(h2-h1)/float(nhi3)
            ibin=int((theta-h1)/dth)+1
            if((ibin.le.0).or.(ibin.gt.nhi3)) then
              w=1.
            else
              th=ibin*dth
              if(theta.GT.0.001*dth) then
                w=PI/(64800.*dth*sin(theta))
              else
                nth0=nth0+1
                if(nth0.le.10) then
                thdeg=theta*facdeg
                write(6,997) thdeg
997    format('theta<0.001*dth (from plot excluded)'
     +                 ,' theta=',e12.4,' deg') 
                endif
                w=0.
              endif
            endif
            theta = facdeg*theta
            call hfill(3,theta,0.,w)
          endif
          call hfill(4,10.*R,0.,1.)
        endif
       endif
      endif
      if((r.le.YZABSO).and.(xo.gt.X1ABSO).and.(x.le.X1ABSO)) then
       if(vect(4).lt.0.) then
        if(ipart.eq.ikine) then
          ref=1.
          call hfill(7,1000.*gekin,0.,1.)
          if((1.+vect(4)).gt.1.e-6) then
            theta=acos(vect(4))-PIBY2
          else
            theta=sqrt(2.*(1.+vect(4)))
          endif
          if(nhi6.gt.0) then
            h1=his6(1)/facdeg
            h2=his6(2)/facdeg
            dth=(h2-h1)/float(nhi6)
            ibin=int((theta-h1)/dth)+1
            if((ibin.le.0).or.(ibin.gt.nhi6)) then
              w=0.
            else
              th=ibin*dth
              if(theta.GT.0.001*dth) then
                w=PI/(64800.*dth*sin(theta))
              else
                thdeg=theta*facdeg
                write(6,997) thdeg
                w=0.
              endif
              theta = facdeg*theta
              call hfill(6,theta,0.,w)
            endif
          endif
        endif
       endif
      endif
         
        IF(NUMED.EQ.1) THEN
          anbstep=anbstep+1.
          EDEP = EDEP+DESTEP
        ENDIF
C
      ITOLD=IEVENT
      XO=X
      RO=R
C
      END
