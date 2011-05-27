*    
*     ********************************************************************
*     * License and Disclaimer                                           *
*     *                                                                  *
*     * The  Geant4 software  is  copyright of the Copyright Holders  of *
*     * the Geant4 Collaboration.  It is provided  under  the terms  and *
*     * conditions of the Geant4 Software License,  included in the file *
*     * LICENSE and available at  http://cern.ch/geant4/license .  These *
*     * include a list of copyright holders.                             *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.  Please see the license in the file  LICENSE  and URL above *
*     * for the full disclaimer and the limitation of liability.         *
*     *                                                                  *
*     * This  code  implementation is the result of  the  scientific and *
*     * technical work of the GEANT4 collaboration.                      *
*     * By using,  copying,  modifying or  distributing the software (or *
*     * any work based  on the software)  you  agree  to acknowledge its *
*     * use  in  resulting  scientific  publications,  and indicate your *
*     * acceptance of all terms of the Geant4 Software license.          *
*     ********************************************************************
*    
      program test

      parameter (nwpawc = 2000000,mtuple=20)
      common/PAWC/PAW(nwpawc)

      real  spec(100)
      real spec1(100)
      real spec3(100)
      real spec4(100)
      real e(3),pd(3)

      call hlimit(nwpawc)
      call hropen(1,'EXAMPLE','testspec.hbook','N',1024,istat)
      if (istat.ne.0) go to 99
      call hbook1(100,'beta0',100,0.,1000.,0.)
      call hbook1(200,'beta1',100,0.,1000.,0.)
      call hbook1(300,'beta2',100,0.,1000.,0.)
      call hbook1(400,'beta3',100,0.,1000.,0.)


          Ap = 24
          Zp = 11
          Xp = -8.4175

          Ad = 24
          Zd = 12
          Xd = -13.9306

          Q = 1.0
          ee = 0.511
          Md = Zd * 938.2796 + (Ad-Zd) * 939.5731 - Xd
          Mp = Md + ee + Q;
c
          do i = 1,10000
 10          rd1 = rndm()
             rd2 = rndm()
             if (rd2 .gt. rd1) then
                rd = rd1
                rd1= rd2
                rd2 = rd
             endif
             
             pmax = 0.0
             psum = 0.0
             
             e(1) = rd2 * Q
             pd(1) = sqrt(e(1)*e(1) + 2.0*e(1)*ee)
             if (pd(1) .gt. pmax) pmax = pd(1)
             psum = psum +pd(1)
             
             e(2) = (rd1-rd2) * Q
             pd(2) = e(2)
             if (pd(2) .gt. pmax) pmax = pd(2)
             psum = psum + pd(2)
             
             e(3) = (1.-rd1) * Q
             pd(3) = sqrt(e(3)*e(3) + 2.0*e(3)*Md)
             if (pd(3) > pmax) pmax = pd(3)
             psum = psum + pd(3)
             if (pmax .gt. (psum-pmax)) goto 10
             k = e(1)*1000./10. + 1
             spec3(k) = spec3(k) +1.
          enddo
          sum1 = 0.
          do i = 1,100
             sum1 = sum1 +spec3(i)
          enddo          
          do i = 1,100
             spec3(i) = spec3(i) /sum1
          enddo


          do i = 1,10000
 20          rd1 = rndm()
             rd2 = rndm()
             
             pmax = 0.0
             psum = 0.0
             
             pd(1) = sqrt(rd2) * sqrt((Q+2.0*ee)*Q)
             e(1)  = sqrt(pd(1)*pd(1) + ee*ee) - ee;
             if (pd(1) .gt. pmax) pmax = pd(1)
             psum = psum +pd(1)

             e(2)  = sqrt(rd1) * Q
             pd(2)= e(2)            
             if (pd(2) .gt. pmax) pmax = pd(2)
             psum = psum + pd(2)
             
             e(3)  = Q - e(1)- e(2)
             if (e(3) .gt. 0.0)  then
                pd(3) = sqrt(e(3)*e(3) + 2.0*e(3)*Md)
                if (pd(3) > pmax) pmax = pd(3)
                psum = psum + pd(3)
             else 
                pmax = Q
                psum = Q
             endif
             
             if (pmax .gt. (psum-pmax)) goto 20
             k = e(1)*1000./10. + 1
             spec4(k) = spec4(k) +1.
          enddo
          sum1 = 0.
          do i = 1,100
             sum1 = sum1 +spec4(i)
          enddo          
          do i = 1,100
             spec4(i) = spec4(i) /sum1
          enddo


          
c
          sum1 = 0.
          sum2 = 0.
          dx = q*1000./100.
          e0 = q/ee + 1.
          do i=1,100
             X=i*dx
             G=X/511+ 1.0
             F=F2(G,zd,ad,e0)
             spec(i)=F*SQRT(G*G-1.0)*(E0-G)*(E0-G)*G
             sum1 = sum1 + spec(i)
             spec1(i) = SQRT(G*G-1.0)*(E0-G)*(E0-G)*G
             sum2 = sum2 + spec1(i)
          enddo
          do i = 1,100
             spec(i) = spec(i) /sum1
             spec1(i) = spec1(i) / sum2
          enddo

          call hpak(100,spec)
          call hpak(200,spec1)
          call hpak(300,spec3)
          call hpak(400,spec4)
          
          call hrout(0,ICYCLE,' ')
c     
          call hrend('EXAMPLE')       
      
*     
 99   continue


          end
      
      FUNCTION F2(G,z0,a0,e0)
      PI=4.0*ATAN(1.0)
      P=SQRT(G*G-1.0) 
      BETA=P/G
      U=Z0/137.0
      S=SQRT(1.0-U*U)
      R=1.48*EXP((1./3.)*ALOG(A0))
      RO=R/386.157
      ETA=U/BETA
      FI=ATAN(S/ETA)
      S2ETA2=S*S+ETA*ETA
      A1=4.0*PI*(1+S)
      A2=GAMMA(2.0*S)
      A2=A2*A2
      A3=EXP(2.0*(S-1.0)*ALOG(2.0*P*RO))
      A4=EXP((S-0.5)*ALOG(S2ETA2))
      A5=2.0*FI*ETA-2.0*S+S/6.0/S2ETA2
      F2=A1/A2*A3*A4*EXP(A5)
      RETURN
      END








