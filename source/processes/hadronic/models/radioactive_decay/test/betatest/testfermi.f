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
      program testfermi

      parameter (nwpawc = 2000000,mtuple=20)
      common/PAWC/PAW(nwpawc)
      
      real  spec(100)
      real spec1(100)
      real spec3(100)
      real spec4(100)
      real e(3),pd(3)
      
      integer Ntotal, A, Z
      real Q
      character*80 filename
      character*1 mode
      do i = 1,100
         spec(i) = 0.
         spec1(i) = 0.
         spec4(i) = 0.
         spec3(i) = 0.
      enddo

      print *, " Parent isotope (A,Z):"
      read (*,*) A,Z
      print *, " beta + or - decay:"
      read(*,*) mode
      print *, " The end point energy (MeV):"
      read(*,*) Q
      print *, " Number of decays to be simulated:"
      read(*,*) Ntotal
      print *, " Name of the Hbook file"
      read(*,*) filename
c
      call hlimit(nwpawc)
      call hropen(1,'EXAMPLE',filename,'N',1024,istat)
      if (istat.ne.0) go to 99
      Q1 = Q*1000
      call hbook1(100,'beta0 ',100,0.,Q1,0.)
      call hbook1(200,'beta1 (No F)',100,0.,Q1,0.)
      call hbook1(300,'beta2 ',100,0.,Q1,0.)
      call hbook1(400,'beta3 (No F)',100,0.,Q1,0.)
      call hbook1(500,'fermi factor',100,0.,Q1,0.)
      

      sum1 = 0.
      sum2 = 0.
      ee = .511
      dx = q/100.
      e0 = q/ee + 1.
      ad = a
      zd = z+1
      if (mode .eq. '+') zd = -(z-1)
c
c now try to work out the normalisation factor for the fermi factor
      fnorm = 0.
      do i=1,100
         X=i*dx
         G=X/ee+ 1.0
         F=F2(G,zd) 
         X = X*1000.
         call hfill(500,X,0.,F)
         if (F.gt.fnorm) fnorm = F
      enddo
      print *, zd, ad
      print *, " Fermi normlization factor = ", fnorm
c ideal spectra with and without fermi correction      
      do i=1,100
         X=i*dx
         G=X/ee+ 1.0
         F=F2(G,zd)
         spec(i)=F*SQRT(G*G-1.0)*(E0-G)*(E0-G)*G
         sum1 = sum1 + spec(i)
         spec1(i) = SQRT(G*G-1.0)*(E0-G)*(E0-G)*G
         sum2 = sum2 + spec1(i)
      enddo
      do i = 1,100
         spec(i) = spec(i) /sum1
         spec1(i) = spec1(i) / sum2
      enddo
      print *, ' '
      print *, " Simulating ..."
      do i = 1,Ntotal
 20      rd1 = rndm()
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
         
         k = e(1)/Q*100 + 1
         spec4(k) = spec4(k) +1.
         G=e(1)/ee + 1.0
         fermi=F2(G,zd)/fnorm
         if (rndm() .le. fermi) then
            spec3(k) = spec3(k) +1.
         endif
         
      enddo


      sum1 = 0.
      sum2 = 0.
      do i = 1,100
         sum1 = sum1 +spec3(i)
         sum2 = sum2 +spec4(i)
      enddo          
      do i = 1,100
         spec3(i) = spec3(i) /sum1
         spec4(i) = spec4(i) /sum2
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
      
      FUNCTION F2(X,z0)
      PI=4.0*ATAN(1.0)
      P=SQRT(X*X-1.0) 
      U=Z0/137.0
      S=SQRT(1.0-U*U) - 1.
      Y = 2*PI*U*SQRT(1+P*P)/P
      A1 = U*U*X*X + P*P/4.
      A2 = abs(Y/(1-exp(-Y)))
      F2=A1**S * A2
      RETURN
      END









