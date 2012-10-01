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
      program testbeta

      parameter (nwpawc = 2000000,mtuple=20)
      common/PAWC/PAW(nwpawc)
      
      real  spec(100)
      real spec1(100)
      
      integer Ntotal, A, Z
      real Q
      character*80 filename
      character*1 mode
      do i = 1,100
         spec(i) = 0.
         spec1(i) = 0.
      enddo

      print *, " Parent isotope (A,Z):"
      read (*,*) A,Z
      print *, " beta + or - decay:"
      read(*,*) mode
      print *, " The end point energy (MeV):"
      read(*,*) Q
      print *, " Name of the Hbook file"
      read(*,*) filename
      print *, " The time scale"
      read(*,*) ti
c
      call hlimit(nwpawc)
      call hropen(1,'EXAMPLE',filename,'N',1024,istat)
      if (istat.ne.0) go to 99
      Q1 = Q*1000
      call hbook1(100,'beta-theory ',100,0.,Q1,0.)
      call hbook1(200,'beta-Geant4 ',100,0.,Q1,0.)
      call hbook1(300,'decay profile', 100, 0.,ti,0.) 
      sum1 = 0.
      sum2 = 0.
      ee = .511
      dx = q/100.
      e0 = q/ee + 1.
      ad = a
      zd = z+1
      if (mode .eq. '+') zd = -(z-1)
c
      do i=1,100
         X=i*dx
         G=X/ee+ 1.0
         F=F2(G,zd)
c         print*, G,F
         spec(i)=F*SQRT(G*G-1.0)*(E0-G)*(E0-G)*G
         sum1 = sum1 + spec(i)
      enddo
      do i = 1,100
         spec(i) = spec(i) /sum1
c         print *, spec(i)
      enddo
      print *, ' '
      print *, " Name of the data file"
      read(*,*) filename
      open(11, file=filename, status='old')
      do while (.true.)
         read(11,*, end=30) e, x, x
         k = int(e/q/10.)
         spec1(k) = spec1(k) +1.
         call hfill (300, x, 0.,1.)
      enddo
 30   sum1 = 0.
      do i = 1,100
         sum1 = sum1 +spec1(i)
      enddo          
      do i = 1,100
         spec1(i) = spec1(i) /sum1
      enddo
      
            
      call hpak(100,spec)
      call hpak(200,spec1)
      
      call hrout(0,ICYCLE,' ')
c     
      call hrend('EXAMPLE')       
      
*     
 99   continue
      
      
      end
      
      FUNCTION F2(X,z0)
      PI=3.14159
      P=SQRT(X*X-1.0) 
      U=Z0/137.0
      S=SQRT(1.0-U*U) - 1.
      Y = 2*PI*U*X/P
      A1 = U*U*X*X + P*P/4.
      A2 = abs(Y/(1-exp(-Y)))
      F2=A1**S * A2
      RETURN
      END









