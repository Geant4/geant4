C     This is an adapted version of subroutine ranecu:
C     A. Padal, J. Sempau Computer Physics Cummunications 175 (2006) 440-450
      function ranecu(dummy)
      
      implicit double precision (a-h,o-z), integer*4 (i-n)
      common/rseed/iseed1,iseed2
      uscale=1.0D0/2.147483563D9

      i1=iseed1/53668
      iseed1=40014*(iseed1-i1*53668)-i1*12211
      if(iseed1.lt.0) iseed1=iseed1+2147483563

      i2=iseed2/52774
      iseed2=40692*(iseed2-i2*52774)-i2*3791
      if(iseed2.lt.0) iseed2=iseed2+2147483399

      iz=iseed1-iseed2
      if(iz.lt.1) iz=iz+2147483562
      ranecu=iz*uscale

      return
      end

C     For debugging
      subroutine dump_avatars()
      DIMENSION IND(20000),JND(20000)                                   P-N00350
      COMMON/BL2/CROIS(19900),K,IND,JND                                 P-N00540

      write(6,*) 'Avatars: (number of avatars = ',k,')'
      do i=1,k
         write(6,*) 'i = ', i
         write(6,*) 'crois(',i,') = ',crois(i)
         write(6,*) 'ind(',i,') = ',ind(i)
         write(6,*) 'jnd(',i,') = ',jnd(i)
      enddo

      return
      end

      subroutine dump_nucleons()
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N00530
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N00550
      integer dumps
      data pf,pf2,pf3 /270.33936,73083.4,19756627./
      save dumps

      IA=IA1+IA2                                                        P-N22500
      write(6,*) 'Nucleons: (number of nucleons = ',IA,')'
      do i=1,IA
         write(6,*) 'i = ', i
         write(6,*) 'x1(',i,') = ',x1(i)
         write(6,*) 'x2(',i,') = ',x2(i)
         write(6,*) 'x3(',i,') = ',x3(i)
      enddo
      do i=1,IA
         write(6,*) 'p1(',i,') = ',p1(i)
         write(6,*) 'p2(',i,') = ',p2(i)
         write(6,*) 'p3(',i,') = ',p3(i)
         write(6,*) 'eps(',i,') = ',eps(i)
         ptot2 = p1(i)**2 + p2(i)**2 + p3(i)**2
         if(ptot2.gt.pf2) then
            write(6,*)'Warning! ptot2/pf2 = ',ptot2/pf2
         endif
      enddo
      dumps = dumps + 1
      return
      end
