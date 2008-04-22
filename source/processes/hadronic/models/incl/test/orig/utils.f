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
