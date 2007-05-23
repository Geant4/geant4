C     Clean up common blocks

      subroutine commoncleanup()

      COMMON/SPL2/ X(100),Y(100),A(100),B(100),C(100),N    
      DIMENSION DTONC(13), DTOND(13)
      COMMON /DTON/DTONC,DTOND,FN

      do i=1,13
         DTONC(i) = 0.0
         DTOND(i) = 0.0
      enddo
      fn = 0.0

      DO i=1,100
         X(i) = 0.0
         Y(i) = 0.0
         A(i) = 0.0
         B(i) = 0.0
         C(i) = 0.0
      ENDDO
      N = 0

      return
      end
