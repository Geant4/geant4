      subroutine hbfinit()

C     Initialize the PAWC common block to a know size and tell HBook
C     about it. Used by the HBookFile class. I would do this in C++, but
C     I don't know how to create the correct style storage for the
C     common block.
C
C     Paul Rensing July 1994

      implicit none

      integer lqpaw, pawc
      PARAMETER (LQPAW = 1000000)
      COMMON /PAWC/ PAWC(LQPAW)

      CALL HLIMIT (LQPAW)
      return 
      end


      subroutine doclose(lun)
      
C     do a fortran close

      close(lun)
      return
      end
