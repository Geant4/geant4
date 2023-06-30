      subroutine flabrw(cOut, lenOut, cMes, lenMes)

      character(*), intent(in) :: cOut
      integer, intent(in)      :: lenOut
      character(lenOut)        :: CHROUT

      character(*), intent(in) :: cMes
      integer, intent(in)      :: lenMes
      character(lenMes)        :: CHRMESS

      CHROUT = cOut
      CHRMESS = cMes


      call flabrt(CHROUT, CHMESS)

      return
      end subroutine flabrw

