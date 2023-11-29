      subroutine oauxfw(cFile, lenFil, IONUMB, cForm, lenFor, IERR)

      character(*), intent(in) :: cFile
      integer, intent(in)      :: lenFil
      character(lenFil)        :: FILE

      character(*), intent(in) :: cForm
      integer, intent(in)      :: lenFor
      character(lenFor)        :: CHSTTS

      FILE = cFile
      CHSTTS = cForm


      call oauxfi(FILE, IONUMB, CHSTTS, IERR)

      return
      end subroutine oauxfw
