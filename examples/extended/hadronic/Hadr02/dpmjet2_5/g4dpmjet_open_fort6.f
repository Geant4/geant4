      subroutine g4dpmjet_open_fort6 (namelen, opened, filename)
      
      character*(*) filename
      integer       namelen
      logical       opened
C
C ------------------------------------------------------------------------------
C
      opened = .TRUE.
      open (unit=6, file=filename(1:namelen), status="UNKNOWN",
     & form="FORMATTED", err=1010)
     
      return
      
 1010 opened = .FALSE.
      return
      
      end
      
