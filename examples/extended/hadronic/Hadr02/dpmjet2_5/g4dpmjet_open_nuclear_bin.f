      subroutine g4dpmjet_open_nuclear_bin (namelen, unit,
     + opened, filename)
      
      character*(*) filename
      integer       namelen
      integer       unit
      logical       opened
C
C ------------------------------------------------------------------------------
C
      opened = .TRUE.
C      write (6,'(A)') filename(1:namelen)
C      close (6)
      open (unit=unit, file=filename(1:namelen), status="OLD", 
     + form="UNFORMATTED", err=1010)
     
      return
      
 1010 opened = .FALSE.
      return
      
      end
      
