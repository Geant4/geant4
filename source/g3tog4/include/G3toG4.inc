*     This code implementation is the intellectual property of
*     the RD44 GEANT4 collaboration.
*    
*     By copying, distributing or modifying the Program (or any work
*     based on the Program) you indicate your acceptance of this statement,
*     and all its terms.
*    
*     $Id: G3toG4.inc,v 1.1 1999-01-07 16:06:44 gunter Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
*    
      integer maxlines
      parameter (maxlines = 750) ! max no. of cmd lines per C++ source file

      integer irot_index
      integer imat_index
      character context*4
      logical doclist      ! Build call list
      logical docode       ! Build C++ code
      logical dogeom       ! Build the Geant geometry as well
      integer lunlist      ! Output logical unit for call list
      integer luncode      ! Output logical unit for code
      integer nfile        ! Current source code file number
      integer nlines       ! Current source code line count
      common /g3tog4_c/ irot_index, imat_index,
     +  doclist, docode, dogeom, lunlist, luncode, context, nfile,
     +  nlines
