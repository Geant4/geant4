*     This code implementation is the intellectual property of
*     the GEANT4 collaboration.
*    
*     By copying, distributing or modifying the Program (or any work
*     based on the Program) you indicate your acceptance of this statement,
*     and all its terms.
*    
*     $Id: rztog4.F,v 1.2 1999-12-05 17:50:17 gcosmo Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
*    
      program rztog4
************************************************************************
*
*       rztog4
*
*       A standalone program to convert a Geant RZ initialization file
*       (containing geometry, materials, etc.) to a Geant4 C++
*       geometry.
*
*       Torre Wenaus, LLNL  6/95
*
*       todo:
*       - attributes
*
*       JSHAPE routine of Jouko Vuoskoski, CERN employed. 
*
*...  History:
*     15-Jan-1997 Lockman : put tog4 in a separate file
************************************************************************
      implicit none
#include "G3toG4.inc"
      integer iargc, lenocc, ln, lnout
      character*256  file, ofile
*
      real h, q
      integer nh, nz
      parameter (nh=1000000, nz=2000000)
      common/pawc/h(nh)
      common/gcbank/q(nz)
C
      COMMON/GCUNIT/LIN,LOUT,NUNITS,LUNITS(5)
      INTEGER LIN,LOUT,NUNITS,LUNITS
      COMMON/GCMAIL/CHMAIL
      CHARACTER*132 CHMAIL
      logical lexist
*
      lout = 6
*
      if( iargc() .gt. 0 ) then
        call getarg(1,file)
        ln = lenocc(file)
      else
        print *,'rztog4 <rzfile> [ <call-list> ]'
        print *,'<rzfile> : name of the ZEBRA rz file'
        print *,'<call-list> : name of output call list file'
        print *,'(def: g3calls.dat)'
        print *,'You need to specify an rz file'
        goto 999
      endif

      if (iargc() .gt. 1 ) then
         call getarg(2,ofile)
      else
         ofile='g3calls.dat'
      end if
      inquire(file=file(:ln),exist=lexist)
      if (.not.lexist) then
         write(6,*) 'rz-file "',file(:ln),'" doesn''t exist.'
         stop 1
      end if
*
***     Geant/Zebra initialization
      call gzebra(nz)
      call gzinit
      call grfile(90,file(:ln),'i')
*
***     conversion initialization
      dogeom = .false.
      context = '----'
      lunlist = 98
      lnout = lenocc(ofile)
      open(unit=lunlist,file=ofile(:lnout),status='unknown')
      print *,'Opened call list file ',ofile(:lnout)
*
      luncode = 99
      luncode = 0
      nfile = 1
      call g3source
*
***     do conversion
      call tog4
*
 999  end

