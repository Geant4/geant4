*    
*     ********************************************************************
*     * License and Disclaimer                                           *
*     *                                                                  *
*     * The  Geant4 software  is  copyright of the Copyright Holders  of *
*     * the Geant4 Collaboration.  It is provided  under  the terms  and *
*     * conditions of the Geant4 Software License,  included in the file *
*     * LICENSE and available at  http://cern.ch/geant4/license .  These *
*     * include a list of copyright holders.                             *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.  Please see the license in the file  LICENSE  and URL above *
*     * for the full disclaimer and the limitation of liability.         *
*     *                                                                  *
*     * This  code  implementation is the result of  the  scientific and *
*     * technical work of the GEANT4 collaboration.                      *
*     * By using,  copying,  modifying or  distributing the software (or *
*     * any work based  on the software)  you  agree  to acknowledge its *
*     * use  in  resulting  scientific  publications,  and indicate your *
*     * acceptance of all terms of the Geant4 Software license.          *
*     ********************************************************************
*    
*    
*     $Id: rztog4.F 67982 2013-03-13 10:36:03Z gcosmo $
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

