*     This code implementation is the intellectual property of
*     the RD44 GEANT4 collaboration.
*    
*     By copying, distributing or modifying the Program (or any work
*     based on the Program) you indicate your acceptance of this statement,
*     and all its terms.
*    
*     $Id: g3tog4.F,v 1.2 1999-05-04 21:17:03 lockman Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
*    
*
*       G3toG4
*
*       Package to convert Geant3 Fortran geometry code to a call list
*       file to be interpreted by the Geant4 geometry call list
*       interpreter, or alternatively, directly to Geant4 code.
*
*       This set of routines is to be linked in front of, so overriding,
*       the standard Geant library.
*
*       It is possible to execute the Geant calls while at the same
*       time building the call list/Geant4 code. In order to do this,
*       these routines must occupy a different name space to that of
*       the real Geant routines. This is provided by the CALL_GEANT
*       cpp flag. If these routines are compiled with this flag,
*       the routine names begin with K rather than G. eg. GSVOLU
*       becomes KSVOLU. Routine names in your source code must be so
*       converted; a perl script is provided to do this. 
*         $$$ provide the script
*       Under normal circumstances it should *not* be necessary to go
*       through this; it is only necessary if during the geometry
*       generation process your code extracts information from Geant
*       about material already generated.
*
*       Torre Wenaus, LLNL  6/95
*
*  To Do
*  - option to divide generated Geant4 code into separate files/routines
*    based on context
*
************************************************************************
*
      subroutine G3toG4(luni,lunc,chopt)
************************************************************************
*
*       G3toG4
*
*       Initialization/setup routine
*
*       luni (call list), lunc (C++ code) logical unit numbers:
*         lun>0: Open output file on unit lun. Filenames used:
*             g3calls.dat      Call list file
*             g4geom.cc        Geant4 C++ geometry code
*         lun<0: File open has been done by the user. Just write to |lun|
*         lun=0: Don't generate this output.
*            ie.  luni=0:  Don't generate the call list
*                 lunc=0:  Don't generate the Geant4 code
*
*       chopt options:
*       'G'  execute the actual Geant calls as well as building the
*            code/call list. In case users use info obtained from Geant
*            during the geometry building process. THIS IS THE DEFAULT
*            at present:
#define CALL_GEANT
*
************************************************************************
      implicit none
      integer luni, lunc
      character chopt*(*)
#include "G3toG4.inc"
*
      print *,'Initializing Geant3 to Geant4 conversion'
#ifdef CALL_GEANT
c      dogeom = index(chopt,'G') + index(chopt,'g') .ne. 0
      dogeom = .true.
#else
      dogeom = .false.
#endif
      context = '----'
      if (luni.eq.0.and.lunc.eq.0) then
        print *,'G3TOG4: No output requested by user. No output'//
     +  ' will be generated.'
      endif
      lunlist = abs(luni)
      luncode = abs(lunc)
      if (lunlist.ne.0) then
        doclist = .true.
      else
        doclist = .false.
      endif
      if (luncode.ne.0) then
        docode = .true.
      else
        docode = .false.
      endif
***     If lun>0, open the file
      if (lunlist.gt.0) then
        open(unit=lunlist,file='g3calls.dat',status='unknown')
      endif
      if (luncode.gt.0) then
        nfile = 1
        call g3source
      endif
*
      end
*
      subroutine g4init
************************************************************************
************************************************************************
      implicit none
#include "G3toG4.inc"
*
      if (luncode.ne.0) then
        write(luncode,
     + '(''//G4GeometryManager* GeoMgr = new G4GeometryManager();'')')
*     call ctocp('void G3G4init();')
      endif
*
      end
*
      subroutine g3header
************************************************************************
*
************************************************************************
      implicit none
      call g4init
      end

      subroutine g3source
************************************************************************
*
************************************************************************
      implicit none
#include "G3toG4.inc"
      character fname*30
      if (luncode.le.0) return
      if (nfile.gt.1) write(luncode,'(''}'')')
      close(luncode)
      write (fname,'(''G3toG4code_'',i2.2,''.cc'')') nfile
      open(unit=luncode,file=fname,status='unknown')
      write(luncode,'(''#include "G3toG4.hh"'')')
      if (nfile.eq.1) call g3header
      write(luncode,'(/''void G3toG4code_'',i2.2,''()'')') nfile
      write(luncode,'(''{'')')
      call ctocp('// init to 0 avoids "unused" warnings')
      call ctocp('G4int nd=0,nh=0,nv=0,imate=0,itmed=0,nmat=0,')
      call ctocp(' isvol=0,ifield=0,nwhi=0,nwdi=0,idtyp=0,ipart=0,')
      call ctocp(' itrtyp=0,nlmat=0,npar=0,ndvmx=0,numed=0,iaxis=0,')
      call ctocp(
     +     ' ndiv=0,irot=0,ival=0,num=0,nmed=0,nbits[100],mode[6];')
      call ctocp('G4String chnam[100];')
      call ctocp('G4String name="",moth="",attr="",only="",shape="";')
      call ctocp('G4String chset="",chdet="",chali="",chpar="";')
      call ctocp('G4double amass=0.,charge=0.,tlife=0.,parval=0.;')
      call ctocp('G4double c0=0.,step=0.,a=0.,dens=0.,radl=0.,x=0.;')
      call ctocp('G4double y=0.,z=0.,theta1=0.,phi1=0.,theta2=0.;')
      call ctocp('G4double phi2=0.,theta3=0.,phi3=0.,fieldm=0.;')
      call ctocp('G4double tmaxfd=0.,stemax=0.,deemax=0.,epsil=0.;')
      call ctocp('G4double stmin=0.,par[100],fact[100],orig[100];')
      call ctocp('G4double bratio[6],aa[100],zz[100],wmat[100];')
      call ctocp('nbits[0]=mode[0]=0;chnam[0]="";par[0]=0.;')
      call ctocp('fact[0]=orig[0]=bratio[0]=aa[0]=zz[0]=wmat[0]=0.;')      
      call ctocp(' ')
      if (nfile.eq.1) then
*     call ctocp('G3G4init();')
         call ctocp(' ')
      endif
      end

      subroutine g3main
************************************************************************
************************************************************************
      implicit none
#include "G3toG4.inc"
      integer i
*
      close(luncode)
      open(unit=luncode,file='G3toG4code.cc',status='unknown')
      do i=1,nfile
        write(luncode,'('' void G3toG4code_'',i2.2,''();'')') i
      enddo
      call ctocp('void G3toG4code()')
      call ctocp('{')
      do i=1,nfile
        write(luncode,'(''   G3toG4code_'',i2.2,''();'')') i
      enddo
      call ctocp('}')
      close(luncode)
      end

      subroutine g3context(cntxt)
************************************************************************
*
*       g3context
*
*       Set the current geometry code context. eg. context can be used
*       to distinguish code for different subdetectors. The Geant4
*       call list interpreter can then execute the code selectively for
*       a particular context only, if desired. Spaces not allowed.
*
************************************************************************
      implicit none
      character*(*) cntxt
#include "G3toG4.inc"
      context = cntxt
      end
*
      subroutine ctocp(string)
************************************************************************
************************************************************************
      implicit none
      character*(*) string
#include "G3toG4.inc"
      write (luncode,*) string
      end
*
      subroutine rtocp(string,x)
************************************************************************
************************************************************************
      implicit none
      character*(*) string
      real x
#include "G3toG4.inc"
      write(luncode,'(4x,a,'' = '',e14.8,'';'')')
     +  string, x
      end
*
      subroutine artocp(string,ax,n)
************************************************************************
************************************************************************
      implicit none
      character*(*) string
      real ax(*)
      integer n,i
#include "G3toG4.inc"
      do i=1,n
        write(luncode,'(''    '',a,''['',i3,''] = '',e14.8,'';'')')
     +    string, i-1, ax(i)
      enddo
      end
*
      subroutine aitocp(string,ai,n)
************************************************************************
************************************************************************
      implicit none
      character*(*) string
      integer ai(*)
      integer n,i
#include "G3toG4.inc"
      do i=1,n
        write(luncode,'(''    '',a,''['',i3,''] = '',i10,'';'')')
     +    string, i-1, ai(i)
      enddo
      end
*
      subroutine astocp(string,ac,n)
************************************************************************
************************************************************************
      implicit none
      character*(*) string, ac(*)
      integer n,i
#include "G3toG4.inc"
c      write(luncode,'(''    G4String '',a,''['',i3,''];'')') string, n
      do i=1,n
        write(luncode,'(''    '',a,''['',i3,''] = "'',a,''";'')')
     +    string, i-1, ac(i)
      enddo
      end
*
      subroutine g3ldpar(par,npar)
************************************************************************
*
*       g3ldpar
*
************************************************************************
      implicit none
*
      integer npar, i
      real par(*)
#include "G3toG4.inc"
*
      if (npar.gt.0) then
        write(luncode,'(''    par['',i4,''] = '',e14.8,'';'')')
     +    (i-1,par(i),i=1,npar)
      endif
      end
*
      subroutine check_lines
************************************************************************
************************************************************************
      implicit none
#include "G3toG4.inc"
      if (luncode.ne.0) then
        nlines = nlines +1
        if (nlines.gt.maxlines) then
          nfile = nfile +1
          call g3source
          nlines = 0
        endif
      endif
      end
