*     This code implementation is the intellectual property of
*     the RD44 GEANT4 collaboration.
*    
*     By copying, distributing or modifying the Program (or any work
*     based on the Program) you indicate your acceptance of this statement,
*     and all its terms.
*    
*     $Id: g3routines.F,v 1.1 1999-01-07 16:06:53 gunter Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
*    
#define CALL_GEANT

#ifndef CALL_GEANT
      subroutine gsvolu(name, shape, nmed, par, npar, ivol)
#else 
      subroutine Ksvolu(name, shape, nmed, par, npar, ivol)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, shape*4, fmt*150
      integer nmed, npar, ivol, k
      real par(npar)
      character rname*6
#include "G3toG4.inc"
      data rname /'GSVOLU'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsvolu(name, shape, nmed, par, npar, ivol)
#endif
      if (npar.ne.0) call checkshape(name, shape, par, npar)
*     
      if (lunlist.ne.0) then
*        write(lunlist,
*     +    '(a4,1x,a6,1x,a4,1x,a4,2i5,<npar>e15.8)')
*     +    context, rname, name, shape, nmed, npar,
*     +    (par(k),k=1,npar)
         write(fmt,'(A,I2,A)')'(a4,1x,a6,1x,a4,1x,a4,2i5,',max(npar,1),
     >        '(1x,e16.8))'
         write(lunlist,fmt) context, rname, name, shape, nmed, npar,
     +        (par(k),k=1,npar)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call g3ldpar(par,npar)
         write(luncode,1000) name, shape, nmed, npar
 1000    format('G4gsvolu(name="',a,'",shape="',a,'",nmed=',i5,
     +        ',par,npar=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gspos(name, num, moth, x, y, z, irot, only)
#else 
      subroutine Kspos(name, num, moth, x, y, z, irot, only)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4, only*4
      integer num, irot
      real x, y, z
      character rname*6
#include "G3toG4.inc"
      data rname /'GSPOS '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gspos(name, num, moth, x, y, z, irot, only)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,i5,1x,a4,3(1x,e16.8),i5,1x,a4)')
     +        context, rname, name, num, moth, x, y, z, irot, only
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('x',x)
         call rtocp('y',y)
         call rtocp('z',z)
         write(luncode,1000) name,num,moth,irot,only
 1000    format('G4gspos(name="',a,'",num=',i5,',moth="',a,
     +        '",x,y,z,irot=',i5,',only="',a,'");')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsposp(name, num, moth, x, y, z, irot, only, par, npar)
#else 
      subroutine Ksposp(name, num, moth, x, y, z, irot, only, par, npar)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4, only*4
      integer num, irot, npar, k
      real x, y, z, par(npar)
      character rname*6, fmt*150
#include "G3toG4.inc"
      data rname /'GSPOSP'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsposp(name, num, moth, x, y, z, irot, only,
     +     par, npar)
#endif
      if (lunlist.ne.0) then
         do k=1,npar
            if (abs(par(k)).gt.1.e10) then
               print *,'Warning: huge junk value in PAR for GSPOS'
               print *,'  zeroed out. Volume ',name
               par(k) = 0.
            endif
         enddo
*        write(lunlist,
*     +    '(a4,1x,a6,1x,a4,i5,1x,a4,3e15.8,i5,1x,a4,
*     +    i5,<npar>e15.8)')
*     +    context, rname, name, num, moth, x, y, z, irot, only,
*     +    npar,
*     +    (par(k),k=1,npar)
         write(fmt,'(A,A,I2,A)')
     >        '(a4,1x,a6,1x,a4,i5,1x,a4,3(1x,e16.8),',       
     +        'i5,1x,a4,i5,',max(npar,1),'(1x,e16.8))'
         write(lunlist,fmt)
     +        context, rname, name, num, moth, x, y, z, irot, only,
     +        npar,
     +        (par(k),k=1,npar)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('x',x)
         call rtocp('y',y)
         call rtocp('z',z)
         call g3ldpar(par,npar)
         write(luncode,1000) name,num,moth,irot,only,npar
 1000    format('G4gsposp(name="',a,'",num=',i5,',moth="',a,
     +        '",x,y,z,irot=',i5,',only="',a,'",par,npar=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsatt(name, attr, ival)
#else 
      subroutine Ksatt(name, attr, ival)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, attr*4
      integer ival
      character rname*6
#include "G3toG4.inc"
      data rname /'GSATT '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsatt(name, attr, ival)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,i12)')
     +        context, rname, name, attr, ival
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         write(luncode,1000) name,attr,ival
 1000    format('G4gsatt(name="',a,'",attr="',a,'",ival=',i10,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsrotm(irot, theta1, phi1, theta2, phi2,
     +     theta3, phi3)
#else 
      subroutine Ksrotm(irot, theta1, phi1, theta2, phi2,
     +     theta3, phi3)
#endif
************************************************************************
************************************************************************
      implicit none
      integer irot
      real theta1, phi1, theta2, phi2, theta3, phi3
      character rname*6
#include "G3toG4.inc"
      data rname /'GSROTM'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsrotm(irot, theta1, phi1, theta2, phi2,
     +     theta3, phi3)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,i5,6f11.5)')
     +        context, rname, irot, theta1, phi1, theta2, phi2,
     +        theta3, phi3
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('theta1',theta1)
         call rtocp('phi1',phi1)
         call rtocp('theta2',theta2)
         call rtocp('phi2',phi2)
         call rtocp('theta3',theta3)
         call rtocp('phi3',phi3)
         write(luncode,1000) irot
 1000    format('G4gsrotm(irot=',i5,
     +        ',theta1,phi1,theta2,phi2,theta3,phi3);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdvn(name, moth, ndiv, iaxis)
#else 
      subroutine Ksdvn(name, moth, ndiv, iaxis)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4
      integer ndiv, iaxis
      character rname*6
#include "G3toG4.inc"
      data rname /'GSDVN '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdvn(name, moth, ndiv, iaxis)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,i5,i3)')
     +        context, rname, name, moth, ndiv, iaxis
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         write(luncode,1000) name, moth, ndiv, iaxis
 1000    format('G4gsdvn(name="',a,'",moth="',a,'",ndiv=',i3,
     +        ',iaxis=',i1,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdvt(name, moth, step, iaxis, numed, ndvmx)
#else 
      subroutine Ksdvt(name, moth, step, iaxis, numed, ndvmx)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4
      real step
      integer iaxis, numed, ndvmx
      character rname*6
#include "G3toG4.inc"
      data rname /'GSDVT '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdvt(name, moth, step, iaxis, numed, ndvmx)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,(1x,e16.8),3i5)')
     +        context, rname, name, moth, step, iaxis, numed, ndvmx
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('step',step)
         write(luncode,1000) name,moth,iaxis,numed,ndvmx
 1000    format('G4gsdvt(name="',a,'",moth="',a,'",step,iaxis=',
     +        i1,',numed=',i4,',ndvmx=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdvx(name, moth, ndiv, iaxis, step, c0, numed, ndvmx)
#else 
      subroutine Ksdvx(name, moth, ndiv, iaxis, step, c0, numed, ndvmx)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4
      integer ndiv, iaxis, numed, ndvmx
      real step, c0
      character rname*6
#include "G3toG4.inc"
      data rname /'GSDVX '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdvx(name, moth, ndiv, iaxis, step, c0, numed,
     +     ndvmx)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,i5,i3,2(1x,e16.8),2i5)')
     +        context, rname, name, moth, ndiv, iaxis,step, c0,
     +        numed, ndvmx
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('step',step)
         call rtocp('c0',c0)
         write(luncode,1000) name,moth,ndiv,iaxis,numed,ndvmx
 1000    format('G4gsdvx(name="',a,'",moth="',a,'",ndiv=',i3,',iaxis=',
     +        i1,',step,c0,numed=',i4,',ndvmx=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdvn2(name, moth, ndiv, iaxis, c0, numed)
#else 
      subroutine Ksdvn2(name, moth, ndiv, iaxis, c0, numed)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4
      integer ndiv, iaxis, numed
      real c0
      character rname*6
#include "G3toG4.inc"
      data rname /'GSDVN2'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdvn2(name, moth, ndiv, iaxis, c0, numed)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,i5,i3,(1x,e16.8),i5)')
     +        context, rname, name, moth, ndiv, iaxis, c0, numed
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('c0',c0)
         write(luncode, 1000) name,moth,ndiv,iaxis,numed
 1000    format('G4gsdvn2(name="',a,'",moth="',a,'",ndiv=',i3,',iaxis=',
     +        i1,',c0,numed=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdvt2(name, moth, step, iaxis, c0, numed, ndvmx)
#else 
      subroutine Ksdvt2(name, moth, step, iaxis, c0, numed, ndvmx)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*4, moth*4
      integer iaxis, numed, ndvmx
      real step, c0
      character rname*6
#include "G3toG4.inc"
      data rname /'GSDVT2'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdvt2(name, moth, step, iaxis, c0, numed, ndvmx)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,(1x,e16.8),i3,(1x,e16.8),2i5)')
     +        context, rname, name, moth, step, iaxis, c0, numed, ndvmx
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('step',step)
         call rtocp('c0',c0)
         write(luncode,1000) name,moth,iaxis,numed,ndvmx
 1000    format('G4gsdvt2(name="',a,'",moth="',a,'",step,iaxis=',
     +        i1,',c0,numed=',i4,',ndvmx=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsmate(imate, name, a, z, dens, radl, absl, ubf, nwbf)
#else 
      subroutine Ksmate(imate, name, a, z, dens, radl, absl, ubf, nwbf)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*(*)
      integer imate, nwbf, k
      real a, z, dens, radl, absl, ubf(nwbf)
      character rname*6, fmt*150
#include "G3toG4.inc"
      data rname /'GSMATE'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsmate
     +     (imate, name, a, z, dens, radl, absl, ubf, nwbf)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,i5,1x,''"'',a,''"'',4e15.8,i3,<nwbf>e15.8)')
*     +        context, rname, imate, name, a, z, dens, radl,
*     +        nwbf, (ubf(k), k=1,nwbf)
         write(fmt,'(A,I3,A)')
     >        '(a4,1x,a6,i5,1x,''"'',a,''"'',4(1x,e16.8),i3,',
     >        max(nwbf,1),'(1x,e16.8))'
         write(lunlist,fmt)
     +        context, rname, imate, name, a, z, dens, radl,
     +        nwbf, (ubf(k), k=1,nwbf)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('a',a)
         call rtocp('z',z)
         call rtocp('dens',dens)
         call rtocp('radl',radl)
         call g3ldpar(ubf,nwbf)
         write(luncode,1000) imate, name, nwbf
 1000    format('G4gsmate(imate=',i4,',name="',a,
     +        '",a,z,dens,radl,npar=',i4,',par);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsmixt(imate, name, a, z, dens, nlmat, wmat)
#else 
      subroutine Ksmixt(imate, name, a, z, dens, nlmat, wmat)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*(*)
      integer imate, nlmat, k, nlmata
      real a(*), z(*), dens, wmat(*)
      character rname*6, fmt*150
#include "G3toG4.inc"
      data rname /'GSMIXT'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsmixt
     +     (imate, name, a, z, dens, nlmat, wmat)
#endif
      if (lunlist.ne.0) then
         nlmata = abs(nlmat)
*         write(lunlist,
*     +        '(a4,1x,a6,i5,1x,''"'',a,''"'',5e15.8,i3,<nlmata>e15.8,
*     +        <nlmata>e15.8,<nlmata>e15.8)')
*     +        context, rname, imate, name, dens,
*     +        nlmat, 
*     +        (a(k), k=1,abs(nlmat)),
*     +        (z(k), k=1,abs(nlmat)),
*     +        (wmat(k), k=1,abs(nlmat))
         write(fmt,'(A,I3,A,I3,A,I3,A)')
     +        '(a4,1x,a6,i5,1x,",a,",5(1x,e16.8),i3,',max(nlmata,1),
     >        '(1x,e16.8),',max(nlmata,1),'(1x,e16.8),',
     >        max(nlmata,1),'(1x,e16.8))'
         write(lunlist,fmt)
     +        context, rname, imate, name, dens,
     +        nlmat, 
     +        (a(k), k=1,abs(nlmat)),
     +        (z(k), k=1,abs(nlmat)),
     +        (wmat(k), k=1,abs(nlmat))
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('dens',dens)
         call artocp('aa',a,abs(nlmat))
         call artocp('zz',z,abs(nlmat))
         call artocp('wmat',wmat,abs(nlmat))
         write(luncode,1000) imate,name,nlmat
 1000    format('G4gsmixt(imate=',i5,',name="',a,
     +        '",aa,zz,dens,nlmat=',i3,',wmat);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gstmed(
     +     itmed, name, nmat, isvol, ifield, fieldm,
     +     tmaxfd, stemax, deemax, epsil, stmin, ubuf, nwbuf)
#else 
      subroutine Kstmed(
     +     itmed, name, nmat, isvol, ifield, fieldm,
     +     tmaxfd, stemax, deemax, epsil, stmin, ubuf, nwbuf)
#endif
************************************************************************
************************************************************************
      implicit none
      character name*(*)
      integer itmed, nmat, isvol, ifield, nwbuf, k
      real fieldm, tmaxfd, stemax, deemax, epsil, stmin, ubuf(nwbuf)
      character rname*6, fmt*150
#include "G3toG4.inc"
      data rname /'GSTMED'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gstmed(
     +     itmed, name, nmat, isvol, ifield, fieldm,
     +     tmaxfd, stemax, deemax, epsil, stmin, ubuf, nwbuf)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,i5,1x,''"'',a,''"'',3i3,6e15.8,i3,<nwbuf>e15.8)')
*     +        context, rname, itmed, name, nmat, isvol, ifield, fieldm,
*     +        tmaxfd, stemax, deemax, epsil, stmin,
*     +        nwbuf, (ubuf(k),k=1,nwbuf)
         write(fmt,'(A,I3,A)') 
     >        '(a4,1x,a6,i5,1x,''"'',a,''"'',3i3,6(1x,e16.8),i3,',
     >        max(nwbuf,1),'(1x,e16.8))'
         write(lunlist,fmt)
     +        context, rname, itmed, name, nmat, isvol, ifield, fieldm,
     +        tmaxfd, stemax, deemax, epsil, stmin,
     +        nwbuf, (ubuf(k),k=1,nwbuf)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('fieldm',fieldm)
         call rtocp('tmaxfd',tmaxfd)
         call rtocp('stemax',stemax)
         call rtocp('deemax',deemax)
         call rtocp('epsil',epsil)
         call rtocp('stmin',stmin)
         call g3ldpar(ubuf,nwbuf)
         write(luncode,1000) itmed,name,nmat,isvol,ifield,nwbuf
 1000    format('G4gstmed(itmed=',i4,',name="',a,'",nmat=',i4,
     +        ',isvol=',i2,',ifield=',i2,',',/
     +        '    fieldm,tmaxfd,stemax,deemax,epsil,stmin,par,npar=',
     +        i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gstpar(itmed, chpar, parval)
#else 
      subroutine Kstpar(itmed, chpar, parval)
#endif
************************************************************************
************************************************************************
      implicit none
      character chpar*(*)
      integer itmed
      real parval
      character rname*6
#include "G3toG4.inc"
      data rname /'GSTPAR'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gstpar (itmed, chpar, parval)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,i5,1x,a4,(1x,e16.8))')
     +        context, rname, itmed, chpar, parval
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         write(luncode,1000) itmed, chpar, parval
 1000    format('G4gstpar(itmed=',i4,',chpar="',a,'",parval=',
     +        (1x,e16.8),');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gspart(
     +     ipart, chpar, itrtyp, amass, charge, tlife, ub, nwb)
#else 
      subroutine Kspart(
     +     ipart, chpar, itrtyp, amass, charge, tlife, ub, nwb)
#endif
************************************************************************
************************************************************************
      implicit none
      character chpar*(*)
      integer ipart, itrtyp, nwb, k
      real amass, charge, tlife, ub(nwb)
      character rname*6, fmt*150
#include "G3toG4.inc"
      data rname /'GSPART'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gspart(
     +     ipart, chpar, itrtyp, amass, charge, tlife, ub, nwb)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,i5,1x,''"'',a,''"'',i3,3e15.8,i3,<nwb>e15.8)')
*     +        context, rname, ipart, chpar, itrtyp, amass, charge, tlife,
*     +        nwb, (ub(k), k=1,nwb)
         write(fmt,'(A,I3,A)') 
     >        '(a4,1x,a6,i5,1x,''"'',a,''"'',i3,3(1x,e16.8),i3,',        
     >        max(nwb,1),'(1x,e16.8))'
         write(lunlist,fmt)
     +        context, rname, ipart, chpar, itrtyp, amass, charge, 
     >        tlife,
     +        nwb, (ub(k), k=1,nwb)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call rtocp('amass',amass)
         call rtocp('charge',charge)
         call rtocp('tlife',tlife)
         call g3ldpar(ub,nwb)
         write(luncode,1000) ipart,chpar,itrtyp,nwb
 1000    format('G4gspart(ipart=',i8,',chpar="',a,'",itrtyp=',i8,
     +        ',amass,charge,'/'         tlife,par,npar=',i4,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdk(ipart, bratio, mode)
#else 
      subroutine Ksdk(ipart, bratio, mode)
#endif
************************************************************************
************************************************************************
      implicit none
      integer ipart, mode(6)
      real bratio(6)
      character rname*6
#include "G3toG4.inc"
      data rname /'GSDK  '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdk(ipart, bratio, mode)
#endif
      if (lunlist.ne.0) then
***   6 is prefixed to the arrays for consistency with other
***   array treatments (count precedes the array)
         write(lunlist,
     +        '(a4,1x,a6,i5,i3,6(1x,e16.8),6i8)')
     +        context, rname, ipart, 6, bratio, mode
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call artocp('bratio',bratio,6)
         call aitocp('mode',mode,6)
         write(luncode,1000) ipart
 1000    format('G4gsdk(ipart=',i8,',bratio,mode);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdet(chset, chdet, nv, chnam, nbits, idtyp, nwhi,
     +     nwdi, iset, idet)
#else 
      subroutine Ksdet(chset, chdet, nv, chnam, nbits, idtyp, nwhi,
     +     nwdi, iset, idet)
#endif
************************************************************************
************************************************************************
      implicit none
      integer nv, nbits(nv), idtyp, nwhi, nwdi, iset, idet, k
      character rname*6, chset*4, chdet*4, chnam(nv)*4, fmt*150
#include "G3toG4.inc"
      data rname /'GSDET '/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdet(chset, chdet, nv, chnam, nbits, idtyp,
     +     nwhi, nwdi, iset, idet)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,1x,a4,1x,a4,i5,<nv>(1x,a4),<nv>i10,i10,2i5)')
*     +        context, rname, chset, chdet, nv, (chnam(k), k=1,nv),
*     +        (nbits(k), k=1,nv), idtyp, nwhi, nwdi
         write(fmt,'(A,I3,A,I3,A)')'(a4,1x,a6,1x,a4,1x,a4,i5,',
     >        max(nv,1),'(1x,a4),',max(nv,1),'i10,i10,2i5)'
         write(lunlist,fmt)
     +        context, rname, chset, chdet, nv, (chnam(k), k=1,nv),
     +        (nbits(k), k=1,nv), idtyp, nwhi, nwdi
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call astocp('chnam',chnam,nv)
         call aitocp('nbits',nbits,nv)
         write(luncode,1000) chset, chdet, nv, idtyp, nwhi, nwdi
 1000    format('G4gsdet(chset="',a,'",chdet="',a,'",nv=',i3,
     +        ',chnam,nbits,idtyp=',i8,','/
     +        '        nwhi=',i8,',nwdi=',i8,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdetv(chset, chdet, idtyp, nwhi, nwdi, iset, idet)
#else 
      subroutine Ksdetv(chset, chdet, idtyp, nwhi, nwdi, iset, idet)
#endif
************************************************************************
************************************************************************
      implicit none
      integer idtyp, nwhi, nwdi, iset, idet
      character rname*6, chset*4, chdet*4
#include "G3toG4.inc"
      data rname /'GSDETV'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdetv(chset, chdet, idtyp,
     +     nwhi, nwdi, iset, idet)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,i10,2i5)')
     +        context, rname, chset, chdet, idtyp, nwhi, nwdi
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         write(luncode,1000) chset, chdet, idtyp, nwhi, nwdi
 1000    format('G4gsdetv(chset="',a,'",chdet="',a,'",idtyp=',i8,
     +        ',nwhi=',i8,',nwdi=',i8,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdeta(chset, chdet, chali, nwhi, nwdi, iali)
#else 
      subroutine Ksdeta(chset, chdet, chali, nwhi, nwdi, iali)
#endif
************************************************************************
************************************************************************
      implicit none
      integer nwhi, nwdi, iali
      character rname*6, chset*4, chdet*4, chali*4
#include "G3toG4.inc"
      data rname /'GSDETA'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdeta(chset, chdet, chali, nwhi, nwdi, iali)
#endif
      if (lunlist.ne.0) then
         write(lunlist,
     +        '(a4,1x,a6,1x,a4,1x,a4,1x,a4,2i5)')
     +        context, rname, chset, chdet, chali, nwhi, nwdi
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         write(luncode,1000) chset, chdet, chali, nwhi, nwdi
 1000    format('G4gsdeta(chset="',a,'",chdet="',a,'",chali="',a,
     +        '",nwhi=',i8,',nwdi=',i8,');')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdeth(chset, chdet, nh, chnam, nbits, orig, fact)
#else 
      subroutine Ksdeth(chset, chdet, nh, chnam, nbits, orig, fact)
#endif
************************************************************************
************************************************************************
      implicit none
      integer nh, nbits(nh), k
      real orig(nh), fact(nh)
      character rname*6, chset*4, chdet*4, chnam(nh)*4, fmt*150
#include "G3toG4.inc"
      data rname /'GSDETH'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdeth(chset, chdet, nh, chnam, nbits,
     +     orig, fact)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,1x,a4,1x,a4,i5,<nh>(1x,a4),<nh>i5,<nh>e15.8,
*     +        <nh>e15.8)')
*     +        context, rname, chset, chdet, nh, (chnam(k), k=1,nh),
*     +        (nbits(k), k=1,nh), (orig(k), k=1,nh), (fact(k), k=1,nh)
         write(fmt,'(A,I3,A,I3,A,I3,A,I3,A)')
     >        '(a4,1x,a6,1x,a4,1x,a4,i5,',max(nh,1),'(1x,a4),',
     >        max(nh,1),'i5,',max(nh,1),'(1x,e16.8),',max(nh,1),
     >        '(1x,e16.8))'
         write(lunlist, fmt)
     +        context, rname, chset, chdet, nh, (chnam(k), k=1,nh),
     +        (nbits(k), k=1,nh), (orig(k), k=1,nh), (fact(k), k=1,nh)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call astocp('chnam',chnam,nh)
         call aitocp('nbits',nbits,nh)
         call artocp('orig',orig,nh)
         call artocp('fact',fact,nh)
         write(luncode,1000) chset,chdet,nh
 1000    format('G4gsdeth(chset="',a,'",chdet="',a,'",nh=',i4,
     +        ',chnam,nbits,orig,fact);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdetd(chset, chdet, nd, chnam, nbits)
#else 
      subroutine Ksdetd(chset, chdet, nd, chnam, nbits)
#endif
************************************************************************
************************************************************************
      implicit none
      integer nd, nbits(nd), k
      character rname*6, chset*4, chdet*4, chnam(nd)*4, fmt*150
#include "G3toG4.inc"
      data rname /'GSDETD'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdetd(chset, chdet, nd, chnam, nbits)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,1x,a4,1x,a4,i5,<nd>(1x,a4),<nd>i5)')
*     +        context, rname, chset, chdet, nd, (chnam(k), k=1,nd),
*     +        (nbits(k), k=1,nd)
         write(fmt,'(A,I3,A,I3,A)')
     +        '(a4,1x,a6,1x,a4,1x,a4,i5,',max(nd,1),'(1x,a4),',
     >        max(nd,1),'i5)'
         write(lunlist,fmt)
     +        context, rname, chset, chdet, nd, (chnam(k), k=1,nd),
     +        (nbits(k), k=1,nd)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call astocp('chnam',chnam,nd)
         call aitocp('nbits',nbits,nd)
         write(luncode,1000) chset, chdet, nd
 1000    format('G4gsdetd(chset="',a,'",chdet="',a,'",nd=',i4,
     +        ',chnam,nbits);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine gsdetu(chset, chdet, nupar, upar)
#else 
      subroutine Ksdetu(chset, chdet, nupar, upar)
#endif
************************************************************************
************************************************************************
      implicit none
      integer nupar, k
      real upar(nupar)
      character rname*6, chset*4, chdet*4, fmt*150
#include "G3toG4.inc"
      data rname /'GSDETU'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call gsdetu(chset, chdet, nupar, upar)
#endif
      if (lunlist.ne.0) then
*         write(lunlist,
*     +        '(a4,1x,a6,1x,a4,1x,a4,i5,<nupar>e15.8)')
*     +        context, rname, chset, chdet, nupar, (upar(k), k=1,nupar)
         write(fmt,'(A,I3,A)')
     +        '(a4,1x,a6,1x,a4,1x,a4,i5,',max(nupar,1),'(1x,e16.8))'
         write(lunlist,fmt)
     +        context, rname, chset, chdet, nupar, (upar(k), k=1,nupar)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''{'')')
         call g3ldpar(upar,nupar)
         write(luncode,1000) chset, chdet, nupar
 1000    format('G4gsdetu(chset="',a,'",chdet="',a,'",npar=',
     +        i4,',par);')
         write(luncode,'(''}'')')
      endif
*     
      end
*     
#ifndef CALL_GEANT
      subroutine ggclos
#else 
      subroutine kgclos
#endif
************************************************************************
************************************************************************
      implicit none
      character rname*6
#include "G3toG4.inc"
      data rname /'GGCLOS'/
*     
      call check_lines
#ifdef CALL_GEANT
      if (dogeom) call ggclos
#endif
      if (lunlist.ne.0) then
         write(lunlist,'(a4,1x,a6)') context, rname
         close(lunlist)
      endif
      if (luncode.ne.0) then
         write(luncode,'(''//GeoMgr->CloseGeometry();'')')
         write(luncode,'(''}'')')
         call g3main
         close(luncode)
      endif
*     
      end

      subroutine checkshape(name, shape, par, npar)
      implicit none
************************************************************************
* convert TRAP, PARA and GTRA to external form
************************************************************************
      character name*4, shape*4
      real ph, par(*), tt, raddeg
      integer npar
      
      raddeg = 180./3.1415926

      if (shape(1:3).eq.'BOX'.and.npar.ne.3) then
         print *,'!! error, BOX with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'TRD1'.and.npar.ne.4) then
         print *,'!! error, TRD1 with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'TRD2'.and.npar.ne.5) then
         print *,'!! error, TRD2 with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'TRAP'.and.npar.ne.35.and.npar.ne.11) then
***   G3 sets 11 to 35. Why?
         print *,'!! error, TRAP with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'TUBE'.and.npar.ne.3) then
         print *,'!! error, TUBE with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'TUBS'.and.npar.ne.5) then
         print *,'!! error, TUBS with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'CONE'.and.npar.ne.5) then
         print *,'!! error, CONE with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'CONS'.and.npar.ne.7) then
         print *,'!! error, CONS with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'SPHE'.and.npar.ne.6) then
         print *,'!! error, SPHE with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'PARA'.and.npar.ne.6) then
         print *,'!! error, PARA with ',npar,' parameters, vol ',name
      endif
      if (shape.eq.'PARA') then
*
*  **    PARA
*
         ph = 0.
         if (par(5).ne.0.) ph = atan2(par(6),par(5))*raddeg
         tt = sqrt(par(5)**2+par(6)**2)
         par(4) = atan(par(4))*raddeg
         if (par(4).gt.90.0) par(4) = par(4)-180.0
         par(5) = atan(tt)*raddeg
         if (ph.lt.0.0) ph = ph + 360.0
         par(6) = PH
      end if
      if (shape.eq.'TRAP') then
*
*  **    TRAP
*
         npar=11
         ph = 0.
         if (par(2).ne.0.) ph = atan2(par(3),par(2))*raddeg
         tt = sqrt(par(2)**2+par(3)**2)
         par(2) = atan(tt)*raddeg
         if (ph.lt.0.0) ph = ph+360.0
         par(3) = ph
         par(7) = atan(par(7))*raddeg
         if (par(7).gt.90.0) par(7) = par(7)-180.0
         par(11)= atan(par(11))*raddeg
         if (par(11).gt.90.0) par(11) = par(11)-180.0

      end if
      end
