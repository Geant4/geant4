*     This code implementation is the intellectual property of
*     the RD44 GEANT4 collaboration.
*    
*     By copying, distributing or modifying the Program (or any work
*     based on the Program) you indicate your acceptance of this statement,
*     and all its terms.
*    
*     $Id: tog4.F,v 1.2 1999-05-04 21:17:14 lockman Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
*    
      subroutine tog4
************************************************************************
*     
*     tog4
*     
*     Perform the translation to G4
*     
************************************************************************
      implicit none
#include "gcbank.inc"
      integer maxdivols
      parameter (maxdivols=20000)
      integer nvol, nrotm, nmate, ntmed, nset, i, jma, nmixt, k, nin,
     >     jdiv, jd, iaxis, ivo, ndiv, numed, npar, natt, ivol, jin,
     >     nparv, npard, nr, irot, konly, nwbuf, isvol, nmat, ifield,
     >     nbits(5000), idtyp, nwhi, nwdi, iset, idet, j, in, jmx,
     >     jdh, jdd, jdu, ndet, nn, nupar, npos, ndvol, ndivols, ii, 
     >     npositioned, iia(10000), imate, smixt

      real c0, step, x, y, z, a, dens, radl, absl, fact(5000), 
     >     fieldm, tmaxfd, stemax, deemax, epsil, stmin, orig(5000), 
     >     upar(5000)

      character shape*4, name*4, dname*4, chonly*4, chmat*20, chtmed*20,
     >     chset*4, chdet*4, chnms(5000)*4, divols(maxdivols)*4
*     
      npositioned = 0
*     
***   count materials and convert
      call bankcnt(jmate,iia, nmate)
      print *,'Materials: ',nmate
      do imate=1,nmate
         ii=iia(imate)
         jma = lq(jmate-ii)
         call uhtoc(iq(jma+1),4,chmat,20)
         a = q(jma+6)
         z = q(jma+7)
         dens = q(jma+8)
         radl = q(jma+9)
         absl = q(jma+10)
         nwbuf = iq(jma-1)-11
         if (jma.gt.0) then
            smixt=q(jma+11)
            nmixt=abs(smixt)
            if (nmixt.le.1) then
               write(6,101) imate, chmat, a, z, dens, radl, absl
               call ksmate(ii, chmat, a, z, dens, radl, absl,
     >              q(jma+12), nwbuf)
            else
               jmx = lq(jma-5)
               write(6,102) imate, chmat, a, z, dens, radl, absl,
     >              (j,q(jmx+j),q(jmx+nmixt+j),q(jmx+2*nmixt+j), 
     >              j=1,nmixt)
               call ksmixt(ii, chmat, q(jmx+1), q(jmx+nmixt+1),
     >              dens, smixt, q(jmx+2*nmixt+1))
            end if
         end if
      enddo
 101  format(1x,i5,1x,A12,f6.2,f5.1,f8.2,2f9.2)
 102  format(1x,i5,1x,A12,f6.2,f5.1,f8.2,2f9.2,1x,i2, f6.2, f5.1, 
     >     f6.2/(57x, i2, f6.2, f5.1, f6.2))
*     
***   count tracking media and convert
      call bankcnt(jtmed,iia, ntmed)
      print *,'Media: ',ntmed
      do i=1,ntmed
         ii=iia(i)
         j = lq(jtmed-ii)
         call uhtoc(iq(j+1),4,chtmed,20)
         nmat = q(j+6)
         isvol = q(j+7)
         ifield = q(j+8)
         fieldm = q(j+9)
         tmaxfd = q(j+10)
         stemax = q(j+11)
         deemax = q(j+12)
         epsil = q(j+13)
         stmin = q(j+14)
         nwbuf = iq(j-1) -14
         call kstmed(ii,chtmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
     +        deemax,epsil,stmin,q(j+15),nwbuf)
      enddo
*     
***   count rotation matrices and convert
      call bankcnt(jrotm,iia, nrotm)
      print *,'Rotms: ',nrotm
      do i=1,nrotm
         ii=iia(i) 
         j = lq(jrotm-ii)
         call ksrotm(ii,q(j+11),q(j+12),q(j+13),q(j+14),q(j+15),q(j+16))
      enddo
*     
***   count volumes
      npos = 0
      call bankcnt(jvolum,iia, nvol)
      print *,'Volumes: ',nvol
***   pull out the names of the volumes which are subvolumes of
***   divided volumes (gsvolu should not be called on these)
      ndivols = 0
      do i=1, nvol
         ii=iia(i)
         j = lq(jvolum-ii)
         nin = q(j+3)
         if (nin.lt.0) then
            jdiv = lq(j-1)
            ivo = q(jdiv+2)
            call uhtoc(iq(jvolum+ivo),4,dname,4)
            ndivols = ndivols +1
            if (ndivols.gt.maxdivols) then
               ndivols = maxdivols
               print *,
     +              '!!!ERROR!!! ndivols array exhausted. ',
     +              'Too many divisions.'
            endif
            divols(ndivols) = dname
         endif
      enddo
***   create the logical volumes (gsvolu)
      ndvol = 0
      do i=1, nvol
         ii=iia(i)
         j = lq(jvolum-ii)
         call uhtoc(iq(jvolum+ii),4,name,4)
         call jshape(q(j+2),shape)
         nin = q(j+3)
         numed = q(j+4)
         npar = q(j+5)
         natt = q(j+6)
         do k=1, ndivols
            if (divols(k).eq.name) then
               ndvol = ndvol +1
c     print *,'Division volume ',name,'; no gsvolu call.'
               goto 11
            endif
         enddo
         call ksvolu(name, shape, numed, q(j+7), npar, ivol)
 11      continue
      enddo
      print *,'Divided volumes: ',ndvol

***   properties of the mother volume
      call uhtoc(iq(jvolum+1),4,name,4)
      j=lq(jvolum-1)
      call jshape(q(j+2),shape)
      print *,'mother volume: ',name,' shape: ',shape
***   convert physical volumes
      do i=1, nvol
         ii=iia(i)
         j = lq(jvolum-ii)
         call uhtoc(iq(jvolum+ii),4,name,4)
         nin = q(j+3)
         numed = q(j+4)
         npar = q(j+5)
         if (nin.lt.0) then
*     ! divided volume
            jdiv = lq(j-1)
            iaxis = q(jdiv+1)
            ivo = q(jdiv+2)
            call uhtoc(iq(jvolum+ivo),4,dname,4)
            jd = lq(jvolum-ivo)
            numed = q(jd+4)
            ndiv = q(jdiv+3)
            c0 = q(jdiv+4)
            step = q(jdiv+5)
            call ksdvn2(dname, name, ndiv, iaxis, c0, numed)
         else if (nin.gt.0) then
*     ! volume not divided. Handle positioning of daughter vols
            do in=1, nin
               jin = lq(j-in)
               ivo = q(jin+2)
               call uhtoc(iq(jvolum+ivo),4,dname,4)
               jd = lq(jvolum-ivo)
               nparv = q(jd+5)  ! NPAR declared in the GSVOLU call
               nr = q(jin+3)
               irot = q(jin+4)
               x = q(jin+5)
               y = q(jin+6)
               z = q(jin+7)
               konly = q(jin+8)
               if (konly.ne.0) then
                  chonly = 'ONLY'
               else
                  chonly = 'MANY'
               endif
               npard = iq(jin-1) -9
               npositioned = npositioned +1
               if (nparv.eq.0) then
*     ! use GSPOSP
                  call ksposp(dname, nr, name, x, y, z, irot, chonly,
     +                 q(jin+10), npard)
               else
*     ! GSPOS
                  call kspos(dname, nr, name, x, y, z, irot, chonly)
               endif
            enddo
         endif
      enddo
*     
***   count sensitive detectors
      call bankcnt(jset,iia, nset)
      print *,'Sets: ',nset
      do i=1,nset
         ii=iia(i)
         j = lq(jset-ii)
         call uhtoc(iq(jset+i),4,chset,4)
         ndet = iq(j-1)
         do k=1,ndet
            jd = lq(j-k)
            call uhtoc(iq(j+k),4,chdet,4)
            call gfdet(chset, chdet, nn, chnms, nbits, idtyp,
     +           nwhi, nwdi, iset, idet)
            call ksdet(chset, chdet, nn, chnms, nbits, idtyp,
     +           nwhi, nwdi, iset, idet)
            jdh = lq(jd-1)
            if (jdh.ne.0) then
               call gfdeth(chset,chdet,nn,chnms,nbits,orig,fact)
               call ksdeth(chset,chdet,nn,chnms,nbits,orig,fact)
            endif
            jdd = lq(jd-2)
            if (jdd.ne.0) then
               call gfdetd(chset,chdet,nn,chnms,nbits)
               call ksdetd(chset,chdet,nn,chnms,nbits)
            endif
            jdu = lq(jd-3)
            if (jdu.ne.0) then
               call gfdetu(chset,chdet,100,nupar,upar)
               call ksdetu(chset,chdet,nupar,upar)
            endif
         enddo
      enddo
      print *,'Positioned volumes (gspos, gsposp):',npositioned
*     
      call kgclos
*     
      end

      subroutine bankcnt(link,iia,nbanks)
************************************************************************
************************************************************************
      implicit none
#include "gcbank.inc"
      integer i, link, nbanks, iia(*)
*     
      nbanks=0
      if (link.eq.0) return
C*      do i=1,9999999
      do i=1,iq(link-2)
C*         if(lq(link-nbanks-1).eq.0.or.iq(link-2).eq.nbanks) goto 10
         if(lq(link-i).ne.0)then
           nbanks = nbanks +1
           iia(nbanks)=i
         endif
      enddo
 10   continue
      end
