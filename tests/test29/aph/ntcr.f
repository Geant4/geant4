*    
*     ********************************************************************
*     * DISCLAIMER                                                       *
*     *                                                                  *
*     * The following disclaimer summarizes all the specific disclaimers *
*     * of contributors to this software. The specific disclaimers,which *
*     * govern, are listed with their locations in:                      *
*     *   http://cern.ch/geant4/license                                  *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.                                                             *
*     *                                                                  *
*     * This  code  implementation is the  intellectual property  of the *
*     * GEANT4 collaboration.                                            *
*     * By copying,  distributing  or modifying the Program (or any work *
*     * based  on  the Program)  you indicate  your  acceptance of  this *
*     * statement, and all its terms.                                    *
*     ********************************************************************
*    
      subroutine ntcr
C
      logical      TRIGEV
      data         TRIGEV/.true./
      logical      TRIGTR
      data         TRIGTR/.true./
	  parameter   (MaxEvt=9999999)
C
      DATA NPRIME/1000/
	  parameter   (NCOLIN=28)
      REAL        XTUPIN(NCOLIN)
      CHARACTER*8 INTAGS(NCOLIN)
      DATA INTAGS/'Nevt','Mtot','Mgam','Mpim',
     &            'Mpip','Mpi0','MKp','MK0','MKm','MaK0','Meta','Metap',
     &            'Mrhom','Mrhop','Mrho0','Momega','Mphi',
     &            'MKS0','MKSC','MaKS0','MaKSC',
     &            'NN','PDG','m','Px','Py','Pz','E'/
C	
	  parameter   (NCOLEM=28)
      REAL        XTUPEM(NCOLEM)
      CHARACTER*8 EMTAGS(NCOLEM)
      DATA EMTAGS/'Nevt','Mtot','Mgam','Mpim',
     &            'Mpip','Mpi0','MKp','MK0','MKm','MaK0','Meta','Metap',
     &            'Mrhom','Mrhop','Mrho0','Momega','Mphi',
     &            'MKS0','MKSC','MaKS0','MaKSC',
     &            'NN','PDG','m','Px','Py','Pz','E'/
C	
  	  parameter   (NCOLEV=21)
      REAL        XTUPEV(NCOLEV)
      CHARACTER*8 EVTAGS(NCOLEV)
      DATA EVTAGS/'Nevt','Mtot','Mgam','Mpim','Mpip','Mpi0',
     &            'MKp','MK0','MKm','MaK0','Meta','Metap',
     &            'Mrhom','Mrhop','Mrho0','Momega','Mphi',
     &            'MKS0','MKSC','MaKS0','MaKSC'/
C
      parameter    (NTRMAX=200)
      INTEGER   itr(NTRMAX)
      INTEGER  ipdg(NTRMAX)
      REAL    vmass(NTRMAX)
      REAL       px(NTRMAX)
      REAL       py(NTRMAX)
      REAL       pz(NTRMAX)
      REAL        e(NTRMAX)
C
      vector temp(1) 
      vector ssg2(1) 
      vector eep(1) 
      vector cut(1) 
      vector gfc(1) 
      vector sqc(1) 
      vector vnevt(1) 
      vector vnp(1) 
      vector d1(1) 
      vector d2(1) 
      vector d3(1) 
      vector d4(1) 
      vector d5(1) 
      vector d6(1) 
      vector d7(1) 
      vector d8(1) 
      vector d9(1) 
      vector d10(1) 
      vector d11(1) 
      vector d12(1) 
C
      OPEN(UNIT=4,STATUS='OLD')
      OPEN(UNIT=44,STATUS='OLD')
C
C             NTUPLE direct access file (unit 35) must be opened in
C             calling program
C
      CALL HBOOK1(  1,'N events generated',1,0.0,2.0,0)
      CALL HBOOKN(20,'CHIPStest Inclusive',NCOLIN,'PAWC',NPRIME,INTAGS)
      CALL HBOOKN(22,'CHIPStest EffMasses',NCOLEM,'PAWC',NPRIME,EMTAGS)
      CALL HBOOKN(25,'CHIPStest Event',NCOLEV,'PAWC',NPRIME,EVTAGS)
C
      read(44,444) rtemp,rssg2,reep,rcut,rgfc,rsqc,ipdgpr,ipdgtg,
     &       nevnt,np,id1,id2,id3,id4,id5,id6,id7,id8,id9,id10,id11,id12
 444  FORMAT(f6.2,5f5.2,3i6,i3,12i2)
      temp(1) = rtemp
      ssg2(1) = rssg2
      eep(1) = reep
      cut(1) = rcut
      gfc(1) = rgfc
      sqc(1) = rsqc
      vnevt(1) = nevnt
      vnp(1) = np
      d1(1) = id1
      d2(1) = id2
      d3(1) = id3
      d4(1) = id4
      d5(1) = id5
      d6(1) = id6
      d7(1) = id7
      d8(1) = id8
      d9(1) = id9
      d10(1) = id10
      d11(1) = id11
      d12(1) = id12
      close(44)
C
	  do i = 1, MaxEvt
         read(4,*,END=900) ntr
C#         print *, ntr
C
         if(ntr.gt.NTRMAX) then
            print *, ' %%%NTCR: Event number = ',i
            print *, ' %%%NTCR: Too many tracks in the ievent! ntr=',ntr
            print *, ' %%%NTCR: ntr truncated at the maximum = ',NTRMAX
         endif
         Mtot = 0
         Mgam = 0
         Mpip = 0
         Mpim = 0
         Mpi0 = 0
         MKp  = 0
         MK0  = 0
         MKm  = 0
         MaK0 = 0
         Meta = 0
         Metap  = 0
         Mrhom  = 0
         Mrhop  = 0
         Mrho0  = 0
         Momega = 0
         Mphi   = 0
         MKS0   = 0
         MKSC   = 0
         MaKS0  = 0
         MaKSC  = 0
C
         do j = 1, ntr
            jm = min(j,NTRMAX)
            read(4,*,END=900) itr(jm),ipdg(jm),
     &                        vmass(jm),px(jm),py(jm),pz(jm),e(jm)
C#          print *, itr(jm),ipdg(jm),vmass(jm),px(jm),py(jm),pz(jm),e(jm)
            if(ipdg(jm).EQ.221 .OR. ipdg(jm).EQ.111) then
               Mpi0 = Mpi0 + 1
            elseif(ipdg(jm).EQ. 211) then
               Mpip = Mpip + 1
            elseif(ipdg(jm).EQ.-211) then
               Mpim = Mpim + 1
            elseif(ipdg(jm).EQ. 311) then
               MK0  = MK0 + 1
            elseif(ipdg(jm).EQ.-321) then
               MKp  = MKp + 1
            elseif(ipdg(jm).EQ.-311) then
               MaK0 = MaK0 + 1
            elseif(ipdg(jm).EQ. 321) then
               MKm  = MKm + 1
            elseif(ipdg(jm).EQ. 22) then
               Mgam = Mgam + 1
            elseif(ipdg(jm).EQ. 331) then
               if(vmass(jm).lt.700.) then
                  Meta  = Meta + 1
               else
                  Metap = Metap + 1
               endif
            elseif(ipdg(jm).EQ. 113) then
               Mrho0  = Mrho0 + 1
            elseif(ipdg(jm).EQ. 213) then
               Mrhop  = Mrhop + 1
            elseif(ipdg(jm).EQ.-213) then
               Mrhom  = Mrhom + 1
            elseif(ipdg(jm).EQ. 223) then
               Momega  = Momega + 1
            elseif(ipdg(jm).EQ. 333) then
               Mphi  = Mphi + 1
            elseif(ipdg(jm).EQ. 313) then
               MKS0  = MKS0 + 1
            elseif(ipdg(jm).EQ. 323) then
               MKSC  = MKSC + 1
            elseif(ipdg(jm).EQ.-313) then
               MaKS0  = MaKS0 + 1
            elseif(ipdg(jm).EQ.-323) then
               MaKSC  = MaKSC + 1
            endif
         enddo
C
         IF(TRIGTR) THEN
            XTUPEV(1)  = Nevt
            XTUPEV(2)  = ntr
            XTUPEV(3)  = Mgam
            XTUPEV(4)  = Mpip
            XTUPEV(5)  = Mpim
            XTUPEV(6)  = Mpi0
            XTUPEV(7)  = MKp
            XTUPEV(8)  = Mk0
            XTUPEV(9)  = MKm
            XTUPEV(10) = MaK0
            XTUPEV(11) = Meta
            XTUPEV(12) = Metap
            XTUPEV(13) = Mrhom
            XTUPEV(14) = Mrhop
            XTUPEV(15) = Mrho0
            XTUPEV(16) = Momega
            XTUPEV(17) = Mphi
            XTUPEV(18) = MKS0
            XTUPEV(19) = MKSC
            XTUPEV(20) = MaKS0
            XTUPEV(21) = MaKSC
            CALL HFN(25,XTUPEV)
            CALL HFILL(1,1.,1.,1.)
         ENDIF
C
         do j = 1, ntr
            jm = min(j,NTRMAX)
            IF(TRIGEV) THEN
C
               XTUPIN(1)  = Nevt
               XTUPIN(2)  = ntr
               XTUPIN(3)  = Mgam
               XTUPIN(4)  = Mpip
               XTUPIN(5)  = Mpim
               XTUPIN(6)  = Mpi0
               XTUPIN(7)  = MKp
               XTUPIN(8)  = Mk0
               XTUPIN(9)  = MKm
               XTUPIN(10) = MaK0
               XTUPIN(11) = Meta
               XTUPIN(12) = Metap
               XTUPIN(13) = Mrhom
               XTUPIN(14) = Mrhop
               XTUPIN(15) = Mrho0
               XTUPIN(16) = Momega
               XTUPIN(17) = Mphi
               XTUPIN(18) = MKS0
               XTUPIN(19) = MKSC
               XTUPIN(20) = MaKS0
               XTUPIN(21) = MaKSC
               XTUPIN(22) = itr(jm) + 1
               XTUPIN(23) = ipdg(jm)
               XTUPIN(24) = vmass(jm)
               XTUPIN(25) = px(jm)
               XTUPIN(26) = py(jm)
               XTUPIN(27) = pz(jm)
               XTUPIN(28) = e(jm)
               CALL HFN(20,XTUPIN)
            ENDIF
         enddo
C
         do j = 1, ntr-1
            jm = min(j,NTRMAX)
            do l = j+1, ntr
               lm = min(l,NTRMAX)
               IF(TRIGEV) THEN
C
                  rpx = px(jm) + px(lm)
                  rpy = py(jm) + py(lm)
                  rpz = pz(jm) + pz(lm)
                  re  =  e(jm) +  e(lm)
                  rmass = re**2-rpx**2-rpy**2-rpz**2
                  if(rmass.lt.0.) then
                     print *, ' %%%NTCR: effmass^2 negative!'
                     rmass = -rmass
                  endif
                  rmass = sqrt(rmass)
                  XTUPEM(1)  = Nevt
                  XTUPEM(2)  = ntr
                  XTUPEM(3)  = Mgam
                  XTUPEM(4)  = Mpip
                  XTUPEM(5)  = Mpim
                  XTUPEM(6)  = Mpi0
                  XTUPEM(7)  = MKp
                  XTUPEM(8)  = Mk0
                  XTUPEM(9)  = MKm
                  XTUPEM(10) = MaK0
                  XTUPEM(11) = Meta
                  XTUPEM(12) = Metap
                  XTUPEM(13) = Mrhom
                  XTUPEM(14) = Mrhop
                  XTUPEM(15) = Mrho0
                  XTUPEM(16) = Momega
                  XTUPEM(17) = Mphi
                  XTUPEM(18) = MKS0
                  XTUPEM(19) = MKSC
                  XTUPEM(20) = MaKS0
                  XTUPEM(21) = MaKSC
                  XTUPEM(22) = 1000.*(itr(jm)+1) + (itr(lm)+1)
                  XTUPEM(23) = 1000.*ipdg(jm) + ipdg(lm)
                  XTUPEM(24) = rmass
                  XTUPEM(25) = rpx
                  XTUPEM(26) = rpy
                  XTUPEM(27) = rpz
                  XTUPEM(28) = re
                  CALL HFN(22,XTUPEM)
               ENDIF
            enddo
         enddo
C
         Nevt = i
      enddo
C
 900  continue
      print *, ' %%%ntcr: ',Nevt,' events read'
C
      CALL HROUT(0,ICYCLE,' ')
C
      close(4)
C
    END
