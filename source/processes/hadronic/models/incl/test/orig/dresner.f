      subroutine dresner(iia,iiz,exc,multen,multep,irun)
      implicit real*8 (a-h,o-z)
      parameter (lpt=18)                                                
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter  (l01=3000)                                                      
      integer blz                                                               
      common /komon/ e(l01), ec(l01), tip(l01), u(l01), v(l01),                 
     1 w(l01), wt(l01), x(l01), xc(l01), y(l01), yc(l01),                       
     2 z(l01), zc(l01), name(l01), nmed(l01), blz(l01)                          
      common /exms/ exma(6), parthr(6)                                  
      common /ecom/ etep(2), se(5), ses(5), sr(5), srs(5), monin(7)             
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /fong/ ievap, ifiss                                                
      common /smacom/ yzero, bzero, yzere, bzere                                
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /cpreq/ ex0, ptest, eps0, rcomp, ucon, rho0(6),                    
     1 nemax0, kalb, nstage, nemax, ipr, mode, modeb                            
c     double precision lamda0, lamda1, lamda2, lamda3, mu0, mu1, mu2,           
c    1 mu3, nu0, nu1, nu2, lambda, mu, nu, muprim, nuprim                       
      parameter (lring=100, lesrc=200, ltsrc=200)                               
      common /qfbl/ loff                                                        
      common /comon/ apr, cosks, cosphi, costh, d, delsig, dkwt, emax,  
     1 emin(lpt), erec, ex, hevsum, hsig, oldwt, sinks, sinphi, sinth,  
     2 umax, uu, zpr, ibert, ityp, lelem, mat,maxbch, maxcas, mxmat, n, 
     3 nabov, namax, nbelo, nbogus, negex, neutno(4), lneutp, ngroup,   
     4 npidk, no, nobch, nocas, nomax, nopart, npart(6), nquit, nbertp  
      parameter (l12=120)                                                       
      parameter (l13=50)                                                        
      logical fisinh
      double precision kefis
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /hevaps/ epart(l12,2), hepart(l13,4), cosevp(3,l12,2),             
     1 coshev(3,l13,4)                                                          
      data itopt /-1/, twit /dp0/, tpeak /dp0/, tsig /dp0/                      
      common /evcm/ y0, b0, t(4,7)                                              
      common /dresm/ exmass(6), exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      logical penbar
      common /emitr/ sos(6), strun(6), zmass(6), q(6), fla(6),                  
     1 flkcou(6), ccoul(6), thresh(6), smalla(6), flz(6), r(6), s(6),           
     2 eye1(6), smom1(6), eps, eye0, gamma, ar, zr, ex1, ex2,                   
     3 xlhs, xrhs, argp, sigma, npart0(6,2), jemiss, penbar                     
      common /holdfb/ index
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /filex/ cpa, wtneut(4), wtnt1(4), wtnt2(4), sprlm1,                
     1 sprlm2, wtpi, wtmu, pmass(lpt), flim0, fiswt, fiswt1, fiswt2,            
     2 swtm, lhist, lnrec, nlrnge, neutnt(4), nhstin, ir0(2),                   
     3 icpt, nobalc, nobale, ifisct,kpart(lpt), krec(8), isrfsr,                
     4 nofis, noskip, jcasc, ifbrk, ipreq, ilvden, jcoul,                       
     5 kneutp, nwrds, neutct, nlimit, nlbuf, ltme, lneut                        
      common /output/ exi(20), vp(3,20), erc(20), npt, jjp(20)
c     dimension ilev(lchan), ill1(lchan), iip1(lchan), vx1(lchan), vy1
c    1 (lchan), vz1(lchan), kip(lchan)
      common /orcm2/ ainfo(20), epkin(80,6), sigke(80,6), delexc(80,6,3)        
     1 , cstiff(120)                                                            
      common /labcos/ coslbp(3), coslbr(3)                                      
      dimension cfong(50)                                                       
      data cfong /2*65.d0,65.d0,86.d0,124.d0,170.d0,186.d0,190.d0,186.d0        
     1 ,170.d0,140.d0,118.d0,100.d0,90.d0,80.d0,73.d0,66.d0,58.d0,52.d0,        
     2 46.d0,42.d0,39.d0,37.d0,37.d0,38.d0,43.d0,62.d0,104.d0,180.d0,222        
     3 .d0,245.d0,250.d0,236.d0,190.d0,140.d0,106.d0,83.d0,73.d0,65.d0,6        
     4 0.d0,55.d0,52.d0,49.d0,47.d0,45.d0,43.d0,42.d0,3*42.d0/                  
      data ifoz /2/, ifong /0/, igraph/1/, ifinal /0/,                          
     1 iscale /1/, iloren /1/, iwt /0/, enrmsa /dph/, eta /0.4d0/,              
     2 smult /15.d0/, ifa /0/, idoub /0/                               
ccccc
      real*4 theta,en,aapr,eal,zzpr,rnpr,multen,multep,tabener
      common/cv/tabener(0:2000)
      real*4 acv,zcv,ecv,tcv,pcv
      common/volant1/nombre,acv(200),zcv(200),ecv(200),
     s               tcv(200),pcv(200)
     
      
      IF(irun.EQ.-1335) write(6,*)'dans dresner',iia,iiz,exc,irun    
      call dres(iia,iiz,exc,npart,hevsum,uu,erec,apr,zpr,wtfis,nocas
     +,irun)
c     write(6,'("npart")')
C npart(), nombre de particules par categorie:n,p,d,t,3He,4He
C apr,zpr=residuel (note AB)
      IF(irun.EQ.-1335) 
     s    write(6,*)'npart',(npart(iloc),iloc=1,6),apr,zpr
      if (ncbal.gt.10.or.apr.le.dp0) go to 120                                   
      if (apr.lt.1000) then                                                     
        casbal=casbal-dnrgy(apr,zpr)                                           
      else                                                                      
        if (atfis(1).eq.1.d3) atfis(1)=dp0                                      
        if (atfis(1).gt.dp0) casbal=casbal-dnrgy(atfis(1),ztfis(1))            
        if (atfis(2).gt.dp0) casbal=casbal-dnrgy(atfis(2),ztfis(2))            
      endif                                                                     
  120 continue                                                                  
      do 130 i=1,6                                                               
      if (i.gt.2) then                                                          
        lmtxx=l13                                                               
      else                                                                      
        lmtxx=l12                                                               
      endif                                                                     
      if (npart(i).gt.lmtxx) then                                               
        nx=i                                                                    
        write (io,*) nx,lmtxx                                                 
        call errorf ('erup')                                                    
      endif                                                                     
  130 continue
  
C      COSPHI=1. 
                                                                       
      do 170 k=1,6                                                               
      np=npart(k)                                                               
      if (np.eq.0) go to 170   
      if (k.gt.2) go to 150   
      do 140 j=1,np    
c      WRITE(6,*) j,k,(cosevp(iab,j,k),iab=1,3)
c      WRITE(6,*) 'Rot:',costh,sinth,cosphi,sinphi                                      
      if (epart(j,k).gt.dp0) casbal=casbal-epart(j,k)-exma(k)                   
      t1=costh*cosevp(1,j,k)+sinth*cosevp(3,j,k)                                
      u1=cosphi*t1-sinphi*cosevp(2,j,k)                                         
      cosevp(2,j,k)=sinphi*t1+cosphi*cosevp(2,j,k)                              
      cosevp(3,j,k)=costh*cosevp(3,j,k)-sinth*cosevp(1,j,k)                     
      cosevp(1,j,k)=u1
  140 continue                                                                  
      go to 170
  150 continue                                                                  
      l=k-2                                                                     
      if (l.eq.2.and.apr.ge.1.d3) np=np-14                                      
      do 160 j=1,np
c      WRITE(6,*) j,l,(coshev(iab,j,l),iab=1,3)
c      WRITE(6,*) 'Rot:',costh,sinth,cosphi,sinphi                                  
      if (hepart(j,l).gt.dp0) casbal=casbal-hepart(j,l)-exma(k)                 
      t1=costh*coshev(1,j,l)+sinth*coshev(3,j,l)                                
      u1=cosphi*t1-sinphi*coshev(2,j,l)                                         
      coshev(2,j,l)=sinphi*t1+cosphi*coshev(2,j,l)                              
      coshev(3,j,l)=costh*coshev(3,j,l)-sinth*coshev(1,j,l)                     
      coshev(1,j,l)=u1                                                          
  160 continue                                                                  
  170 continue                                                                  
      do 220 j=1,6 
      np=npart(j)                                                               
      if (np) 250,220,190                                                       
  190 if (j-2) 200,200,210                                                      
  200 continue                                                                  
c     write (100,*) sngl(wt(no)),(sngl(epart(k,j)),k=1,np),((sngl(cosevp        
c    1 (i,k,j)),i=1,3),k=1,np)                                                  
c     write (100,*) (sngl(epart(k,j)),k=1,np),((sngl(cosevp        
c    1 (i,k,j)),i=1,3),k=1,np)                                                  
      nwrds=1+nwrds+4*np+1                                                      
      go to 220                                                                 
  210 l=j-2                                                                     
      np1=np                                                                    
      if (apr.ge.1000.and.l.eq.2) np1=np1-14                                    
c     write (100,*) sngl(wt(no)),(sngl(hepart(k,l)),k=1,np),((sngl              
c    1 (coshev(i,k,l)),i=1,3),k=1,np1)                                          
c     write (100,*) (sngl(hepart(k,l)),k=1,np),((sngl              
c    1 (coshev(i,k,l)),i=1,3),k=1,np1)                                          
      nwrds=1+nwrds+4*np+1                                                      
  220 continue                                                                  
      go to 251
  250 call errorf ('analyz.4')
c neutron
251	continue
	nombre=0
      np=npart(1)                                                               
      do 252,k=1,np
      if (sngl(epart(k,1)).gt.0) then
c     write (101,*) (sngl(epart(k,1))),(sngl(cosevp        
c    1 (3,k,1))) 
c      WRITE(6,*) k,(sngl(cosevp(iab,k,1)),iab=1,3)                                       
      en=sngl(epart(k,1))
      theta=acos(sngl(cosevp(3,k,1)))*180/3.1415
      phicv=atan2(sngl(cosevp(2,k,1)),sngl(cosevp(1,k,1)))
     s      *180/3.1415
      nombre=nombre+1
      acv(nombre)=1.
      zcv(nombre)=0.
      ecv(nombre)=en
      tcv(nombre)=theta
      pcv(nombre)=phicv
      IF(irun.EQ.-1335)
     sWRITE(6,*) 'Here 1',k,'cosevp(*,',k,',1)',acv(nombre),nombre
call hfill(11,en,0.,1.)
call hfill(511,en,0.,1.)
      ijk=en
      if(ijk.ge.2000)ijk=2000
      tabener(ijk)=tabener(ijk)+1
c      write(6,*)' evap',ijk,tabener(ijk)        
      endif
  252 continue
c probleme triton fission ou residuel
      kkk1=0
      kkk2=0
      if (npart(4).lt.16) then
c residuel
      aapr=apr
      zzpr=zpr
      rnpr=apr-zpr
call hfill(24,aapr,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
c      WRITE(6,*) 'Here 2'
      nombre=nombre+1
      acv(nombre)=aapr
      zcv(nombre)=zzpr
      ecv(nombre)=erec
      tcv(nombre)=0.
      pcv(nombre)=0.
      IF(irun.EQ.-1335) WRITE(6,*) 'Here 2',acv(nombre),nombre
C angle ?
      endif
      if (npart(4).ge.16) then
c fission ?
      kkk1=2
      kkk2=14
      np=npart(4)-kkk2
      aprr1=sngl(hepart(11+np,2))
c      aprr2=sngl(hepart(12+np,2))
      zprr1=sngl(hepart(13+np,2))
      eal1=sngl(hepart(2+np,2))
      theta1=acos(sngl(coshev(3,2+np,2)))*180/3.1415
      phicv1=atan2(sngl(coshev(2,2+np,2)),sngl(coshev(1,2+np,2)))
     s      *180/3.1415
c      WRITE(6,*) 'Here 3',2+np,'coshev(*,',2+np,',1)'
      nombre=nombre+1
      acv(nombre)=aprr1
      zcv(nombre)=zprr1
      ecv(nombre)=eal1
      tcv(nombre)=theta1
      pcv(nombre)=phicv1
      IF(irun.EQ.-1335) 
     sWRITE(6,*) 'Here 3',2+np,'coshev(*,',2+np,',1)',acv(nombre),nombre
      aprr2=sngl(hepart(12+np,2))
      zprr2=sngl(hepart(14+np,2))
      eal2=sngl(hepart(6+np,2))
      theta2=acos(sngl(coshev(3,6+np,2)))*180/3.1415
      phicv2=atan2(sngl(coshev(2,6+np,2)),sngl(coshev(1,6+np,2)))
     s      *180/3.1415
c      WRITE(6,*) 'Here 4',6+np,'coshev(*,',6+np,',1)'
      nombre=nombre+1
      acv(nombre)=aprr2
      zcv(nombre)=zprr2
      ecv(nombre)=eal2
      tcv(nombre)=theta2
      pcv(nombre)=phicv2
      IF(irun.EQ.-1335) 
     sWRITE(6,*) 'Here 4',6+np,'coshev(*,',6+np,',1)',acv(nombre),nombre
cfb      write(6,'("toto")')
cfb      write(6,*) aprr1,aprr2,zprr1,zprr2
      aapr=aprr1
      zzpr=zprr1
cfb      write(6,*) aapr,zzpr
      rnpr=aapr-zzpr
call hfill(24,aapr,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
      aapr=aprr2
      zzpr=zprr2
      rnpr=aapr-zzpr
call hfill(24,aapr,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
      endif
c deuton
      np=npart(3)
      do 253,k=1,np
      aapr=2.
      zzpr=1.
      rnpr=1.
      if (sngl(hepart(k,1)).gt.0) then
      eal=sngl(hepart(k,1))
call hfill(24,aapr,0.,1.)
call hfill(103,eal,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
c      WRITE(6,*) 'Here 5',k,'coshev(*,',k,',1)'
      nombre=nombre+1
      acv(nombre)=aapr
      zcv(nombre)=zzpr
      ecv(nombre)=eal
      theta=acos(sngl(coshev(3,k,1)))*180/3.1415
      tcv(nombre)=theta      
      phicv=atan2(sngl(coshev(2,k,1)),sngl(coshev(1,k,1)))
     s      *180/3.1415
      pcv(nombre)=phicv
      IF(irun.EQ.-1335)
     sWRITE(6,*) 'Here 5',k,'coshev(*,',k,',1)',acv(nombre),nombre
      endif
  253 continue
c triton
      np=npart(4)-kkk2
      do 254,k=1,np
      aapr=3.
      zzpr=1.
      rnpr=2.
      if (sngl(hepart(k,2)).gt.0) then
      eal=sngl(hepart(k,2))
call hfill(24,aapr,0.,1.)
call hfill(104,eal,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
c      WRITE(6,*) 'Here 6',k,'coshev(*,',k,',1)'
      nombre=nombre+1
      acv(nombre)=aapr
      zcv(nombre)=zzpr
      ecv(nombre)=eal
      theta=acos(sngl(coshev(3,k,2)))*180/3.1415
      tcv(nombre)=theta   
      phicv=atan2(sngl(coshev(2,k,2)),sngl(coshev(1,k,2)))
     s      *180/3.1415
      pcv(nombre)=phicv
      IF(irun.EQ.-1335)
     sWRITE(6,*) 'Here 6',k,'coshev(*,',k,',1)',acv(nombre),nombre
      endif
  254 continue
c 3he
      np=npart(5)
      do 255,k=1,np
      aapr=3.
      zzpr=2.
      rnpr=1.
      if (sngl(hepart(k,3)).gt.0) then
      eal=sngl(hepart(k,3))
call hfill(24,aapr,0.,1.)
call hfill(105,eal,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
c      WRITE(6,*) 'Here 7',k,'coshev(*,',k,',3)'
      nombre=nombre+1
      acv(nombre)=aapr
      zcv(nombre)=zzpr
      ecv(nombre)=eal
      theta=acos(sngl(coshev(3,k,3)))*180/3.1415
      tcv(nombre)=theta   
      phicv=atan2(sngl(coshev(2,k,3)),sngl(coshev(1,k,3)))
     s      *180/3.1415
      pcv(nombre)=phicv
      IF(irun.EQ.-1335)
     sWRITE(6,*) 'Here 7',k,'coshev(*,',k,',3)',acv(nombre),nombre
      endif
  255 continue
c 4he
      np=npart(6)
      do 256,k=1,np
      aapr=4.
      zzpr=2.
      rnpr=2.
      if (sngl(hepart(k,4)).gt.0) then
      eal=sngl(hepart(k,4))
call hfill(24,aapr,0.,1.)
call hfill(106,eal,0.,1.)
call hfill(247,zzpr,0.,1.)
call hfill(248,rnpr,zzpr,1.)
c      WRITE(6,*) 'Here 8',k,'coshev(*,',k,',4)'
      nombre=nombre+1
      acv(nombre)=aapr
      zcv(nombre)=zzpr
      ecv(nombre)=eal
      theta=acos(sngl(coshev(3,k,4)))*180/3.1415
      tcv(nombre)=theta   
      phicv=atan2(sngl(coshev(2,k,4)),sngl(coshev(1,k,4)))
     s      *180/3.1415
      pcv(nombre)=phicv
      IF(irun.EQ.-1335)
     sWRITE(6,*) 'Here 8',k,'coshev(*,',k,',4)',acv(nombre),nombre
      endif
  256 continue
c proton
      np=npart(2)
      do 257,k=1,np
      if (sngl(epart(k,2)).gt.0) then
      eal=sngl(epart(k,2))
call hfill(102,eal,0.,1.)
call hfill(851,eal,0.,1.)
c      WRITE(6,*) 'Here 9',k,'cosevp(*,',k,',2)'
      nombre=nombre+1
      acv(nombre)=1.
      zcv(nombre)=1.
      ecv(nombre)=eal
      theta=acos(sngl(cosevp(3,k,2)))*180/3.1415
      tcv(nombre)=theta         
      phicv=atan2(sngl(cosevp(2,k,2)),sngl(cosevp(1,k,2)))
     s      *180/3.1415
      IF(irun.EQ.-1335)
     sWRITE(6,*) 'Here 9',k,'cosevp(*,',k,',2)',acv(nombre),nombre
      endif
  257 continue
      multen=npart(1)
      multep=npart(2)
  260 rewind(8)
      rewind(9)
      rewind(10)
      end

      subroutine dres (m2,m3,t1,npart,sochpe,u,erec,apr,zpr,wtfis,nocas
     +,irun)        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c     ==================================================================        
c     dres     calls the following subroutines and functions                    
c        DOST       ENERGY     ERROR      EXP        EXPRNF     FSINFO          
c        GAMR       GTISO      QNRG       SQRT                                  
c     ==================================================================        
c***  npart0(j,1) is the particle counter for the j-th type, for any            
c***  evaporation of particles from the original nucleus before fissioni        
c***  a zero value is inserted in energy tables after evap. products            
c***  from original nucleus                                                     
c***  npart0(j,2) is the particle counter for the j-th type for                 
c***  evaporation of the first fragment                                         
c***  a zero value will be inserted between evaporation product                 
c***  energies of first and second fragments                                    
c***  consequently npart(j) will all be a value of 2 greater than true          
c***  value since zero also inserted after original nucleus evaporation         
      logical fising                                                            
      common /inout/ lu1, lu2, itty, iscrt                                      
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      parameter (lpt=18)                                                        
      common /filex/ cpa, wtneut(4), wtnt1(4), wtnt2(4), sprlm1,                
     1 sprlm2, wtpi, wtmu, pmass(lpt), flim0, fiswt, fiswt1, fiswt2,            
     2 swtm, lhist, lnrec, nlrnge, neutnt(4), nhstin, ir0(2),                   
     3 icpt, nobalc, nobale, ifisct,kpart(lpt), krec(8), isrfsr,                
     4 nofis, noskip, jcasc, ifbrk, ipreq, ilvden, jcoul,                       
     5 kneutp, nwrds, neutct, nlimit, nlbuf, ltme, lneut                        
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      parameter (l04=35+lpt)                                                    
      common /comon/ dimdim(l04), nbertp                                        
      equivalence (nbt,nbertp)                                                  
      common /forcn/ fkey                                                       
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      parameter (l12=120)                                                       
      parameter (l13=50)                                                        
      common /hevaps/ epart(l12,2), hepart(l13,4), cosevp(3,l12,2),             
     1 coshev(3,l13,4)                                                          
      common /exms/ exma(6), parthr(6)                                          
      common /dresm/ exmass(6), exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      common /output/ exi(20), vp(3,20), erc(20), npt, jjp(20)                  
      logical fisinh                                                            
      double precision kefis                                                    
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /fong/ ievap, ifiss                                                
      common /smacom/ yzero, bzero, yzere, bzere                                
      common /labcos/ coslbp(3), coslbr(3)                                      
      logical penbar                                                            
      common /emitr/ sos(6), strun(6), zmass(6), q(6), fla(6),                  
     1 flkcou(6), ccoul(6), thresh(6), smalla(6), flz(6), r(6), s(6),           
     2 eye1(6), smom1(6), eps, eye0, gamma, ar, zr, ex1, ex2,                   
     3 xlhs, xrhs, argp, sigma, npart0(6,2), jemiss, penbar                     
      dimension npart(6)                                                        
      data nwr1 /0/, nwr2 /0/, ballim /0.1d0/, um /931.20793d0/                 
      if (parthr(1).le.0.) then                                                
        do 10 i=1,6                                                             
        fla(i)=ia(i)                                                            
        flz(i)=iz(i)                                                            
        zmass(i)=exmass(i)+um*fla(i)                                            
   10   exma(i)=exmass(i)                                                       
        exbe=dnrgy(8.d0,4.d0)-2.d0*exmass(6)                                   
        parthr(1)=1.0d+30                                                       
        parthr(2)=1.0d+30                                                       
        parthr(3)=exmass(1)+exmass(2)-exmass(3)                                 
        parthr(4)=exmass(1)+exmass(3)-exmass(4)                                 
        parthr(5)=exmass(2)+exmass(3)-exmass(5)                                 
        parthr(6)=exmass(2)+exmass(4)-exmass(6)                                 
      endif                                                                     
      fkey=dp0                                                                  
      do 20 i=1,6                                                               
      npart(i)=0                                                                
   20 smom1(i)=dp0                                                              
      ja=m2                                                                     
      jz=m3                                                                     
      u=t1                                                                      
      nfiss=0                                                                   
      fisinh=.false.                                                            
      penbar=.false.                                                            
      fising=.false.                                                            
c.......................................................................        
c///// start of loop /////                                                      
   30 continue                                                                  
      fja=ja                                                                    
      fjz=jz                                                                    
      ex2=dnrgy(fja,fjz)                                                       
      rnmass=ex2+um*ja                                                          
      rnmass=rnmass+u                                                           
   40 continue                                                                  
c.......................................................................        
c///// Fermi breakup routines /////                                             
      if (ja.gt.20.or.(ja.gt.13.and.u.gt.44.0)) go to 140                       
      if (ifbrk.eq.0.and.ja.ge.6) go to 140 
c      WRITE(6,*) 'CALL FERMI BREAKUP'                                    
      call moment (jz,ja,coslbr(1),coslbr(2),coslbr(3),u,erec,ierr)             
      if (ierr) 70,50,140                                                       
   50 continue                                                                  
      u=dp0                                                                     
      erec=dp0                                                                  
      nhev=0                                                                    
      do 60 i=1,npt                                                             
      ipt=0                                                                     
      if (jjp(i).le.5) ipt=jjp(i)                                               
      if (jjp(i).eq.7) ipt=6                                                    
      if (ipt.gt.0) then                                                        
        if (exi(i).ne.dp0) call errorf ('dres.1')                               
        ja=ja-ia(ipt)                                                           
        jz=jz-iz(ipt)                                                           
        npart(ipt)=npart(ipt)+1                                                 
        smom1(ipt)=smom1(ipt)+erc(i)                                            
        if (ipt.le.2) then                                                      
          epart(npart(ipt),ipt)=erc(i)                                          
          cosevp(1,npart(ipt),ipt)=vp(1,i)                                      
          cosevp(2,npart(ipt),ipt)=vp(2,i)                                      
          cosevp(3,npart(ipt),ipt)=vp(3,i)
c	  WRITE(6,*) 'cosevp(*,',npart(ipt),',',ipt,',)',
c     s     (cosevp(iloc,npart(ipt),ipt),iloc=1,3)    
        else                                                                    
          hepart(npart(ipt),ipt-2)=erc(i)                                       
          coshev(1,npart(ipt),ipt-2)=vp(1,i)                                    
          coshev(2,npart(ipt),ipt-2)=vp(2,i)                                    
          coshev(3,npart(ipt),ipt-2)=vp(3,i)                                    
c	  WRITE(6,*) 'coshev(*,',npart(ipt),',',ipt-2,',)',
c     s     (coshev(iloc,npart(ipt),ipt-2),iloc=1,3)    
        endif                                                                   
      else                                                                      
        nhev=nhev+1                                                             
        if (nhev.gt.1) call errorf ('dres.2')                                   
        u=exi(i)                                                                
        erec=erc(i)                                                             
        coslbr(1)=vp(1,i)                                                       
        coslbr(2)=vp(2,i)                                                       
        coslbr(3)=vp(3,i)                                                       
c	  WRITE(6,*) 'coslbr(*)',
c     s     (coslbr(iloc),iloc=1,3)    
      endif                                                                     
   60 continue                                                                  
      go to 230                                                                 
   70 continue                                                                  
      do 80 j=1,6                                                               
      if (ja.eq.ia(j).and.jz.eq.iz(j)) go to 90                                 
   80 continue                                                                  
      go to 230                                                                 
   90 continue                                                                  
      eps=u+erec                                                                
      npart(j)=npart(j)+1                                                       
      smom1(j)=smom1(j)+eps                                                     
      if (j.gt.2) go to 110                                                     
      epart(npart(j),j)=eps                                                     
      do 100 i=1,3                                                              
  100 cosevp(i,npart(j),j)=coslbr(i)                                            
      go to 130                                                                 
  110 kemiss=j-2                                                                
      hepart(npart(j),kemiss)=eps                                               
      do 120 i=1,3                                                              
  120 coshev(i,npart(j),kemiss)=coslbr(i)                                       
  130 erec=dp0                                                                  
      u=dp0                                                                     
      ja=0                                                                      
      jz=0                                                                      
      go to 230                                                                 
c.......................................................................        
c///// calculate channel probabilities /////                                    
  140 continue                                                                  
      jrn=0                                                                     
      gamma=dp1+erec/rnmass                                                     
      a=ja                                                                      
      z=jz                                                                      
      ex1=ex2                                                                   
      xlhs=u+erec+ex1                                                           
      if (ja.eq.8.and.jz.eq.4) then                                             
        jrn=1                                                                   
        k=6                                                                     
        q(6)=-exbe                                                              
        go to 160                                                               
      endif                                                                     
      call dlect (a,z,ja,jz,u,kfiss)                                            
      if (kfiss.eq.1) go to 240                                                 
      if (sigma.le.dp0) then                                                    
        if (ievap.eq.0) go to 180                                               
        if (ievap.ne.0) go to 200                                               
      endif                                                                     
c.......................................................................        
c///// choose an emitted particle /////                                         
      call nombre_generer(rndm)
      uran=rndm*sigma                                                   
      sum=dp0                                                                   
      do 150 j=1,6                                                              
      k=j                                                                       
c     write(6,'("r(j)")')     
c     write(6,*) r(1),r(2),r(3),r(4),r(5),r(6)
      sum=r(j)+sum                                                              
      if (sum.ge.uran) go to 160                                                
  150 continue                                                                  
      call errorf ('dres.3')                                                    
  160 continue                                                                  
      jemiss=k                                                                  
ccc   modiff  
c     jemiss=6
c     write(6,*) a,z,u,jemiss
      ar=a-fla(jemiss)                                                          
      zr=z-flz(jemiss)                                                          
      if (nobale.eq.0) then                                                     
        ex2=q(jemiss)-exmass(jemiss)+ex1                                        
      else                                                                      
        ex2=dnrgy(ar,zr)                                                       
      endif                                                                     
c.......................................................................        
c///// emit a particle /////                                                    
  170 continue                                                                  
      call emit (a,z,ja,jz,u,erec,kfiss,npart,jrn)                              
      if (kfiss.ne.0) go to 260                                                 
      if (fkey.ne.dp0) fkey=dp0                                                 
      go to 40                                                                  
c.......................................................................        
c///// normal emission blocked - RAL /////                                      
  180 continue                                                                  
c                                                                               
c      check for proton barrier penetration                                     
c      we allow for a barrier transmission probability                          
c      of down to 10**-6.we use a gammow type relation.                         
c                                                                               
      if (argp.lt.dp0) go to 200                                                
      mm=int(a-dp1)                                                             
c**    coulomb barrier                                                          
      vc=1.1076d0*(z-dp1)/rmass(mm)                                             
      eta=argp/vc                                                               
      if (eta.gt.dp1) go to 190                                                 
      arge=0.60072d0*sqrt(z*rmass(mm)/eta)*(acos(eta)-eta*sqrt(dp1-             
     1 eta*eta))                                                                
      if (arge.gt.13.8d0) go to 200                                             
c                                                                               
c      we predict barrier penetration                                           
c      pick ke uniformly in range 0.0 to argp                                   
c                                                                               
  190 continue                                                                  
      call nombre_generer(rndm)
      eps=argp*rndm                                                             
      penbar=.true.                                                             
      jemiss=2                                                                  
      call errorf ('dres.4')                                                    
      go to 170                                                                 
c.......................................................................        
c///// normal emission blocked - ORNL /////                                     
  200 continue                                                                  
c///// switch of FKEY - to be eliminated in future /////                        
      if (fkey.ne.dp1) then                                                     
        fkey=dp1                                                                
        go to 140                                                               
      endif                                                                     
c.......................................................................        
c///// no further evaporation /////                                             
      if (ja.gt.0) go to 230                                                    
      if (ja.eq.5) go to 40                                                     
      if (ifbrk.ne.0) call errorf ('dres.5')                                    
      do 210 j=1,6                                                              
      if (jz.eq.iz(j).and.ja.eq.ia(j)) go to 220                                
  210 continue                                                                  
      go to 40                                                                  
c.......................................................................        
c///// store - residual nuc is of emitted particle type /////                   
  220 continue                                                                  
      if (u.gt.parthr(j)) go to 40                                              
      jemiss=j                                                                  
      eps=u+erec                                                                
      npart(jemiss)=npart(jemiss)+1                                             
      smom1(jemiss)=smom1(jemiss)+eps                                           
      ja=0                                                                      
      jz=0                                                                      
      erec=dp0                                                                  
      u=dp0                                                                     
      if (jemiss.le.2) then                                                     
        epart(npart(jemiss),jemiss)=eps                                         
        cosevp(1,npart(jemiss),jemiss)=coslbr(1)                                
        cosevp(2,npart(jemiss),jemiss)=coslbr(2)                                
        cosevp(3,npart(jemiss),jemiss)=coslbr(3)                                
      else                                                                      
        kemiss=jemiss-2                                                         
        hepart(npart(jemiss),kemiss)=eps                                        
        coshev(1,npart(jemiss),kemiss)=coslbr(1)                                
        coshev(2,npart(jemiss),kemiss)=coslbr(2)                                
        coshev(3,npart(jemiss),kemiss)=coslbr(3)                                
      endif                                                                     
c.......................................................................        
  230 continue                                                                  
      if (ja.gt.0.and.(jz.le.1.or.(ja-jz).le.1.or.ja.le.5)) write (1,         
     1 360) nocas,ja,jz,ex1,u,q(1),q(2),fkey                                    
      if (fisinh) then                                                          
c.......................................................................        
c///// fission termination /////                                                
        ztfis(nfiss)=real(jz)                                                   
        atfis(nfiss)=real(ja)                                                   
        recfis(nfiss)=erec                                                      
        utfis(nfiss)=u                                                          
        if (fising) go to 280                                                   
        hevs2=smom1(3)+smom1(5)+smom1(6)+smom1(4)                               
        ifisct=ifisct+1                                                         
        fiswt=fiswt+wtfis                                                       
        if (atfis(1).eq.dp0) atfis(1)=1.d3                                      
        apr=atfis(1)*1000+atfis(2)                                              
        zpr=ztfis(1)*1000+ztfis(2)                                              
        index=npart(4)+1                                                        
        hepart(index,2)=dp0                                                     
        hepart(index+1,2)=recfis(1)                                             
        hepart(index+2,2)=ufis(1)                                               
        hepart(index+3,2)=hevs1-hevs0                                           
        hepart(index+4,2)=utfis(1)                                              
        hepart(index+5,2)=recfis(2)                                             
        hepart(index+6,2)=ufis(2)                                               
        hepart(index+7,2)=hevs2-hevs1                                           
        hepart(index+8,2)=utfis(2)                                              
        hepart(index+9,2)=hevs0                                                 
        hepart(index+10,2)=atfis(1)                                             
        hepart(index+11,2)=atfis(2)                                             
        hepart(index+12,2)=ztfis(1)                                             
        hepart(index+13,2)=ztfis(2)                                             
        npart(4)=npart(4)+14                                                    
        erec=recfis(1)+recfis(2)                                                
        u=utfis(1)+utfis(2)                                                     
      else                                                                      
c.......................................................................        
c///// non-fission termination /////                                            
        apr=real(ja)                                                            
        zpr=real(jz)                                                            
      endif                                                                     
c.......................................................................        
      sochpe=smom1(3)+smom1(4)+smom1(5)+smom1(6)                                
      return                                                                    
c                                                                               
c***********************************************************************        
c                                                                               
c          fission                                                              
c                                                                               
c***********************************************************************        
c                                                                               
c.......................................................................        
c///// ORNL fission routine /////                                               
  240 continue                                                                  
      kfiss=0                                                                   
      af=a                                                                      
      zf=z                                                                      
      ef=u                                                                      
      rfmass=ex1                                                                
      ernf=erec                                                                 
      do 250 i=1,3                                                              
      coslf0(i)=coslbr(i)                                                       
  250 continue                                                                  
      fising=.true.                                                             
      call fsinfo                                                               
      if (ncall.eq.1) then                                                      
c---- fission failed--------                                                    
        fising=.false.                                                          
        fisinh=.false.                                                          
        nfiss=0                                                                 
        kfiss=1                                                                 
        go to 140                                                               
      endif                                                                     
      xrhs=r1mass+ufis(1)+kefis(1)+r2mass+ufis(2)+kefis(2)                      
      if (abs(xlhs-xrhs).ge.ballim) then                                        
        nwr1=nwr1+1                                                             
        if (nwr1.le.10) write (1,370) xlhs,xrhs,ja,jz,jemiss                  
      endif                                                                     
      go to 280                                                                 
c.......................................................................        
c///// RAL fission routine /////                                                
  260 continue                                                                  
      do 270 i=1,3                                                              
      coslf0(i)=coslbr(i)                                                       
  270 continue                                                                  
      fising=.true.                                                             
      call fissed (ja,jz,u,erec)                                                
      if (.not.fisinh) then                                                     
c---- fission failed--------                                                    
        fising=.false.                                                          
        nfiss=0                                                                 
        penbar=.false.                                                          
        go to 30                                                                
      endif                                                                     
c.......................................................................        
c///// common fission code /////                                                
  280 continue                                                                  
c******* pick up the appropriate fragment                                       
      nfiss=nfiss+1                                                             
      ja=afis(nfiss)                                                            
      jz=zfis(nfiss)                                                            
      u=ufis(nfiss)                                                             
      erec=kefis(nfiss)                                                         
      fkey=dp0                                                                  
      penbar=.false.                                                            
      if (nfiss.eq.2) then                                                      
        fising=.false.                                                          
        go to 310                                                               
      endif                                                                     
      hevs0=smom1(3)+smom1(5)+smom1(6)+smom1(4)                                 
      do 290 i=1,6                                                              
      numpar=npart(i)+1                                                         
      npart(i)=numpar                                                           
  290 continue                                                                  
      do 300 j=1,3                                                              
  300 coslbr(j)=coslf1(j)                                                       
      go to 340                                                                 
  310 continue                                                                  
      hevs1=smom1(3)+smom1(5)+smom1(6)+smom1(4)                                 
      do 320 i=1,6                                                              
      numpar=npart(i)+1                                                         
      npart(i)=numpar                                                           
  320 continue                                                                  
      do 330 j=1,3                                                              
  330 coslbr(j)=coslf2(j)                                                       
  340 continue                                                                  
      do 350 i=1,6                                                              
      numpar=npart(i)                                                           
      if (i.lt.3) epart(numpar,i)=dp0                                           
      if (i.gt.2) hepart(numpar,i-2)=dp0                                        
      do 350 j=1,3                                                              
      if (i.lt.3) cosevp(j,numpar,i)=dp0                                        
      if (i.gt.2) coshev(j,numpar,i-2)=dp0                                      
  350 continue                                                                  
      go to 30                                                                  
c.......................................................................        
c                                                                               
  360 format (/5x,'nocas=',i7,' a=',i2,' z=',i2,' exm=',g11.4,' u=',g11.        
     1 4,' q1(n)=',g11.4,' q(p)=',g11.4,' fkey=',f5.1)                          
  370 format (' ***fission xlhs=',g10.3,' xrhs=',g10.3,' ja=',i3,' jz='         
     1 ,i3,' jemiss=',i2)                                                       
      end                                                                       
      subroutine redmas (ixedt,l1,l2)                                           
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  REDMAS controls all the initiation phase for the Fermi breakup model.        
c  It read the data library, FBRDATA, and establishes the mass data             
c  arrays. It calls MASEDT to initialize the branching ratio arrays.            
c  Using IXEDT = 0 eliminates the long printout of the data.                    
c                                                                               
c    REFERENCES.......  FBRK                                                    
c    EXTERNALS........  EXITA     MASEDT    BRDTAB                              
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      dimension istat(lstat)                                                    
      data rho0 /1.414d0/                                                       
      data ierr /0/                                                             
c     L1 = output file; L2 = data file                                          
      lout=l1                                                                   
      lmass=l2                                                                  
      itop=17                                                                   
      atop=itop                                                                 
      unstab=dp1                                                                
      loff=0                                                                    
      do 10 i=1,11                                                              
      do 10 j=1,13                                                              
   10 ip(i,j)=0                                                                 
      do 20 i=1,lstat                                                           
      wtl(i)=dp1                                                                
   20 continue                                                                  
      read (lmass,*) nr,nr1,nr2                                                 
      nrm=lpart                                                                 
c                                                                               
c       read particle stable nuclei                                             
c                                                                               
      i=0                                                                       
      do 40 i1=1,nr                                                             
      read (lmass,*) zin,ain,xmin,nsin                                          
      if (ain.gt.atop) then                                                     
        call fixm (zin,ain,xmin)                                                
        read (lmass) idmm                                                       
        go to 40                                                                
      endif                                                                     
      i=i+1                                                                     
      a(i)=ain                                                                  
      z(i)=zin                                                                  
      xm(i)=xmin                                                                
      ns(i)=nsin                                                                
      jz(i)=z(i)                                                                
      ja(i)=a(i)                                                                
      axm(i)=xm(i)+a(i)*931.478d0                                               
      xminv(i)=dp1/axm(i)                                                       
      iz1=z(i)+1.01d0                                                           
      in1=a(i)-z(i)+1.01d0                                                      
      ip(iz1,in1)=i                                                             
      if (i.ge.8.or.i.eq.6) call fixm (zin,ain,xmin)                            
      nss=ns(i)                                                                 
      jns(i)=loff                                                               
      tmin=dph*(ja(i)-2*jz(i))                                                  
      tmin=abs(tmin)                                                            
      lo=loff+1                                                                 
      loff=loff+nss                                                             
      if (loff.gt.lstat) go to 110                                              
      read (lmass,*) (e(j),istat(j),tiso(j),j=lo,loff)                          
      do 30 j=lo,loff                                                           
      tiso(j)=max(tiso(j),tmin)                                                 
   30 continue                                                                  
   40 continue                                                                  
      nr=i                                                                      
c     write (lout,180) nr,loff                                                  
      do 50 i=1,loff                                                            
      parity(i)=dp1                                                             
      if (istat(i).lt.0) parity(i)=-parity(i)                                   
      istat(i)=iabs(istat(i))                                                   
      aj(i)=dph*(istat(i)-1)                                                    
   50 continue                                                                  
      rho(1)=dp0                                                                
      do 60 i=2,nr                                                              
      if (a(i).gt.a(i-1)) go to 60                                              
      if (a(i).eq.a(i-1).and.z(i).gt.z(i-1)) go to 60                           
c     write (lout,150) a(i),z(i)                                                
      ierr=ierr+1                                                               
   60 continue                                                                  
      if (ierr.ne.0) go to 120                                                  
c                                                                               
c       read particle unstable res.nuc.                                         
c                                                                               
      i3=nr+1                                                                   
      i4=nr+nr1                                                                 
      i=0                                                                       
      do 70 j=i3,i4                                                             
      read (lmass,*) zin,ain,xmin                                               
      call fixm (zin,ain,xmin)                                                  
      if (ain.gt.atop) go to 70                                                 
      i=i+1                                                                     
      i2=i+nr                                                                   
      a(i2)=ain                                                                 
      z(i2)=zin                                                                 
      xm(i2)=xmin                                                               
      jz(i2)=z(i2)                                                              
      ns(i2)=0                                                                  
      ja(i2)=a(i2)                                                              
      axm(i2)=xm(i2)+a(i2)*931.478d0                                            
      iz1=z(i2)+dp1                                                             
      in1=a(i2)-z(i2)+dp1                                                       
      ip(iz1,in1)=i2                                                            
   70 continue                                                                  
      nr1=i                                                                     
c     write (lout,190) nr1                                                      
      i4=nr+nr1                                                                 
      do 80 i=1,i4                                                              
      ia=ja(i)                                                                  
      if (ia.lt.10) then                                                        
        if (ia.eq.9) rho(i)=3.25d0                                              
        if (ia.eq.8) rho(i)=2.83d0                                              
        if (ia.eq.7) rho(i)=2.42d0                                              
        if (ia.eq.6) rho(i)=2.02d0                                              
        if (ia.eq.5) rho(i)=2.02d0                                              
        if (ia.le.4) rho(i)=1.2d0                                               
      else                                                                      
        rho(i)=rho0*a(i)**dpth+dp1                                              
      endif                                                                     
   80 continue                                                                  
      call masedt (ixedt)                                                       
      n=0                                                                       
      do 90 i=1,lmas                                                            
      if (i.gt.nr) go to 90                                                     
      if (nklev(jns(i)+1).gt.0) go to 90                                        
      n=n+1                                                                     
   90 continue                                                                  
      if (n.gt.lmas1) call errorf ('redmas.1')                                  
      do 100 i=1,nr2                                                            
      read (lmass,*) zin,ain,xmin                                               
      call fixm (zin,ain,xmin)                                                  
  100 continue                                                                  
      return                                                                    
  110 continue                                                                  
c     write (lout,160)                                                          
      go to 130                                                                 
  120 continue                                                                  
c     write (lout,170) ierr                                                     
  130 call errorf ('redmas.2')                                                  
c                                                                               
  140 format (/3x,' no. of particle stable nuclei =',i2/' no. of particl        
     1e unstable res.nuc.=',i2/' in FBRDATA library')                           
  150 format (/' stable mass a=',i3,' z=',i3,' is out of order')                
  160 format (' array length LSTAT must be increased')                          
  170 format (2x,i4,' errors in reading FBRDATA library')                       
  180 format (//' new nr =',i4,' required array length lstat=',i6)              
  190 format (//' new nr1=',i4)                                                 
      end                                                                       
      subroutine fixm (z,a,xm)                                                  
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      dimension jprime(250)                                                     
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      data jprime /19*0,10*1,8*2,7*3,7*4,6*5,5*6,6*7,5*8,5*9,5*10,5*11,4        
     1 *12,5*13,4*14,4*15,4*16,4*17,4*18,4*19,4*20,4*21,4*22,3*23,4*24,4        
     2 *25,3*26,4*27,3*28,4*29,3*30,3*31,4*32,3*33,3*34,4*35,3*36,3*37,3        
     3 *38,3*39,3*40,3*41,3*42,3*43,3*44,3*45,3*46,3*47,3*48,3*49,3*50,3        
     4 *51,3*52,3*53,3*54,55,55,3*56,3*57,3*58,2*59,2*60/                       
      data delamu /0.296064d0/                                                  
      ia=a                                                                      
      iz=z                                                                      
      j=ia-2*iz-jprime(ia)+10                                                   
      if (j.ge.1.and.j.le.20) waps(ia,j)=xm+delamu*ia                           
      return                                                                    
      end                                                                       
      function amass (iz,ia,iok)                                                
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  AMASS returns the total mass as the function value for arguements            
c  IZ and IA.  IOK is the index of the total mass in the AXM array and          
c  of the mass excess in the XM array. IOK is returned as 0 when the            
c  mass is not in the library (in which case the breakup routines do            
c  not apply to the mass.)                                                      
c                                                                               
c    REFERENCES.......  LEVBRK    CBREAK    FBRK                                
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      iz1=iz+1                                                                  
      if (iz1.le.0.or.iz1.gt.11) go to 10                                       
      in1=ia-iz+1                                                               
      if (in1.le.0.or.in1.gt.13) go to 10                                       
      ipt=ip(iz1,in1)                                                           
      if (ipt.le.0.or.ipt.gt.lmas) go to 10                                     
      iok=ipt                                                                   
      amass=axm(ipt)                                                            
      return                                                                    
   10 iok=0                                                                     
      return                                                                    
      end                                                                       
      function igetd (k)                                                        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (lpt=18)                                                        
      common /comon/ apr, cosks, cosphi, costh, d, delsig, dkwt, emax,          
     1 emin(lpt), erec, ex, hevsum, hsig, oldwt, sinks, sinphi, sinth,          
     2 umax, uu, zpr, ibert, ityp, lelem, mat,maxbch, maxcas, mxmat, n,         
     3 nabov, namax, nbelo, nbogus, negex, neutno(4), lneutp, ngroup,           
     4 npidk, no, nobch, nocas, nomax, nopart, npart(6), nquit, nbertp          
      igetd=-777777                                                             
      if (k.eq.1) igetd=nocas                                                   
      return                                                                    
      end                                                                       
      function fltrnf (x)                                                       
      implicit doubleprecision(a-h,o-z)                                         


      parameter (lpt=18)
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c       return the pseudo-random number 0 < r < +1                             
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      common /comon/ apr, cosks, cosphi, costh, d, delsig, dkwt, emax,          
     1 emin(lpt), erec, ex, hevsum, hsig, oldwt, sinks, sinphi, sinth,          
     2 umax, uu, zpr, ibert, ityp, lelem, mat,maxbch, maxcas, mxmat, n,         
     3 nabov, namax, nbelo, nbogus, negex, neutno(4), lneutp, ngroup,           
     4 npidk, no, nobch, nocas, nomax, nopart, npart(6), nquit, nbertp          
      parameter (p=2.**24, q=2.**(-24), r=2.**(-48), gb=1136868, gs=6328        
     1 637)                                                                     
c                                                                               
c     a=gs*rans                                                                 
c        this expression for b is valid only if gb+gs<2**24                     
c     b=gb*rans+gs*ranb+aint(a*q)                                               
c     rans=a-aint(a*q)*p                                                        
c     ranb=b-aint(b*q)*p                                                        
c     fltrnf=(ranb*p+rans)*r                                                    
c     if (nocas.ge.151) then
c     write(1,*) fltrnf
c     endif
c     read(1,*) fltrnf
c     nrand=nrand+1                                                             
      call nombre_generer(rndm)
      fltrnf=rndm
      return                                                                    
      end                                                                       
      subroutine azirn (sin,cos)                                                
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
   10 continue                                                                  
      r1=sflraf(dumm)                                                           
      r1sq=r1*r1                                                                
      r2=fltrnf(dumm)                                                           
      r2sq=r2*r2                                                                
      rsq=r1sq+r2sq                                                             
      if (dp1-rsq) 10,20,20                                                     
   20 continue                                                                  
      rsq=dph*rsq                                                               
      cos=(rsq-r1sq)/rsq                                                        
      sin=r1*r2/rsq                                                             
      return                                                                    
      end                                                                       
      block data bd00                                                           
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (lpt=18)                                                        
      parameter (lmt=50)                                                        
      parameter (lmt1=lmt+1)                                                    
      parameter (lel=50)                                                        
      common /ioutsd/ iouts                                                     
      character*8 mpart                                                         
      common /aas/ ncolc(13)                                                    
      common /aasch/ mpart(lpt)                                                 
      dimension ncolx(8)                                                        
      equivalence (ncolx(1),ncolc(6))                                           
      common /chin/ charge(lpt), rpart(lpt), radfm(260),                        
     1 rcoul(260), instab(lpt)                                                  
      common /filex/ cpa, wtneut(4), wtnt1(4), wtnt2(4), sprlm1,                
     1 sprlm2, wtpi, wtmu, pmass(lpt), flim0, fiswt, fiswt1, fiswt2,            
     2 swtm, lhist, lnrec, nlrnge, neutnt(4), nhstin, ir0(2),                   
     3 icpt, nobalc, nobale, ifisct,kpart(lpt), krec(8), isrfsr,                
     4 nofis, noskip, jcasc, ifbrk, ipreq, ilvden, jcoul,                       
     5 kneutp, nwrds, neutct, nlimit, nlbuf, ltme, lneut                        
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      common /inout/ in, io, itty, iscrt                                        
      parameter (l3=lel*lmt)                                                    
      common /az/ exlow(5), lowaz, negaz, nrnorm(4)                             
      common /rmsx/ rms(lpt), rmsinv(lpt)                                       
      common /fong/ ievap, ifiss                                                
      common /smacom/ yzero, bzero, yzere, bzere                                
      common /pseudo/ nseudo                                                    
      common /spred/ nspred, nwsprd                                             
      common /n1cl/ n1col                                                       
      common /histp/ nhstp                                                      
      common /comon/ apr, cosks, cosphi, costh, d, delsig, dkwt, emax,          
     1 emin(lpt), erec, ex, hevsum, hsig, oldwt, sinks, sinphi, sinth,          
     2 umax, uu, zpr, ibert, ityp, lelem, mat,maxbch, maxcas, mxmat, n,         
     3 nabov, namax, nbelo, nbogus, negex, neutno(4), lneutp, ngroup,           
     4 npidk, no, nobch, nocas, nomax, nopart, npart(6), nquit, nbertp          
      parameter (lels=121,npoits=20*(lels-1))                                   
      common /elstic/ elas, fone, totels, sgels(lel), sige(npoits),             
     1 ef(npoits), es(npoits), f1(npoits), locf1(lels), locsig(lels),           
     2 noel(lmt), id(lel,lmt), nelstp, nledit, noelas                           
      common /hie/ ehin, ehipi, ihie, ncmosc                                    
      common /joint/ ibertp                                                     
      logical bfhist                                                            
      common /buffer/ bfhist, mxhist, iqhist                                    
      common /elo/ elow(2), rngelo(lpt,lmt), sigelo(lpt,lmt)                    
      common /forisa/ ichoic, iexisa                                            
      parameter (liso=60, liso5=liso-5)                                         
      common /proba/ zrev(liso), arev(liso), rtarg(liso,2),                     
     1 sigtrg(6,liso), sigis(lel,lmt), sigisa(lel,lmt), nainp,                  
     2 iazdex(lel,lmt)                                                          
      common /dkc/ dkcons(lpt)                                                  
      common /dresm/ exmass(6), exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /inpu/ andit, ctofe, ctofen                                        
c record 4                                                                      
      data maxcas /100/                                                         
      data maxbch /1/                                                           
      data mxmat /1/                                                            
      data noelas /0/                                                           
      data nlimit /400000/                                                      
      data neutct /0/                                                           
      data lhist /23068672/                                                     
      data nlbuf /0/                                                            
      data jcasc /1/                                                            
      data iexisa /0/                                                           
      data ichoic /23/                                                          
      data ifbrk /1/                                                            
      data ipreq /0/                                                            
      data ilvden /0/                                                           
      data jcoul /1/                                                            
c record 5                                                                      
      data ibertp /1/                                                           
      data nelstp /0/                                                           
      data lneutp /1/                                                           
      data nhstp /1/                                                            
      data mxhist /27/                                                          
      data nwsprd /0/                                                           
      data nofis /1/                                                            
      data noskip /10/                                                          
c record 5a                                                                     
      data krec /4*0,-1,2*0,1/                                                  
c record 5b                                                                     
      data kpart /lpt*0/                                                        
c record 6                                                                      
      data nspred /1/                                                           
      data nbogus /1/                                                           
      data npidk /0/                                                            
      data n1col /0/                                                            
      data ltme /0/                                                             
      data nledit /0/                                                           
      data nlrnge /0/                                                           
      data isrfsr /0/                                                           
      data icpt /0/                                                             
c record 7                                                                      
      data emax /3495.d0/                                                       
      data emin /dp1,2.d1,.14875d0,1.0d+30,.14875d0,2*.11261d0,4*dp1,           
     1 1.0d+30,.52614d0,2*1.0d-6,.52614d0,dp1,1.0d-6/                           
      data ctofe /-1.d0/                                                        
      data elas /50.d0/                                                         
      data totels /dp0/                                                         
      data id /l3*0/                                                            
      data wtmu /0.1d0/                                                         
      data wtpi /0.1d0/                                                         
      data flim0 /-1.d0/                                                        
      data swtm /2.5d-1/                                                        
c     the use of swtm=0.25 will prevent error messages in hmcnp                 
c define unit numbers                                                           
      data in /1/                                                               
      data io /2/                                                               
      data iscrt /3/                                                            
      data nhstin /9/                                                           
c initialize random number generator                                            
      data nrand /0/                                                            
c other data statements                                                         
      data wtneut /4*dp0/                                                       
      data neutnt /4*0/                                                         
      data wtnt2 /4*dp0/                                                        
      data wtnt1 /4*dp0/                                                        
      data ifisct, fiswt, fiswt1, fiswt2 /0,3*dp0/                              
      data sprlm1 /0.1d0/, sprlm2 /dp1/                                         
      data exlow /dp0,1.0d+30,dp0,-1.0d+30,dp0/, nrnorm /4*0/                   
      data negaz /0/                                                            
      data lowaz /0/, neutno /4*0/, nobch /0/, negex /0/, nocas /1/,            
     1 nbelo /0/, nabov /0/                                                     
      data rms /dp1,dp0,6.7228d0,dp0,6.7228d0,8.8803d0,8.8803d0,                
     1 .50025d0,.33397d0,1.32723d0,dp1,dp0,1.9006d0,2*dp0,1.9006d0,dp1,         
     2 dp0/                                                                     
      data rmsinv /dp1,dp0,.14875d0,dp0,.14875d0,.11261d0,.11261d0,             
     1 1.999d0,2.9943d0,.75345d0,dp1,dp0,.52614d0,2*dp0,.52614d0,dp1,           
     2 dp0/                                                                     
      data charge /dp1,dp0,dp1,dp0,-1.d0,dp1,-1.d0,2*dp1,dp2,dp2,dp0,           
     1 dp1,dp0,dp0,-1.d0,-1.d0,dp0/                                             
      data instab /0,0,-1,1,-1,1,1,0,0,0,0,0,1,1,1,1,0,0/                       
      data rpart /lpt*dp0/                                                      
      data pmass /938.280d0,939.573d0,139.567d0,134.96d0,139.567d0,             
     1 2*105.66d0,1875.627d0,2808.951d0,2808.421d0,3727.418d0,dp0,              
     2 493.67d0,2*497.67d0,493.67d0,938.280d0,939.573d0/                        
c     data exmass /8.3651d0,7.5831d0,13.7222d0,15.8383d0,15.8193d0,             
c    1 3.6084d0/                                                                
      data dkcons /dp0,dp0,0.001281d0,4.0d+05,0.001281d0,1.518d-5,1.518d        
     1 -5,5*dp0,0.002696d0,0.0006435d0,0.3738d0,0.002696d0,2*dp0/               
      data iouts /0/                                                            
      data mpart /'proton','neutron','pi+','pi0','pi-','mu+','mu-','h-2'        
     1 ,'h-3','he-3','he-4','photon','K+','K0long','K0shrt','K-','p-bar'        
     2 ,'n-bar'/                                                                
      data ncolc /13*0/                                                         
      data ehin /3495.d0/                                                       
      data ehipi /2495.d0/                                                      
      data arev /dp1,dp2,2*dp3,dp4,liso5*dp0/                                   
      data zrev /3*dp1,dp2,dp2,liso5*dp0/                                       
      data nainp /5/                                                            
      end                                                                       
      block data dbd01                                                           
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (lpt=18)                                                        
      logical isz, isn                                                          
      parameter (inn=150, izz=98)                                               
      common /dcook/ sz(izz), sn(inn), con(2), amean(240), pz(izz),              
     1 pn(inn), isz(izz), isn(inn)                                              
      parameter (l01=3000)                                                      
      integer blz                                                               
      common /komon/ e(l01), ec(l01), tip(l01), u(l01), v(l01),                 
     1 w(l01), wt(l01), x(l01), xc(l01), y(l01), yc(l01),                       
     2 z(l01), zc(l01), name(l01), nmed(l01), blz(l01)                          
      parameter (l02=150)                                                       
      parameter (l05=60)                                                        
      parameter (lmt=50)                                                        
      parameter (lmt1=lmt+1)                                                    
      parameter (lel=50)                                                        
      common /comona/ ea(l02), eb(l02), tipa(l02), tipb(l02), ua(l02),          
     1 ub(l02), va(l02), vb(l02), wa(l02), wb(l02), wta(l02),                   
     2 wtb(l02), arg(lmt), denh(lmt), hsigg(lpt,lmt), sigmx(lpt,lmt1),          
     3 a(lel,lmt), den(lel,lmt), eion(lel,lmt), sigg(6,lel,lmt),                
     4 zz(lel,lmt), alpha(l05), beta(l05), gam(l05),ep(l05),                    
     5 wtfas(l05), namea(l02), nel(lmt), kind(l05)                              
      common /save1/ wts, es, xs, ys, zs, nsav                                  
      common /cwrtm/ icount(2)                                                  
      parameter (l1=lpt*lmt, l2=lpt*lmt1, l3=lel*lmt, l4=6*l3)                  
      data icount /2*0/                                                         
      data nsav /0/                                                             
      data a, zz, sigg /l3*dp0,l3*dp0,l4*dp0/                                   
      data nel /lmt*0/                                                          
      data hsigg, sigmx /l1*dp0,l2*dp0/                                         
      data blz /l01*0/                                                          
c                                                                               
c     data tables of dcook et. al.  aaec/tm392, supplimented by g+c.             
c                                                                               
      data sz /8*dp0,-.11d0,-.81d0,-2.91d0,-4.17d0,-5.72d0,-7.8d0,-8.97d        
     1 0,-9.7d0,-10.1d0,-10.7d0,-11.38d0,-12.07d0,-12.55d0,-13.24d0,-13.        
     2 93d0,-14.71d0,-15.53d0,-16.37d0,-17.36d0,-18.6d0,-18.7d0,-18.01d0        
     3 ,-17.87d0,-17.08d0,-16.6d0,-16.75d0,-16.5d0,-16.35d0,-16.22d0,-16        
     4 .41d0,-16.89d0,-16.43d0,-16.68d0,-16.73d0,-17.45d0,-17.29d0,-17.4        
     5 4d0,-17.82d0,-18.62d0,-18.27d0,-19.39d0,-19.91d0,-19.14d0,-18.26d        
     6 0,-17.4d0,-16.42d0,-15.77d0,-14.37d0,-13.91d0,-13.1d0,-13.11d0,-1        
     7 1.43d0,-10.89d0,-10.75d0,-10.62d0,-10.41d0,-10.21d0,-9.85d0,-9.47        
     8 d0,-9.03d0,-8.61d0,-8.13d0,-7.46d0,-7.48d0,-7.2d0,-7.13d0,-7.06d0        
     9 ,-6.78d0,-6.64d0,-6.64d0,-7.68d0,-7.89d0,-8.41d0,-8.49d0,-7.88d0,        
     $ -6.3d0,-5.47d0,-4.78d0,-4.37d0,-4.17d0,-4.13d0,-4.32d0,-4.55d0,-5        
     $ .04d0,-5.28d0,-6.06d0,-6.28d0,-6.87d0,-7.20d0,-7.74d0/                   
      data (sn(i),i=1,110) /8*dp0,10.3d0,5.66d0,6.8d0,7.53d0,7.55d0,7.21        
     1 d0,7.44d0,8.07d0,8.94d0,9.81d0,10.6d0,11.39d0,12.54d0,13.68d0,14.        
     2 34d0,14.19d0,13.83d0,13.5d0,13.d0,12.13d0,12.6d0,13.26d0,14.13d0,        
     3 14.92d0,15.52d0,16.38d0,17.16d0,17.55d0,18.03d0,17.59d0,19.03d0,1        
     4 8.71d0,18.8d0,18.99d0,18.46d0,18.25d0,17.76d0,17.38d0,16.72d0,15.        
     5 62d0,14.38d0,12.88d0,13.23d0,13.81d0,14.9d0,14.86d0,15.76d0,16.2d        
     6 0,17.62d0,17.73d0,18.16d0,18.67d0,19.69d0,19.51d0,20.17d0,19.48d0        
     7 ,19.98d0,19.83d0,20.2d0,19.72d0,19.87d0,19.24d0,18.44d0,17.61d0,1        
     8 7.1d0,16.16d0,15.9d0,15.33d0,14.76d0,13.54d0,12.63d0,10.65d0,10.1        
     9 d0,8.89d0,10.25d0,9.79d0,11.39d0,11.72d0,12.43d0,12.96d0,13.43d0,        
     $ 13.37d0,12.96d0,12.11d0,11.92d0,11.d0,10.8d0,10.42d0,10.39d0,9.69        
     $ d0,9.27d0,8.93d0,8.57d0,8.02d0,7.59d0,7.33d0,7.23d0,7.05d0,7.42d0        
     $ ,6.75d0,6.6d0,6.38d0/                                                    
      data (sn(i),i=111,150) /6.36d0,6.49d0,6.25d0,5.85d0,5.48d0,4.53d0,        
     1 4.3d0,3.39d0,2.35d0,1.66d0,.81d0,0.46d0,-.96d0,-1.69d0,-2.53d0,-3        
     2 .16d0,-1.87d0,-.41d0,.71d0,1.66d0,2.62d0,3.22d0,3.76d0,4.1d0,4.46        
     3 d0,4.83d0,5.09d0,5.18d0,5.17d0,5.1d0,5.01d0,4.97d0,5.09d0,5.03d0,        
     4 4.93d0,5.28d0,5.49d0,5.50d0,5.37d0,5.30d0/                               
      data pz /dp0,5.44d0,dp0,2.76d0,dp0,3.34d0,dp0,2.7d0,dp0,2.5d0,dp0,        
     1 2.46d0,dp0,2.09d0,dp0,1.62d0,dp0,1.62d0,dp0,1.83d0,dp0,1.73d0,dp0        
     2 ,1.35d0,dp0,1.54d0,dp0,1.28d0,0.26d0,0.88d0,0.19d0,1.35d0,-.05d0,        
     3 1.52d0,-.09d0,1.17d0,.04d0,1.24d0,0.29d0,1.09d0,.26d0,1.17d0,.23d        
     4 0,1.15d0,-.08d0,1.35d0,0.34d0,1.05d0,.28d0,1.27d0,dp0,1.05d0,dp0,        
     5 1.d0,.09d0,1.2d0,.2d0,1.4d0,.93d0,1.d0,-.2d0,1.19d0,.09d0,.97d0          
     6 ,dp0,.92d0,.11d0,.68d0,.05d0,.68d0,-.22d0,.79d0,.09d0,.69d0,.01d0        
     7 ,.72d0,dp0,.4d0,.16d0,.73d0,dp0,.46d0,.17d0,.89d0,dp0,.79d0,dp0,.        
     8 89d0,dp0,.81d0,-.06d0,.69d0,-.2d0,.71d0,-.12d0,.72d0,dp0,.77d0/          
      data (pn(i),i=1,125) /dp0,5.98d0,dp0,2.77d0,dp0,3.16d0,dp0,3.01d0         
     1 ,dp0,2.5d0,dp0,2.67d0,dp0,1.8d0,dp0,1.67d0,dp0,1.86d0,dp0,2.04d0         
     2 ,dp0,1.64d0,dp0,1.44d0,dp0,1.54d0,dp0,1.3d0,dp0,1.27d0,dp0,1.29d0        
     3 ,.08d0,1.41d0,-.08d0,1.5d0,-.05d0,2.24d0,-.47d0,1.43d0,-.15d0,1.4        
     4 4d0,.06d0,1.56d0,.25d0,1.57d0,-.16d0,1.46d0,dp0,.93d0,.01d0,.62d0        
     5 ,-.5d0,1.42d0,.13d0,1.52d0,-.65d0,.8d0,-.08d0,1.29d0,-.47d0,1.25d        
     6 0,-.44d0,.97d0,.08d0,1.65d0,-.11d0,1.26d0,-.46d0,1.06d0,0.22d0,1.        
     7 55d0,-.07d0,1.37d0,0.1d0,1.2d0,-.27d0,.92d0,-.35d0,1.19d0,dp0,1.0        
     8 5d0,-.25d0,1.61d0,-.21d0,.9d0,-.21d0,.74d0,-.38d0,.72d0,-.34d0,.9        
     9 2d0,-.26d0,.94d0,.01d0,.65d0,-.36d0,.83d0,.11d0,.67d0,.05d0,1.d0,        
     $ .51d0,1.04d0,.33d0,.68d0,-.27d0,.81d0,.09d0,.75d0,.17d0,.86d0,.14        
     $ d0,1.1d0,-.22d0,.84d0,-.47d0,.48d0,.02d0,.88d0,.24d0,.52d0,.27d0,        
     $ .41d0,-.05/                                                              
      data (pn(i),i=126,150) /0.38d0,.15d0,.67d0,dp0,.61d0,dp0,.78d0,dp0        
     1 ,.67d0,dp0,.67d0,dp0,.79d0,dp0,.6d0,.04d0,.64d0,-.06d0,.45d0,.05d        
     2 0,.26d0,-.22d0,.39d0,dp0,.39d0/                                          
c                                                                               
c Julich mean value level density parameters                                    
c                                                                               
      data (amean(i),i=1,100) /.125d0,.25d0,.375d0,.5d0,.625d0,.75d0,.87        
     1 5d0,1.d0,1.125d0,1.25d0,1.375d0,1.5d0,1.625d0,1.75d0,1.875d0,2.d0        
     2 ,2.125d0,2.25d0,2.375d0,3.94d0,2.63d0,2.75d0,2.88d0,3.55d0,4.35d0        
     3 ,3.25d0,3.38d0,3.96d0,3.63d0,3.75d0,3.88d0,4.82d0,4.44d0,4.43d0,4        
     4 .43d0,4.42d0,4.63d0,5.66d0,5.81d0,5.95d0,5.49d0,6.18d0,7.11d0,6.9        
     5 6d0,7.2d0,7.73d0,6.41d0,6.85d0,6.77d0,6.91d0,7.26d0,7.2d0,6.86d0,        
     6 8.06d0,7.81d0,7.82d0,8.41d0,8.13d0,7.19d0,8.35d0,8.13d0,8.02d0,8.        
     7 93d0,8.9d0,9.69d0,9.65d0,10.55d0,9.38d0,9.72d0,10.66d0,11.98d0,12        
     8 .76d0,12.1d0,12.86d0,13.03d0,12.81d0,12.54d0,12.65d0,12.d0,12.69d        
     9 0,14.05d0,13.33d0,13.28d0,13.23d0,13.17d0,8.66d0,11.09d0,10.4d0,1        
     $ 3.47d0,10.17d0,12.22d0,11.62d0,12.95d0,13.15d0,13.57d0,12.87d0,16        
     $ .16d0,14.71d0,15.69d0,14.09d0/                                           
      data (amean(i),i=101,200) /18.56d0,16.22d0,16.67d0,17.13d0,17.d0,1        
     1 6.86d0,15.33d0,15.61d0,16.77d0,17.93d0,17.45d0,16.97d0,17.88d0,17        
     2 .58d0,15.78d0,16.83d0,17.49d0,16.03d0,15.08d0,16.74d0,17.74d0,17.        
     3 43d0,18.14d0,17.06d0,19.01d0,17.02d0,17.02d0,17.02d0,18.51d0,17.2        
     4 d0,16.75d0,16.97d0,16.94d0,16.91d0,17.69d0,15.55d0,14.56d0,14.35d        
     5 0,16.55d0,18.29d0,17.8d0,17.05d0,21.31d0,19.15d0,19.51d0,19.87d0,        
     6 20.39d0,20.9d0,21.85d0,22.89d0,25.68d0,24.64d0,24.91d0,23.24d0,22        
     7 .85d0,22.46d0,21.98d0,21.64d0,21.75d0,21.85d0,21.77d0,21.69d0,23.        
     8 74d0,21.35d0,23.03d0,20.66d0,21.81d0,20.77d0,22.18d0,22.58d0,22.5        
     9 5d0,21.45d0,21.16d0,21.02d0,20.87d0,22.09d0,22.d0,21.28d0,23.05d0        
     $ ,21.70d0,21.45d0,22.28d0,23.d0,22.11d0,23.56d0,22.83d0,24.88d0,22        
     $ .64d0,23.27d0,23.89d0,23.92d0,23.94d0,21.16d0,22.3d0,21.75d0,21.1        
     $ 9d0,20.72d0,20.24d0,21.34d0,19.d0/                                       
      data (amean(i),i=201,240) /17.93d0,17.85d0,15.7d0,13.54d0,11.78d0,        
     1 10.02d0,10.98d0,10.28d0,11.72d0,13.81d0,14.46d0,15.3d0,16.15d0,16        
     2 .99d0,17.84d0,18.68d0,19.53d0,20.37d0,21.22d0,22.06d0,22.91d0,23.        
     3 75d0,24.6d0,25.44d0,26.29d0,27.13d0,27.98d0,28.82d0,29.67d0,30.71        
     4 d0,30.53d0,31.45d0,29.63d0,30.15d0,30.65d0,30.27d0,29.52d0,30.08d        
     5 0,29.8d0,29.87/                                                          
c                                                                               
c flags for deformed nuclei                                                     
c                                                                               
      data isz /53*.false.,25*.true.,7*.false.,13*.true./                       
      data isn /85*.false.,37*.true.,7*.false.,21*.true./                       
      data con /.142d0,.12d0/                                                   
c                                                                               
      end                                                                       
      block data bd02                                                           
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /bertl/ izero(2509)                                                
      common /blnk/ sf, andit, amasno, zee, einc, ctofe, casesn                 
     1 , prtin, trsym, pnms, sqnm, dncms, rcpmv,poms,ipec(13),isw(13)        
     2 , iout(6), nor, inpt, no, nmas, in, iv, ip, med, inc, ifca, ifcc,        
     3 not, it, ifc, knot, i1, i2, i6, ik, i5, i4, i3, ka, itote, itoti,        
     4 itot2, itot3, nwds, nrt, ln                                              
      common/sblnk/ pt(48), pnbc(5), space(178), hvn(3), hvp(3), awd(3),        
     1 fvnp(3), vnvp(3), pmac(3), ppan(3), thpn(3), ffptfn(3), tffn(3),         
     2 tffp(3), cfepn(6), fmpn(6), fmax(7), crdt(25), s(37), xi(3), dcos        
     3 (3), d(6), ce(21), wkrpn(6), pm(4), e(4), pxyz(12), c(3), eco(2),        
     4 col(24), cc(12), pnidk(23), out(40), fcn, fcp, pgcnt, pacnt,             
     5 pecnt, value2, pppda, ppmda, ppnda, ppnna, clcfe, value1, ans,           
     6 begru, frand, ex, sign, clsm, efrp, efrn, aabs, strkp, rlke,             
     7 p1oe1, polc, pols, sopc, sops, p2, any, sn, absec, com, snt, cst,        
     8 com2, value3, univ, unive, univer, com1, ftr, com4, com3, a,             
     9 ctofen                                                                   
      common /blnk05/ curr(11)                                                  
      common /blnkl/ esps(481), plvc(961), pgvc(441)                            
      common /hold1/ hold                                                       
      common /holdfb/ index                                                     
      common /exms/ exma(6), parthr(6)                                          
      common /alnfct/ f(0:40)                                                   
      common /masval/ xmasss(3,3)                                               
      data xmasss /-1.000000d+00,9.382800d+02,1.876560d+03,9.395730d+02,        
     1 1.875627d+03,2.808421d+03,1.879146d+03,2.808951d+03,3.727418d+03/        
      data f /41*dp0/                                                           
      data parthr /6*-1.0d+30/                                                  
      data hold /dp0/                                                           
      data begru /dp0/, pacnt, pgcnt, pecnt /3*0.0/                             
      data izero /2509*0/                                                       
      end                                                                       
      block data bd03                                                           
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /ecom/ etep(2), se(5), ses(5), sr(5), srs(5), monin(7)             
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /fong/ ievap, ifiss                                                
      common /smacom/ yzero, bzero, yzere, bzere                                
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      common /cpreq/ ex0, ptest, eps0, rcomp, ucon, rho0(6),                    
     1 nemax0, kalb, nstage, nemax, ipr, mode, modeb                            
      double precision lamda0, lamda1, lamda2, lamda3, mu0, mu1, mu2,           
     1 mu3, nu0, nu1, nu2, lambda, mu, nu, muprim, nuprim                       
      common /cmgcom/ p0(6), p1(6), p2(6), lamda0(6), lamda1(6), lamda2(        
     1 6), lamda3(6), mu0(6), mu1(6), mu2(6), mu3(6), nu0(6), nu1(6),           
     2 nu2(6), itt(6), ematch, lambda, mu, nu, muprim, nuprim, ecoul,           
     3 ezero, acrit, rad0, rcompa, rcompb, deltab, rfix1, rfix2, delta(6        
     4 ), p, q, r, emx                                                          
      common /emita/ zmass(6), fla(6), flz(6), rr(6), s(6), eps, gamma          
      parameter (lring=100, lesrc=200, ltsrc=200)                               
      common /srscom/ tip0, x0, y0, z0, a0, b0, cuts, smu, e0, rsmu,            
     1 xcut, twit, tpeak, tsig, w0, sring1(lring), sring2(lring),               
     2 sprob(lring), esrc(lesrc), eprb(lesrc), tsrc(ltsrc), tprb(ltsrc)         
     3 , isopt, neprb, itopt, ietab, nring, ntprb, jbr, itip0                   
      common /qfbl/ loff                                                        
      data sring1, sring2, sprob /lring*dp0,lring*dp0,lring*dp0/                
      data esrc, eprb /lesrc*dp0,lesrc*dp0/                                     
      data tsrc, tprb /ltsrc*dp0,ltsrc*dp0/                                     
      data tip0, x0, y0, z0, a0, b0, cuts, smu /8*dp0/, isopt /0/, e0 /8        
     1 .d2/                                                                     
      data itopt /-1/, twit /dp0/, tpeak /dp0/, tsig /dp0/                      
      data rho0 /dp0,dp0,3*.8d0,1.2d0/, rcomp /1.23d0/, mode /1/, modeb         
     1 /0/, nemax0 /0/, ptest /1.d-3/, itt /4,3,3,2,2,1/, kalb /0/,             
     2 delta /2*dp0,3*0.8d0,1.2d0/                                              
      data zmass /9.395730d+02,9.387910d+02,1.876138d+03,2.809462d+03,2.        
     1 809443d+03,3.728440d+03/, fla /dp1,dp1,dp2,dp3,dp3,dp4/, flz /dp0        
     2 ,3*dp1,2*dp2/, emx /2.d2/                                                
      data nklev /lstat*0/                                                      
c     data z91 /91.d0/, ucut /dp4/, ievap /0/, qfactr /dp10/, ismal /2/,        
c    1 iversn /1/, yzero /1.5d0/, bzero /dp10/, yzere /1.5d0/, bzere /8         
c    2 .d0/, ro /1.2d0/, kemax /6/                                              
c     data monin, se, ses, sr, srs /7*0,20*dp0/                                 
      end                                                                       
      function betad (n)                                                        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  BETAD samples the beta distribution given by                                 
c           p(v)=v**dph*(1-v)**(3*n/2-4)/beta(3/2,3*n/2-3)                      
c  for n .ge. 3.                                                                
c  For use in HETC, only the n=3 case is retained                               
c                                                                               
c    REFERENCES.......  MOMENT                                                  
c    EXTERNALS........  RANF      ALOG      MOD                                 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      common /hold1/ hold                                                       
      if (hold.ne.dp0) then                                                     
        tau2=hold                                                               
        hold=dp0                                                                
      else                                                                      
   10   continue                                                                
        r1=sflraf(dumm)                                                         
        sum=r1**2                                                               
        r2=sflraf(dumm)                                                         
        sum=sum+r2**2                                                           
        if (sum.ge.1) go to 10                                                  
        fac=-log(sum)/sum                                                       
        tau2=fac*r1**2                                                          
        hold=fac*r2**2                                                          
      endif                                                                     
      x=tau2-log(fltrnf(dumm))                                                  
      if (hold.ne.0) then                                                       
        tau2=hold                                                               
        hold=dp0                                                                
      else                                                                      
   20   continue                                                                
        r1=sflraf(dumm)                                                         
        sum=r1**2                                                               
        r2=sflraf(dumm)                                                         
        sum=sum+r2**2                                                           
        if (sum.ge.1) go to 20                                                  
        fac=-log(sum)/sum                                                       
        tau2=fac*r1**2                                                          
        hold=fac*r2**2                                                          
      endif                                                                     
      y=tau2-log(fltrnf(dumm))                                                  
      betad=x/(x+y)                                                             
      return                                                                    
      end                                                                       
      subroutine cbreak (ipt,ex,n,kip,ilev,t,ki)                                
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c///////////////////////////////////////////////////////////////////////        
c  CBREAK implements the energy-dependent factors in the breakup channel        
c  probabilites for an arbitrary excited nucleus and selects a breakup          
c  channel.                                                                     
c///////////////////////////////////////////////////////////////////////        
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
      dimension kip(lchan), ilev(lchan)                                         
      kk=ipdat(ipt)                                                             
      nct=npdat(ipt)                                                            
      if (nct.eq.0) go to 80                                                    
      if (nct.eq.1) then                                                        
        if ((ex+q0(kk)).le.dp0) go to 80                                        
        itr(1)=kk                                                               
        kk=1                                                                    
        ki=1                                                                    
        kj=2                                                                    
        wtt=dp1                                                                 
        wff(1)=dp1                                                              
        go to 60                                                                
      endif                                                                     
      kk=kk-1                                                                   
      k=kk+nct                                                                  
      kj=1                                                                      
      wtt=dp0                                                                   
   10 continue                                                                  
      kk=kk+1                                                                   
      if (kk.gt.k) go to 40                                                     
      t0=ex+q0(kk)                                                              
      if (t0.le.dp0) go to 10                                                   
      ijpp=mod(izaps(kk),1000000)/10000                                         
      n0k=2                                                                     
      if (ijpp.gt.0) n0k=3                                                      
      if (n0k.gt.2) texp=1.5d0*n0k-2.5d0                                        
      itr(kj)=kk                                                                
      w00=w0(kk)                                                                
      coulx=coulch(kk)                                                          
      if (coulx.eq.dp0) go to 30                                                
c     ICMB=5                                                                    
      if (n0k.gt.2) go to 20                                                    
      couly=coulck(kk)                                                          
      tem2=sqrt(t0)                                                             
      eta=coulx/tem2                                                            
      rhox=couly*tem2                                                           
      w00=w00*max(1.0d-10,dclmb1(rhox,eta,ml))                                   
      go to 30                                                                  
   20 continue                                                                  
c     ICMB=3                                                                    
      tem2=coulx/t0                                                             
      if (tem2.gt.100.d0) tem2=100.d0                                           
      w00=w00*exp(-tem2)                                                        
   30 continue                                                                  
      if (n0k.eq.2) wtt=wtt+w00*sqrt(t0)                                        
      if (n0k.ne.2) wtt=wtt+w00*exp(texp*log(t0))                               
      wff(kj)=wtt                                                               
      kj=kj+1                                                                   
      go to 10                                                                  
   40 continue                                                                  
      if (kj.eq.1) go to 80                                                     
      fran=wtt*fltrnf(dumm)                                                     
c                                                                               
c     binary search                                                             
c                                                                               
      ki=kj-1                                                                   
      lo=1                                                                      
      lhi=ki                                                                    
      if (fran.le.wff(1)) lhi=1                                                 
   50 continue                                                                  
      if ((lhi-lo).gt.1) then                                                   
        ll=(lhi+lo)/2                                                           
        if (fran.gt.wff(ll)) then                                               
          lo=ll                                                                 
        else                                                                    
          lhi=ll                                                                
        endif                                                                   
        go to 50                                                                
      endif                                                                     
      kk=lhi                                                                    
c                                                                               
c     decay channel selected                                                    
c                                                                               
   60 continue                                                                  
      i=itr(kk)                                                                 
      q=q0(i)                                                                   
      t=q+ex                                                                    
      iizap=izaps(i)                                                            
      jjix=iizap/1000000                                                        
      jj=jjix/3+1                                                               
      ix=mod(jjix,3)+1                                                          
      kip(3)=mod(iizap,1000000)/10000                                           
      n=2                                                                       
      if (kip(3).gt.0) n=3                                                      
      kip(2)=mod(iizap,10000)/100                                               
      kip(1)=mod(iizap,100)                                                     
      do 70 j=1,n                                                               
      ilev(j)=1                                                                 
   70 continue                                                                  
      ilev(ix)=jj                                                               
      iok1=1                                                                    
      return                                                                    
c                                                                               
c       kj=1; res.nuc. does not have enough ex.en. to breakup                   
c                                                                               
   80 continue                                                                  
      n=1                                                                       
      ki=0                                                                      
      kip(1)=ipt                                                                
      ilev(1)=-1                                                                
      t=dp0                                                                     
      iok1=1                                                                    
      return                                                                    
      end                                                                       
      subroutine cbrset (izr,iar,ex,nmax)                                       
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  CBRSET calculates all the excitation energy-independent factors that         
c  enter the breakup channel probabilities for an arbitrary excited             
c  state.                                                                       
c                                                                               
c    REFERENCES.......  FBRK                                                    
c    ENTRY POINTS.....  CBRAK1                                                  
c    EXTERNALS........  AMASS     CH2       CH3       CH4       CH5             
c                       CH6       CH7       EXITA     SQRT      CLMB            
c                       EXP       ALOG      RANF                                
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
c      dimension kip(lchan), s(lchan), nss(lchan), inss(lchan), ilev             
c     1 (lchan), jlev(lchan), jpp(lchan), e0(lchan), fact(lchan), npp  
      dimension  s(lchan), nss(lchan), inss(lchan),              
     1jlev(lchan), jpp(lchan), e0(lchan), fact(lchan), npp            
     2 (lchan), ict(lchan), ef(lchan), wtl1(lchan)                              
      common /holdfb/ index                                                     
      data fact /dp1,dp2,6.d0/                                                  
      data rho0 /1.414d0/                                                       
      data hbc3 /7.683545d+06/, con1 /0.72d0/, con2 /.00516d0/, con3            
     1 /.007168d0/                                                              
      amr=amass(izr,iar,iok)                                                    
      if (iok.ne.0) go to 10                                                    
      iok1=0                                                                    
      return                                                                    
   10 continue                                                                  
      iok1=1                                                                    
      ipt=iok                                                                   
      esmax=1.001d0*ex                                                          
      fzr=izr                                                                   
      k=index-1                                                                 
      nct=0                                                                     
      ioff=0                                                                    
      fac=v(iar)/hbc3                                                           
c                                                                               
c       first loop over no. of particles in breakup                             
c                                                                               
      nop=nmax                                                                  
      w2=dp1                                                                    
      nop=min(iar,nmax)                                                         
      do 210 n=2,nop                                                            
      w2=w2*fac                                                                 
      w235=dph*w2*gam(n)                                                        
c                                                                               
c       find all channels with thes number of particles                         
c                                                                               
      if (n-3) 30,40,20                                                         
   20 call errorf ('cbrset.1')                                                  
   30 call ch2 (nchan)                                                          
C assign modified for gfortran : AB 5/2009
C      assign 110 to jbr2 
       jbr2=110                                                       
      go to 50                                                                  
   40 call ch3 (nchan)                                                          
c      assign 90 to jbr2                                                         
       jbr2=90                                                       
   50 continue                                                                  
      if (nchan.le.0) go to 210                                                 
c                                                                               
c       loop over channels                                                      
c                                                                               
      do 200 i=1,nchan                                                          
c                                                                               
c       loop over no. of particles in channel                                   
c                                                                               
      q=amr                                                                     
      sum=dp0                                                                   
      prod=dp1                                                                  
      do 60 np=1,n                                                              
      iip=izap(np,i)                                                            
      jpp(np)=iip                                                               
      nss(np)=ns(iip)                                                           
      inss(np)=jns(iip)                                                         
      amm=axm(iip)                                                              
      q=q-amm                                                                   
      sum=sum+amm                                                               
      prod=prod*amm                                                             
   60 continue                                                                  
      tem4=esmax+q                                                              
      if (tem4.le.dp0) go to 200                                                
      r=prod/sum                                                                
      facm=r*sqrt(r)                                                            
      ww=w235*facm                                                              
c                                                                               
      coulx=dp0                                                                 
      if (jpp(n-1).eq.1) go to 80                                               
      if (n.eq.2) then                                                          
        iip1=jpp(1)                                                             
        iip2=jpp(2)                                                             
        tem1=dp1/sqrt(xminv(iip1)+xminv(iip2))                                  
        coulx=con2*z(iip1)*z(iip2)*tem1                                         
        couly=con3*(rho(iip1)+rho(iip2))*tem1                                   
      else                                                                      
        tem=dp0                                                                 
        do 70 np=1,n                                                            
        iip=jpp(np)                                                             
        if (iip.eq.1) go to 70                                                  
        iz1=izr-jz(iip)+1                                                       
        if (iz1.le.0) go to 70                                                  
        in1=iar-ja(iip)-iz1+2                                                   
        iip1=ip(iz1,in1)                                                        
        if (iip1.le.0) then                                                     
          aa=(iz1+in1-2)                                                        
          rho2=rho(iip)+rho0*aa**dpth+dp1                                       
        else                                                                    
          rho2=rho(iip)+rho(iip1)                                               
        endif                                                                   
        couly=dp0                                                               
        zz=z(iip)                                                               
        tem=tem+zz*(fzr-zz)/rho2                                                
   70   continue                                                                
        coulx=tem*con1                                                          
      endif                                                                     
   80 continue                                                                  
c                                                                               
c     loop over excited states for a given breakup_                             
c                                                                               
      ns1=nss(1)                                                                
      lof1=inss(1)                                                              
      lof2=inss(2)                                                              
      do 190 j1=1,ns1                                                           
      if (jpp(1).eq.jpp(2).and.j1.gt.1) go to 190                               
      jlev(1)=j1                                                                
      e0(1)=e(lof1+j1)                                                          
      ajj=aj(lof1+j1)                                                           
      s(1)=dp2*ajj+dp1                                                          
      wtl1(1)=wtl(lof1+j1)                                                      
      lo2=1                                                                     
      ns2=1                                                                     
      if (j1.eq.1) ns2=nss(2)                                                   
      do 180 j2=lo2,ns2                                                         
      e0(2)=e(lof2+j2)                                                          
      ef(2)=e0(1)+e0(2)                                                         
      if (ef(2).ge.tem4) go to 180                                              
      jlev(2)=j2                                                                
      ajj=aj(lof2+j2)                                                           
      s(2)=s(1)*(ajj+ajj+dp1)                                                   
      wtl1(2)=wtl(lof2+j2)                                                      
C GO TO modified for gfortran : AB 5/2009
C      go to jbr2, (110,90)
      IF(jbr2.EQ.110) GO TO 110
      IF(jbr2.EQ.90) GO TO 90                                                      
   90 continue                                                                  
      j3=1                                                                      
      if (j1.eq.1.and.j2.eq.1) j3=nss(3)                                        
      if (jpp(3).eq.jpp(2).and.j2.gt.1) go to 180                               
      lo3=1                                                                     
      lof3=inss(3)                                                              
  100 continue                                                                  
      e0(3)=e(lof3+j3)                                                          
      ef(3)=ef(2)+e0(3)                                                         
      if (ef(3).ge.tem4) go to 170                                              
      jlev(3)=j3                                                                
      ajj=aj(lof3+j3)                                                           
      s(3)=s(2)*(ajj+ajj+dp1)                                                   
      wtl1(3)=wtl(lof3+j3)                                                      
  110 continue                                                                  
      tem3=q-ef(n)                                                              
      k=k+1                                                                     
      if (k.gt.mch1.or.nct.gt.mch6) call errorf ('cbrset.2')                    
      nct=nct+1                                                                 
      q0(k)=tem3                                                                
      gg=dp1                                                                    
      if (n.eq.2) then                                                          
        jj=jlev(1)+jlev(2)-2                                                    
        ix=2                                                                    
        if (jlev(1).gt.1) ix=1                                                  
        jjix=3*jj+ix-1                                                          
        izaps(k)=jpp(1)+100*(jpp(2)+10000*jjix)                                 
        wtl0=wtl1(1)*wtl1(2)                                                    
        if (jpp(1).eq.jpp(2).and.e0(1).eq.e0(2)) gg=dp2                         
        go to 160                                                               
      endif                                                                     
      jj=jlev(1)+jlev(2)+jlev(3)-3                                              
      ix=3                                                                      
      if (jlev(1).gt.1) ix=1                                                    
      if (jlev(2).gt.1) ix=2                                                    
      jjix=3*jj+ix-1                                                            
      izaps(k)=jpp(1)+100*(jpp(2)+100*(jpp(3)+100*jjix))                        
      nkind=0                                                                   
      wtl0=dp1                                                                  
      do 120 iii=1,n                                                            
      ict(iii)=jpp(iii)                                                         
      wtl0=wtl0*wtl1(iii)                                                       
  120 continue                                                                  
      do 140 iii=1,n                                                            
c                                                                               
c     check that this particle has not already been counted                     
c                                                                               
      if (ict(iii).eq.0) go to 140                                              
      nkind=nkind+1                                                             
      npp(nkind)=1                                                              
      if (iii.eq.n) go to 140                                                   
      et=e0(iii)                                                                
      ica=ict(iii)                                                              
c                                                                               
c     now scan the rest of the particles for a match                            
c                                                                               
      i1=iii+1                                                                  
      do 130 ii=i1,n                                                            
      if (ict(ii).gt.ica) go to 140                                             
      if (ict(ii).eq.ica.and.e0(ii).eq.et) then                                 
        npp(nkind)=npp(nkind)+1                                                 
        ict(ii)=0                                                               
      endif                                                                     
  130 continue                                                                  
  140 continue                                                                  
c                                                                               
c       compute g factor                                                        
c                                                                               
      if (n-nkind.gt.1) then                                                    
        do 150 iii=1,nkind                                                      
        npi=npp(iii)                                                            
  150   gg=gg*fact(npi)                                                         
      elseif (n-nkind.eq.1) then                                                
        gg=dp2                                                                  
      endif                                                                     
  160 continue                                                                  
      ioff=ioff+n                                                               
      w0(k)=ww*s(n)/gg                                                          
      w0(k)=w0(k)*wtl0                                                          
      coulch(k)=coulx                                                           
      coulck(k)=couly                                                           
      if (n.eq.2) go to 180                                                     
  170 j3=j3-1                                                                   
      if (j3.ge.lo3) go to 100                                                  
  180 continue                                                                  
  190 continue                                                                  
  200 continue                                                                  
  210 continue                                                                  
      ipdat(ipt)=index                                                          
      npdat(ipt)=nct                                                            
      index=index+nct                                                           
      return                                                                    
      end                                                                       
      subroutine ch2 (l)                                                        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  CH2 determines all particle combinations allowed for 2-body breakup          
c  subject to the condition that the first particle must be one                 
c  of the first LPART in the list IALL.                                         
c                                                                               
c    REFERENCES.......  LEVBRK    CBREAK                                        
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
c  allowed light particle breakup products                                      
      dimension iall(10)                                                        
      data iall /1,2,3,4,5,6,7,8,9,10/                                          
c               n,H1,H2,H3,He3,H4,He4,Li4,He5,Li5                               
      data n /2/, nm1 /1/, fn /dp2/                                             
      i=0                                                                       
      fizr=irz                                                                  
      ztest=max(fizr,dp1)                                                       
      fiar=ira                                                                  
      do 40 j1=1,nrm                                                            
      i1=iall(j1)                                                               
      if (fn*a(i1).gt.fiar) go to 50                                            
      if (z(i1).ge.ztest) go to 40                                              
      fa(1)=a(i1)                                                               
      fla(1)=fa(1)                                                              
      fz(1)=z(i1)                                                               
      flz(1)=fz(1)                                                              
c                                                                               
c       check to see if this is a valid breakup mode                            
c                                                                               
      fz(n)=fizr-flz(nm1)                                                       
      fa(n)=fiar-fla(nm1)                                                       
      do 20 i2=i1,nr                                                            
      if (a(i2)-fa(n)) 20,10,40                                                 
   10 continue                                                                  
      if (z(i2).eq.fz(n)) go to 30                                              
   20 continue                                                                  
      go to 40                                                                  
   30 i=i+1                                                                     
      izap(1,i)=i1                                                              
      izap(2,i)=i2                                                              
   40 continue                                                                  
   50 continue                                                                  
      l=i                                                                       
      return                                                                    
      end                                                                       
      subroutine ch3 (l)                                                        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  CH3 determines all particle combinations allowed for 3-body breakup          
c  subject to the condition that the 1st and 2nd particle must be               
c  one of the first LPART in the list IALL.                                     
c                                                                               
c    REFERENCES.......  LEVBRK    CBREAK                                        
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
c  allowed light particle breakup products                                      
      dimension iall(10)                                                        
      data iall /1,2,3,4,5,6,7,8,9,10/                                          
c               n,H1,H2,H3,He3,H4,He4,Li4,He5,Li5                               
      data n /3/, nm1 /2/, fn /dp3/                                             
      i=0                                                                       
      fizr=irz                                                                  
      ztest=max(fizr,dp1)                                                       
      fiar=ira                                                                  
      do 50 j1=1,nrm                                                            
      i1=iall(j1)                                                               
      if (fn*a(i1).gt.fiar) go to 60                                            
      if (z(i1).ge.ztest) go to 50                                              
      fa(1)=a(i1)                                                               
      fla(1)=fa(1)                                                              
      fz(1)=z(i1)                                                               
      flz(1)=fz(1)                                                              
      do 40 j2=j1,nrm                                                           
      i2=iall(j2)                                                               
      if (fla(1)+(fn-dp1)*a(i2).gt.fiar) go to 50                               
      flz(2)=flz(1)+z(i2)                                                       
      if (flz(2).ge.ztest) go to 40                                             
      fa(2)=a(i2)                                                               
      fla(2)=fla(1)+fa(2)                                                       
      fz(2)=z(i2)                                                               
c                                                                               
c       check to see if this is a valid breakup mode                            
c                                                                               
      fz(n)=fizr-flz(nm1)                                                       
      fa(n)=fiar-fla(nm1)                                                       
      do 20 i3=i2,nr                                                            
      if (a(i3)-fa(n)) 20,10,40                                                 
   10 continue                                                                  
      if (z(i3).eq.fz(n)) go to 30                                              
   20 continue                                                                  
      go to 40                                                                  
   30 i=i+1                                                                     
      izap(1,i)=i1                                                              
      izap(2,i)=i2                                                              
      izap(3,i)=i3                                                              
   40 continue                                                                  
   50 continue                                                                  
   60 continue                                                                  
      l=i                                                                       
      return                                                                    
      end                                                                       
      subroutine dlect (a,z,ja,jz,u,kfiss)                                      
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      parameter (lpt=18)                                                        
      common /filex/ cpa, wtneut(4), wtnt1(4), wtnt2(4), sprlm1,                
     1 sprlm2, wtpi, wtmu, pmass(lpt), flim0, fiswt, fiswt1, fiswt2,            
     2 swtm, lhist, lnrec, nlrnge, neutnt(4), nhstin, ir0(2),                   
     3 icpt, nobalc, nobale, ifisct,kpart(lpt), krec(8), isrfsr,                
     4 nofis, noskip, jcasc, ifbrk, ipreq, ilvden, jcoul,                       
     5 kneutp, nwrds, neutct, nlimit, nlbuf, ltme, lneut                        
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      common /forcn/ fkey                                                       
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      common /dresm/ exmass(6), exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      logical fisinh                                                            
      double precision kefis                                                    
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /fong/ ievap, ifiss                                                
      logical penbar                                                            
      common /emitr/ sos(6), strun(6), zmass(6), q(6), fla(6),                  
     1 flkcou(6), ccoul(6), thresh(6), smalla(6), flz(6), r(6), s(6),           
     2 eye1(6), smom1(6), eps, eye0, gamma, ar, zr, ex1, ex2,                   
     3 xlhs, xrhs, argp, sigma, npart0(6,2), jemiss, penbar                     
      data nwr1 /0/, nwr2 /0/, ballim /0.1d0/, um /931.20793d0/                 
      if (fkey.eq.dp1) go to 10                                                 
c.......................................................................        
c///// ORNL fission test /////                                                  
      if (ievap.ne.0.and.kfiss.eq.0.and..not.fisinh.and.nofis.ne.0) then        
        if (z.ge.z91.and.(a-z).ge.z91.and.u.ge.ucut) then                       
          evapp=dp1/(dp1+dp1/gamr(a,z,u))                                       
          fissp=dp1-evapp                                                       
          call nombre_generer(rndm)
          rfiss=rndm                                                            
          if (rfiss.le.fissp) then                                              
            kfiss=1                                                             
            fisinh=.true.                                                       
            go to 90                                                            
          endif                                                                 
        endif                                                                   
      endif                                                                     
c.......................................................................        
      flkcou(1)=dp0                                                             
      flkcou(2)=ddost(1,z-flz(2))                                                
      flkcou(3)=flkcou(2)+6.d-02                                                
      flkcou(4)=flkcou(2)+1.2d-01                                               
      flkcou(6)=ddost(2,z-flz(6))                                                
      flkcou(5)=flkcou(6)-6.d-02                                                
      ccoul(1)=dp1                                                              
      ccou2=ddost(3,z-flz(2))                                                    
      ccoul(2)=ccou2+dp1                                                        
      ccoul(3)=ccou2*1.5d0+dp3                                                  
      ccoul(4)=ccou2+dp3                                                        
      ccoul(6)=ddost(4,z-flz(6))*dp2+dp2                                         
      ccoul(5)=dp2*ccoul(6)-dp1                                                 
   10 continue                                                                  
      kfiss=0                                                                   
c                                                                               
c*****************************************************************              
c                                                                               
c       main loop to compute the emission probabilities                         
c       for the 6 clusters.1=n,2=p,3=d,4=t,5=he3,6=he4                          
c                                                                               
c       alternate coding for ORNL and RAL models                                
c                                                                               
c*****************************************************************              
c                                                                               
      if (ievap.eq.0) go to 60                                                  
c.......................................................................        
c///// ORNL channel probability section /////                                   
      sigma=dp0                                                                 
      do 20 j=1,6                                                               
      s(j)=dp0                                                                  
      sos(j)=dp0                                                                
      zz=z-flz(j)                                                               
      aa=a-fla(j)                                                               
      q(j)=1.0d+10                                                              
      if (aa.lt.fla(j).or.zz.lt.flz(j).or.aa.le.zz) go to 20                    
      mm=ja-ia(j)                                                               
      q(j)=dnrgy(aa,zz)-ex1+exmass(j)                                          
      thresh(j)=q(j)+.88235d0*flkcou(j)*flz(j)*zz/(rmass(mm)+rho(j))            
      nn=aa-zz                                                                  
      izz=zz                                                                    
      smalla(j)=dgeta(u-thresh(j),izz,nn,ilvden,isdum)                           
      corr=dp0                                                                  
      if (fkey.eq.dp0) corr=cam4(izz)+cam5(nn)                                  
      arg=u-thresh(j)-corr                                                      
      if (arg.ge.dp0) then                                                      
        s(j)=sqrt(smalla(j)*arg)*dp2                                            
        sos(j)=dp10*s(j)                                                        
      endif                                                                     
   20 continue                                                                  
      n1=1                                                                      
      do 30 j=1,6                                                               
      if (sos(j).ge.1250.d0) then                                               
        ses=max(s(1),s(2),s(3),s(4),s(5),s(6))                                  
        go to 40                                                                
      endif                                                                     
   30 continue                                                                  
      n1=2                                                                      
   40 continue                                                                  
      do 50 j=1,6                                                               
      r(j)=dp0                                                                  
      if (s(j).eq.dp0) go to 50                                                 
      js=sos(j)+dp1                                                             
      mm=ja-ia(j)                                                               
      if (n1.ne.1.and.js.lt.1000) then                                          
        fjs=js                                                                  
        strun(j)=fjs-dp1                                                        
        eye1(j)=(pp1(js)+(pp1(js+1)-pp1(js))*(sos(j)-strun(j)))/smalla(j        
     1  )**2                                                                    
      else                                                                      
        if (n1.ne.1) then                                                       
          sas=exp(s(j)-50.d0)                                                   
        else                                                                    
          sas=exp(s(j)-ses)                                                     
        endif                                                                   
        eye1(j)=(s(j)**2-dp3*s(j)+dp3)*sas/(dp4*smalla(j)**2)                   
        fjs=js                                                                  
        strun(j)=fjs-dp1                                                        
      endif                                                                     
      if (j.gt.1) then                                                          
        r(j)=ccoul(j)*rmass(mm)**2*eye1(j)                                      
      else                                                                      
        if (n1.ne.1.and.js.lt.1000) then                                        
          eye0=(pp0(js)+(pp0(js+1)-pp0(js))*(sos(j)-strun(j)))/smalla(j)        
        else                                                                    
          eye0=(s(j)-dp1)*dph*sas/smalla(j)                                     
        endif                                                                   
        r(j)=rmass(mm)**2*alph(mm)*(eye1(j)+bet(mm)*eye0)                       
      endif                                                                     
      r(j)=max(r(j),dp0)                                                        
      sigma=sigma+r(j)                                                          
   50 continue                                                                  
      go to 90                                                                  
c.......................................................................        
c///// RAL channel probability section /////                                    
   60 continue                                                                  
      do 70 j=1,6                                                               
      zz=z-flz(j)                                                               
      aa=a-fla(j)                                                               
      s(j)=dp0                                                                  
      r(j)=dp0                                                                  
      argp=-dp1                                                                 
c this defeats PENBAR ?????                                                     
      q(j)=1.0d+10                                                              
      if (aa.lt.fla(j).or.zz.lt.flz(j).or.aa.le.zz) go to 70                    
      q(j)=dnrgy(aa,zz)-ex1+exmass(j)                                          
      temp=dp2*zz/aa                                                            
      mm=ja-ia(j)                                                               
      coulef=0.846927d0*flkcou(j)*flz(j)*zz/(rmass(mm)+rho(j))                  
      if (j.gt.1) coulef=coulef/(dp1+0.005d0*u/flz(j))                          
ccc   modif
c     rmass(mm)=1.22*(mm/1.)**(0.3333)
c     rho(j)=1.22*(j/1.)**(0.3333)
c     coulef=1.44d0*flz(j)*zz/((rmass(mm)+rho(j))+2.)                  
      thresh(j)=q(j)+coulef                                                     
      nn=int(aa-zz)                                                             
      izz=int(zz)                                                               
      smalla(j)=dgeta(u-thresh(j),izz,nn,ilvden,isdum)                           
      corr=cam4(izz)+cam5(nn)                                                   
      if (fkey.eq.dp1) corr=dp0                                                 
      arg=u-thresh(j)-corr                                                      
c                                                                               
c     collect available s.e. for proton for possible use                        
c     in barrier penetration when below neutron threshhold.                     
c                                                                               
      if (j.eq.2) argp=u-q(2)-corr                                              
      if (arg.lt.dp0) go to 70                                                  
      s(j)=sqrt(smalla(j)*arg)*dp2                                              
      call drein1 (j,s(j),smalla(j),eye1(j),eye0)                               
      if (j.gt.1) then                                                          
        r(j)=rmass(mm)*rmass(mm)*ccoul(j)*eye1(j)                               
      else                                                                      
        r(j)=rmass(mm)*rmass(mm)*alph(mm)*(eye1(j)+bet(mm)*eye0)                
      endif                                                                     
      r(j)=max(r(j),dp0)                                                        
   70 continue                                                                  
      smax=max(s(1),s(2),s(3),s(4),s(5),s(6))                                   
      sigma=dp0                                                                 
      do 80 j=1,6                                                               
      r(j)=r(j)*exp(s(j)-smax)                                                  
      sigma=sigma+r(j)                                                          
   80 continue                                                                  
   90 continue                                                                  
      return                                                                    
      end                                                                       
      function dclmb1 (rho,eta,ml)                                               
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c                                                                               
c     ETA = C2*Z1*Z2*SQRT(M/E)                                                  
c     RHO = C3*(R1+R2)*SQRT(ME)                                                 
c     M = reduced mass in MeV                                                   
c     E = kinetic energy in C of M                                              
c     C2 = .00516 = fine structure constant / sqrt(2)                           
c     C3 = .007165 = sqrt(2)/(hbar-c)                                           
c                                                                               
      parameter (ln0=81, lt0=21)                                                
      parameter (ln1=61, lt1=61)                                                
      common /trancm/ psi0(ln0), trans0(lt0,ln0), x0(lt0), f0(ln0),             
     1 delp0, delx0, psi1(ln1), trans1(lt1,ln1), x1(lt1), f1(ln1),              
     2 delp1, delx1                                                             
      data pi /3.14159d0/                                                       
      data c0 /.11225d0/, c1 /dph/, gamma /0.5772157d0/, s3 /0.2020569d0        
     1 /, s4/0.08232323d0/                                                      
      y=dp2*eta                                                                 
      psi=rho*y                                                                 
      if (rho.gt.y) then                                                        
        if (psi.gt.dp4.and.psi.lt.50.d0) then                                   
          prob=dclmb2(rho,eta,dumm)                                              
        else                                                                    
          x=exp(log(eta)/6.d0)                                                  
          x3=x*x*x                                                              
          temp=c0+c1*x3                                                         
          temp=temp+rho*x                                                       
          arg=dp1-y*x/temp                                                      
          prob=sqrt(arg)                                                        
        endif                                                                   
        ml=0                                                                    
      else                                                                      
        x=rho/y                                                                 
        if (psi.le.psi0(1)) then                                                
          t=min(pi*y,69.06d0)                                                   
          cx=t/(exp(t)-dp1)                                                     
          t1=cos(rho)*(dp1-.75d0*psi**2+dp5*x*psi**2)-dph*psi*rho*              
     1    sin(rho)                                                              
          t2=dp1+dph*psi*(dp1-x/6.d0)                                           
          if (eta.gt.dp1) then                                                  
            t3=log(psi)+dp2*gamma-dp1+dp1/(12.d0*eta**2)+dp1/(12.d1*            
     1        eta**4)                                                           
          else                                                                  
            t3=log(dp2*rho)+gamma-dp1/(dp1+eta**2)+s3*eta**2+s4*eta**4          
          endif                                                                 
          g=t1+psi*t2*t3                                                        
          f=cx*rho*t2                                                           
          prob=cx/(g**2+f**2)                                                   
          ml=3                                                                  
        elseif (psi.le.psi0(ln0)) then                                          
          if (x.le.x0(1)) then                                                  
            temp=log(psi/psi0(1))                                               
            j0=1+int(temp/delp0)                                                
            j0=min(max(j0,1),ln0-1)                                             
            temp=temp-(j0-1)*delp0                                              
            t=f0(j0)+(f0(j0+1)-f0(j0))*temp/delp0                               
            xk=x*sqrt(psi)                                                      
            prob=(dp1+3.33d-1*x+3.d-1*xk+1.d-1*xk**2)*exp(t)                    
            t=min(pi*y,69.06d0)                                                 
            cx=t/(exp(t)-dp1)                                                   
            prob=cx/prob**2                                                     
            ml=1                                                                
          else                                                                  
            temp1=log(x/x0(1))                                                  
            i0=1+int(temp1/delx0)                                               
            i0=min(max(i0,1),lt0-1)                                             
            temp1=temp1-(i0-1)*delx0                                            
            temp2=log(psi/psi0(1))                                              
            j0=1+int(temp2/delp0)                                               
            j0=min(max(j0,1),ln0-1)                                             
            temp2=temp2-(j0-1)*delp0                                            
            t1=trans0(i0,j0)+(trans0(i0+1,j0)-trans0(i0,j0))*temp1/delx0        
            t2=trans0(i0,j0+1)+(trans0(i0+1,j0+1)-trans0(i0,j0+1))*temp1        
     1      /delx0                                                              
            t=t1+(t2-t1)*temp2/delp0                                            
            ml=2                                                                
            prob=exp(t)                                                         
          endif                                                                 
        elseif (psi.le.psi1(ln1)) then                                          
          if (x.le.x1(1)) then                                                  
            temp=log(psi/psi1(1))                                               
            j0=1+int(temp/delp1)                                                
            j0=min(max(j0,1),ln1-1)                                             
            temp=temp-(j0-1)*delp1                                              
            t=f1(j0)+(f1(j0+1)-f1(j0))*temp/delp1                               
            xk=x*sqrt(psi)                                                      
            prob=(dp1+3.33d-1*x+3.d-1*xk+1.d-1*xk**2)*exp(t)                    
            t=min(pi*y,69.06d0)                                                 
            cx=t/(exp(t)-dp1)                                                   
            prob=cx/prob**2                                                     
            ml=1                                                                
          else                                                                  
            temp1=log(x/x1(1))                                                  
            i0=1+int(temp1/delx1)                                               
            i0=min(max(i0,1),lt1-1)                                             
            temp1=temp1-(i0-1)*delx1                                            
            temp2=log(psi/psi1(1))                                              
            j0=1+int(temp2/delp1)                                               
            j0=min(max(j0,1),ln1-1)                                             
            temp2=temp2-(j0-1)*delp1                                            
            t1=trans1(i0,j0)+(trans1(i0+1,j0)-trans1(i0,j0))*temp1/delx1        
            t2=trans1(i0,j0+1)+(trans1(i0+1,j0+1)-trans1(i0,j0+1))*temp1        
     1      /delx1                                                              
            t=t1+(t2-t1)*temp2/delp1                                            
            ml=2                                                                
            prob=exp(t)                                                         
          endif                                                                 
        else                                                                    
          prob=dclmb2(rho,eta,dumm)                                              
          ml=4                                                                  
        endif                                                                   
      endif                                                                     
      dclmb1=prob                                                                
      return                                                                    
      end                                                                       
      function dclmb2 (rho,eta,t1)                                               
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (ln=101)                                                        
      dimension t0(ln)                                                          
      data t0 /dp0,.1083d0,.1369d0,.1572d0,.1736d0,.1876d0,.2d0,.2113d0,        
     1 .2216d0,.2312d0,.2403d0,.2489d0,.2571d0,.265d0,.2725d0,.2798d0,.2        
     2 869d0,.2938d0,.3006d0,.3071d0,.3136d0,.3199d0,.3261d0,.3322d0,.33        
     3 82d0,.3442d0,.3499d0,.3557d0,.3615d0,.3672d0,.3729d0,.3785d0,.384        
     4 1d0,.3897d0,.3952d0,.4008d0,.4063d0,.4118d0,.4173d0,.4228d0,.4283        
     5 d0,.4338d0,.4393d0,.4448d0,.4504d0,.4559d0,.4615d0,.4671d0,.4728d        
     6 0,.4784d0,.4841d0,.4899d0,.4957d0,.5015d0,.5074d0,.5133d0,.5193d0        
     7 ,.5253d0,.5315d0,.5376d0,.5439d0,.5503d0,.5567d0,.5632d0,.5698d0,        
     8 .5765d0,.5833d0,.5903d0,.5973d0,.6045d0,.6118d0,.6193d0,.6269d0,.        
     9 6346d0,.6426d0,.6507d0,.659d0,.6675d0,.6763d0,.6853d0,.6945d0,.70        
     $ 4d0,.7139d0,.724d0,.7345d0,.7453d0,.7566d0,.7683d0,.7805d0,.7932d        
     $ 0,.8065d0,.8205d0,.8352d0,.8508d0,.8673d0,.8849d0,.9038d0,.9243d0        
     $ ,.9467d0,.9715d0,dp1/                                                    
      data x1 /1.d-2/, xi /1.d2/                                                
      x=dp1/(dp1+sqrt(dph*rho*eta))                                             
      if (x.lt.x1) then                                                         
        temp=t0(2)*(x/x1)**(dpth)                                               
      else                                                                      
        i=xi*x                                                                  
        i=i+1                                                                   
        i=max(min(i,ln-1),2)                                                    
        temp=t0(i)+(t0(i+1)-t0(i))*(x-i*x1)/x1                                  
      endif                                                                     
      t1=dp1-temp                                                               
      prob=dp1-dp2*t1*eta/rho                                                   
      dclmb2=max(prob,dp0)                                                       
      return                                                                    
      end                                                                       
      function ddost (i,z)                                                       
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      if (z-70.d0) 30,10,10                                                     
   10 ddost=t(i,7)                                                               
   20 return                                                                    
   30 if (z-10.d0) 40,40,50                                                     
   40 ddost=t(i,1)                                                               
      go to 20                                                                  
   50 n=.1d0*z+dp1                                                              
      x=10*n                                                                    
      x=(x-z)*.1d0                                                              
      ddost=x*t(i,n-1)+(dp1-x)*t(i,n)                                            
      go to 20                                                                  
      end                                                                       
                                                 
      subroutine emit (a,z,ja,jz,u,erec,kfiss,npart,jrn)                        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (lpt=18)                                                        
      common /filex/ cpa, wtneut(4), wtnt1(4), wtnt2(4), sprlm1,                
     1 sprlm2, wtpi, wtmu, pmass(lpt), flim0, fiswt, fiswt1, fiswt2,            
     2 swtm, lhist, lnrec, nlrnge, neutnt(4), nhstin, ir0(2),                   
     3 icpt, nobalc, nobale, ifisct,kpart(lpt), krec(8), isrfsr,                
     4 nofis, noskip, jcasc, ifbrk, ipreq, ilvden, jcoul,                       
     5 kneutp, nwrds, neutct, nlimit, nlbuf, ltme, lneut                        
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      common /forcn/ fkey                                                       
      parameter (l12=120)                                                       
      parameter (l13=50)                                                        
      common /hevaps/ epart(l12,2), hepart(l13,4), cosevp(3,l12,2),             
     1 coshev(3,l13,4)                                                          
      common /exms/ exma(6), parthr(6)                                          
      common /dresm/ exmass(6), exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      logical fisinh                                                            
      double precision kefis                                                    
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /fong/ ievap, ifiss                                                
      common /labcos/ coslbp(3), coslbr(3)                                      
      logical penbar                                                            
      common /emitr/ sos(6), strun(6), zmass(6), q(6), fla(6),                  
     1 flkcou(6), ccoul(6), thresh(6), smalla(6), flz(6), r(6), s(6),           
     2 eye1(6), smom1(6), eps, eye0, gamma, ar, zr, ex1, ex2,                   
     3 xlhs, xrhs, argp, sigma, npart0(6,2), jemiss, penbar                     
      common /inout/ lu1, lu2, itty, iscrt                                      
      dimension npart(6), c(3)                                                  
      data nwr1 /0/, nwr2 /0/, ballim /0.1d0/, um /931.20793d0/                 
      kfiss=0                                                                   
      k=jemiss                                                                  
      if (jrn.ne.0) then                                                        
        eps=u-q(k)                                                              
        go to 40                                                                
      endif                                                                     
      if (penbar) go to 40                                                      
c.......................................................................        
c///// ORNL routine for cluster energy                                          
      if (ievap.ne.0) then                                                      
        js=sos(k)+dp1                                                           
        if (js.ge.1000) then                                                    
          ratio2=(s(k)**3-6.d0*s(k)**2+15.d0*s(k)-15.d0)/((dp2*s(k)**2-         
     1    6.d0*s(k)+6.d0)*smalla(k))                                            
        else                                                                    
          ratio2=(pp2(js)+(pp2(js+1)-pp2(js))*(sos(k)-strun(k)))/smalla         
     1    (k)                                                                   
        endif                                                                   
        epsav=dp2*ratio2                                                        
c.......................................................................        
c///// RAL routine for cluster energy /////                                     
      else                                                                      
        call drein2 (s(k),smalla(k),eye2)                                       
        epsav=eye2/eye1(k)                                                      
      endif                                                                     
c.......................................................................        
      if (k.eq.1) then                                                          
        mm=ja-ia(k)                                                             
        epsav=(epsav+bet(mm))/(dp1+bet(mm)*eye0/eye1(k))                        
      endif                                                                     
c.......................................................................        
c///// RAL fission test /////                                                   
      if (ievap.eq.0.and.jz.ge.71.and.k.eq.1.and..not.fisinh.and.nofis          
     1 .ne.0) then                                                              
        agoes=(u-7.d0)/(epsav+7.d0)                                             
        if (agoes.lt.dp1.or.z.gt.88.d0) agoes=dp1                               
        probf=dfprob(z,a,u)/agoes                                                
        call nombre_generer(rndm)
        if (probf.gt.rndm) then                                         
          fisinh=.true.                                                         
          kfiss=1                                                               
          go to 90                                                              
        endif                                                                   
      endif                                                                     
c.......................................................................        
      eps0=dph*epsav                                                            
      smax=u-thresh(jemiss)                                                     
      if (smax.lt.eps0) go to 20                                                
   10 continue                                                                  
      f1=exprnf(dumm)                                                           
      f2=exprnf(dumm)                                                           
      ss=eps0*(f1+f2)                                                           
      if (ss.ge.smax) go to 10                                                  
      go to 30                                                                  
   20 continue                                                                  
      call nombre_generer(rndm)
      r1=rndm                                                                   
      test=r1*exp(smax*(dp1-r1)/eps0)                                           
      call nombre_generer(rndm)
      r2=rndm                                                                   
      if (r2.ge.test) go to 20                                                  
      ss=smax*r1                                                                
   30 continue                                                                  
      eps=ss+thresh(jemiss)-q(jemiss)                                           
c     eps=0.1+thresh(jemiss)-q(jemiss)                                           
   40 continue                                                                  
      rnmass=ex2+ar*um                                                          
      unew=(u-q(jemiss))-eps                                                    
      if (unew.lt.dp0) then                                                     
        test1=unew+u*1.0d-04                                                    
        if (test1.lt.dp0) write (15,100) u,unew,ja,jz,jemiss                   
        eps=max((u-q(jemiss)),dp0)                                              
        unew=dp0                                                                
      endif                                                                     
      rnmass=rnmass+unew                                                        
c     xkp=dph*eps*(dp2*rnmass+eps)/(rnmass+eps+zmass(jemiss))                   
      xkp=eps*(rnmass)/(rnmass+zmass(jemiss))                   
      xkt=eps-xkp                                                               
      parmas=exmass(jemiss)+um*fla(jemiss)                                      
      epp=xkp+parmas                                                            
      etp=xkt+rnmass                                                            
      aaa=xkp*(xkp+dp2*parmas)
c erreur trouvee par D.Dore decembre 98
	bbb=aaa*(gamma-dp1)*(gamma+dp1)
	if(bbb.ge.0)then                                                        
      bbb=sqrt(aaa*(gamma-dp1)*(gamma+dp1))
      else
      bbb=0.
      endif                                     
      call gtiso (c(1),c(2),c(3))                                               
      coscm=coslbr(1)*c(1)+coslbr(2)*c(2)+coslbr(3)*c(3)                        
      eps=gamma*epp+bbb*coscm-parmas                                            
      erec=gamma*etp-bbb*coscm-rnmass                                           
      aa1=aaa**2*(dp1-coscm**2)                                                 
      dnp=sqrt(aa1+(bbb*epp+gamma*aaa*coscm)**2)                                
      dnt=sqrt(aa1+(bbb*etp-gamma*aaa*coscm)**2)                                
      c1=bbb*epp+aaa*(gamma-dp1)*coscm                                          
      c2=bbb*etp-aaa*(gamma-dp1)*coscm                                          
      do 50 i=1,3                                                               
      coslbp(i)=(c1*coslbr(i)+aaa*c(i))/dnp                                     
      coslbr(i)=(c2*coslbr(i)-aaa*c(i))/dnt                                     
   50 continue                                                                  
      u=unew                                                                    
      npart(jemiss)=npart(jemiss)+1                                             
      ja=ja-ia(jemiss)                                                          
      jz=jz-iz(jemiss)                                                          
      xrhs=u+eps+erec+exmass(jemiss)+ex2                                        
      if (abs(xlhs-xrhs).ge.ballim) then                                        
        nwr2=nwr2+1                                                             
        if (nwr2.lt.10) write (15,110) xlhs,xrhs,ja,jz,jemiss                  
      endif                                                                     
c*****store,end of normal cycle                                                 
      smom1(jemiss)=smom1(jemiss)+eps                                           
      if (npart(jemiss).le.0) call errorf ('emit')                              
      if (jemiss.gt.2) go to 70                                                 
      epart(npart(jemiss),jemiss)=eps                                           
      do 60 i=1,3                                                               
   60 cosevp(i,npart(jemiss),jemiss)=coslbp(i)                                  
      go to 90                                                                  
   70 kemiss=jemiss-2                                                           
      hepart(npart(jemiss),kemiss)=eps                                          
      do 80 i=1,3                                                               
   80 coshev(i,npart(jemiss),kemiss)=coslbp(i)                                  
   90 continue                                                                  
c.......................................................................        
      return                                                                    
c                                                                               
  100 format (' *** u=',g10.3,' unew=',g10.3,' ja=',i3,' jz=',i3,' jemis        
     1s=',i2)                                                                   
  110 format (' *** xlhs=',g10.3,' xrhs=',g10.3,' ja=',i3,' jz=',i3,' je        
     1miss=',i2)                                                                
      end                                                                       
      function dnrgy (a,z)                                                     
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      dimension jprime(250)                                                     
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      cam(a,z,b)=8.36755d0*a-0.78244d0*z-17.0354d0*a*(dp1-1.84619d0*(a          
     1 -dp2*z)**2/a**2)+25.8357d0*b**2*(dp1-1.71219d0*(a-dp2*z)**2/a**2)        
     2 *(dp1-0.62025d0/b**2)**2+0.779d0*z*(z-dp1)*(dp1-1.5849d0/b**2+1.2        
     3 273d0/a+1.5772d0/(a*b))/b-0.4323d0*exp(1.3333333333d0*log(z))*           
     4 (dp1-0.57811d0/b-0.14518d0/b**2+0.49597d0/a)/b                           
      data jprime /19*0,10*1,8*2,7*3,7*4,6*5,5*6,6*7,5*8,5*9,5*10,5*11,4        
     1 *12,5*13,4*14,4*15,4*16,4*17,4*18,4*19,4*20,4*21,4*22,3*23,4*24,4        
     2 *25,3*26,4*27,3*28,4*29,3*30,3*31,4*32,3*33,3*34,4*35,3*36,3*37,3        
     3 *38,3*39,3*40,3*41,3*42,3*43,3*44,3*45,3*46,3*47,3*48,3*49,3*50,3        
     4 *51,3*52,3*53,3*54,55,55,3*56,3*57,3*58,2*59,2*60/                       
c......................................................................         
c      deltas (a,b)=(a-1.0)/(1.0+(124.0/b**2))                                  
c      i=a                                                                      
c      a3=exp(log(a)/dp3)                                                       
c      jprime(i)=deltas(a,a3)                                                   
c......................................................................         
      i=a                                                                       
      kz=z                                                                      
      j=i-2*kz-jprime(i)+10                                                     
      if (j.gt.0.and.j.le.20.and.i.le.250) then                                 
        dnrgy=waps(i,j)                                                        
        if (dnrgy.ne.dp0) return                                               
      endif                                                                     
      n=a-z                                                                     
      a3=exp(0.3333333333d0*log(a))                                             
      dnrgy=cam(a,z,a3)+cam2(kz)+cam3(n)                                       
      return                                                                    
      end                                                                       
      subroutine errorf (loc)                                                   
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      character*8 loc                                                           
      parameter (lpt=18)                                                        
      common /inout/ in, io, itty, iscrt                                        
      nocas=igetd(1)                                                            
      write (io,10) loc,nocas                                                   
   10 format (//' error exit in subroutine ',a8,' nocas=',i8)                   
      write (itty,10) loc,nocas                                                 
      stop                                                                      
      end                                                                       
      function exprnf (a)                                                       
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      ri=dp0                                                                    
   10 continue                                                                  
      call nombre_generer(rndm)
      x=rndm                                                                    
      z=x                                                                       
   20 continue                                                                  
      call nombre_generer(rndm)
      y=rndm                                                                    
      if (z-y) 50,50,30                                                         
   30 continue                                                                  
      call nombre_generer(rndm)
      z=rndm                                                                    
      if (z-y) 20,40,40                                                         
   40 ri=ri+dp1                                                                 
      go to 10                                                                  
   50 exprnf=x+ri                                                               
      return                                                                    
      end                                                                       
      subroutine fisdis                                                         
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      logical fisinh                                                            
      double precision kefis                                                    
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /labcos/ coslbp(3), coslbr(3)                                      
      dimension g(3)                                                            
cccc                                                                            
c     select angles in center of momentum system isotropically                  
cccc                                                                            
      call gtiso (g(1),g(2),g(3))                                               
      coscm=coslf0(1)*g(1)+coslf0(2)*g(2)+coslf0(3)*g(3)                        
cccc                                                                            
c     compute fragment kinetic energies and velocities relativistically         
cccc                                                                            
      psq=dp2*amdiff*amc1*amc2/amcf                                             
      pscorr=(dp1-amdiff/(dp2*amcf))*(dp1+amdiff/(dp2*amc2))*(dp1+amdiff        
     1 /(dp2*amc1))                                                             
      psq=psq*pscorr                                                            
      gamma=dp1+ernff/amcf
c erreur trouvee par D.Dore decembre 98
	bbb=aaa*(gamma-dp1)*(gamma+dp1)
	if(bbb.ge.0)then                                                        
      bbb=sqrt(aaa*(gamma-dp1)*(gamma+dp1))
      else
      bbb=0.
      endif                                                                   
                                    
      ecm1=(amcf**2+amc1**2-amc2**2)/(dp2*amcf)                                 
      ecm2=(amcf**2-amc1**2+amc2**2)/(dp2*amcf)                                 
      kefis(1)=gamma*ecm1+bbb*coscm-amc1                                        
      kefis(2)=gamma*ecm2-bbb*coscm-amc2                                        
      aa0=psq**2*(dp1-coscm**2)                                                 
      dn1=sqrt(aa0+(bbb*ecm1+gamma*psq*coscm)**2)                               
      dn2=sqrt(aa0+(bbb*ecm2-gamma*psq*coscm)**2)                               
      c1=bbb*ecm1+psq*(gamma-dp1)*coscm                                         
      c2=bbb*ecm2-psq*(gamma-dp1)*coscm                                         
      do 10 i=1,3                                                               
      coslf1(i)=(c1*coslf0(i)+psq*g(i))/dn1                                     
      coslf2(i)=(c2*coslbr(i)-psq*g(i))/dn2                                     
   10 continue                                                                  
      return                                                                    
      end                                                                       
      subroutine fissed (ja,jz,u,erec)                                          
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
c                                                                               
c      subroutine to pick post fission parameters for nucleus                   
c      jz,ja excited to u and recoiling with erec.                              
c                                                                               
      logical fisinh                                                            
      double precision kefis                                                    
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /labcos/ coslbp(3), coslbr(3)                                      
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      common /fbarr/ ef                                                         
      common /inout/ lu1, lu2, itty, iscrt                                      
c      dimension g(3), evodba(4), evodbs(4)   
      dimension  evodba(4), evodbs(4)                                       
c     following data for assymetric vs symetric picking                         
      data evodba, evodbs/18.8d0,18.1d0,18.1d0,18.5d0,20.43d0,20.43d0,          
     1 20.43d0,20.712d0/                                                        
      data sigmaa, y00, aamean /6.5d0,0.1176d0,140.d0/                          
c     asym gauss width,1.0/level density parameter and high mass mean.          
      data afact, zfact, enmass /8.3674d0,0.7826d0,939.5124d0/                  
c     mass diff from 1 amu for neutron,diff p & n masses and n mass             
      if (.not.fisinh) then                                                     
c***********************************************************************        
c                                                                               
c  logic error fissed called with fisinh unset.                                 
c                                                                               
c***********************************************************************        
        write (lu2,60) ja,jz,u,erec                                             
        call errorf ('fissed')                                                  
      endif                                                                     
c***********************************************************************        
c                                                                               
c       pick the masses                                                         
c                                                                               
c***********************************************************************        
      z=real(jz)                                                                
      a=real(ja)                                                                
      nck2=0                                                                    
   10 continue                                                                  
      nck2=nck2+1                                                               
      if (nck2.gt.10) then                                                      
c***********************************************************************        
c                                                                               
c     fission failure                                                           
c                                                                               
c***********************************************************************        
        write (lu2,50) ja,jz,u,erec,zfis(1),afis(1),zfis(2),afis(2)             
        fisinh=.false.                                                          
        return                                                                  
      endif                                                                     
      temp=z*z/a                                                                
      if (temp.le.35.d0) then                                                   
        if (jz.le.88) go to 20                                                  
      elseif (u.le.62.d0) then                                                  
c***********************************************************************        
c                                                                               
c   high z fission mass distribution. competition for sym vs assym              
c  simple symmetric to assymetric data fit                                      
c                                                                               
c***********************************************************************        
        arg=-0.36d0*u                                                           
        arg=4.87d+03*exp(arg)                                                   
        proba=arg/(dp1+arg)                                                     
        if (fltrnf(dumm).le.proba) then                                         
c***********************************************************************        
c                                                                               
c      assymetric fission                                                       
c                                                                               
c***********************************************************************        
          a1=dgaussn(aamean,sigmaa)                                              
          go to 30                                                              
        endif                                                                   
      endif                                                                     
c***********************************************************************        
c                                                                               
c     find assymetric barrier for width computation                             
c     assymetric barrier from seaborg and vandenbosch                           
c     phys.rev. 88,507 (1952) phys.rev. 110,507 (1958)                          
c     find if ee eo oe or oo nucleus                                            
c     na=1 odd-odd,2 even-odd,3 odd-even,4 even-even                            
c                                                                               
c***********************************************************************        
      in=ja-jz                                                                  
      na=1                                                                      
      if (jz.eq.2*(jz/2)) na=na+1                                               
      if (in.eq.2*(in/2)) na=na+2                                               
      temp=z*z/a                                                                
      ef=evodba(na)-0.36d0*temp                                                 
   20 continue                                                                  
      upr=u-ef                                                                  
      if (upr.gt.1.d+02) upr=1.d+02                                             
      sigmas=0.425d0*(upr*(dp1-0.005d0*upr)+9.35d0)                             
c     parametre utilise par Atchison
c     sigmas=0.425d0*(upr*(dp1-0.05d0*upr)+0.935d0)                             
c***********************************************************************        
c                                                                               
c     sigmas is symmetrin fission mass width.taken from systematic of           
c     neuzil & fairhall phys.rev. 129,2705,(1963)                               
c                                                                               
c      low z fission is always symmetric                                        
c                                                                               
c      high z fission is sometimes.loop back here from 1 loop if                
c      symmetric fission predicted.                                             
c                                                                               
c***********************************************************************        
      amean=dph*a                                                               
      a1=dgaussn(amean,sigmas)                                                   
   30 continue                                                                  
c***********************************************************************        
c                                                                               
c      1 loop for assymmetrin fission returns to here.                          
c                                                                               
c***********************************************************************        
      afis(1)=aint(a1)                                                          
c***********************************************************************        
c                                                                               
c    check for low final a                                                      
c                                                                               
c***********************************************************************        
      if (afis(1).lt.dp5) afis(1)=dp5                                           
      if (a-afis(1).lt.dp5) afis(1)=a-dp5                                       
      afis(2)=a-afis(1)                                                         
c***********************************************************************        
c                                                                               
c        pick the charge                                                        
c                                                                               
c***********************************************************************        
      z1=65.5d0*afis(1)/(131.d0+afis(1)**0.666667d0)                            
      z2=65.5d0*afis(2)/(131.d0+afis(2)**0.666667d0)                            
      z1=z1+dph*(z-z1-z2)                                                       
c***********************************************************************        
c                                                                               
c       we use constant charge density with a 2 unit gaussian smearing          
c                                                                               
c***********************************************************************        
      z1=dgaussn(z1,dp2)                                                         
      zfis(1)=aint(z1)                                                          
      zfis(2)=z-zfis(1)                                                         
c***********************************************************************        
c                                                                               
c   check for reasonable z a combinations.                                      
c                                                                               
c***********************************************************************        
      if (zfis(1).ge.afis(1)) go to 10                                          
      if (zfis(2).ge.afis(2)) go to 10                                          
      if (zfis(1).lt.dp1) go to 10                                              
      if (zfis(2).lt.dp1) go to 10                                              
c***********************************************************************        
c                                                                               
c      compute binding energy and actual masses of fragments                    
c                                                                               
c***********************************************************************        
      be0=afact*a-zfact*z-dnrgy(a,z)                                           
      rm0=enmass*a-zfact*z-be0                                                  
      be1=afact*afis(1)-zfact*zfis(1)-dnrgy(afis(1),zfis(1))                   
      rm1=enmass*afis(1)-zfact*zfis(1)-be1                                      
      be2=afact*afis(2)-zfact*zfis(2)-dnrgy(afis(2),zfis(2))                   
      rm2=enmass*afis(2)-zfact*zfis(2)-be2                                      
c***********************************************************************        
c                                                                               
c      pick recoil kinetic energy.use systematic of ...........                 
c      unik et.al. proc.3rd iaea symp.on phy.& chem. fision,rochester.vo        
c                                                                               
c***********************************************************************        
      totkm=0.13323d0*z*z/a**0.33333333d0-11.4d0                                
c     parametre utilise par Atchison
c     totkm=0.1065d0*z*z/a**0.333333d0+20.1d0
c***********************************************************************        
c                                                                               
c      use a width of 15% value at half height.                                 
c                                                                               
c***********************************************************************        
      sigmak=0.084d0*totkm                                                      
c***********************************************************************        
c                                                                               
c   check event is energetically possible                                       
c                                                                               
c***********************************************************************        
      temp2=u+be1+be2-be0                                                       
      nck=0                                                                     
   40 continue                                                                  
      totke=dgaussn(totkm,sigmak)                                                
      if (nck.gt.10) go to 10                                                   
      nck=nck+1                                                                 
      if (totke.gt.temp2) go to 40                                              
c***********************************************************************        
c                                                                               
c      pick excitation from equidistribution of original plus energy bal        
c                                                                               
c***********************************************************************        
      temp=(temp2-totke)/a                                                      
      ufis(1)=afis(1)*temp                                                      
      ufis(2)=afis(2)*temp                                                      
cccc                                                                            
c     find total masses, including excitation energies, at evap time            
cccc                                                                            
      amcf=rm0+u                                                                
      amc1=rm1+ufis(1)                                                          
      amc2=rm2+ufis(2)                                                          
      amdiff=amcf-amc1-amc2                                                     
      if (amdiff.lt.dp0) go to 40                                               
c***********************************************************************        
c                                                                               
c     amdiff= ekin should be satisfied                                          
c                                                                               
c***********************************************************************        
      ernff=erec                                                                
      call fisdis                                                               
      return                                                                    
c                                                                               
   50 format ('---> fission failed: ja=',i5,'  jz=',i5/'-       u=',1pe1        
     1 0.3,'      erec=',e10.3/'    zfis1=',e10.3,'     afis1=',e10.3/'         
     2   zfis2=',e10.3,'     afis2=',e10.3)                                     
   60 format (//'  logic error in fissed.called with fisinh flag',' unse        
     1t.',2i10,2f10.5)                                                          
      end                                                                       
      function dfprob (z,a,u)                                                    
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c                                                                               
c   function to compute the fission probability.                                
c     for z<90 uses ..........                                                  
c    uses statistical model fits.                                               
c                                                                               
      common /fbarr/ ef                                                         
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      parameter (lnfis=18)                                                      
      dimension slope(lnfis), anort(lnfis)                                      
      double precision i0, i1                                                   
      data a1, a2, a3, const /0.2185024d0,16.70314d0,321.175d0,                 
     1 0.3518099d0/                                                             
      data c1, c2, c3 /1.089257d0,0.01097896d0,31.08551d0/                      
c     forz> 89 and <100 ........                                                
c     use the systematics of vandenbosch & huizenga with a ball park            
c     observation that fission probability drops for most nuclei at 6 m         
      data slope /0.23d0,0.233d0,0.12225d0,0.14727d0,0.13559d0,0.15735d0        
     1 ,0.16597d0,0.17589d0,0.18018d0,0.19568d0,0.16313d0,0.17123d0,            
     2 6*0.17d0/                                                                
c     data anort /221.6,226.9,229.75,234.04,238.88,241.34,243.04,245.52,        
c    1 246.84,250.18,254.0,257.8/                                               
      data anort /219.4d0,226.9d0,229.75d0,234.04d0,238.88d0,241.34d0,          
     1 243.04d0,245.52d0,246.84d0,250.18d0,254.0d0,257.8d0,261.3d0,             
     2 264.8d0,268.3d0,271.8d0,275.3d0,278.8d0/                                 
      dfprob=dp0                                                                 
      iz=int(z)                                                                 
      if (iz.gt.88) then                                                        
c                                                                               
c     high z fission probability                                                
c                                                                               
        iz=iz-88                                                                
      if (iz.gt.lnfis.or.u.lt.6.0) go to 20                                     
        ja=int(a)                                                               
        gamnf=slope(iz)*(real(ja)-anort(iz))                                    
        gamnf=dp10**gamnf                                                       
        go to 10                                                                
      endif                                                                     
      in=int(a-z)                                                               
      se=dnrgy((a-dp1),z)-dnrgy(a,z)+8.3674d0-cam2(iz)-cam3(in-1)             
      if ((2*(iz/2).ne.iz).or.(2*(in/2).ne.in)) se=se-cam4(iz)-cam5(in-1        
     1 )                                                                        
      if ((2*(iz/2).ne.iz).and.(2*(in/2).ne.in)) se=se-cam4(iz)-cam5(in-        
     1 1)                                                                       
      x=z*z/a                                                                   
      ef=x*(a1*x-a2)+a3+se                                                      
      if (ef.gt.u) go to 20                                                     
      an=0.125d0*(a-dp1)                                                        
      an1=dph/an                                                                
      an2=0.25d0*an1/an                                                         
      af=x-c3                                                                   
      af=an*(c1+c2*af*af)                                                       
      a1thrd=a**0.33333333d0                                                    
      s=dp2*sqrt(an*(u-se))                                                     
      if (s.gt.dp10) then                                                       
        i0=an1*(s-dp1)                                                          
        i1=s*an2*((s+s-6.d0)+6.d0)                                              
        gamnf=const*(((0.76d0*i1-5.d-02*i0)*a1thrd+1.93d0*i1)*a1thrd+           
     1  1.66d0*i0)                                                              
        s2=dp2*sqrt(af*(u-ef))                                                  
        exps=dp0                                                                
        if ((s-s2).gt.-150.d0) exps=exp(s-s2)                                   
        gamnf=gamnf*exps*af/(s2-dp1)                                            
      else                                                                      
        e1=exp(s)                                                               
        i0=((s-dp1)*e1+dp1)*an1                                                 
        i1=((6.d0+s*(s+s-6.d0))*e1+s*s-6.d0)*an2                                
        gamnf=const*(((0.76d0*i1-5.d-02*i0)*a1thrd+1.93d0*i1)*a1thrd+           
     1  1.66d0*i0)                                                              
        s=dp2*sqrt(af*(u-ef))                                                   
        e2=((s-dp1)*exp(s)+dp1)/af                                              
        gamnf=gamnf/e2                                                          
      endif                                                                     
   10 continue                                                                  
      dfprob=dp1/(dp1+gamnf)                                                     
   20 continue                                                                  
      return                                                                    
      end                                                                       
      subroutine fsinfo                                                         
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      common /labcos/ coslbp(3), coslbr(3)                                      
      common /smacom/ yzero, bzero, yzere, bzere                                
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      logical fisinh                                                            
      double precision kefis                                                    
      common /fishun/ afis(10), zfis(10), ufis(10), kefis(10), amdiff,          
     1 atfis(10), ztfis(10), utfis(10), recfis(10), coslf0(3), rnmass,          
     2 rfmass, coslf1(3), coslf2(3), ernff, amcf, amc1, amc2, fisinh            
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /orcm2/ ainfo(20), epkin(80,6), sigke(80,6), delexc(80,6,3)        
     1 , cstiff(120)                                                            
      common /inout/ lu1, lu2, itty, iscrt                                      
c      dimension g(3)                                                            
      data cnver /931.143d0/, pcks /2.d-02/, ekcalm /dp10/,                     
     1 kecalm /50/, zcalm /8.d1/, ecalm /40.d1/                                 
      interr=0                                                                  
      ncall=0                                                                   
      d0=dp0                                                                    
cccc                                                                            
c     ro    = 1.5d0                                                             
c     ro passed in common from block data                                       
cccc                                                                            
      rocc=ro                                                                   
cccc                                                                            
c     the mass differences are computed from actual wapstra mass tables         
c     set the scale factors: these are quite arbitrary at present               
c     a more complicated scaling function for a could be used                   
cccc                                                                            
      scala=af/afnp                                                             
      scalke=dp1                                                                
      omegao=sqrt(dp4*ef*afnp/b0)                                               
      rfmass=dnrgy(af,zf)                                                      
cccc                                                                            
c     calculate the value of delef, for incoming ef                             
cccc                                                                            
      a1m=afnp-a2m                                                              
      call interp (a1m)                                                         
cccc                                                                            
c     for neptunium, calculate qamax, the value of the probability funct        
c     for a2, at a2 = a2m                                                       
cccc                                                                            
      call qmaxx (a2m,qamax,q)                                                  
cccc                                                                            
c     multiply the calculated value by qfactr to account for variations         
c     in a2m as ef changes                                                      
cccc                                                                            
      qamax=qamax*qfactr                                                        
cccc                                                                            
c     the value of a2, heavy fragment, will be chosen by rejection              
c     double interpolate in ef and a for neptunium                              
cccc                                                                            
      ciamin=iamin                                                              
      acall=dp0                                                                 
   10 continue                                                                  
      acall=acall+1                                                             
      if (acall.gt.acalm) go to 110                                             
      ekcall=dp0                                                                
      zcall=dp0                                                                 
      ecall=dp0                                                                 
      kecall=0                                                                  
      call nombre_generer(rndm)
      aprime=afnp/dp2+rndm*(afnp/dp2-ciamin)                            
      a2np=aprime                                                               
      a2=aprime*scala                                                           
      iint=1                                                                    
   20 continue                                                                  
      ia2np=a2np                                                                
      cia2np=ia2np                                                              
      dela=a2np-cia2np                                                          
      ia2=ia2np                                                                 
      delkeh=dp0                                                                
      delech=dp0                                                                
      k2=ke+1                                                                   
      delkel=sigke(ia2-100,ke)*(dp1-dela)+sigke(ia2-99,ke)*dela                 
      delecl=delexc(ia2-100,ke,ismal)*(dp1-dela)+delexc(ia2-99,ke,ismal)        
     1 *dela                                                                    
      if (ke.eq.kemax) go to 30                                                 
      delkeh=sigke(ia2-100,k2)*(dp1-dela)+sigke(ia2-99,k2)*dela                 
      delech=delexc(ia2-100,k2,ismal)*(dp1-dela)+delexc(ia2-99,k2,ismal)        
     1 *dela                                                                    
   30 continue                                                                  
      delken=delkel*(dp1-delef)+delkeh*delef                                    
      delecn=delecl*(dp1-delef)+delech*delef                                    
      emaxn=delecn+ef                                                           
cccc                                                                            
c     compute smasum for neptunium and assume arg is invariant                  
cccc                                                                            
      call keave (a2np,a2,smasmm)                                               
      arg=smasmm*emaxn                                                          
      z1pnp=z1p                                                                 
      z2pnp=z2p                                                                 
      if (iint.eq.0) go to 40                                                   
      call qmaxx (a2np,qa,q)                                                    
      qaa=qa/qamax                                                              
      if (qaa.gt.dp1) write (6,120) qfactr,qaa,qa,qamax                       
      call nombre_generer(rndm)
      raa=rndm                                                                  
      if (raa.gt.qaa) go to 10                                                  
cccc                                                                            
c     choose nearest integer value of a2 for true af                            
cccc                                                                            
      ia2=a2+dph                                                                
      a2=ia2                                                                    
      a1=af-a2                                                                  
      iah1=a1-58                                                                
      if (iah1.lt.1) then                                                       
        iah1=1                                                                  
        a1=59.d0                                                                
        a2=af-a1                                                                
      endif                                                                     
      iah2=a2-58                                                                
      if (iah2.lt.1) then                                                       
        iah2=1                                                                  
        a2=59.d0                                                                
        a1=af-a2                                                                
        iah1=a1-58                                                              
      endif                                                                     
cccc                                                                            
c     assume the integer values of a2 and a1 have been chosen by process        
c     use integer value of a2 from now on to establish corresponding a2n        
c     for further interpolation to get ave rel ke and z2p,etc                   
cccc                                                                            
      a2np=a2/scala                                                             
      iint=0                                                                    
      go to 20                                                                  
   40 continue                                                                  
      rocc=ro                                                                   
      ekenpm=ekinr4                                                             
cccc                                                                            
c     set value of ekso                                                         
cccc                                                                            
      ekso=dp0                                                                  
      coulnp=ekenpm-ekso                                                        
cccc                                                                            
c     insert a theory to caculate drat here, if desired.                        
cccc                                                                            
      drat=dp1                                                                  
cccc                                                                            
c     set ratio of def1/def2 eq to lowest possible value                        
c     the ratio of stiffness coefficients is assumed equal to that given        
c     using fong's coefficients, from block data                                
cccc                                                                            
      drat=cstiff(iah2)/cstiff(iah1)                                            
      drat=drat*(a1/a2)**0.66666667d0                                           
cccc                                                                            
      coulm=scalke*coulnp                                                       
      ekeafm=coulm+ekso                                                         
cccc                                                                            
c     compute interpolated value of rmcbar                                      
cccc                                                                            
      ekcall=ekcall+1                                                           
      if (ekcall.gt.ekcalm) go to 50                                            
      interr=1                                                                  
      ekinr4=coulm                                                              
      z2p=zprobf(af,zf,a2)                                                      
      z1p=zf-z2p                                                                
      izz2=z2p+5.d-04                                                           
      cizz2=izz2                                                                
      izz1=z1p+5.d-04                                                           
      cizz1=izz1                                                                
      delz1p=z1p-cizz1                                                          
      delz2p=z2p-cizz2                                                          
      r1mass=dnrgy(a1,cizz1)*(dp1-delz1p)+dnrgy(a1,cizz1+dp1)*delz1p          
      r2mass=dnrgy(a2,cizz2)*(dp1-delz2p)+dnrgy(a2,cizz2+dp1)*delz2p          
      pairl=pairr(a1,izz1+1,a2,izz2,corr1h,corr2l)                              
      pairh=pairr(a1,izz1,a2,izz2+1,corr1l,corr2h)                              
      pairsm=pairl*(dp1-delz2p)+pairh*delz2p                                    
      rmcbar=rfmass+ef-r1mass-r2mass-pairsm                                     
cccc                                                                            
c     use invariance of arg to find ebar for af                                 
cccc                                                                            
      c=(a2-z2p*dp2)/a2                                                         
      smala2=a2*(dp1+yzero*c*c)/bzero                                           
      c=(a1-z1p*dp2)/a1                                                         
      smala1=a1*(dp1+yzero*c*c)/bzero                                           
      smasum=smala1+smala2                                                      
      ebar=arg/smasum                                                           
cccc                                                                            
c     fix scalke here, if desired                                               
c     scalke = z1pnp*z2pnp/(z1p*z2p)                                            
c     find mean relative kinetic energy by scaling from neptunium values        
cccc                                                                            
      delke=scalke*delken                                                       
      coulm=scalke*coulnp                                                       
      ekeafm=coulm+ekso                                                         
cccc                                                                            
c     compute the average value of deformation energy                           
cccc                                                                            
      dbar=rmcbar-ekeafm-ebar                                                   
      if (dbar.lt.dp0) then                                                     
        dbar=dp0                                                                
        ebarnw=rmcbar-ekeafm                                                    
        if (ebarnw.le.dp0) ekeafm=rmcbar-ebar                                   
        if (ekeafm.le.dp0) go to 110                                            
      endif                                                                     
   50 continue                                                                  
      if (iversn.eq.1) go to 70                                                 
      dlowke=dbar                                                               
cccc                                                                            
c     compute d0 and ro                                                         
cccc                                                                            
      vs=coulm**2/(coulm-dp2*dbar)                                              
      c12v=c12s*rocc                                                            
      v=c12v*z1p*z2p                                                            
      rocc=v/vs                                                                 
c     should iterate here to get new z2p, or use correct solution of cub        
      interr=interr+1                                                           
      c12sn=c12v/rocc                                                           
      c12s=c12sn                                                                
      couls=c12s*z1p*z2p                                                        
      ckeav=coulm/couls                                                         
      if (dbar.ne.dp0) d0=dbar/((dp1/ckeav)-dp1)**2                             
cccc                                                                            
c     choose the value of the relative kinetic energy of fission fragmen        
c     choose relative kinetic energy for af from a truncated qaussian           
cccc                                                                            
   60 continue                                                                  
      kecall=kecall+1                                                           
      if (kecall.gt.kecalm) go to 10                                            
      ekin=gaurn(ekin)*delke+ekeafm                                             
      coulm=ekin-ekso                                                           
      if (coulm.le.dp0) go to 60                                                
      dlowke=d0*(vs/coulm-dp1)**2                                               
      em=rmcbar-ekin-dlowke                                                     
      if (em.lt.1.0d-8) go to 60                                                
cccc                                                                            
c     assume ekin has been chosen by process above                              
c     choose z1 and z2 for the two fission fragments                            
c     find sigmaz of gaussian z2 distribution                                   
c     should find z2p(ekin) here                                                
cccc                                                                            
      arg=smasum*(rmcbar-ekin-dlowke)                                           
   70 continue                                                                  
      fact=sqrt(arg)                                                            
      sigmaz=(dp1-fm/fact)*sqrt(fact/(dp2*baz*smasum))                          
cccc                                                                            
c     determine minimum and maximum z2                                          
c     assign arbitrary values to zlow and zhigh, initially                      
cccc                                                                            
      zl=z2p-izdel                                                              
      zh=z2p+izdel                                                              
cccc                                                                            
c     choose z2 from a truncated gaussian distribution                          
cccc                                                                            
   80 continue                                                                  
      zcall=zcall+1                                                             
      if (zcall.gt.zcalm) go to 10                                              
      z=sigmaz*gaurn(z)+z2p                                                     
      if (z.le.zl.or.z.ge.zh) go to 80                                          
cccc                                                                            
c     assume z has been chosen by process above                                 
cccc                                                                            
      z2=z                                                                      
      z1=zf-z2                                                                  
cccc                                                                            
c     select nearest integers                                                   
cccc                                                                            
      izz1=z1+dph                                                               
      izz2=z2+dph                                                               
      z1=izz1                                                                   
      z2=izz2                                                                   
cccc                                                                            
c     compute difference between initial energy and masses of fragments         
c     at the characteristic level;i.e.,rmc                                      
cccc                                                                            
      r1mass=dnrgy(a1,z1)                                                      
      r2mass=dnrgy(a2,z2)                                                      
      pairsm=pairr(a1,izz1,a2,izz2,corr1,corr2)                                 
      rmc=rfmass+ef-r1mass-r2mass-pairsm                                        
      if (iversn.eq.1) ekin=c12*z1*z2/(dp1-pcks)                                
cccc                                                                            
c     set minimum deformation energy and determine maximum excitation           
c     energy of the fission fragments                                           
cccc                                                                            
      def=dbar                                                                  
      if (iversn.ne.1) then                                                     
        cke=(ekin-ekso)/(c12s*z1*z2)                                            
        if (cke.gt.dp0) then                                                    
          dlowke=d0*(dp1/cke-dp1)**2                                            
          def=dlowke                                                            
        endif                                                                   
      endif                                                                     
      e1maxp=rmc-ekin-def                                                       
cccc                                                                            
c     now, do check on energy conservation                                      
c     select new z1,z2 pair                                                     
cccc                                                                            
      if (e1maxp.le.1.0d-8) go to 80                                            
      e=e1maxp                                                                  
cccc                                                                            
c     select scission point excitation energies for chosen fragments            
cccc                                                                            
      c=(a2-z2*dp2)/a2                                                          
      smala2=a2*(dp1+yzero*c*c)/bzero                                           
      c=(a1-z1*dp2)/a1                                                          
      smala1=a1*(dp1+yzero*c*c)/bzero                                           
   90 continue                                                                  
      ecall=ecall+1                                                             
      if (ecall.gt.ecalm) go to 10                                              
      smasum=smala1+smala2                                                      
      call nombre_generer(rndm)
      e1c=e*rndm                                                                
      qearg=(dp2*sqrt(smala1*e1c)+dp2*sqrt(smala2*(e-e1c))-dp2*sqrt             
     1 (smasum*e1maxp))                                                         
      qe=exp(qearg)                                                             
      qe=qe*(e1c*(e-e1c)/e1maxp**2)**1.5d0                                      
      call nombre_generer(rndm)
      re=rndm                                                                   
      if (re.gt.qe) go to 90                                                    
      e2c=e-e1c                                                                 
cccc                                                                            
c     assume the scission point excitation energies have been chosen;           
c     find excitation energies at evap time                                     
cccc                                                                            
      def2=def/(dp1+drat)                                                       
      def1=drat*def2                                                            
cccc                                                                            
c     if magic no corrections have been added, they should be in corr1 a        
cccc                                                                            
      e1=e1c+def1+corr1                                                         
      e2=e2c+def2+corr2                                                         
cccc                                                                            
c     find total masses, including excitation energies, at evap time            
cccc                                                                            
      am1=a1*cnver+r1mass+e1                                                    
      am2=a2*cnver+r2mass+e2                                                    
      amcf=rfmass+ef+af*cnver                                                   
      amc1=am1                                                                  
      amc2=am2                                                                  
      amdiff=rfmass+ef-r1mass-r2mass-e1-e2                                      
cccc                                                                            
c     amdiff= ekin should be satisfied                                          
cccc                                                                            
      ernff=ernf                                                                
      call fisdis                                                               
      eke1=kefis(1)                                                             
      eke2=kefis(2)                                                             
      afis(1)=a1                                                                
      afis(2)=a2                                                                
      zfis(1)=z1                                                                
      zfis(2)=z2                                                                
      ufis(1)=e1                                                                
      ufis(2)=e2                                                                
  100 continue                                                                  
      return                                                                    
cccc                                                                            
  110 continue                                                                  
      ncall=ncall+1                                                             
      write (1,130) acall,zcall,ecall,kecall,af,zf                            
      go to 100                                                                 
c                                                                               
  120 format (1h ,'fsinfo ',2x,'qfactr =',1pe12.5,2x,'qaa =',e12.5,2x,'q        
     1a =',e12.5,2x,'qamax =',e12.5)                                            
  130 format (1h ,'fission failed: acall=',f8.2,2x,'zcall=',f8.2,2x,'eca        
     1ll=',f8.2,2x,'kecall=',i5,2x,'af=',f6.1,2x,'zf=',f6.1)                    
      end                                                                       
      function gammar (x)                                                       
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  GAMMAR computes the gamma funtion of x.                                      
c                                                                               
c    REFERENCES.......  GMS                                                     
c    EXTERNALS........  EXP                                                     
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      if (x.gt.1) then                                                          
        z=x-dp1                                                                 
      else                                                                      
        z=x                                                                     
   10   continue                                                                
        if (z.le.dp0) then                                                      
          z=z+dp1                                                               
          go to 10                                                              
        endif                                                                   
      endif                                                                     
      t1=z+dph                                                                  
      tzg=t1+dp5                                                                
      t1=tzg**t1                                                                
      t1=exp(-tzg)*t1*2.50662827465d0                                           
      gamz=t1*(dp1+76.18009173d0/(z+dp1)-86.50532033d0/(z+dp2)+                 
     1 24.01409822d0/(z+dp3)-1.231739516d0/(z+dp4)+.120858003d-2/(z+dp5)        
     2 -.536382d-5/(z+6.d0))                                                    
      if (x.le.1) then                                                          
   20   continue                                                                
        gamz=gamz/z                                                             
        adz=abs(z-x)                                                            
        if (adz.ge.dph) then                                                    
          z=z-dp1                                                               
          go to 20                                                              
        endif                                                                   
      endif                                                                     
      gammar=gamz                                                               
      return                                                                    
      end                                                                       
      function gamr (a,z,unew)                                                  
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      add=dp0                                                                   
      jz=z                                                                      
      if ((jz/2)*2.ne.jz) add=.12d0                                             
      if (a-z.le.153.d0) then                                                   
        gamr=dp10**(0.14d0*(a-z)-0.276d0*z+5.460+add)                           
      else                                                                      
        gamr=dp10**(5.d-02*(a-z)-0.276d0*z+19.23d0+add)                         
      endif                                                                     
      return                                                                    
      end                                                                       
      function gaurn (x)                                                        
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
   10 continue                                                                  
      qi=dp0                                                                    
   20 continue                                                                  
      call nombre_generer(rndm)
      q1=rndm                                                                   
      q3=q1                                                                     
   30 continue                                                                  
      call nombre_generer(rndm)
      q2=rndm                                                                   
      if (q3-q2) 60,60,40                                                       
   40 continue                                                                  
      call nombre_generer(rndm)
      q3=rndm                                                                   
      if (q3-q2) 30,50,50                                                       
   50 qi=qi+dp1                                                                 
      go to 20                                                                  
   60 y=q1+qi                                                                   
      qi=dp0                                                                    
   70 continue                                                                  
      call nombre_generer(rndm)
      q1=rndm                                                                   
      q3=q1                                                                     
   80 continue                                                                  
      call nombre_generer(rndm)
      q2=rndm                                                                   
      if (q3-q2) 110,110,90                                                     
   90 continue                                                                  
      call nombre_generer(rndm)
      q3=rndm                                                                   
      if (q3-q2) 80,100,100                                                     
  100 qi=qi+dp1                                                                 
      go to 70                                                                  
  110 z=q1+qi                                                                   
      test=(y-dp1)**2/dp2                                                       
      if (test-z) 120,120,10                                                    
  120 continue                                                                  
      call nombre_generer(rndm)
      fl=rndm                                                                   
      if (fl.ge.dph) go to 130                                                  
      y=-y                                                                      
  130 gaurn=y                                                                   
      return                                                                    
      end                                                                       
       function dgaussn (xmean,sd)                                               
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
c                                                                               
c****************************************************                           
c*                                                  *                           
c*  compute random gaussian number for given        *                           
c*  mean and s.d.                                   *                           
c*  uses mean of sum of 12 uniform r.n"s            *                           
c*                                                  *                           
c****************************************************                           
      a=dp0                                                                     
      n=0                                                                       
   10 continue                                                                  
      a=a+fltrnf(dumm)                                                          
      n=n+1                                                                     
      if (n.lt.12) go to 10                                                     
      dgaussn=(a-6.d0)*sd+xmean                                                  
      return                                                                    
      end                                                                       
      SUBROUTINE generateur_init(IJ,KL)                                         
      INTEGER*2 I97, J97                                                        
      REAL*8 u(97), c, cd, cm                                                   
      LOGICAL*1 test                                                            
                                                       
      COMMON /raset1/ u, c, cd, cm, i97, j97, test                              
      DATA test /.FALSE./                                                                                
      IF( ij .LT. 0  .OR.  ij .GT. 31328  .OR.                                  
     *    kl .lt. 0  .or.  kl .GT. 30081 ) then                                 
          print '(A)', ' The first random number seed must have a value         
     *between 0 and 31328'                                                      
          print '(A)',' The second seed must have a value between 0 and         
     *30081'                                                                    
            STOP                                                                
      END IF                                                                    
                                                                                
      i = MOD(ij/177, 177) + 2                                                  
      j = MOD(ij    , 177) + 2                                                  
      k = MOD(kl/169, 178) + 1                                                  
      l = MOD(kl,     169)                                                      
                                                                                
      DO ii = 1, 97                                                             
         s = 0.0                                                                
         t = 0.5                                                                
         DO jj = 1, 24                                                          
            m = MOD(MOD(i*j, 179)*k, 179)                                       
            i = j                                                               
            j = k                                                               
            k = m                                                               
            l = MOD(53*l+1, 169)                                                
            IF (MOD(l*m, 64) .GE. 32) s = s + t                                 
            t = 0.5 * t                                                         
         END DO                                                                 
         u(ii) = s                                                              
      END DO                                                                    
                                                                                
      c = 362436.0 / 16777216.0                                                 
      cd = 7654321.0 / 16777216.0                                               
      cm = 16777213.0 /16777216.0                                               
                                                                                
      i97 = 97                                                                  
      j97 = 33                                                                  
                                                                                
      test = .TRUE.                                                             
      RETURN                                                                    
      END                                                                       
      function dgeta (e,iz,in,mode,is)                                           
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      logical isz, isn                                                          
      parameter (inn=150, izz=98)                                               
      common /dcook/ sz(izz), sn(inn), con(2), amean(240), pz(izz),              
     1 pn(inn), isz(izz), isn(inn)                                              
      ia=in+iz                                                                  
      aa=ia                                                                     
      if (mode) 10,20,30                                                        
   10 continue                                                                  
ccc   original HETC parametrisation pour alevel
      zz=iz                                                                     
      temp=aa*(dp1+y0*(aa-dp2*zz)**2/aa**2)/b0                                  
      go to 40                                                                  
   20 continue                                                                  
ccc   default- Gilbert-Cameron-Ignatyuk
      if (iz.gt.izz.or.in.gt.inn) go to 30                                      
      if (in.ge.9.and.iz.ge.9) then                                             
        if (isz(iz).or.isn(in)) then                                            
          is=2                                                                  
        else                                                                    
          is=1                                                                  
        endif                                                                   
        st=sz(iz)+sn(in)                                                        
        fa=(9.17d-3*st+con(is))                                                 
      else                                                                      
        if (ia.le.240) then                                                     
          fa=amean(ia)/aa                                                       
        else                                                                    
          fa=(amean(240)+(2.5d0-0.1d0*amean(240))*(real(ia)-240.d0))/aa         
        endif                                                                   
      endif                                                                     
      u=max(5.d-02*(e-pz(iz)-pn(in)),1.d-05)                                    
      gu=(dp1-exp(-u))/u                                                        
      temp=aa*(fa*gu+(dp1-gu)*(0.1375d0-8.36d-05*aa))                           
      go to 40                                                                  
   30 continue                                                                  
ccc   Julich
      if (ia.le.240) then                                                       
         temp=amean(ia)                                                         
      else                                                                      
         temp=(amean(240)+(2.5d0-1.d-1*amean(240))*(real(ia)-240.d0))           
      endif                                                                     
   40 continue                                                                  
      dgeta=temp                                                                 
      return                                                                    
      end                                                                       
      subroutine gms                                                            
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  GMS computes the interaction radius RH and interaction volume                
c     V=(4PI/3)*RH**3 for masses A=1,...,22.                                    
c  It also calculates the statistical factor                                    
c     GAM(N)=(2PI)**(-X)/GAMMA(X), X=(3N-3)/2                                   
c  for N=2,...,LCHAN.                                                           
c                                                                               
c    REFERENCES.......  MASEDT                                                  
c    EXTERNALS........  ALOG      EXP       GAMMAR                              
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      data con /4.18879d0/, pi2 /6.2831852d0/                                   
      f=log(pi2)                                                                
      do 10 n=2,lchan                                                           
      x=n                                                                       
      tem=-dp3*(x-1)/dp2                                                        
      w3=exp(f*tem)                                                             
      w5=gammar(dp3*(x-1)/dp2)                                                  
      gam(n)=w3/w5                                                              
   10 continue                                                                  
      rm=1.08d0*12.d0**dpth+0.73d0                                              
      do 20 i=1,22                                                              
      am=i                                                                      
      if (i.le.4) then                                                          
        r=2.08d0                                                                
      elseif (i.le.11) then                                                     
        r=1.85d0*am**dpth                                                       
        r=max(r,rm)                                                             
      else                                                                      
        r=1.08d0*am**dpth+0.73d0                                                
      endif                                                                     
      rh(i)=r                                                                   
   20 continue                                                                  
      do 30 i=1,22                                                              
      v(i)=con*rh(i)**3                                                         
   30 continue                                                                  
      return                                                                    
      end                                                                       
      subroutine gtiso (u,v,w)                                                  
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
   10 continue                                                                  
      call nombre_generer(rndm)
      z=rndm                                                                    
      x=0.687368d0*sflraf(dumm)                                                 
      y=0.687368d0*sflraf(dumm)                                                 
      xsq=x*x                                                                   
      ysq=y*y                                                                   
      zsq=z*z                                                                   
      d=xsq+ysq+zsq                                                             
      if (d*d-z) 20,20,10                                                       
   20 u=dp2*x*z/d                                                               
      v=dp2*y*z/d                                                               
      w=(zsq-xsq-ysq)/d                                                         
      return                                                                    
      end                                                                       
      subroutine interp (a1)                                                    
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c     the purpose of this subroutine is to interpolate in epperson's            
c     data tables for relative kinetic energy and total excitation energ        
c     scission time - excitation energy of fissioning nucleus                   
c     the latter are supplied by  r.g.alsmiller's fits                          
c                                                                               
c     specifically, delef and ke are computed for a given ef                    
c                                                                               
      common /smacom/ yzero, bzero, yzere, bzere                                
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /orcm2/ ainfo(20), epkin(80,6), sigke(80,6), delexc(80,6,3)        
     1 , cstiff(120)                                                            
cccc                                                                            
c     if kemax.gt.6,change dimension s in exitt                                 
cccc                                                                            
      a2=afnp-a1                                                                
      ia1=a1+dph                                                                
      iaf=afnp+dph                                                              
      ia2=iaf-ia1                                                               
      iahm=afnp/dp2                                                             
      iahm=iaf-iahm                                                             
      if (ia1.ge.iahm) then                                                     
        iatm=ia2                                                                
        ia2=ia1                                                                 
        ia1=iatm                                                                
      endif                                                                     
      iaa=ia2-100                                                               
      delef=dp0                                                                 
      j1=kemax+4                                                                
      j2=j1+kemax-2                                                             
      do 10 j=j1,j2                                                             
      if (ef.gt.ainfo(j)) go to 10                                              
      ke=j-kemax-3                                                              
      go to 20                                                                  
   10 continue                                                                  
      ke=kemax                                                                  
      delef=dp0                                                                 
      delkex=sigke(iaa,ke)                                                      
      epkinx=epkin(iaa,ke)                                                      
      e1maxp=delexc(iaa,ke,ismal)                                               
      ern=ainfo(kemax+3)                                                        
      go to 30                                                                  
   20 continue                                                                  
      delef=ainfo(j+1)-ainfo(j)                                                 
      delef=(ef-ainfo(j))/delef                                                 
      if (delef.lt.dp0) delef=dp0                                               
      delkex=(dp1-delef)*sigke(iaa,ke)+delef*sigke(iaa,ke+1)                    
      epkinx=(dp1-delef)*epkin(iaa,ke)+delef*epkin(iaa,ke+1)                    
      e1maxp=delexc(iaa,ke,ismal)*(dp1-delef)+delexc(iaa,ke+1,ismal)            
     1 *delef                                                                   
cccc                                                                            
c     put kinetic energies of neptunium(236,93) compound nucleus in ainf        
c     ainfo(kemax+3)                                                            
cccc                                                                            
      ern=ainfo(ke+3)*(dp1-delef)+ainfo(ke+4)*delef                             
   30 continue                                                                  
      delken=delkex                                                             
      e1maxp=e1maxp+ef                                                          
      emaxn=e1maxp                                                              
      delke=scalke*delken                                                       
      ekinr4=epkinx-ern                                                         
      z2p=zprobf(afnp,zfnp,a2)                                                  
      z1p=zfnp-z2p                                                              
      c=(a1-z1p*dp2)/a1                                                         
      smala1=a1*(dp1+yzero*c*c)/bzero                                           
      c=(a2-z2p*dp2)/a2                                                         
      smala2=a2*(dp1+yzero*c*c)/bzero                                           
      smasum=smala1+smala2                                                      
      arg=smasum*emaxn                                                          
cccc                                                                            
c     return delken and arg in common to qmaxx later                            
cccc                                                                            
      return                                                                    
      end                                                                       
      subroutine keave (a2np,a2,smasum)                                         
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /smacom/ yzero, bzero, yzere, bzere                                
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /orcm2/ ainfo(20), epkin(80,6), sigke(80,6), delexc(80,6,3)        
     1 , cstiff(120)                                                            
cccc                                                                            
c     double interpolate in ef and a for neptunium                              
cccc                                                                            
      ia2=a2np+dph                                                              
      ekinh=dp0                                                                 
      k2=ke+1                                                                   
      ekinl=epkin(ia2-100,ke)*(dp1-dela)+epkin(ia2-99,ke)*dela                  
      if (ke.ne.kemax) ekinh=epkin(ia2-100,k2)*(dp1-dela)+epkin(ia2-99          
     1 ,k2)*dela                                                                
      ekinr4=ekinl*(dp1-delef)+ekinh*delef-ern                                  
      z2pnp=zprobf(afnp,zfnp,a2np)                                              
      z1pnp=zfnp-z2pnp                                                          
      a1np=afnp-a2np                                                            
      c=(a1np-z1pnp*dp2)/a1np                                                   
      smala1=a1np*(dp1+yzero*c*c)/bzero                                         
      c=(a2np-z2pnp*dp2)/a2np                                                   
      smala2=a2np*(dp1+yzero*c*c)/bzero                                         
      smasum=smala1+smala2                                                      
      z1p=z1pnp                                                                 
      z2p=z2pnp                                                                 
      return                                                                    
      end                                                                       
      subroutine levbrk (ipt,llev,n,kip,ilev,t)                                 
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c///////////////////////////////////////////////////////////////////////        
c  FBRLEV selects a breakup channel for a specific level.                       
c///////////////////////////////////////////////////////////////////////        
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      parameter (nmx=3)                                                         
      common /fbrl/ q1(mch5), wf(mch5), il0(lstat), n1(mch5),                   
     1 klev1(mch5), izap1(nmx,mch5)                                             
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
      dimension kip(nmx), ilev(nmx)                                             
      l=il0(jns(ipt)+llev)                                                      
      ki=nklev(jns(ipt)+llev)                                                   
      if (ki.le.0) then                                                         
        n=1                                                                     
        kip(1)=ipt                                                              
        ilev(1)=-1                                                              
        t=dp0                                                                   
        return                                                                  
      endif                                                                     
      nlo=l+1                                                                   
      nup=nlo+ki-1                                                              
      wtot=wf(nup)                                                              
      call nombre_generer(rndm)
      fran=wtot*rndm                                                            
      do 10 kk=nlo,nup                                                          
      if (wf(kk).gt.fran) go to 20                                              
   10 continue                                                                  
      kk=nup                                                                    
c                                                                               
c       decay channel selected                                                  
c                                                                               
   20 continue                                                                  
      t=q1(kk)                                                                  
      n=n1(kk)                                                                  
      do 30 j=1,n                                                               
      kip(j)=izap1(j,kk)                                                        
      ilev(j)=1                                                                 
   30 continue                                                                  
      jjix=klev1(kk)                                                            
      ix=mod(jjix,3)+1                                                          
      ilev(ix)=jjix/3+1                                                         
      return                                                                    
      end                                                                       
      subroutine levset (ipt,llev)                                              
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c  LEVSET calculated the branching probabilities for each nuclear               
c  level in the data library.                                                   
c                                                                               
c    REFERENCES.......  MOMENT    MASEDT                                        
c    EXTERNALS........  AMASS     CH2       CH3       CH4       EXITA           
c                       SQRT      VL        EXP       CLMB      ALOG            
c                       RANF                                                    
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      parameter (nmx=3)                                                         
      common /fbrl/ q1(mch5), wf(mch5), il0(lstat), n1(mch5),                   
     1 klev1(mch5), izap1(nmx,mch5)                                             
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
c      dimension kip(nmx), s(nmx), nss(nmx), inss(nmx), ilev(nmx), jlev  
      dimension  s(nmx), nss(nmx), inss(nmx),  jlev         
     1 (nmx), jpp(nmx), e0(nmx), fact(lchan), npp(nmx), ict(nmx), ef(nmx        
     2 ), wtl1(lchan)                                                           
      common /qfbl/ loff                                                        
      data lo3 /1/                                                              
      data fact /dp1,dp2,6.d0/                                                  
      data hbc3 /7.683545d+06/, con1 /0.72d0/, con2 /.00516d0/, con3            
     1 /.007168d0/, thibit /dp0/                                                
c                                                                               
c scan levels for instability                                                   
c                                                                               
      izr=jz(ipt)                                                               
      iar=ja(ipt)                                                               
      estar=e(jns(ipt)+llev)                                                    
      ti0=tiso(jns(ipt)+llev)                                                   
      pi=parity(jns(ipt)+llev)                                                  
      ajj0=aj(jns(ipt)+llev)                                                    
      emax=estar+1.d-2                                                          
      fzr=izr                                                                   
      k=0                                                                       
      fac=v(iar)/hbc3                                                           
      amr=amass(izr,iar,iok)                                                    
      if (iok.eq.0) then                                                        
        iok1=0                                                                  
        return                                                                  
      endif                                                                     
c                                                                               
c       first loop over no. of particles in breakup                             
c                                                                               
      nop=nmx                                                                   
      w2=dp1                                                                    
      if (nop.gt.iar) nop=iar                                                   
      wtot=dp0                                                                  
      kj=1                                                                      
      do 200 n=2,nop                                                            
      w2=w2*fac                                                                 
      w235=dph*w2*gam(n)                                                        
c                                                                               
c       find all channels with thes number of particles                         
c                                                                               
      if (n.le.2) then                                                          
        call ch2 (nchan)                                                        
C assign modified for gfortran compilator : AB 5/2009
C        assign 160 to jbr                                                       
C        assign 60 to jbr2
        jbr=160
	jbr2=60                                                       
      else                                                                      
        call ch3 (nchan)                                                        
c        assign 160 to jbr                                                       
c        assign 40 to jbr2                                                       
        jbr=160
	jbr2=40                                                       
      endif                                                                     
      if (nchan.le.0) go to 200                                                 
c                                                                               
c       loop over channels                                                      
c                                                                               
      do 190 i=1,nchan                                                          
c                                                                               
c       loop over no. of particles in channel                                   
c                                                                               
      q=amr                                                                     
      sum=dp0                                                                   
      prod=dp1                                                                  
      do 10 np=1,n                                                              
      iip=izap(np,i)                                                            
      jpp(np)=iip                                                               
      nss(np)=ns(iip)                                                           
      inss(np)=jns(iip)                                                         
      amm=axm(iip)                                                              
      q=q-amm                                                                   
      sum=sum+amm                                                               
      prod=prod*amm                                                             
   10 continue                                                                  
      tem4=emax+q                                                               
      if (tem4.le.dp0) go to 190                                                
      r=prod/sum                                                                
      facm=r*sqrt(r)                                                            
      ww=w235*facm                                                              
      couly=dp0                                                                 
      if (n.eq.2) then                                                          
        tem0=xminv(jpp(1))+xminv(jpp(2))                                        
        tem1=dp1/sqrt(tem0)                                                     
        iip1=jpp(1)                                                             
        iip2=jpp(2)                                                             
        couly=con3*(rho(iip1)+rho(iip2))*tem1                                   
      endif                                                                     
c                                                                               
      coulx=dp0                                                                 
      if (jpp(n-1).eq.1) go to 30                                               
      if (n.eq.2) then                                                          
        coulx=con2*z(iip1)*z(iip2)*tem1                                         
      else                                                                      
        tem=dp0                                                                 
        do 20 np=1,n                                                            
        iip=jpp(np)                                                             
        if (iip.gt.1) then                                                      
          zz=z(iip)                                                             
          iz1=izr-jz(iip)+1                                                     
          if (iz1.gt.0) then                                                    
            in1=iar-ja(iip)-iz1+2                                               
            iip1=ip(iz1,in1)                                                    
            rho2=rho(iip)+rho(iip1)                                             
            tem=tem+zz*(fzr-zz)/rho2                                            
          endif                                                                 
        endif                                                                   
   20   continue                                                                
        coulx=tem*con1                                                          
      endif                                                                     
   30 continue                                                                  
c                                                                               
c       loop over excited states for a given breakup_                           
c                                                                               
      ns1=nss(1)                                                                
      ns2=nss(2)                                                                
      lof1=inss(1)                                                              
      lof2=inss(2)                                                              
      do 180 j1=1,ns1                                                           
      if (jpp(1).eq.jpp(2).and.j1.gt.1) go to 180                               
      jlev(1)=j1                                                                
      e0(1)=e(lof1+j1)                                                          
      ajj1=aj(lof1+j1)                                                          
      ti1=tiso(lof1+j1)                                                         
      pi1=pi*parity(lof1+j1)                                                    
      s(1)=dp2*ajj1+dp1                                                         
      wtl1(1)=wtl(lof1+j1)                                                      
      lo2=1                                                                     
      ns2=1                                                                     
      if (j1.eq.1) ns2=nss(2)                                                   
      do 170 j2=lo2,ns2                                                         
      jlev(3)=1                                                                 
      e0(2)=e(lof2+j2)                                                          
      ef(2)=e0(1)+e0(2)                                                         
      if (ef(2).ge.tem4) go to 170                                              
      jlev(2)=j2                                                                
      ajj2=aj(lof2+j2)                                                          
      ti2=tiso(lof2+j2)                                                         
      pi2=pi1*parity(lof2+j2)                                                   
      s(2)=s(1)*(ajj2+ajj2+dp1)                                                 
      wtl1(2)=wtl(lof2+j2)                                                      
      j3=1                                                                      
C GO TO modified for gfortran : AB 5/2009
C      go to jbr2, (60,40)
      IF(jbr2.EQ.60) GO TO 60                                                       
      IF(jbr2.EQ.40) GO TO 40                                                       
   40 continue                                                                  
      if (j1.eq.1.and.j2.eq.1) j3=nss(3)                                        
      if (jpp(3).eq.jpp(2).and.j2.gt.1) go to 170                               
      lof3=inss(3)                                                              
   50 continue                                                                  
      e0(3)=e(lof3+j3)                                                          
      ef(3)=ef(2)+e0(3)                                                         
      if (ef(3).ge.tem4) go to 160                                              
      jlev(3)=j3                                                                
      ajj3=aj(lof3+j3)                                                          
      ti3=tiso(lof3+j3)                                                         
      pi3=pi2*parity(lof3+j3)                                                   
      s(3)=s(2)*(ajj3+ajj3+dp1)                                                 
      wtl1(3)=wtl(lof3+j3)                                                      
   60 continue                                                                  
      tem3=q-ef(n)                                                              
      t0=estar+tem3                                                             
      if (t0.lt.dp0) go to 150                                                  
      k=k+1                                                                     
      q01=tem3                                                                  
      gg=dp1                                                                    
      ww1=ww                                                                    
      if (n.lt.3) then                                                          
        wtl0=wtl1(1)*wtl1(2)                                                    
        tmax=ti1+ti2                                                            
        if (ti0.gt.tmax) ww1=ww1*thibit                                         
        tmin=abs(ti1-ti2)                                                       
        if (ti0.lt.tmin) ww1=ww1*thibit                                         
        lmin=0                                                                  
        ajmax=ajj1+ajj2                                                         
        if (ajj0.gt.ajmax) then                                                 
          lmin=ajj0-ajmax+1.d-3                                                 
        else                                                                    
          ajmin=abs(ajj1-ajj2)                                                  
          if (ajj0.lt.ajmin) lmin=ajmin-ajj0+1.d-3                              
        endif                                                                   
        if (lmin.gt.0) pi2=pi2*(-dp1)**lmin                                     
        if (pi2.lt.dp0) lmin=lmin+1                                             
        lmax=ajj0+ajmax                                                         
        if (lmin.gt.lmax) then                                                  
          k=k-1                                                                 
          go to 150                                                             
        endif                                                                   
        ww1=ww1*vl(t0*couly**2,lmin)                                            
        if (jpp(1).eq.jpp(2).and.e0(1).eq.e0(2)) gg=dp2                         
        go to 110                                                               
      endif                                                                     
      ajmax=ajj1+ajj2+ajj3                                                      
      lmax=ajj0+ajmax                                                           
      lmin=0                                                                    
      if (ajj0.ge.ajmax) then                                                   
        lmin=ajj0-ajmax                                                         
      else                                                                      
        if (ajj1-ajj2-ajj3.ge.dp0) then                                         
          ajmin=ajj1-ajj2-ajj3                                                  
        elseif (abs(ajj2-ajj3)-ajj1.ge.dp0) then                                
          ajmin=abs(ajj2-ajj3)-ajj1                                             
        else                                                                    
          ajmin=mod(ajmax,dp1)                                                  
        endif                                                                   
        if (ajj0.le.ajmin) lmin=ajmin-ajj0                                      
      endif                                                                     
      if (lmin.gt.0) pi3=pi3*(-dp1)**lmin                                       
      if (pi3.lt.dp0) lmin=lmin+1                                               
      if (lmin.gt.lmax) then                                                    
        k=k-1                                                                   
        go to 150                                                               
      endif                                                                     
      tmax=ti1+ti2+ti3                                                          
      if (ti0.gt.tmax) ww1=ww1*thibit                                           
      t1m=max(ti1,ti2,ti3)                                                      
      tmin=max(dp0,dp2*t1m-tmax)                                                
      if (ti0.lt.tmin) ww1=ww1*thibit                                           
      nkind=0                                                                   
      wtl0=dp1                                                                  
      do 70 iii=1,n                                                             
      ict(iii)=jpp(iii)                                                         
      wtl0=wtl0*wtl1(iii)                                                       
   70 continue                                                                  
      do 90 iii=1,n                                                             
c                                                                               
c       check that this particle has not                                        
c       already been counted                                                    
c                                                                               
      if (ict(iii).eq.0) go to 90                                               
      nkind=nkind+1                                                             
      npp(nkind)=1                                                              
      if (iii.eq.n) go to 90                                                    
      et=e0(iii)                                                                
      ica=ict(iii)                                                              
c                                                                               
c       now scan the rest of the particles for a match                          
c                                                                               
      i1=iii+1                                                                  
      do 80 ii=i1,n                                                             
      if (ict(ii).gt.ica) go to 90                                              
      if (ict(ii).eq.ica.and.e0(ii).eq.et) then                                 
        npp(nkind)=npp(nkind)+1                                                 
        ict(ii)=0                                                               
      endif                                                                     
   80 continue                                                                  
   90 continue                                                                  
c                                                                               
c       compute g factor                                                        
c                                                                               
      if (n-nkind.gt.1) then                                                    
        do 100 iii=1,nkind                                                      
        npi=npp(iii)                                                            
  100   gg=gg*fact(npi)                                                         
      elseif (n-nkind.eq.1) then                                                
        gg=dp2                                                                  
      endif                                                                     
  110 continue                                                                  
      if (ww1.le.dp0) then                                                      
        k=k-1                                                                   
        go to 150                                                               
      endif                                                                     
      w00=ww1*s(n)*wtl0/gg                                                      
      m0=n                                                                      
      if (m0.gt.2) texp=1.5d0*m0-2.5d0                                          
      if (coulx.eq.dp0) go to 130                                               
      if (m0.eq.2) go to 120                                                    
      tem2=coulx/t0                                                             
      if (tem2.gt.100.d0) tem2=100.d0                                           
      w00=w00*exp(-tem2)                                                        
      go to 130                                                                 
  120 continue                                                                  
      tem2=sqrt(t0)                                                             
      eta=coulx/tem2                                                            
      rhox=couly*tem2                                                           
      w00=w00*max(1.0d-10,dclmb1(rhox,eta,ml))                                   
  130 continue                                                                  
      if (m0.eq.2) then                                                         
        wtot=wtot+w00*sqrt(t0)                                                  
      else                                                                      
        wtot=wtot+w00*exp(texp*log(t0))                                         
      endif                                                                     
      jj=jlev(1)+jlev(2)+jlev(3)-3                                              
      ix=1                                                                      
      wf(loff+kj)=wtot                                                          
      q1(loff+kj)=q01+estar                                                     
      n1(loff+kj)=m0                                                            
      do 140 l=1,m0                                                             
      izap1(l,loff+kj)=jpp(l)                                                   
      if (jlev(l).gt.1) ix=l                                                    
  140 continue                                                                  
      jjix=3*jj+ix-1                                                            
      klev1(loff+kj)=jjix                                                       
      kj=kj+1                                                                   
  150 continue                                                                  
  160 j3=j3-1                                                                   
      if (j3.ge.lo3) go to 50                                                   
  170 continue                                                                  
  180 continue                                                                  
  190 continue                                                                  
  200 continue                                                                  
      if (k.eq.0.or.kj.eq.1) then                                               
c                                                                               
c       kj=1; res.nuc. does not have enough ex.en. to breakup                   
c                                                                               
        return                                                                  
      endif                                                                     
      if (k.eq.0) kj=1                                                          
      ki=kj-1                                                                   
      if ((loff+ki).gt.mch5) go to 220                                          
      il0(jns(ipt)+llev)=loff                                                   
      nklev(jns(ipt)+llev)=ki                                                   
      if (ki.le.0) go to 210                                                    
      loff=loff+ki                                                              
  210 continue                                                                  
      return                                                                    
  220 continue                                                                  
      call errorf ('levset')                                                    
      end                                                                       
      subroutine masedt (ixedt)                                                 
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c.......................................................................        
c  MASEDT initializes the channel probability arrays by calling                 
c  CBRSET and LEVSET.  When called with a non-zero arguement, it will           
c  print out the nuclear data library with branching ratios for each            
c  level.                                                                       
c.......................................................................        
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c    REFERENCES.......  REDMAS                                                  
c    EXTERNALS........  GMS       LEVBRK    EXITA                               
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      parameter (nmx=3)                                                         
      common /fbrl/ q1(mch5), wf(mch5), il0(lstat), n1(mch5),                   
     1 klev1(mch5), izap1(nmx,mch5)                                             
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      dimension ja3(nmx), il3(nmx), jz3(nmx)                                    
      dimension est0(20)                                                        
c.......................................................................        
      data est0 /13*100,7*40/                                                   
c The definition of est0 allows all possible channels for masses le. 13,        
c but only those open at 40 MeV for masses gt. 13.                              
c.......................................................................        
      call gms                                                                  
      nmax=lchan                                                                
      if (ixedt.ne.0) write (lout,100)                                          
      do 60 ia=1,itop                                                           
      iup=min(ia+1,11)                                                          
      do 50 iz1=1,iup                                                           
      iz=iz1-1                                                                  
      in=ia-iz                                                                  
      in1=in+1                                                                  
      if (in1.gt.13.or.in1.lt.1) go to 50                                       
      ipt=ip(iz1,in1)                                                           
      if (ipt.eq.0) then                                                        
        if (ixedt.ne.0) write (lout,90) iz,ia                                   
        go to 50                                                                
      endif                                                                     
      est=est0(ia)                                                              
      ira=ia                                                                    
      irz=iz                                                                    
      call cbrset (iz,ia,est,nmax)                                              
      nss=ns(ipt)                                                               
      if (ixedt.ne.0) write (lout,70) iz,ia,xm(ipt),nss,rho(ipt),rh(ia)         
      if (nss.eq.0.or.ia.eq.itop) go to 50                                      
      lo=jns(ipt)+1                                                             
      lhi=lo+nss-1                                                              
      do 40 i=lo,lhi                                                            
      if (nklev(i).ne.0) call errorf ('masedt')                                 
      kk=i-lo+1                                                                 
      call levset (ipt,kk)                                                      
      nk=nklev(i)                                                               
      if (nk.le.0) go to 30                                                     
      nlo=il0(i)+1                                                              
      nup=nlo+nk-1                                                              
      wtot=wf(nup)                                                              
      wtl(i)=unstab                                                             
      if (ixedt.eq.0) go to 40                                                  
      write (lout,110) e(i),aj(i),parity(i),tiso(i),wtl(i),nk                   
      if (ixedt.ne.2) go to 40                                                  
      prob=dp1                                                                  
      do 20 k=nlo,nup                                                           
      if (nk.gt.1) then                                                         
        if (k.eq.nlo) prob=wf(k)/wtot                                           
        if (k.gt.nlo) prob=(wf(k)-wf(k-1))/wtot                                 
      endif                                                                     
      q=q1(k)                                                                   
      n=n1(k)                                                                   
      do 10 l=1,n                                                               
      iiip=izap1(l,k)                                                           
      ja3(l)=ja(iiip)                                                           
      jz3(l)=jz(iiip)                                                           
      il3(l)=1                                                                  
   10 continue                                                                  
      jjix=klev1(k)                                                             
      ix=mod(jjix,3)+1                                                          
      il3(ix)=jjix/3+1                                                          
      write (lout,120) prob,q,(jz3(l),ja3(l),il3(l),l=1,n)                      
   20 continue                                                                  
      go to 40                                                                  
   30 continue                                                                  
      wtl(i)=dp1                                                                
      if (ixedt.ne.0) write (lout,80) e(i),aj(i),parity(i),tiso(i),wtl(i        
     1 )                                                                        
   40 continue                                                                  
   50 continue                                                                  
   60 continue                                                                  
      return                                                                    
c                                                                               
   70 format (/5x,'z =',i3,'    a =',i3,'    mass excess =',f10.5,'    n        
     1umber of states =',i3,'    r(coul) =',f6.2,'    r(int) =',f6.2)           
   80 format (10x,'ex =',f10.5,'  j =',f4.1,'  pi = ',f4.1,'  t =',f4.1,        
     1 '  wtl =',g11.5)                                                         
   90 format (/5x,'z =',i3,'    a =',i3,'    not tabulated')                    
  100 format (1h1,25x,'edit of mass data input'//)                              
  110 format (10x,'ex =',f10.5,'  j =',f4.1,'  pi = ',f4.1,'  t =',f4.1,        
     1 '  wtl =',g11.5,'   unstable -',i4,' channels')                          
  120 format (15x,'prob=',g12.5,'  q=',g12.5,4(' -- ',3i3))                     
      end                                                                       
      subroutine moment (izr,iar,uu,vv,ww,ex,erec,ierr)                         
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c.......................................................................        
c These modifications imply the use of real masses in MOMENT.                   
c The arguements to MOMENT are Z and A of the residual nucleus, the             
c direction cosines U, V, and W of its recoil momentum, the excitation          
c energy EX, and the recoil energy EREC.                                        
c                                                                               
c The flag IERR is 0 for a normal execution, +1 if the mass of the              
c residual is not in the data tables, and is -1 if the residual is              
c particle  stable.                                                             
c.......................................................................        
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c    REFERENCES.......  FBRK                                                    
c    EXTERNALS........  RANF      BETAD     SQRT      SIN       COS             
c                       ACOS      SPECTR    LEVBRK    EXITA                     
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      parameter (rparm1=1.53960072)                                             
      real *8rani, ranj, rijk, rijk0, ranb, rans                                
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      parameter (mch1=13370, mch6=690, lpart=10, lchan=3, mch5=1520,            
     1 mch2=lpart*(lpart+1)/2, lmas=119, lmas1=70, lstat=605)                   
      common /fbrc1/ coulch(mch1), coulck(mch1), q0(mch1), w0(mch1),            
     1 fz(lchan), fa(lchan), flz(lchan), fla(lchan), gam(lchan),                
     2 wff(mch6), ipdat(lmas), npdat(lmas), izaps(mch1),                        
     3 itr(mch6), izap(lchan,mch2)                                              
      common /fbrc0/ atop, z(lmas), a(lmas), e(lstat), aj(lstat),               
     1 wtl(lstat), tiso(lstat), parity(lstat), unstab, xm(lmas),                
     2 axm(lmas), rho(lmas), xminv(lmas), v(22), rh(22), nklev(lstat),          
     3 jns(lmas1), ip(11,13), lout, lmass, iok1, itop, nr,                      
     4 nrm, irz, ira, jz(lmas), ja(lmas), ns(lmas)                              
      common /output/ exi(20), vp(3,20), erc(20), npt, jjp(20)                  
      dimension ilev(lchan), ill1(lchan), iip1(lchan), vx1(lchan), vy1          
     1 (lchan), vz1(lchan), kip(lchan)                                          
      ierr=0                                                                    
      aaaa=amass(izr,iar,iok)                                                   
      if (iok.gt.0) go to 10                                                    
      ierr=1                                                                    
      return                                                                    
   10 continue                                                                  
      ipt0=iok                                                                  
c      call cbreak (ipt0,ex,n,kip,ilev,t,ki,iok1)                                
      call cbreak (ipt0,ex,n,kip,ilev,t,ki)
      if (n.gt.1) go to 20                                                      
      ierr=-1                                                                   
      return                                                                    
   20 continue                                                                  
      vs=sqrt(dp2*erec/aaaa)                                                    
c.......................................................................        
c this modification is intended to conserve energy in a primary breakup         
c.......................................................................        
      amt=dp0                                                                   
      do 30 i=1,n                                                               
   30 amt=amt+axm(kip(i))                                                       
      vs=vs*sqrt(aaaa/amt)                                                      
c.......................................................................        
      vx=vs*uu                                                                  
      vy=vs*vv                                                                  
      vz=vs*ww                                                                  
      ntot=n                                                                    
      nextra=0                                                                  
      nt=1                                                                      
      npt=0                                                                     
      izr1=izr                                                                  
      iar1=iar                                                                  
   40 continue                                                                  
      ar=axm(kip(1))+axm(kip(2))                                                
      if (n.gt.2) ar=ar+axm(kip(3))                                             
      m1=n-1                                                                    
      do 50 n2=1,m1                                                             
      ipt=kip(n2)                                                               
      iap=ja(ipt)                                                               
      izp=jz(ipt)                                                               
      iar=iar-iap                                                               
      izr=izr-izp                                                               
      ap=axm(ipt)                                                               
      ar=ar-ap                                                                  
c                                                                               
c       iar,izr are the mass and charge of resid.nuc.                           
c       after one of the particles emitted in breakup                           
c                                                                               
      rn2=fltrnf(dumm)                                                          
      call azirn (sinp,cosp)
                          
      costh=dp1-dp2*rn2                                                         
      sox=dp2*ap*t/(ar*(ar+ap))                                                 
      np=n+1-n2                                                                 
      xsq=dp1                                                                   
      if (np.gt.2) xsq=betad(np)                                                
      sox=sox*xsq                                                               
      t=t*(dp1-xsq)                                                             
      vt=sqrt(sox)                                                              
c                                                                               
c       vt is the recoil vel. of resid.nuc.                                     
c       t is the total k.e. of res.nuc.in c.m. system                           
c                                                                               
      sinsq=dp1-costh*costh                                                     
      sinth=sqrt(sinsq)                                                         
      vzse=vt*costh                                                             
      vyse=vt*sinth*sinp                                                        
      vxse=vt*sinth*cosp 
      
c      WRITE(6,*) 'costh,sinth,cosp,sinp:',costh,sinth,cosp,sinp                                                       
c                                                                               
c     vzse,vyse,vxse are recoil vel. components of resid.nuc.                   
c                                                                               
      vz=vz+vzse                                                                
      vy=vy+vyse                                                                
      vx=vx+vxse                                                                
c                                                                               
c     vz,vy,vx are vel.comp. of resid.nuc.                                      
c                                                                               
      v1=ar/ap+dp1                                                              
      vzp=vz-vzse*v1                                                            
      vyp=vy-vyse*v1                                                            
      vxp=vx-vxse*v1                                                            
c                                                                               
c     vzp,vyp,vxp are the vel.comp. of emitted particle in breakup              
c                                                                               
      llev=ilev(n2)                                                             
      if (nklev(jns(ipt)+llev).ne.0) then                                       
        nextra=nextra+1                                                         
        nt=nt+1                                                                 
        iip1(nextra)=ipt                                                        
        ill1(nextra)=llev                                                       
        vx1(nextra)=vxp                                                         
        vy1(nextra)=vyp                                                         
        vz1(nextra)=vzp                                                         
      else                                                                      
        npt=npt+1                                                               
        jjp(npt)=ipt                                                            
        exi(npt)=e(jns(ipt)+llev)                                               
        vp(1,npt)=vxp                                                           
        vp(2,npt)=vyp                                                           
        vp(3,npt)=vzp                                                           
      endif 
c      WRITE(6,*) 'in moment1: npt,vp:',npt,(vp(iloc,npt),iloc=1,3)
   50 continue                                                                  
      ipt=kip(n)                                                                
      if (n.gt.1) then                                                           
        llev=ilev(n)                                                            
        if (nklev(jns(ipt)+llev).ne.0) then                                     
          nextra=nextra+1                                                       
          nt=nt+1                                                               
          iip1(nextra)=ipt                                                      
          ill1(nextra)=llev                                                     
          vx1(nextra)=vx                                                        
          vy1(nextra)=vy                                                        
          vz1(nextra)=vz                                                        
          go to 60                                                              
        endif                                                                   
        npt=npt+1                                                               
        jjp(npt)=ipt                                                            
        exi(npt)=e(jns(ipt)+llev)                                               
        vp(1,npt)=vx                                                            
        vp(2,npt)=vy                                                            
        vp(3,npt)=vz                                                            
      endif                                                                     
c      WRITE(6,*) 'in moment2: npt,vp:',npt,(vp(iloc,npt),iloc=1,3)
   60 continue                                                                  
      call levbrk (ipt,llev,n,kip,ilev,t)                                       
      if (nextra.ne.0) then                                                     
        nk=nextra                                                               
        nextra=nextra-1                                                         
        ipt=iip1(nk)                                                            
        llev=ill1(nk)                                                           
        vx=vx1(nk)                                                              
        vy=vy1(nk)                                                              
        vz=vz1(nk)                                                              
        iar=ja(ipt)                                                             
        izr=jz(ipt)                                                             
c        call levbrk (ipt,llev,n,kip,ilev,t,0)                                   
        call levbrk (ipt,llev,n,kip,ilev,t)
        if (n.lt.2) call errorf ('moment.1')                                    
c.......................................................................        
c this modification conserves energy in a secondary breakup                     
c.......................................................................        
        a0=axm(ipt)                                                             
        amt=axm(kip(1))+axm(kip(2))                                             
        if (n.gt.2) amt=amt+axm(kip(3))                                         
        fact=sqrt(a0/amt)                                                       
        vx=vx*fact                                                              
        vy=vy*fact                                                              
        vz=vz*fact                                                              
c.......................................................................        
        ntot=ntot-1+n                                                           
        go to 40                                                                
      endif                                                                     
      if (ntot.gt.1.and.npt.ne.ntot) call errorf ('moment.2')                   
      izr=izr1                                                                  
      iar=iar1                                                                  
      do 80 i=1,npt                                                             
      vsq=vp(1,i)**2+vp(2,i)**2+vp(3,i)**2                                      
      erc(i)=dph*axm(jjp(i))*vsq                                                
      vsq=sqrt(vsq)                                                             
      do 70 j=1,3                                                               
   70 vp(j,i)=vp(j,i)/vsq                                                       
   80 continue                                                                  
      return                                                                    
c                                                                               
      end                                                                       
      SUBROUTINE nombre_generer(nombre)                                         
      REAL*8 nombre                                                             
      REAL*8 u(97), c, cd, cm                                                   
      INTEGER*2 i97, j97                                                        
      LOGICAL*1 test                                                            
      common /raset1/ u, c, cd, cm, i97, j97, test                              
                                                                                
       nombre = u(I97) - u(J97)                                                 
       IF( nombre .LT. 0.0 ) nombre = nombre + 1.0                              
                                                                                
       U(I97) = nombre                                                          
       I97 = I97 - 1                                                            
       IF(I97 .EQ. 0) I97 = 97                                                  
                                                                                
       J97 = J97 - 1                                                            
       IF(J97 .EQ. 0) J97 = 97                                                  
                                                                                
        C = C - CD                                                              
        IF( C .LT. 0.0 )  C = C + CM                                            
                                                                                
        nombre = nombre - C                                                     
        IF( nombre .LT. 0.0 ) nombre = nombre + 1.0                             
      return                                                                    
      END                                                                       
      function pairr (a1,izz1,a2,izz2,corr1,corr2)                              
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
cccc                                                                            
c     choose nearest integers                                                   
cccc                                                                            
      iaf=af+dph                                                                
      ia2=a2+dph                                                                
      ia1=iaf-ia2                                                               
      inn1=ia1-izz1                                                             
      inn2=ia2-izz2                                                             
      corr1=cam4(izz1)+cam5(inn1)                                               
      corr2=cam4(izz2)+cam5(inn2)                                               
      pairr=corr1+corr2                                                         
      return                                                                    
      end                                                                       
      subroutine qmaxx (a,qa,q)                                                 
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /evcml/ cam2(130), cam3(200), waps(250,20)                         
      common /evcm/ y0, b0, t(4,7)                                              
      common /orcm0/ af, zf, ef, ernf, e1, e2, r1mass, r2mass, z91,             
     1 ucut, ncall                                                              
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /inout/ lu1, lu2, itty, iscrt                                      
cccc                                                                            
c     omegao = sqrt(dp4*ef*af/b0)                                               
cccc                                                                            
      delke=delken                                                              
      qa=dp0                                                                    
      q=dp1                                                                     
      a2=a                                                                      
cccc                                                                            
c     a1=af-a2                                                                  
cccc                                                                            
      a1=afnp-a2                                                                
      qnorm=dp1                                                                 
      fact=sqrt(arg)                                                            
      smasum=smala1+smala2                                                      
      if (iversn.ne.2) then                                                     
        fn=2.375d0                                                              
        delke=sqrt(z1p*z2p/smasum)                                              
        delke=delke*sqrt(fact)                                                  
        q=q/((a1+a2)*(a1**1.66666667d0+a2**1.66666667d0)**1.5d0)                
      else                                                                      
        fn=3.625d0                                                              
cccc                                                                            
c        this change was put in to give new  relativistic wke                   
cccc                                                                            
        fn=fm                                                                   
        q=q*sqrt(ekinr4*arg)/smasum                                             
        q=q/(a1+a2)                                                             
      endif                                                                     
      q3=exp(dp2*fact-omegao)                                                   
      crit=dp1-fn/fact                                                          
      if (crit.le.dp0) then                                                     
        write (lu2,10) crit,fn,fact,iversn,arg,qnorm,delke                      
        return                                                                  
      endif                                                                     
      q2=crit*arg**2.5d0                                                        
      q=q/sqrt(baz*(a1+a2))                                                     
      q5=qnorm*(a1*a2)**4*sqrt(smala1*smala2)/smasum**3                         
      q=q*q5                                                                    
      q=q/smasum**2.5d0                                                         
      q5=q*q2*q3                                                                
      qa=q5*delke                                                               
      if (qa.le.dp0) write (lu2,20) qa,q5,delke,q,q2,q3,qnorm,baz,a1,a2         
     1 ,smala1,smala2,arg,crit,fact,fn,omegao,ekinr4,delken,z1p,z2p,af          
     2 ,zf,ef,afnp,zfnp,y0,b0,e1maxp                                            
      return                                                                    
c                                                                               
   10 format (1h ,'crit=',1pe15.6,2x,1p2e15.6,2x,i10,3e15.6)                    
   20 format (1h1,'qmaxx has qamax le. 0.0',2x//(1x,1p8e15.6,2x))               
      end                                                                       
      function sflraf (x)                                                       
      implicit doubleprecision(a-h,o-z)                                         
      parameter (lpt=18)
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c        return the pseudo-random number -1 < r < +1                            
      parameter (rparm1=1.53960072)                                             
      real*8 rani, ranj, rijk, rijk0, ranb, rans                                

      common /comon/ apr, cosks, cosphi, costh, d, delsig, dkwt, emax,          
     1 emin(lpt), erec, ex, hevsum, hsig, oldwt, sinks, sinphi, sinth,          
     2 umax, uu, zpr, ibert, ityp, lelem, mat,maxbch, maxcas, mxmat, n,         
     3 nabov, namax, nbelo, nbogus, negex, neutno(4), lneutp, ngroup,           
     4 npidk, no, nobch, nocas, nomax, nopart, npart(6), nquit, nbertp          
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand                 
      parameter (two=dp2, half=dph)                                             
      parameter (p=2.**24, q=2.**(-24), r=2.**(-48), gb=1136868, gs=6328        
     1 637)                                                                     
c                                                                               
c     a=gs*rans                                                                 
c        this expression for b is valid only if gb+gs<2**24                     
c     b=gb*rans+gs*ranb+aint(a*q)                                               
c     rans=a-aint(a*q)*p                                                        
c     ranb=b-aint(b*q)*p                                                        
c     sflraf=two*((ranb*p+rans)*r-half)                                         
c     if (nocas.ge.151) then
c     write(11,*) sflraf 
c     endif
c     read(11,*) sflraf 
c     nrand=nrand+1                                                             
      call nombre_generer(rndm)
      sflraf=2.*rndm-1.
      return                                                                    
      end                                                                       
      function vl (rk2,l)                                                       
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
c Function VL computes the angular momentum barrier factor for angular          
c momentum L.                                                                   
c                                                                               
c    REFERENCES.......  LEVBRK                                                  
c    EXTERNALS........  EXP                                                     
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\c         
      if (l.eq.0) then                                                          
        vl=dp1                                                                  
        go to 60                                                                
      elseif (l.gt.5) then                                                      
        fl=l                                                                    
        arg=fl*(fl+dp1)/rk2                                                     
        vl=exp(-arg)                                                            
        go to 60                                                                
      endif                                                                     
      go to (10,20,30,40,50), l                                                 
   10 continue                                                                  
      vl=rk2/(rk2+1)                                                            
      go to 60                                                                  
   20 continue                                                                  
      rk4=rk2**2                                                                
      vl=rk4/(rk4+dp3*rk2+9.d0)                                                 
      go to 60                                                                  
   30 continue                                                                  
      rk4=rk2*rk2                                                               
      rk6=rk4*rk2                                                               
      vl=rk6/(rk6+6.d0*rk4+45.d0*rk2+225.d0)                                    
      go to 60                                                                  
   40 continue                                                                  
      rk4=rk2*rk2                                                               
      y=(rk4-45.d0*rk2+105.d0)**2+rk2*(dp10*rk2-105.d0)**2                      
      vl=rk4**2/y                                                               
      go to 60                                                                  
   50 continue                                                                  
      rk4=rk2*rk2                                                               
      y=rk2*(rk4-105.d0*rk2+945.d0)**2+(15.d0*rk4-42.d1*rk2+945.d0)**2          
      vl=rk4*rk4*rk2/y                                                          
   60 continue                                                                  
      return                                                                    
      end                                                                       
      function zprobf (af,zf,a2)                                                
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /fong/ ievap, ifiss                                                
      common /orcm1/ omegao, smala1, smala2, rocc, ern, c12s, ekinr4,           
     1 z1p, z2p, baz, c12, delef, e1maxp, ro, acalm, qfactr, a2m,               
     2 delken, scalke, arg, dela, delke, afnp, zfnp, fm, ke,                    
     3 kemax, izdel, iamin, ismal, iversn                                       
      common /orcm2/ ainfo(20), epkin(80,6), sigke(80,6), delexc(80,6,3)        
     1 , cstiff(120)                                                            
      data cnver /931.143d0/                                                    
cccc                                                                            
c     table lookup of fong's dma(in mev) and dza for mass corrections           
cccc                                                                            
      ro=rocc                                                                   
      a1=af-a2                                                                  
      ia1=a1+dph                                                                
      ia2=a2+dph                                                                
      iaa=ia1                                                                   
      k=iaa-59                                                                  
      if (k.lt.1) return                                                        
      kdel=ia2-iaa                                                              
      k2=k+kdel                                                                 
      if (k2.lt.1) return                                                       
      baz=dp0                                                                   
      c12=dp0                                                                   
      c12s=dp1/(a1**0.333333333d0+a2**0.333333333d0)                            
      c12s=c12s*1.44d0/ro                                                       
      c12=c12s                                                                  
cccc                                                                            
c     change ba and za                                                          
c     from fong,wing,1964,p.151,fong's book                                     
cccc                                                                            
      za1=(a1+3.d-03*a1**2)/(dp2+1.d-02*a1)                                     
      za2=(a2+3.d-03*a2**2)/(dp2+1.d-02*a2)                                     
cccc                                                                            
c     ba below is already in mev                                                
cccc                                                                            
      ba1=3.258d0-(60.22d0/a1**dph)+(431.6d0/a1)                                
      ba2=3.258d0-(60.22d0/a2**dph)+(431.6d0/a2)                                
      pac=11.51d0*(dp1/a1**dph+dp1/a2**dph)                                     
      ama1=0.0089794d0*a1**2-2.20717d0*a1+33.448d0+a1*(931.481d0-cnver)         
      ama2=0.0089794d0*a2**2-2.20717d0*a2+33.448d0+a2*(931.481d0-cnver)         
      cz=-ba1*(zf-za1)**2/dp2-ba2*za2**2/dp2                                    
      cz=cz-pac-ama1-ama2                                                       
      az=-(ba1+ba2)/dp2                                                         
      bz=(ba1*(zf-za1)+ba2*za2)                                                 
      itrmax=2                                                                  
      do 10 itr=1,itrmax                                                        
      z2p=(bz-c12*zf)/(-dp2*(az+c12))                                           
      z1p=zf-z2p                                                                
      baz=-az-c12                                                               
      c12=ekinr4/(z1p*z2p)                                                      
   10 continue                                                                  
      zprobf=z2p                                                                
      return                                                                    
      end                                                                       
                                                                    
ccc
      subroutine insig 
      implicit doubleprecision(a-h,o-z)                                         
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.        
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,          
     2 dp2th=dp2/dp3)                                                           
      common /dr/ cam4(130), cam5(200),rmass(300)
      common /dre/ exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /ev/ cam2(130), cam3(200)
      common /evc/ t(4,7)                                              
      open(unit=1,file='finbcd') 
      read (1,*) (ia(j),j=1,6),(iz(j),j=1,6)                                 
      read (1,*) (rho(j),j=1,6),(omega(j),j=1,6)                             
      read (1,*) (exmm(j),j=1,6)                                             
      read (1,*) (cam2(j),j=1,130)                                           
      read (1,*) (cam3(j),j=1,200)                                           
      read (1,*) (cam4(j),j=1,130)                                           
      read (1,*) (cam5(j),j=1,200)                                           
      read (1,*) ((t(i,j),j=1,7),i=1,3)                                      
      read (1,*) (rmass(j),j=1,297)                                          
      do 50 k=1,7                                                               
   50 t(4,k)=dp0                                                                
      rewind(1)
      close(1)
      return                                                                    
      end                                                                       




      subroutine setidt                                                     
      implicit doubleprecision(a-h,o-z)                                         
c        put the machine designator, date, and time into hi                     
c        in the form   x mm/dd/yy hh:mm:ss                                      
      dimension jd(3), jt(3)                                                    
c                                                                               
c      call idate (jd)                                                           
c      call itime (jt)                                                           
      write (6,10) jd(2),jd(1),mod(jd(3),100),jt                               
      write (6,10) jd(2),jd(1),mod(jd(3),100),jt                               
   10 format(2x,i2,1h/,i2,1h/,i2,1x,i2,1h:,i2,1h:,i2)                           
      return                                                                    
      end                                                                       
C Autre effacement collis et decay CV le 17/9/98                                        

      
      subroutine angid(rlke,cst)                                                
      real*8 r

      if (rlke.lt.5.d+02) go to 60                                              
      tesiso=0.75d0                                                             
      call nombre_generer(r)                                                    
      if (r.le.tesiso) go to 60                                                 
c**** backward/forward                                                          
   50 continue                                                                  
      call nombre_generer(r)                                                    
c****  test for direction                                                       
      cst=0.9999995d0                                                           
      snt=0.003162d0                                                            
      if (r.le.0.5) go to 70                                                    
      cst=-0.9999995d0                                                          
      go to 70                                                                  
c****  isotropic                                                                
   60 call polrn (cst,snt)                                                      
   70 return                                                                    
      end                                                                       

      subroutine polrn (cs,si)                                                  
      real*8 r1,r2
      parameter (rparm1=1.53960072)                                             
   10 continue                                                                  
      call nombre_generer(r1)                                                    
      call nombre_generer(r2)                                                    
      r1s=r1*r1                                                                 
      r2s=r2*r2                                                                 
      e=r1s+r2s                                                                 
      p=2.*r1*r2/e                                                             
      if (e-rparm1*p) 20,20,10                                                  
   20 si=p                                                                      
      cs=(r2s-r1s)/e                                                            
      return                                                                    
      end                                                                       

      subroutine angli1 (coscm)  
      real*8 rant1,eta2
   10 continue                  
      call nombre_generer(rant1)                                                    
      call nombre_generer(eta2)                                                    
      coscm=2.*rant1-1.      
      p=.25+.75*coscm*coscm  
      if (eta2.gt.p) goto 10 
      return                 
      end                   

      subroutine onepin (ecm,rmisob)  
      real*8 eta1,eta2
      rm=938.2796                    
      rmimin=1080.
      rmimax=1573.
      rmmax=ecm-rm                    
      if (rmmax.gt.rmimax) rmmax=rmimax
      rangem=rmmax-rmimin             
      pmax=phase(ecm,rm,rmimin)*200.   
c                                                                       onepin  
c        fix to avoid zero-divides of marginal isobars                  onepin  
c        ie in call to "phase/pmax", when pmax =0  [mrc 12 june 83]     onepin  
c                                                                       onepin  
      rmisob=rmimin            
      if (pmax.le.0.) return    
c                                                                       onepin  
c the max of p(mi) is less than or equal (max.sig33)*(max.phase)        onepin  
c for a given ecm.sig33 is the pi+ + proton cross section.              onepin  
c now select the isobar mass by the von neuman rejection method.        onepin  
   70 continue                                                          onepin  
      call nombre_generer(eta1)
      call nombre_generer(eta2)
      rmisob=rmimin+rangem*eta1               
      tcm=rmisob-1077.811                                               onepin  
      izero=0                                                           onepin  
      ione=1                                                            onepin  
c     call ptotal (tcm,ione,izero,izero,sig33,absor)                    onepin  
      call sspn(tcm,sipn)
c     write(2,*) tcm,sig33,sipn
c     pm=sig33*phase(ecm,rm,rmisob)/pmax 
      pm=sipn*phase(ecm,rm,rmisob)/pmax 
      if (eta2.gt.pm) goto 70                                           onepin  
      return                                                            onepin  
      end                                                               onepin  

      function phase (ecm,m3,mi)                                        phase   
c     ==================================================================phase   
c     phase calculates the two body phase factor                        phase   
c                                                                       phase   
c     phase    calls the following subroutines and functions            phase   
c     ==================================================================phase   
c                                                                       phase   
      real*4 m3,mi
      e3=(ecm*ecm+m3*m3-mi*mi)/(2.*ecm)                                 phase   
      p3sq=e3*e3-m3*m3                                                  phase   
      e4=ecm-e3                                                         phase   
      if ((e3.le.m3).or.(e4.le.mi)) goto 10                             phase   
      phase=sqrt(p3sq)*e3*e4/ecm                                        phase   
      return                                                            phase   
   10 continue                                                          phase   
      phase=0.                                                          phase   
      return                                                            phase   
      end                                                               phase   

      subroutine ptotal (e,ind9,ind10,ind11,scatt,absor)                ptotal  
c     ==================================================================ptotal  
c     ptotal, given the kinetic energy in the c.m.routine will return   ptotal  
c     the total i-i or i-j pi-nucleon cross section for that energy     ptotal  
c     and the absorption cross section as well as the pi0-nucleon       ptotal  
c     cross section.  ind9=1 means an i-j collision,ind10=1 means an    ptotal  
c     i-i collision, and ind11=1 means an i-0 collisision               ptotal  
c                                                                       ptotal  
c     ptotal   calls the following subroutines and functions            ptotal  
c     ==================================================================ptotal  
c                                                                       ptotal  

      real scaii(67), scaij(67), scao(67)                               ptotal  
      data scaii /0.0,2.,4.,6.,12.,19.,29.,42.,58.,77.,99.,129.,159.,179ptotal  
     1 .,191.,195.,185.,162.,143.,124.,109.,95.,83.,73.,64.,56.,50.,44.,ptotal  
     2 40.,36.,33.,31.,29.,27.,25.,24.,22.,21.,20.,19.,18.,17.,17.,16.,1ptotal  
     3 6.,16.,16.,16.,17.,17.,18.,18.,19.,20.,21.,22.,23.,24.,24.,25.,25ptotal  
     4 .,26.,26.,26.,27.,27.,27./                                       ptotal  
      data scaij /0.0,2.,4.,6.,8.,11.,14.,18.,22.,29.,37.,47.,55.,62.,66ptotal  
     1 .,67.,65.,59.,52.,44.,39.,36.,33.,31.,29.,28.,27.,27.,27.,27.,27.ptotal  
     2 ,27.,27.,28.,29.,30.,31.,33.,34.,36.,38.,42.,45.,47.,48.,47.,45.,ptotal  
     3 42.,40.,38.,38.,38.,39.,40.,42.,45.,49.,52.,56.,58.,59.,58.,57.,5ptotal  
     4 4.,50.,47.,45./                                                  ptotal  
      data scao /0.0,2.,4.,6.,10.,15.,22.,30.,40.,53.,68.,88.,107.,121.,ptotal  
     1 129.,131.,125.,111.,98.,84.,74.,66.,58.,52.,47.,42.,39.,36.,34.,3ptotal  
     2 2.,30.,29.,28.,28.,27.,27.,27.,27.,27.,28.,28.,30.,31.,32.,32.,32ptotal  
     3 .,31.,29.,29.,28.,28.,28.,29.,30.,32.,34.,36.,38.,40.,42.,42.,42.ptotal  
     4 ,42.,40.,39.,37.,36./                                            ptotal  
c above tables are based on cern/hera 69-1 report                       ptotal  
c above tables are good for c.m.k.e. =660mev or lab.pi k.e.=990mev      ptotal  
c    or lab pi momentum = 1.12 gev/c                                    ptotal  
      if (e.lt.660.) goto 10                                            ptotal  
      write (6,40) e                                                    ptotal  
      e=659.99                                                          ptotal  
   10 absor=0.0                                                         ptotal  
c the following is based on fact that above data are at equally         ptotal  
c spaced energy intervals of 10 mev                                     ptotal  
      rie=.1*e                                                          ptotal  
      ie=rie                                                            ptotal  
      ie=ie+1                                                           ptotal  
      de=rie-float(ie-1)                                                ptotal  
      ae=1.-de                                                          ptotal  
      if (ind9.eq.1) goto 20                                            ptotal  
      if (ind10.eq.1) goto 30                                           ptotal  
c zero ind(9),ind(10),and ind(11) before returning                      ptotal  
      ind9=0                                                            ptotal  
      ind10=0                                                           ptotal  
      ind11=0                                                           ptotal  
      scatt=ae*scao(ie)+de*scao(ie+1)                                   ptotal  
      return                                                            ptotal  
   20 continue                                                          ptotal  
c zero ind(9),ind(10),and ind(11) before returning                      ptotal  
      ind9=0                                                            ptotal  
      ind10=0                                                           ptotal  
      ind11=0                                                           ptotal  
      scatt=ae*scaii(ie)+de*scaii(ie+1)                                 ptotal  
      return                                                            ptotal  
   30 continue                                                          ptotal  
c zero ind(9),ind(10),and ind(11) before returning                      ptotal  
      ind9=0                                                            ptotal  
      ind10=0                                                           ptotal  
      ind11=0                                                           ptotal  
      scatt=ae*scaij(ie)+de*scaij(ie+1)                                 ptotal  
      return                                                            ptotal  
c                                                                       ptotal  
   40 format ('0energy=',f10.4,' too large (>659.9) for ptotal')        ptotal  
      end                                                               ptotal  
      subroutine sspn(tcm,sipn) 
cc
cc    attention donne sigma a partir de ecin lab!!
cc
      rmpi=138.
      rm=938.27
      x=tcm+rmpi+rm
      Y=X*X                                                   
      Q2=(Y-1076.**2)*(Y-800.**2)/Y/4.                       
      Q3=(SQRT(Q2))**3                                   
      F3=Q3/(Q3+180.0**3)                               
      sipn=326.5*f3/(4.*((X-1215.)/110.)**2+1.)     
   10 continue
      return
      END                                             


