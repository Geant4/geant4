      subroutine init_dresner(RACINE)
      implicit real*8 (a-h,o-z)
      CHARACTER*80 FILENAME,RACINE,FILEDAT  
      parameter (max=500)
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
      common /fong/ iEVAp, ifiss                                                
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
                                 
      common /evcm/ y0, b0, t(4,7)                                              
      common /dresm/ exmass(6), exmm(6), rho(6), omega(6), ia(6), iz(6)         
      common /dresc/ pp0(1001), pp1(1001), pp2(1001), cam4(130), cam5(20        
     1 0), rmass(300), alph(300), bet(300)                                      
      logical penbar
      common /demitr/ sos(6), strun(6), zmass(6), q(6), fla(6),                  
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
      data index  / 1 /
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
      data exmass /8.3651d0,7.5831d0,13.7222d0,15.8383d0,15.8193d0,             
     1 3.6084d0/                                                                
      real*4 zcn_proden,acn_proden,excn_proden,tabener
      common /inp_gemden/ zcn_proden,acn_proden,excn_proden,jiji
	common/cv/tabener(0:2000)
      common/part/apart(500),zpart(500),exci(500),pxp(500)
     +,pyp(500),pzp(500),xlive(500),xp(500),yp(500),zp(500)
     +,temiss(500),aemis(500),zemis(500),tzero,erel(500)
     +,exint(2,500),ether(500),spinx(500),spiny(500),spinz(500)
     +,erot(500),proven(500),tcas(500),mult
      common/serie/aev(5000),zev(5000),iccmax
      common/pres/iexist(0:300,0:300)
      common/estart/edebut
      common/dec/alevel,aclust,tlimite,npartic
      common/info1/xemax(20,20),deg(20,20,20),exm(20,20,20)
     #,ase(20,20,20),zse(20,20,20)
      common/info2/amaxn,zmaxn
      common/info3/nniv(20,20),nmo(20,20,20)

      real*4 zzz(300),ccc(300),aijk

      data loff /0/                                                             
      data tip0, x0, y0, z0, a0, b0, cuts, smu /8*dp0/, isopt /0/, e0 /8        
     1 .d2/                                                                     
      data itopt /-1/, twit /dp0/, tpeak /dp0/, tsig /dp0/                      
      data zmass /9.395730d+02,9.387910d+02,1.876138d+03,2.809462d+03,2.        
     1 809443d+03,3.728440d+03/, fla /dp1,dp1,dp2,dp3,dp3,dp4/, flz /dp0        
     2 ,3*dp1,2*dp2/, emx /2.d2/                                                
      data z91 /91.d0/, ucut /dp4/, ievap /0/, qfactr /dp10/, ismal /2/,        
     1 iversn /1/, yzero /1.5d0/, bzero /dp10/, yzere /1.5d0/, bzere /8         
     2 .d0/, ro /1.2d0/, kemax /6/                                              
      data monin, se, ses, sr, srs /7*0,20*dp0/

C ....... semences pour verification ......                                     
      IJ = 1802                                                                 
      KL = 9373                                                                 
      call generateur_init(ij,kl)                                               
      do 40 k=1,250                                                             
      do 40 j=1,20                                                              
   40 waps(k,j)=0.                                                              
      nbertp=8
c      filename='bertini2'
c      REWIND(nbertp)
      
      FILEDAT='/bertini2'
      long=0
      i=0
      DO WHILE(long.EQ.0)
      	i=i+1
	IF(RACINE(i:i).EQ.' ') long=i-1
      END DO
      FILENAME=RACINE(1:long)//FILEDAT           
      open(unit=nbertp,file=filename,status='old')
      read (nbertp,*) (pp0(j),pp1(j),pp2(j),j=1,1001)                           
      read (nbertp,*) (ia(j),j=1,6),(iz(j),j=1,6)                               
      read (nbertp,*) (rho(j),j=1,6),(omega(j),j=1,6)                           
      read (nbertp,*) (exmm(j),j=1,6)                                           
      read (nbertp,*) (cam2(j),j=1,130)                                         
      read (nbertp,*) (cam3(j),j=1,200)                                         
      read (nbertp,*) (cam4(j),j=1,130)                                         
      read (nbertp,*) (cam5(j),j=1,200)                                         
      read (nbertp,*) ((t(i,j),j=1,7),i=1,3)                                    
      read (nbertp,*) (rmass(j),j=1,297)                                        
      read (nbertp,*) (alph(j),j=1,297)                                         
      read (nbertp,*) (bet(j),j=1,297)                                          
      read (nbertp,*) ((waps(i,j),i=1,250),j=1,20)                              
      do 50,k=1,7
   50 t(4,k)=0.
      waps(1,11)=exmass(1)                                                      
      waps(1,9)=exmass(2)                                                       
      waps(2,10)=exmass(3)                                                      
      waps(3,11)=exmass(4)                                                      
      waps(3,9)=exmass(5)                                                       
      waps(4,10)=exmass(6)                                                      
      waps(4,8)=exmass(5)+exmass(2)+2.91d0                                      
      waps(4,12)=exmass(4)+exmass(1)+2.9d0                                      
      waps(5,7)=exmass(5)+dp2*exmass(2)+2.7d0                                   
      waps(5,13)=exmass(4)+dp2*exmass(1)+2.69d0                                 
      io=6
      nbertp=9
      ndisk=10
ccccc
      ifbrk=1
      ilvden=0 
c     ilvden=-1
c     yzere=0.
c     yzero=0.
c     bzere=8.
c     bzero=8.
      ievap=0
      IF(ievap.EQ.0) WRITE(6,*) 'RAL fission in Dresner'
      nofis=1
      costh=1.
      sinth=0.
      
C Added by AB (8/2008); Was NOT initialised before....
      cosphi=1.
      sinphi=0.
      
      coslbr(1)=0.
      coslbr(2)=0.
      coslbr(3)=1.
ccccc
c      filename='bertini3'
      FILEDAT = '/bertini3'                                                   
      FILENAME=RACINE(1:long)//FILEDAT
c      REWIND(11)
c      REWIND(ndisk)
                                                           
      OPEN(UNIT =11,STATUS='UNKNOWN',FILE=FILENAME)
      close(11)
      open(unit=ndisk,file=filename,status='old')
      if (ievap.ne.1) then                                                      
        read (ndisk,*) ainfom                                                      
        read (ndisk,*) epkinm,sigkem                                                
        read (ndisk,*) delexcm                                                     
        go to 70                                                                
      endif                                                                     
      read (ndisk,*) (ainfo(j),j=1,20)                                           
      itrigg=ainfo(16)                                                          
      if (itrigg.ne.0) iversn=itrigg                                            
      afnp=ainfo(1)                                                             
      zfnp=ainfo(2)                                                             
      read (ndisk,*) ((epkin(i,j),i=1,80),j=1,6)
     +,((sigke(i,j),i=1,80),j=1,6)                             
      read (ndisk,*) (((delexc(i,j,k),i=1,80),j=1,6),k=1,3)                        
cccc                                                                            
c     the following excitation energies are the computed excitations for        
c     epperson's incident proton energies on u235,92                            
c     ef1 = ainfo(10)                                                           
c     ef2 = ainfo(11)                                                           
c     ef3 = ainfo(12)                                                           
c     ef4 = ainfo(13)                                                           
c     ef5 = ainfo(14)                                                           
c     ef6 = ainfo(15)                                                           
cccc                                                                            
      b0=bzero                                                                  
      y0=yzero                                                                  
      if (ismal.gt.2) go to 30                                                  
      bo1=dp1/ainfo(19)                                                         
      bo2=dp1/ainfo(17)                                                         
      bc=dp1/b0                                                                 
      bfactr=(bc-bo1)/(bo2-bo1)                                                 
      do 20 j=1,80                                                              
      do 20 ke=1,6                                                              
      delex1=delexc(j,ke,1)                                                     
      deldel=delexc(j,ke,2)-delex1                                              
      delexc(j,ke,ismal)=delex1+bfactr*deldel                                   
   20 continue                                                                  
   30 continue                                                                  
      if (iversn.eq.2) fm=.125d0                                                
      if (iversn.eq.1) fm=1.25d0                                                
cccc                                                                            
c     vary these initial values to determine qamax                              
cccc                                                                            
      a2m=133.d0                                                                
cccc                                                                            
c     iamin is the minimum a-value for a fission fragment                       
cccc                                                                            
      izdel=7                                                                   
      iamin=70                                                                  
cccc                                                                            
c     acalm is the total number of tries for a given fission history            
cccc                                                                            
      acalm=4.d2                                                                
      if=0                                                                      
      ib=1                                                                      
      do 80 ih=69,165,2                                                         
      if=if+1                                                                   
      iah1=ih-58                                                                
      cstiff(iah1)=cfong(if)                                                    
      cstiff(iah1+1)=cfong(if)+dph*(cfong(if+1)-cfong(if))                      
   80 continue                                                                  
      do 90 ih=59,68                                                            
      iah1=ih-58                                                                
      cstiff(iah1)=cfong(ib)                                                    
   90 continue                                                                  
      do 60 ih=167,178                                                          
      iah1=ih-58                                                                
      cstiff(iah1)=cfong(if)                                                    
   60 continue                                                                  

   70 continue
c     filename='bertini1'
      FILEDAT = '/bertini1'                                                   
      FILENAME=RACINE(1:long)//FILEDAT 
c      REWIND(nbertp)
      open(unit=nbertp,file=filename,status='old')
      call redmas(0,io,nbertp)
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cfb
