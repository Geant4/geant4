C*********This routine is the same as the drein1 in LAHET code********** 
      subroutine drein1 (j,s,a,eye1,eye0)                               
      implicit doubleprecision(a-h,o-z)                                 
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
c                                                                       
c     compute statistical theory emission integrals                     
c     for s<dph use a series expansion.                                 
c     for s>dph the explicit relationship                               
c     return for neutrons (and compute eye0)                            
c     return 1 for all others.                                          
c                                                                       
c     coeficients for series expansions                                 
c                                                                       
c     correction of c1 by r. e. prael                                   
      dimension c0(7), c1(7)                                            
      data c0 /0.66666667d0,0.25d0,0.06666667d0,0.01388889d0,           
     1 0.00238095d0,0.00034722d0,0.00004409d0/                          
      data c1 /0.53333333d0,0.16666667d0,0.03809524d0,0.00694444d0,     
     1 0.0010582d0,0.00013889d0,0.0000160d0/                            
c                                                                       
      exps=dp0                                                          
      if (s.lt.1.d+02) exps=exp(-s)                                     
      if (s.lt.dph) go to 10                                            
c///// explicit relation                                                
      b=dph/a                                                           
      eye1=b*b*(dp3+s*(s-dp3)+exps*(dph*s*s-dp3))                       
      if (j.eq.1) eye0=b*(s-dp1+exps)                                   
      return                                                            
c///// small s series expansion                                         
   10 continue                                                          
c                                                                       
c  eye1=(1.0/(8*a*a))*(s**4/4)*(sum n=0 to 7:8*s**n/(n!*(n+2)*(n+4))    
c                                                                       
      eye1=dp1                                                          
      b=dp1                                                             
      do 20 n=1,7                                                       
      b=b*s                                                             
      c=b*c1(n)                                                         
      if (c.lt.1.0d-7) go to 30                                         
      eye1=eye1+c                                                       
   20 continue                                                          
   30 continue                                                          
      b=s*s/a                                                           
      eye1=eye1*exps*b*b*0.03125d0                                      
      if (j.gt.1) return                                                
c eye0 (neutrons only)=(.5/a)*s**2/2*(sum n=0 to 7:2*s**n/(n!*(n+2))    
      eye0=dp1                                                          
      b=dp1                                                             
      do 40 n=1,7                                                       
      b=b*s                                                             
      c=b*c0(n)                                                         
      if (c.lt.1.0d-7) go to 50                                         
      eye0=eye0+c                                                       
   40 continue                                                          
   50 continue                                                          
      eye0=exps*eye0*s*s*0.25d0/a                                       
      return                                                            
      end                                                               
