C*********This routine is the same as the drein2 in LAHET code**********
      subroutine drein2 (s,a,eye2)                                      
      implicit doubleprecision(a-h,o-z)                                 
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
c                                                                       
c     compute statistical theory emission integrals                     
c                                                                       
c     compute third integral                                            
c                                                                       
c     for s<dph use a series expansion.                                 
c     for s>dph the explicit relationship                               
c                                                                       
c     coeficients for series expansions                                 
c                                                                       
      dimension c2(7)                                                   
      data c2 /0.45714286d0,0.125d0,0.02539683d0,0.00416667d0,          
     1 0.0005772d0,0.00006944d0,0.0000074d0/                            
      exps=dp0                                                          
      if (s.lt.1.d+02) exps=exp(-s)                                     
      if (s.lt.dph) go to 10                                            
c///// explicit relation                                                
      b=s*s                                                             
      eye2=0.25d0*(s*(15.d0-s*(6.d0-s))-15.d0+(15.d0+0.125d0*b*(b-12.d0)
     1 )*exps)/(a*a*a)                                                  
      return                                                            
   10 continue                                                          
c                                                                       
c///// series expansion                                                 
c    eye2=(1/(32a**3))*s**6/6*(sum n=0 to 7:48*s**n/(n!(n+2)(n+4)(n+6)) 
c                                                                       
      eye2=dp1                                                          
      b=dp1                                                             
      do 20 n=1,7                                                       
      b=b*s                                                             
      c=b*c2(n)                                                         
      if (c.lt.1.0d-7) go to 30                                         
      eye2=eye2+c                                                       
   20 continue                                                          
   30 continue                                                          
      b=0.25d0*s*s/a                                                    
      eye2=eye2*b*b*b*exps*0.33333333d0                                 
      return                                                            
      end                                                               
