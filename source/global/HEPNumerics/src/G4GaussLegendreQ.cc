// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GaussLegendreQ.cc,v 1.1 1999-01-07 16:08:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4GaussLegendreQ.hh"


G4GaussLegendreQ::G4GaussLegendreQ( function pFunction )
   : G4VGaussianQuadrature(pFunction)
{
   ;
}



// ----------------------------------------------------------------------------
//
// Constructor for GaussLegendre quadrature method. The value nLegendre set the
// accuracy required, i.e the number of points where the function pFunction will
// be evaluated during integration. The constructor creates the arrays for 
// abscissas and weights that used in Gauss-Legendre quadrature method. 
// The values a and b are the limits of integration of the pFunction.
// nLegendre MUST BE EVEN !!!

G4GaussLegendreQ::G4GaussLegendreQ( function pFunction,
				    G4int nLegendre           )
   : G4VGaussianQuadrature(pFunction)
{
   const G4double tolerance = 1.6e-10 ;
   G4int i, j,   k = nLegendre ;
   fNumber = (nLegendre + 1)/2 ;
   if(2*fNumber != k)
   {
      G4Exception("Invalid nLegendre in G4GaussLegendreQ::G4GaussLegendreQ") ;
   }
   G4double newton, newton1, temp1, temp2, temp3, temp ;

   fAbscissa = new G4double[fNumber] ;
   fWeight   = new G4double[fNumber] ;
      
   for(i=1;i<=fNumber;i++)      // Loop over the desired roots
   {
      newton = cos(pi*(i - 0.25)/(k + 0.5)) ;  // Initial root approximation
      do
      {               // loop of Newton's method               
	 temp1 = 1.0 ;
	 temp2 = 0.0 ;
	 for(j=1;j<=k;j++)
	 {
	    temp3 = temp2 ;
	    temp2 = temp1 ;
	    temp1 = ((2.0*j - 1.0)*newton*temp2 - (j - 1.0)*temp3)/j ;
	 }
	 temp = k*(newton*temp1 - temp2)/(newton*newton - 1.0) ;
	 newton1 = newton ;
	 newton  = newton1 - temp1/temp ;       // Newton's method
      }
      while(fabs(newton - newton1) > tolerance) ;
	 
      fAbscissa[fNumber-i] =  newton ;
      fWeight[fNumber-i] = 2.0/((1.0 - newton*newton)*temp*temp) ;
   }
}


// -------------------------------------------------------------------------------
//
// Returns the integral of the function to be pointed by fFunction between a and b,
// by 2*fNumber point Gauss-Legendre integration: the function is evaluated exactly
// 2*fNumber Times at interior points in the range of integration. Since the weights
// and abscissas are, in this case, symmetric around the midpoint of the range of
// integration, there are actually only fNumber distinct values of each.

G4double 
G4GaussLegendreQ::Integral(G4double a, G4double b) const 
{
   G4int i ;
   G4double xDiff, xMean, dx, integral ;
   
   xMean = 0.5*(a + b) ;
   xDiff = 0.5*(b - a) ;
   integral = 0.0 ;
   for(i=0;i<fNumber;i++)
   {
      dx = xDiff*fAbscissa[i] ;
      integral += fWeight[i]*(fFunction(xMean + dx) + fFunction(xMean - dx)) ;
   }
   return integral *= xDiff ;
}

// -------------------------------------------------------------------------------
//
// Returns the integral of the function to be pointed by fFunction between a and b,
// by ten point Gauss-Legendre integration: the function is evaluated exactly
// ten Times at interior points in the range of integration. Since the weights
// and abscissas are, in this case, symmetric around the midpoint of the range of
// integration, there are actually only five distinct values of each

G4double 
   G4GaussLegendreQ::QuickIntegral(G4double a, G4double b) const 
{
   G4int i ;
   G4double xDiff, xMean, dx, integral ;
   
   // From Abramowitz M., Stegan I.A. 1964 , Handbook of Math... , p. 916
   
   static G4double abscissa[] = { 0.148874338981631, 0.433395394129247,
                                  0.679409568299024, 0.865063366688985,
				  0.973906528517172                      } ;
   
   static G4double weight[] =   { 0.295524224714753, 0.269266719309996, 
                                  0.219086362515982, 0.149451349150581,
				  0.066671344308688                      } ;
   xMean = 0.5*(a + b) ;
   xDiff = 0.5*(b - a) ;
   integral = 0.0 ;
   for(i=0;i<5;i++)
   {
      dx = xDiff*abscissa[i] ;
      integral += weight[i]*(fFunction(xMean + dx) + fFunction(xMean - dx)) ;
   }
   return integral *= xDiff ;
}


// -------------------------------------------------------------------------
//
// Returns the integral of the function to be pointed by fFunction between a and b,
// by 96 point Gauss-Legendre integration: the function is evaluated exactly
// ten Times at interior points in the range of integration. Since the weights
// and abscissas are, in this case, symmetric around the midpoint of the range of
// integration, there are actually only five distinct values of each

G4double 
   G4GaussLegendreQ::AccurateIntegral(G4double a, G4double b) const 
{
   G4int i ;
   G4double xDiff, xMean, dx, integral ;
   
   // From Abramowitz M., Stegan I.A. 1964 , Handbook of Math... , p. 919
   
   static 
   G4double abscissa[] = { 
                           0.016276744849602969579, 0.048812985136049731112,
                           0.081297495464425558994, 0.113695850110665920911,
                           0.145973714654896941989, 0.178096882367618602759,  // 6
                           
			   0.210031310460567203603, 0.241743156163840012328,
			   0.273198812591049141487, 0.304364944354496353024,
			   0.335208522892625422616, 0.365696861472313635031,  // 12
			   
			   0.395797649828908603285, 0.425478988407300545365,
			   0.454709422167743008636, 0.483457973920596359768,
			   0.511694177154667673586, 0.539388108324357436227,  // 18
			   
			   0.566510418561397168404, 0.593032364777572080684,
			   0.618925840125468570386, 0.644163403784967106798,
			   0.668718310043916153953, 0.692564536642171561344,  // 24
			   
			   0.715676812348967626225, 0.738030643744400132851,
			   0.759602341176647498703, 0.780369043867433217604,
			   0.800308744139140817229, 0.819400310737931675539,  // 30
			   
			   0.837623511228187121494, 0.854959033434601455463,
			   0.871388505909296502874, 0.886894517402420416057,
			   0.901460635315852341319, 0.915071423120898074206,  // 36
			   
			   0.927712456722308690965, 0.939370339752755216932,
			   0.950032717784437635756, 0.959688291448742539300,
			   0.968326828463264212174, 0.975939174585136466453,  // 42
			   
			   0.982517263563014677447, 0.988054126329623799481,
			   0.992543900323762624572, 0.995981842987209290650,
			   0.998364375863181677724, 0.999689503883230766828   // 48
                                                                            } ;
   
   static 
   G4double weight[] =   {  
                           0.032550614492363166242, 0.032516118713868835987,
                           0.032447163714064269364, 0.032343822568575928429,
			   0.032206204794030250669, 0.032034456231992663218,  // 6
			   
			   0.031828758894411006535, 0.031589330770727168558,
			   0.031316425596862355813, 0.031010332586313837423,
			   0.030671376123669149014, 0.030299915420827593794,  // 12
			   
			   0.029896344136328385984, 0.029461089958167905970,
			   0.028994614150555236543, 0.028497411065085385646,
			   0.027970007616848334440, 0.027412962726029242823,  // 18
			   
			   0.026826866725591762198, 0.026212340735672413913,
			   0.025570036005349361499, 0.024900633222483610288,
			   0.024204841792364691282, 0.023483399085926219842,  // 24
			   
			   0.022737069658329374001, 0.021966644438744349195,
			   0.021172939892191298988, 0.020356797154333324595,
			   0.019519081140145022410, 0.018660679627411467385,  // 30
			   
			   0.017782502316045260838, 0.016885479864245172450,
			   0.015970562902562291381, 0.015038721026994938006,
			   0.014090941772314860916, 0.013128229566961572637,  // 36
			   
			   0.012151604671088319635, 0.011162102099838498591,
			   0.010160770535008415758, 0.009148671230783386633,
			   0.008126876925698759217, 0.007096470791153865269,  // 42
			   
			   0.006058545504235961683, 0.005014202742927517693,
			   0.003964554338444686674, 0.002910731817934946408,
			   0.001853960788946921732, 0.000796792065552012429   // 48
                                                                            } ;
   xMean = 0.5*(a + b) ;
   xDiff = 0.5*(b - a) ;
   integral = 0.0 ;
   for(i=0;i<48;i++)
   {
      dx = xDiff*abscissa[i] ;
      integral += weight[i]*(fFunction(xMean + dx) + fFunction(xMean - dx)) ;
   }
   return integral *= xDiff ;
}

