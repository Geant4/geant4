//
//  This is the standard right-hand side for equation of motion.
//
//    The only case another is required is when using a moving reference
//     frame ... or extending the class to include additional Forces,
//     eg an electric field
//
//           10.11.98   V.Grichine
//
#include "G4EqMagElectricField.hh"

void  
G4EqMagElectricField::
SetChargeMomentumMass( const G4double particleCharge, // e+ units
			            const G4double MomentumXc,
                                    const G4double particleMass)
{
   fElectroMagCof =  eplus*c_squared ;
   fElectroMagCof *= particleCharge/particleMass; 
  
}



void
G4EqMagElectricField::EvaluateRhsGivenB( const G4double y[],
			                 const G4double Field[],
				               G4double dydx[] ) const
{

   // Components of y:
   //    0-2 dr/ds, 
   //    3-5 dv/ds  
   //  FCof() = charge/mass

   G4double vSquared = y[3]*y[3] + y[4]*y[4] + y[5]*y[5] ;

   G4double vModule = sqrt(vSquared) ;

   G4double vDotE = y[3]*Field[3] + y[4]*Field[4] + y[5]*Field[5] ;

   G4double gammaReci = sqrt(1 - vSquared/c_squared) ;

   dydx[0] = y[3]/vModule ;                         
   dydx[1] = y[4]/vModule ;                         
   dydx[2] = y[5]/vModule ;                        

   dydx[3] = fElectroMagCof*(Field[3] + (y[4]*Field[2] - y[5]*Field[1])/c_light -
                        y[3]*vDotE/c_light/c_light )*gammaReci/vModule ;
   
   dydx[4] = fElectroMagCof*(Field[4] + (y[5]*Field[0] - y[3]*Field[2])/c_light -
                        y[4]*vDotE/c_light/c_light )*gammaReci/vModule ; 
 
   dydx[5] = fElectroMagCof*(Field[5] + (y[3]*Field[1] - y[4]*Field[0])/c_light -
                        y[5]*vDotE/c_light/c_light)*gammaReci/vModule ;  
   return ;
}
