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
G4EqMagElectricField::SetChargeMomentumMass(G4double particleCharge, // e+ units
		                            G4double MomentumXc,
                                            G4double particleMass)
{
   fElectroMagCof =  eplus*particleCharge ;
   fMassCof = particleMass*particleMass/c_light/c_light; 
}



void
G4EqMagElectricField::EvaluateRhsGivenB(const G4double y[],
			                const G4double Field[],
				              G4double dydx[] ) const
{

   // Components of y:
   //    0-2 dr/ds, 
   //    3-5 dp/ds - momentum derivatives 

   G4double pSquared = y[3]*y[3] + y[4]*y[4] + y[5]*y[5] ;

   G4double cof2     = sqrt( pSquared + fMassCof )/c_light ;

   G4double pModuleInverse  = 1.0/sqrt(pSquared) ;

   G4double cof1     = fElectroMagCof*pModuleInverse ;

   //  G4double vDotE = y[3]*Field[3] + y[4]*Field[4] + y[5]*Field[5] ;


   dydx[0] = y[3]*pModuleInverse ;                         
   dydx[1] = y[4]*pModuleInverse ;                         
   dydx[2] = y[5]*pModuleInverse ;                        

   dydx[3] = cof1*(cof2*Field[3] + (y[4]*Field[2] - y[5]*Field[1])) ;
   
   dydx[4] = cof1*(cof2*Field[4] + (y[5]*Field[0] - y[3]*Field[2])) ; 
 
   dydx[5] = cof1*(cof2*Field[5] + (y[3]*Field[1] - y[4]*Field[0])) ;  
   return ;
}
