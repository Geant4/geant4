//
//
//  This is the standard right-hand side for equation of motion.
//        This version of the right-hand side includes
//        the three components of the particle's spin.
//
//            J. Apostolakis, February 8th, 1999
//            P. Gumplinger,  February 8th, 1999
//
#ifndef G4MAG_SPIN_EQRHS
#define G4MAG_SPIN_EQRHS

#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"

class G4Mag_SpinEqRhs: public G4Mag_EqRhs{

   public:

     G4Mag_SpinEqRhs( G4MagneticField* MagField ) :
		       G4Mag_EqRhs( MagField ) {};
    ~G4Mag_SpinEqRhs() {} ; 

     void SetChargeMomentumMass(const G4double particleCharge, // in e+ units
                                const G4double MomentumXc,
                                const G4double mass); 

     //  Given the value of the magnetic field B, this function 
     //   calculates the value of the derivative dydx.
     //
     void EvaluateRhsGivenB( const  G4double y[],
			     const  G4double B[3],
			     G4double dydx[] ) const;

   private:

     G4double omegac;
     G4double anomaly;
     G4double ParticleCharge;

};

#endif /* G4MAG_SPIN_EQRHS */
