// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Mag_SpinEqRhs.hh,v 1.4 2000-11-01 15:15:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4Mag_SpinEqRhs
//
// Class description:
//
// This is the standard right-hand side for equation of motion.
// This version of the right-hand side includes the three components
// of the particle's spin.

// History:
// - Created: J.Apostolakis, P.Gumplinger - February 8th, 1999.

#ifndef G4MAG_SPIN_EQRHS
#define G4MAG_SPIN_EQRHS

#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"

class G4Mag_SpinEqRhs : public G4Mag_EqRhs
{
   public:  // with description

     G4Mag_SpinEqRhs( G4MagneticField* MagField )
       : G4Mag_EqRhs( MagField ) {;}
    ~G4Mag_SpinEqRhs() {;}
       // Constructor and destructor. No actions.

     void SetChargeMomentumMass(G4double particleCharge, // in e+ units
                                G4double MomentumXc,
                                G4double mass); 

     void EvaluateRhsGivenB( const  G4double y[],
			     const  G4double B[3],
			     G4double dydx[] ) const;
       // Given the value of the magnetic field B, this function 
       // calculates the value of the derivative dydx.

   private:

     G4double omegac;
     G4double anomaly;
     G4double ParticleCharge;

};

#endif /* G4MAG_SPIN_EQRHS */
