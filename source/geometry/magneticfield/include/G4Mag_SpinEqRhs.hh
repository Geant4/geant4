//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Mag_SpinEqRhs.hh,v 1.7 2001-07-17 02:13:43 gum Exp $
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
// - Modified: D. Cote-Ahern, P.Gumplinger - April 11th, 2001.

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

     G4double E;
     G4double gamma;
     G4double beta;

};

#endif /* G4MAG_SPIN_EQRHS */
