//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4Mag_UsualEqRhs.hh,v 1.7 2006-06-29 18:23:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4Mag_UsualEqRhs
//
// Class description:
//
// This is the standard right-hand side for equation of motion.
// The only case another is required is when using a moving reference
// frame ... or extending the class to include additional Forces,
// eg an electric field

// History:
// - Created: J. Apostolakis, January 13th 1997.
// --------------------------------------------------------------------

#ifndef G4MAG_USUAL_EQRHS
#define G4MAG_USUAL_EQRHS

#include "G4Mag_EqRhs.hh"

class G4MagneticField;

class G4Mag_UsualEqRhs : public G4Mag_EqRhs
{
   public:  // with description

     G4Mag_UsualEqRhs( G4MagneticField* MagField );
    ~G4Mag_UsualEqRhs();
       // Constructor and destructor. No actions.

     void EvaluateRhsGivenB( const G4double y[],
                             const G4double B[3],
                                   G4double dydx[] ) const;
       // Given the value of the magnetic field B, this function 
       // calculates the value of the derivative dydx.

     virtual void SetChargeMomentumMass( G4double particleCharge, // in e+ units
                                         G4double MomentumXc,
                                         G4double mass);
     
  private:

    G4double  fInvCurrentMomentumXc;   // This extra state enables us 
                                    // to save a square root in a
                                    // critical method.
};

#endif /* G4MAG_USUAL_EQRHS */
