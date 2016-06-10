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
// $Id: G4Mag_SpinEqRhs.hh 69699 2013-05-13 08:50:30Z gcosmo $
//
//
// class G4Mag_SpinEqRhs
//
// Class description:
//
// This is the equation of motion for a particle with spin in a pure
// magnetic field. The three components of the particle's spin are
// treated utilising BMT equation.

// History:
// - Created: J.Apostolakis, P.Gumplinger - February 8th, 1999.
// --------------------------------------------------------------------

#ifndef G4MAG_SPIN_EQRHS
#define G4MAG_SPIN_EQRHS

#include "G4Types.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ChargeState.hh"

class G4MagneticField;

class G4Mag_SpinEqRhs : public G4Mag_EqRhs
{
   public:  // with description

     G4Mag_SpinEqRhs( G4MagneticField* MagField );
    ~G4Mag_SpinEqRhs();
       // Constructor and destructor. No actions.

     void SetChargeMomentumMass(G4ChargeState particleCharge,
                                G4double MomentumXc,
                                G4double mass); 

     void EvaluateRhsGivenB( const  G4double y[],
                             const  G4double B[3],
                                    G4double dydx[] ) const;
       // Given the value of the magnetic field B, this function 
       // calculates the value of the derivative dydx.

     inline void SetAnomaly(G4double a) { anomaly = a; }
     inline G4double GetAnomaly() const { return anomaly; }
       // set/get magnetic anomaly

   private:

     G4double charge, mass, magMoment, spin;

     G4double omegac, anomaly;

     G4double beta, gamma;

};

#endif /* G4MAG_SPIN_EQRHS */
