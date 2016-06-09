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
// $Id: G4Mag_SpinEqRhs.hh,v 1.12 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4MagneticField;

class G4Mag_SpinEqRhs : public G4Mag_EqRhs
{
   public:  // with description

     G4Mag_SpinEqRhs( G4MagneticField* MagField );
    ~G4Mag_SpinEqRhs();
       // Constructor and destructor. No actions.

     void SetChargeMomentumMass(G4double particleCharge, // in e+ units
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

     G4double omegac;
     G4double anomaly;
     G4double pcharge;
     G4double E;
     G4double gamma;
     G4double beta;
};

#endif /* G4MAG_SPIN_EQRHS */
