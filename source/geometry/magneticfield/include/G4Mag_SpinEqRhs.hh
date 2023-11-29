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
// G4Mag_SpinEqRhs
//
// Class description:
//
// This is the equation of motion for a particle with spin in a pure
// magnetic field. The three components of the particle's spin are
// treated utilising BMT equation.

// Created: J.Apostolakis, P.Gumplinger - 08.02.1999
// --------------------------------------------------------------------
#ifndef G4MAG_SPIN_EQRHS_HH
#define G4MAG_SPIN_EQRHS_HH

#include "G4Types.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ChargeState.hh"

class G4MagneticField;

class G4Mag_SpinEqRhs : public G4Mag_EqRhs
{
   public:

     G4Mag_SpinEqRhs( G4MagneticField* MagField );
    ~G4Mag_SpinEqRhs() override;
       // Constructor and destructor. No actions.

     void SetChargeMomentumMass(G4ChargeState particleCharge,
                                G4double MomentumXc,
                                G4double mass) override; 

     void EvaluateRhsGivenB( const  G4double y[],
                             const  G4double B[3],
                                    G4double dydx[] ) const override;
       // Given the value of the magnetic field B, this function 
       // calculates the value of the derivative dydx.

     inline void SetAnomaly(G4double a) { anomaly = a; }
     inline G4double GetAnomaly() const { return anomaly; }
       // set/get magnetic anomaly

   private:

     G4double charge=0.0, mass=0.0, magMoment=0.0, spin=0.0;
     G4double omegac=0.0, anomaly=0.0011659208;
     G4double beta=0.0, gamma=0.0;
};

#endif
