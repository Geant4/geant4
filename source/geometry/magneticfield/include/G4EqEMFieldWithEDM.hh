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
// G4EqEMFieldWithEDM
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field, with spin tracking for both MDM and
// EDM terms.

// Created: Kevin Lynch, 19.02.2009 - Based on G4EqEMFieldWithSpin
// -------------------------------------------------------------------
#ifndef G4EQEMFIELDWITHEDM_HH
#define G4EQEMFIELDWITHEDM_HH

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"

class G4ElectroMagneticField;

class G4EqEMFieldWithEDM : public G4EquationOfMotion
{
  public:

    G4EqEMFieldWithEDM(G4ElectroMagneticField* emField );

   ~G4EqEMFieldWithEDM() override;
  
    void  SetChargeMomentumMass(G4ChargeState particleCharge, // in e+ units
                                G4double MomentumXc,
                                G4double mass) override;

    void EvaluateRhsGivenB(const G4double y[],
                           const G4double Field[],
                                 G4double dydx[] ) const override;
      // Given the value of the electromagnetic field, this function 
      // calculates the value of the derivative dydx.

    inline void SetAnomaly(G4double a) { anomaly = a; }
    inline G4double GetAnomaly() const { return anomaly; }
      // set/get magnetic anomaly

    inline void SetEta(G4double n) { eta = n; }
    inline G4double GetEta() const { return eta; }
      // set/get EDM eta parameter

  private:

    G4double charge{0.}, mass{0.}, magMoment{0.}, spin{0.};

    G4double fElectroMagCof{0.} ;
    G4double fMassCof{0.};

    G4double omegac{0.}, anomaly{0.0011659208}, eta{0.};
    G4double beta{0.}, gamma{0.};
};

#endif
