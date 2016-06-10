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
// $Id: G4EqEMFieldWithSpin.hh 69699 2013-05-13 08:50:30Z gcosmo $
//
//
// class G4EqEMFieldWithSpin
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field.

// History:
// - Created. Chris Gong, P.Gumplinger, 30.08.2007
// -------------------------------------------------------------------

#ifndef G4EQEMFIELDWITHSPIN_hh
#define G4EQEMFIELDWITHSPIN_hh

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"

class G4ElectroMagneticField;

class G4EqEMFieldWithSpin : public G4EquationOfMotion
{
  public:  // with description

    G4EqEMFieldWithSpin(G4ElectroMagneticField *emField );

    ~G4EqEMFieldWithSpin();

    void  SetChargeMomentumMass(G4ChargeState particleCharge, // in e+ units
                                G4double MomentumXc,
                                G4double mass);

    void EvaluateRhsGivenB(const G4double y[],
                           const G4double Field[],
                                 G4double dydx[] ) const;
      // Given the value of the electromagnetic field, this function 
      // calculates the value of the derivative dydx.

    inline void SetAnomaly(G4double a) { anomaly = a; }
    inline G4double GetAnomaly() const { return anomaly; }
      // set/get magnetic anomaly

  private:

    G4double charge, mass, magMoment, spin;

    G4double fElectroMagCof ;
    G4double fMassCof;

    G4double omegac, anomaly;
    G4double beta, gamma;

};

#endif /* G4EQEMFIELDWITHSPIN */
