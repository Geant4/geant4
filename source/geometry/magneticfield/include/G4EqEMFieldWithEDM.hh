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
// $Id: G4EqEMFieldWithEDM.hh,v 1.2 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4EqEMFieldWithEDM
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field, with spin tracking for both MDM and
// EDM terms.

// History:
// - Created. Kevin Lynch, 19.02.2009, based on G4EqEMFieldWithSpin
// -------------------------------------------------------------------

#ifndef G4EQEMFIELDWITHEDM_hh
#define G4EQEMFIELDWITHEDM_hh

#include "G4EquationOfMotion.hh"

class G4ElectroMagneticField;

class G4EqEMFieldWithEDM : public G4EquationOfMotion
{
  public:  // with description

    G4EqEMFieldWithEDM(G4ElectroMagneticField *emField );

    ~G4EqEMFieldWithEDM();
  
    void  SetChargeMomentumMass(G4double particleCharge, // in e+ units
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

    inline void SetEta(G4double n) { eta = n; }
    inline G4double GetEta() const { return eta; }
      // set/get EDM eta parameter

  private:

    G4double fElectroMagCof ;
    G4double fMassCof;
    G4double omegac;
    G4double anomaly;
    G4double eta;
    G4double pcharge;
    G4double E;
    G4double gamma;
    G4double beta;

};

#endif /* G4EQEMFIELDWITHEDM */
