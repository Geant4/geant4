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
//---------------------------------------------------------------------------
//
// ClassName:   G4ChargeExchangePhysics
//
// Author: 19 November 2008 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4ChargeExchangePhysics_h
#define G4ChargeExchangePhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class G4ChargeExchangePhysics : public G4VPhysicsConstructor
{
public: 
  explicit G4ChargeExchangePhysics(G4int ver = 1);
  ~G4ChargeExchangePhysics() override = default;

  void ConstructParticle() override;
 
  void ConstructProcess() override;

  void SetLowEnergyLimit(G4double val) { fLowEnergyLimit = val; }

  void SetCrossSectionFactor(G4double val) { fXSFactor = val; }

  G4ChargeExchangePhysics& operator=
  (const G4ChargeExchangePhysics& right) = delete;
  G4ChargeExchangePhysics(const G4ChargeExchangePhysics&) = delete;

private:

  G4double fLowEnergyLimit;
  G4double fXSFactor{1.0};
};


#endif








