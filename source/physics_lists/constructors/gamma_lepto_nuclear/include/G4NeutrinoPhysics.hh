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
// ClassName:   G4NeutrinoPhysics
//
// Author: 2023 V. Ivanchenko extracted from G4EmExtraPhysics
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4NeutrinoPhysics_h
#define G4NeutrinoPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4NeutrinoPhysicsMessenger.hh"

class G4NeutrinoPhysics : public G4VPhysicsConstructor
{
public:

  G4NeutrinoPhysics(G4int ver = 1);

  ~G4NeutrinoPhysics() override;

  void ConstructParticle() override;
  void ConstructProcess() override;

  void NuETotXscActivated(G4bool val);
  void SetNuOscillation(G4bool val);
  void SetNuEleCcBias(G4double bf);
  void SetNuEleNcBias(G4double bf);
  void SetNuNucleusBias(G4double bf);
  void SetNuOscDistanceBias(G4double bf);
  void SetNuDetectorName(const G4String& dn);
  void SetNuOscDistanceName(const G4String& dn);

  G4NeutrinoPhysics& operator=(const G4NeutrinoPhysics& right) = delete;
  G4NeutrinoPhysics(const G4NeutrinoPhysics&) = delete;

private:

  G4bool fNuETotXscActivated = false;
  G4bool fNuOscillation = true;

  G4double fNuEleCcBias = 1.0;
  G4double fNuEleNcBias = 1.0;
  G4double fNuNucleusBias = 1.0;
  G4double fNuOscDistanceBias = 1.0;

  G4String fNuDetectorName = "0";
  G4String fNuOscDistanceName = "0";

  G4NeutrinoPhysicsMessenger* theMessenger;
  G4int verbose;
};

#endif





