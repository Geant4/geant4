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

// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file lArTPCSD.hh
/// \brief Definition of the CaTS::lArTPCSD class

#pragma once

#include "lArTPCHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ScintillationTrackInformation.hh"
#include <G4MaterialPropertyVector.hh>
#include <G4String.hh>
#include <G4Types.hh>
class G4Material;
class G4PhysicsTable;
class G4TouchableHistory;
class G4Step;
class G4HCofThisEvent;
class G4MaterialPropertiesTable;
class G4PhysicsOrderedFreeVector;

class lArTPCSD : public G4VSensitiveDetector
{
 public:
  lArTPCSD(G4String name);
  virtual ~lArTPCSD() = default;
  void Initialize(G4HCofThisEvent* hitCollection) final;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) final;
  void EndOfEvent(G4HCofThisEvent* hitCollection) final;

 private:
  G4int materialIndex;
  const G4Material* aMaterial;
  G4MaterialPropertiesTable* aMaterialPropertiesTable;
  //
  // properties related to Scintillation
  //
  G4MaterialPropertyVector* Fast_Intensity{ nullptr };
  G4MaterialPropertyVector* Slow_Intensity{ nullptr };
  G4double YieldRatio{ 0 };        // slowerRatio,
  G4double FastTimeConstant{ 0 };  // TimeConstant,
  G4double SlowTimeConstant{ 0 };  // slowerTimeConstant,
  G4ScintillationType ScintillationType;
  //
  // properties related to Cerenkov
  //
  G4MaterialPropertyVector* Rindex{ nullptr };
  G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals{ nullptr };
  const G4PhysicsTable* thePhysicsTable{ 0 };
  G4double Pmin{ 0 };
  G4double Pmax{ 0 };
  G4double dp{ 0 };
  G4double nMax{ 0 };
  G4bool first{ false };
  G4bool verbose{ false };
  G4int tCphotons{ 0 };
  G4int tSphotons{ 0 };
  G4double NumElectrons(G4double e, G4double ds);
  lArTPCHitsCollection* flArTPCHitsCollection{ nullptr };
  G4int fHCID{ 0 };
};
