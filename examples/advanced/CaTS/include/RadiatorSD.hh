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
/// \file RadiatorSD.hh
/// \brief Definition of the CaTS::RadiatorSD class

#pragma once

#include "G4VSensitiveDetector.hh"
#include "G4ScintillationTrackInformation.hh"
#include <G4MaterialPropertyVector.hh>
#include <G4String.hh>
#include <G4Types.hh>
class G4Step;
class G4Material;
class G4HCofThisEvent;
class G4MaterialPropertiesTable;
class G4PhysicsOrderedFreeVector;
class G4PhysicsTable;
class G4TouchableHistory;

class RadiatorSD : public G4VSensitiveDetector
{
 public:
  RadiatorSD(G4String name);
  virtual ~RadiatorSD() = default;
  // methods from base class
  void Initialize(G4HCofThisEvent* hitCollection) final;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) final;
  void EndOfEvent(G4HCofThisEvent* hitCollection) final;

 private:
  const G4Material* aMaterial;
  G4MaterialPropertiesTable* aMaterialPropertiesTable;
  //
  // properties related to Scintillation
  //
  G4MaterialPropertyVector* Fast_Intensity;
  G4MaterialPropertyVector* Slow_Intensity;
  G4double YieldRatio;        // slowerRatio,
  G4double FastTimeConstant;  // TimeConstant,
  G4double SlowTimeConstant;  // slowerTimeConstant,
  G4ScintillationType ScintillationType;
  //
  // properties related to Cerenkov
  //
  G4MaterialPropertyVector* Rindex;
  G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals;
  const G4PhysicsTable* thePhysicsTable;
  G4double Pmin{ 0 };
  G4double Pmax{ 0 };
  G4double dp{ 0 };
  G4double nMax{ 0 };
  G4bool first{ true };
  G4bool verbose{ false };
  G4int tCphotons{ 0 };
  G4int tSphotons{ 0 };
#ifdef WITH_G4OPTICKS
  //
  // info needed for generating Cerenkov photons on the GPU;
  //
  G4double maxCos{ 0.0 };
  G4double maxSin2{ 0.0 };
  G4double beta{ 0.0 };
  G4double beta1{ 0.0 };
  G4double beta2{ 0.0 };
  G4double BetaInverse{ 0.0 };
  G4double MeanNumberOfPhotons1{ 0.0 };
  G4double MeanNumberOfPhotons2{ 0.0 };
  G4int Sphotons{ 0 };  // number of scintillation photons this step
  G4int Cphotons{ 0 };  // number of Cerenkov photons this step
  const G4double ScintillationTime{ 0.0 };
  const G4int scntId{ 1 };
#endif
};
