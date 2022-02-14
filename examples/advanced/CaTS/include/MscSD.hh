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
/// \file MscSD.hh
/// \brief Definition of the CaTS::MscSD class

#pragma once

#include <G4String.hh>
#include <G4Types.hh>
#include "G4VSensitiveDetector.hh"
#include "MscHit.hh"
class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class MscSD : public G4VSensitiveDetector
{
 public:
  MscSD(G4String);
  ~MscSD() = default;
  void Initialize(G4HCofThisEvent*) final;
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) final;
  void EndOfEvent(G4HCofThisEvent* hitCollection) final;

 private:
  MscHitsCollection* fMscHitsCollection{ nullptr };
  G4int fHCID{ 0 };
  G4bool verbose{ false };
  G4bool done{ false };
};
