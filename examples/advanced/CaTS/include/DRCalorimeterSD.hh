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
/// \file DRCalorimeterSD.hh
/// \brief Definition of the CaTS::DRCalorimeterSD class

#pragma once

#include "G4VSensitiveDetector.hh"
#include "DRCalorimeterHit.hh"
#include <G4String.hh>
#include <G4Types.hh>
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class DRCalorimeterSD : public G4VSensitiveDetector
{
 public:
  DRCalorimeterSD(G4String);
  ~DRCalorimeterSD();
  void Initialize(G4HCofThisEvent*) final;
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) final;
  void EndOfEvent(G4HCofThisEvent* hitCollection) final;

 private:
  DRCalorimeterHitsCollection* fDRCalorimeterHitsCollection;
  G4int fHCID{ 0 };
  G4bool verbose{ false };
};
