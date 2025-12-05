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
// CaTS (Calorimetry and Tracking Simulation)
//
// Authors: Hans Wenzel and Soon Yung Jun
//          (Fermi National Accelerator Laboratory)
//
// History: October 18th, 2021 : first implementation
//
// ********************************************************************

/// \file DetectorConstruction.hh
/// \brief Definition of the CaTS::DetectorConstruction class

#pragma once

#include "G4String.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VUserDetectorConstruction.hh"

#include "TrackerSD.hh"
#include "MscSD.hh"
#include "lArTPCSD.hh"
#include "CalorimeterSD.hh"
#include "DRCalorimeterSD.hh"
#include "RadiatorSD.hh"
#include "PhotonSD.hh"
#include "InteractionSD.hh"

#include <memory>
#include <string>
#include <unordered_map>

class ColorReader;
class G4GDMLParser;
class G4VPhysicalVolume;

class DetectorConstruction final : public G4VUserDetectorConstruction
{
 public:
  explicit DetectorConstruction(G4String fname);
  ~DetectorConstruction() final {};
  DetectorConstruction& operator=(const DetectorConstruction& right) = delete;
  DetectorConstruction(const DetectorConstruction&) = delete;

  G4VPhysicalVolume* Construct() final;
  void ConstructSDandField() final;

 private:
  void ReadGDML();
  void UpdateGeometry();

  // The map of sensitive detectors
  using SDMap = std::unordered_map<std::string,
    std::function<std::unique_ptr<G4VSensitiveDetector>(const G4String&)>>;

  inline const SDMap& GetSDMap()
  {
    static const SDMap sdMap = {
      {"PhotonDetector", [](const G4String& name)
        { return std::make_unique<PhotonSD>(name); }},
      {"Target", [](const G4String& name)
        { return std::make_unique<InteractionSD>(name); }},
      {"Tracker", [](const G4String& name)
        { return std::make_unique<TrackerSD>(name); }},
      {"Msc", [](const G4String& name)
        { return std::make_unique<MscSD>(name); }},
      {"lArTPC", [](const G4String& name)
        { return std::make_unique<lArTPCSD>(name); }},
      {"Radiator", [](const G4String& name)
        { return std::make_unique<RadiatorSD>(name); }},
      {"Calorimeter", [](const G4String& name)
        { return std::make_unique<CalorimeterSD>(name); }},
      {"DRCalorimeter", [](const G4String& name)
        { return std::make_unique<DRCalorimeterSD>(name); }}
    };
    return sdMap;
  }

private:
  G4bool verbose{ false };
  G4String gdmlFile;
  G4GDMLParser* parser{ nullptr };
  ColorReader* fReader{ nullptr };
};
