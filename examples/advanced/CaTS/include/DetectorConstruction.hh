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
/// \file DetectorConstruction.hh
/// \brief Definition of the CaTS::DetectorConstruction class

#pragma once

#include "G4VUserDetectorConstruction.hh"
#include <G4String.hh>
class G4VPhysicalVolume;
class ColorReader;
class G4GDMLParser;

class DetectorConstruction final : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction(G4String fname);
  ~DetectorConstruction() final;
  DetectorConstruction& operator=(const DetectorConstruction& right) = delete;
  DetectorConstruction(const DetectorConstruction&)                  = delete;
  void ReadGDML();
  G4VPhysicalVolume* Construct() final;
  void ConstructSDandField() final;
  void UpdateGeometry();

 private:
  G4String gdmlFile;
  G4GDMLParser* parser{ nullptr };
  ColorReader* fReader{ nullptr };
  G4bool verbose{ false };
};
