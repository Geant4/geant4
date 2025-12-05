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
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the CaTS::DetectorConstruction class

// Geant4 headers
#include "G4RunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"
#include "G4GDMLParser.hh"
#include "G4SDManager.hh"
// project headers
#include "ConfigurationManager.hh"
#include "DetectorConstruction.hh"
#include "ColorReader.hh"
// c++ headers
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction(G4String fname)
  : G4VUserDetectorConstruction()
  , gdmlFile(fname)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Import the application geometry from a GDML file
  ReadGDML();

  const G4GDMLAuxMapType* auxmap = parser->GetAuxMap();
  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  if(verbose)
  {
    G4cout << "Found " << auxmap->size()
	   << " volume(s) with auxiliary information." << G4endl;
  }
  for(G4GDMLAuxMapType::const_iterator iter = auxmap->begin();
      iter != auxmap->end(); iter++)
  {
    if(verbose)
    {
      G4cout << "Volume " << ((*iter).first)->GetName()
             << " has the following list of auxiliary information: " << G4endl;
    }
    for(G4GDMLAuxListType::const_iterator vit = (*iter).second.begin();
        vit != (*iter).second.end(); vit++)
    {
      if(verbose)
      {
        G4cout << "Type: " << (*vit).type << " Value: " << (*vit).value
	       << G4endl;
      }
      if((*vit).type == "StepLimit")
      {
        G4UserLimits* fStepLimit = new G4UserLimits(atof((*vit).value));
        ((*iter).first)->SetUserLimits(fStepLimit);
      }
    }
  }
  G4VPhysicalVolume* worldPhysVol = parser->GetWorldVolume();
  if(ConfigurationManager::getInstance()->isDumpgdml())
  {
    std::ifstream ifile;
    G4String fileName = ConfigurationManager::getInstance()->getGDMLFileName();
    ifile.open(fileName);
    if(ifile)
    {
      G4cout << "**************************************************" << G4endl;
      G4cout << fileName << " already exists!!!" << G4endl;
      G4cout << "No new gdml dump created!!!" << G4endl;
      G4cout << "**************************************************" << G4endl;
    }
    else
    {
      G4cout << "Writing: " << fileName << G4endl;
      parser->Write(fileName, worldPhysVol);
    }
  }
  return worldPhysVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  const G4GDMLAuxMapType* auxmap = parser->GetAuxMap();
  if(verbose)
  {
    G4cout << "Found " << auxmap->size()
           << " volume(s) with auxiliary information." << G4endl;
  }

  // Add a sensitive detector into the SD manager
  for(G4GDMLAuxMapType::const_iterator iter = auxmap->begin();
      iter != auxmap->end(); iter++)
  {
    G4LogicalVolume* volume = (*iter).first;
    G4String volName = volume->GetName();

    if(verbose)
    {
      G4cout << "Volume " << volName
             << " has the following list of auxiliary information: " << G4endl;
    }
    for(G4GDMLAuxListType::const_iterator vit = (*iter).second.begin();
        vit != (*iter).second.end(); vit++)
    {
      G4String auxType = (*vit).type;
      G4String auxValue = (*vit).value;

      if(verbose)
      {
        G4cout << "Aux Type: " << auxType << " Value: " << auxValue << G4endl;
      }
      if(auxType == "SensDet")
      {
        if(verbose)
        {
          G4cout << "Found sensitive Detector: " << auxValue << G4endl;
          G4cout << "Attaching it to Volume:  " << volName << G4endl;
        }

        const auto& sdMap = GetSDMap();

        if (auto it = sdMap.find(auxValue); it != sdMap.end())
        {
	  // construct SD with given name
          auto sd = it->second(volume->GetName() + "_" + auxValue);
          G4SDManager::GetSDMpointer()->AddNewDetector(sd.get());
          volume->SetSensitiveDetector(sd.release());
        }
        else
        {
          G4cerr << "Unknown sensitive detector type: " << auxValue << G4endl;
        }
      }
      else if(auxType == "Solid" && auxValue == "True")
      {
        G4VisAttributes* visibility = new G4VisAttributes();
        visibility->SetForceSolid(true);
        G4VisAttributes* visatt = new G4VisAttributes(
          volume->GetVisAttributes()->GetColour());
        visatt->SetVisibility(true);
        visatt->SetForceSolid(true);
        visatt->SetForceAuxEdgeVisible(true);
        volume->SetVisAttributes(visatt);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ReadGDML()
{
  fReader = new ColorReader;
  parser  = new G4GDMLParser(fReader);
  parser->Read(gdmlFile, false);

  G4VPhysicalVolume* World = parser->GetWorldVolume();
  //----- GDML parser makes world invisible, this is a hack to make it
  // visible again...
  G4LogicalVolume* pWorldLogical = World->GetLogicalVolume();
  pWorldLogical->SetVisAttributes(0);
  G4cout << World->GetTranslation() << G4endl << G4endl;

  if(verbose)
  {
    G4cout << "Found World:  " << World->GetName() << G4endl;
    G4cout << "World LV:  " << World->GetLogicalVolume()->GetName() << G4endl;
    G4cout << "Found " << G4LogicalVolumeStore::GetInstance()->size()
	   << " logical volumes." << G4endl;
    G4cout << "Found " << G4PhysicalVolumeStore::GetInstance()->size()
	   << " physical volumes." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
