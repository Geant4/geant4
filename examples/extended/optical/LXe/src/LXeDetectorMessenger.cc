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
/// \file optical/LXe/src/LXeDetectorMessenger.cc
/// \brief Implementation of the LXeDetectorMessenger class
//
//
#include "LXeDetectorMessenger.hh"
#include "LXeDetectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Scintillation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorMessenger::LXeDetectorMessenger(LXeDetectorConstruction* detector)
 : fLXeDetector(detector)
{
  //Setup a command directory for detector controls with guidance
  fDetectorDir = new G4UIdirectory("/LXe/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fVolumesDir = new G4UIdirectory("/LXe/detector/volumes/");
  fVolumesDir->SetGuidance("Enable/disable volumes");
 
  //Various commands for modifying detector geometry
  fDimensionsCmd =
    new G4UIcmdWith3VectorAndUnit("/LXe/detector/dimensions",this);
  fDimensionsCmd->SetGuidance("Set the dimensions of the detector volume.");
  fDimensionsCmd->SetParameterName("scint_x","scint_y","scint_z",false);
  fDimensionsCmd->SetDefaultUnit("cm");

  fHousingThicknessCmd = new G4UIcmdWithADoubleAndUnit
    ("/LXe/detector/housingThickness",this);
  fHousingThicknessCmd->SetGuidance("Set the thickness of the housing.");
  fHousingThicknessCmd->SetParameterName("d_mtl",false);
  fHousingThicknessCmd->SetDefaultUnit("cm");

  fPmtRadiusCmd = new G4UIcmdWithADoubleAndUnit
    ("/LXe/detector/pmtRadius",this);
  fPmtRadiusCmd->SetGuidance("Set the radius of the PMTs.");
  fPmtRadiusCmd->SetParameterName("radius",false);
  fPmtRadiusCmd->SetDefaultUnit("cm");

  fNxCmd = new G4UIcmdWithAnInteger("/LXe/detector/nx",this);
  fNxCmd->SetGuidance("Set the number of PMTs along the x-dimension.");
  fNxCmd->SetParameterName("nx",false);

  fNyCmd = new G4UIcmdWithAnInteger("/LXe/detector/ny",this);
  fNyCmd->SetGuidance("Set the number of PMTs along the y-dimension.");
  fNyCmd->SetParameterName("ny",false);

  fNzCmd = new G4UIcmdWithAnInteger("/LXe/detector/nz",this);
  fNzCmd->SetGuidance("Set the number of PMTs along the z-dimension.");
  fNzCmd->SetParameterName("nz",false);

  fSphereCmd = new G4UIcmdWithABool("/LXe/detector/volumes/sphere",this);
  fSphereCmd->SetGuidance("Enable/Disable the sphere.");

  fReflectivityCmd = new G4UIcmdWithADouble("/LXe/detector/reflectivity",this);
  fReflectivityCmd->SetGuidance("Set the reflectivity of the housing.");

  fWlsCmd = new G4UIcmdWithABool("/LXe/detector/volumes/wls",this);
  fWlsCmd->SetGuidance("Enable/Disable the WLS slab");

  fLxeCmd = new G4UIcmdWithABool("/LXe/detector/volumes/lxe",this);
  fLxeCmd->SetGuidance("Enable/Disable the main detector volume.");

  fNFibersCmd = new G4UIcmdWithAnInteger("/LXe/detector/nfibers",this);
  fNFibersCmd->SetGuidance("Set the number of WLS fibers in the WLS slab.");

  fUpdateCmd = new G4UIcommand("/LXe/detector/update",this);
  fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
  fUpdateCmd->SetGuidance
    ("Must be run before beamOn if detector has been changed.");

  fDefaultsCmd = new G4UIcommand("/LXe/detector/defaults",this);
  fDefaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  fDefaultsCmd->SetGuidance("(Update still required)");

  fMainScintYield=new G4UIcmdWithADouble("/LXe/detector/MainScintYield",this);
  fMainScintYield->SetGuidance("Set scinitillation yield of main volume.");
  fMainScintYield->SetGuidance("Specified in photons/MeV");

  fWLSScintYield = new G4UIcmdWithADouble("/LXe/detector/WLSScintYield",this);
  fWLSScintYield->SetGuidance("Set scintillation yield of WLS Slab");
  fWLSScintYield->SetGuidance("Specified in photons/MeV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorMessenger::~LXeDetectorMessenger()
{
  delete fDimensionsCmd;
  delete fHousingThicknessCmd;
  delete fPmtRadiusCmd;
  delete fNxCmd;
  delete fNyCmd;
  delete fNzCmd;
  delete fUpdateCmd;
  delete fDetectorDir;
  delete fVolumesDir;
  delete fDefaultsCmd;
  delete fSphereCmd;
  delete fWlsCmd;
  delete fLxeCmd;
  delete fNFibersCmd;
  delete fReflectivityCmd;
  delete fMainScintYield;
  delete fWLSScintYield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fDimensionsCmd ){
    fLXeDetector->SetDimensions(fDimensionsCmd->GetNew3VectorValue(newValue));
  }
  else if (command == fHousingThicknessCmd){
    fLXeDetector->SetHousingThickness(fHousingThicknessCmd
                                     ->GetNewDoubleValue(newValue));
  }
  else if (command == fPmtRadiusCmd){
    fLXeDetector->SetPMTRadius(fPmtRadiusCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fNxCmd){
    fLXeDetector->SetNX(fNxCmd->GetNewIntValue(newValue));
  }
  else if (command == fNyCmd){
    fLXeDetector->SetNY(fNyCmd->GetNewIntValue(newValue));
  }
  else if (command == fNzCmd){
    fLXeDetector->SetNZ(fNzCmd->GetNewIntValue(newValue));
  }
  else if (command == fUpdateCmd){
    fLXeDetector->UpdateGeometry();
  }
  else if (command == fDefaultsCmd){
    fLXeDetector->SetDefaults();
  }
  else if (command == fSphereCmd){
    fLXeDetector->SetSphereOn(fSphereCmd->GetNewBoolValue(newValue));
  }
  else if (command == fReflectivityCmd){
    fLXeDetector
      ->SetHousingReflectivity(fReflectivityCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWlsCmd){
    fLXeDetector->SetWLSSlabOn(fWlsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fLxeCmd){
    fLXeDetector->SetMainVolumeOn(fLxeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fNFibersCmd){
    fLXeDetector->SetNFibers(fNFibersCmd->GetNewIntValue(newValue));
  }
  else if (command == fMainScintYield){
   fLXeDetector->SetMainScintYield(fMainScintYield->GetNewDoubleValue(newValue));
  }
  else if (command == fWLSScintYield){
    fLXeDetector->SetWLSScintYield(fWLSScintYield->GetNewDoubleValue(newValue));
  }
}
