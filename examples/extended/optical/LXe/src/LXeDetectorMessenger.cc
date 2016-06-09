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

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeDetectorMessenger::LXeDetectorMessenger(LXeDetectorConstruction* LXeDetect)
:LXeDetector(LXeDetect)
{
  //Setup a command directory for detector controls with guidance
  detectorDir = new G4UIdirectory("/LXe/detector/");
  detectorDir->SetGuidance("Detector geometry control");

  volumesDir = new G4UIdirectory("/LXe/detector/volumes/");
  volumesDir->SetGuidance("Enable/disable volumes");
   
  //Various commands for modifying detector geometry
  dimensionsCmd = 
    new G4UIcmdWith3VectorAndUnit("/LXe/detector/dimensions",this);
  dimensionsCmd->SetGuidance("Set the dimensions of the detector volume.");
  dimensionsCmd->SetParameterName("scint_x","scint_y","scint_z",false);
  dimensionsCmd->SetDefaultUnit("cm");

  housingThicknessCmd = new G4UIcmdWithADoubleAndUnit
    ("/LXe/detector/housingThickness",this);
  housingThicknessCmd->SetGuidance("Set the thickness of the housing.");
  housingThicknessCmd->SetParameterName("d_mtl",false);
  housingThicknessCmd->SetDefaultUnit("cm");

  pmtRadiusCmd = new G4UIcmdWithADoubleAndUnit
    ("/LXe/detector/pmtRadius",this);
  pmtRadiusCmd->SetGuidance("Set the radius of the PMTs.");
  pmtRadiusCmd->SetParameterName("radius",false);
  pmtRadiusCmd->SetDefaultUnit("cm");

  nxCmd = new G4UIcmdWithAnInteger("/LXe/detector/nx",this);
  nxCmd->SetGuidance("Set the number of PMTs along the x-dimension.");
  nxCmd->SetParameterName("nx",false);

  nyCmd = new G4UIcmdWithAnInteger("/LXe/detector/ny",this);
  nyCmd->SetGuidance("Set the number of PMTs along the y-dimension.");
  nyCmd->SetParameterName("ny",false);
  
  nzCmd = new G4UIcmdWithAnInteger("/LXe/detector/nz",this);
  nzCmd->SetGuidance("Set the number of PMTs along the z-dimension.");
  nzCmd->SetParameterName("nz",false);

  sphereCmd = new G4UIcmdWithABool("/LXe/detector/volumes/sphere",this);
  sphereCmd->SetGuidance("Enable/Disable the sphere."); 

  reflectivityCmd = new G4UIcmdWithADouble("/LXe/detector/reflectivity",this);
  reflectivityCmd->SetGuidance("Set the reflectivity of the housing.");

  wlsCmd = new G4UIcmdWithABool("/LXe/detector/volumes/wls",this);
  wlsCmd->SetGuidance("Enable/Disable the WLS slab");

  lxeCmd = new G4UIcmdWithABool("/LXe/detector/volumes/lxe",this);
  wlsCmd->SetGuidance("Enable/Disable the main detector volume.");

  nFibersCmd = new G4UIcmdWithAnInteger("/LXe/detector/nfibers",this);
  nFibersCmd->SetGuidance("Set the number of WLS fibers in the WLS slab.");

  updateCmd = new G4UIcommand("/LXe/detector/update",this);
  updateCmd->SetGuidance("Update the detector geometry with changed values.");
  updateCmd->SetGuidance
    ("Must be run before beamOn if detector has been changed.");
  
  defaultsCmd = new G4UIcommand("/LXe/detector/defaults",this);
  defaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  defaultsCmd->SetGuidance("(Update still required)");

  MainScintYield=new G4UIcmdWithADouble("/LXe/detector/MainScintYield",this);
  MainScintYield->SetGuidance("Set scinitillation yield of main volume.");
  MainScintYield->SetGuidance("Specified in photons/MeV");

  WLSScintYield = new G4UIcmdWithADouble("/LXe/detector/WLSScintYield",this);
  WLSScintYield->SetGuidance("Set scintillation yield of WLS Slab");
  WLSScintYield->SetGuidance("Specified in photons/MeV");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeDetectorMessenger::~LXeDetectorMessenger()
{
  delete dimensionsCmd;
  delete housingThicknessCmd;
  delete pmtRadiusCmd;
  delete nxCmd;
  delete nyCmd;
  delete nzCmd;
  delete updateCmd;
  delete detectorDir;
  delete volumesDir;
  delete defaultsCmd;
  delete sphereCmd;
  delete wlsCmd;
  delete lxeCmd;
  delete nFibersCmd;
  delete reflectivityCmd;
  delete MainScintYield;
  delete WLSScintYield;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == dimensionsCmd ){ 
    LXeDetector->SetDimensions(dimensionsCmd->GetNew3VectorValue(newValue));
  }
  else if (command == housingThicknessCmd){
    LXeDetector->SetHousingThickness(housingThicknessCmd
				     ->GetNewDoubleValue(newValue));
  }
  else if (command == pmtRadiusCmd){
    LXeDetector->SetPMTRadius(pmtRadiusCmd->GetNewDoubleValue(newValue));
  }
  else if (command == nxCmd){
    LXeDetector->SetNX(nxCmd->GetNewIntValue(newValue));
  }
  else if (command == nyCmd){
    LXeDetector->SetNY(nyCmd->GetNewIntValue(newValue));
  }
  else if (command == nzCmd){
    LXeDetector->SetNZ(nzCmd->GetNewIntValue(newValue));
  }
  else if (command == updateCmd){
    LXeDetector->UpdateGeometry();
  }
  else if (command == defaultsCmd){
    LXeDetector->SetDefaults();
  }
  else if (command == sphereCmd){
    LXeDetector->SetSphereOn(sphereCmd->GetNewBoolValue(newValue));
  }
  else if (command == reflectivityCmd){
    LXeDetector
      ->SetHousingReflectivity(reflectivityCmd->GetNewDoubleValue(newValue));
  }
  else if (command == wlsCmd){
    LXeDetector->SetWLSSlabOn(wlsCmd->GetNewBoolValue(newValue));
  }
  else if (command == lxeCmd){
    LXeDetector->SetMainVolumeOn(lxeCmd->GetNewBoolValue(newValue));
  }
  else if (command == nFibersCmd){
    LXeDetector->SetNFibers(nFibersCmd->GetNewIntValue(newValue));
  }
  else if (command == MainScintYield){
   LXeDetector->SetMainScintYield(MainScintYield->GetNewDoubleValue(newValue));
  }
  else if (command == WLSScintYield){
    LXeDetector->SetWLSScintYield(WLSScintYield->GetNewDoubleValue(newValue));
  }
}



