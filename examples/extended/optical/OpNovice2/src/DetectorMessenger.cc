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
/// \file optical/OpNovice2/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 77288 2013-11-22 10:52:58Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include <sstream>
#include <iostream>

#include "G4OpticalSurface.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det)
{ 
  fOpticalDir = new G4UIdirectory("/opnovice2/");
  fOpticalDir->SetGuidance("Parameters for optical simulation.");

  fSurfaceTypeCmd = new G4UIcmdWithAString("/opnovice2/surfaceType", this);
  fSurfaceTypeCmd->SetGuidance("Surface type.");
  fSurfaceTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceTypeCmd->SetToBeBroadcasted(false);
  
  fSurfaceFinishCmd = new G4UIcmdWithAString("/opnovice2/surfaceFinish", this);
  fSurfaceFinishCmd->SetGuidance("Surface finish.");
  fSurfaceFinishCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceFinishCmd->SetToBeBroadcasted(false);
  
  fSurfaceModelCmd =
    new G4UIcmdWithAString("/opnovice2/surfaceModel", this);
  fSurfaceModelCmd->SetGuidance("surface model.");
  fSurfaceModelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceModelCmd->SetToBeBroadcasted(false);
 
  fSurfaceSigmaAlphaCmd =
    new G4UIcmdWithADouble("/opnovice2/surfaceSigmaAlpha", this);
  fSurfaceSigmaAlphaCmd->SetGuidance("surface sigma alpha");
  fSurfaceSigmaAlphaCmd->SetGuidance(" parameter.");
  fSurfaceSigmaAlphaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceSigmaAlphaCmd->SetToBeBroadcasted(false);

  fSurfaceMatPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/surfaceProperty", this);
  fSurfaceMatPropVectorCmd->SetGuidance("Set material property vector");
  fSurfaceMatPropVectorCmd->SetGuidance(" for the surface.");
  fSurfaceMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceMatPropVectorCmd->SetToBeBroadcasted(false);

  fBoxMatPropVectorCmd = new G4UIcmdWithAString("/opnovice2/boxProperty", this);
  fBoxMatPropVectorCmd->SetGuidance("Set material property vector for ");
  fBoxMatPropVectorCmd->SetGuidance("the box.");
  fBoxMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBoxMatPropVectorCmd->SetToBeBroadcasted(false);

  fBoxMatConstPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/boxConstProperty", this);
  fBoxMatConstPropVectorCmd->SetGuidance("Set material constant property ");
  fBoxMatConstPropVectorCmd->SetGuidance("for the box.");
  fBoxMatConstPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBoxMatConstPropVectorCmd->SetToBeBroadcasted(false);

  fWorldMatPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/worldProperty", this);
  fWorldMatPropVectorCmd->SetGuidance("Set material property vector ");
  fWorldMatPropVectorCmd->SetGuidance("for the world.");
  fWorldMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWorldMatPropVectorCmd->SetToBeBroadcasted(false);

  fWorldMatConstPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/worldConstProperty", this);
  fWorldMatConstPropVectorCmd->SetGuidance("Set material constant property");
  fWorldMatConstPropVectorCmd->SetGuidance(" for the world.");
  fWorldMatConstPropVectorCmd->
    AvailableForStates(G4State_PreInit, G4State_Idle);
  fWorldMatConstPropVectorCmd->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fSurfaceFinishCmd;
  delete fSurfaceTypeCmd;
  delete fSurfaceModelCmd;
  delete fSurfaceSigmaAlphaCmd;
  delete fSurfaceMatPropVectorCmd;
  delete fBoxMatPropVectorCmd;
  delete fBoxMatConstPropVectorCmd;
  delete fWorldMatPropVectorCmd;
  delete fWorldMatConstPropVectorCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{    
  //    FINISH              
  if (command == fSurfaceFinishCmd) {
    if (newValue == "polished") {
      fDetector->SetSurfaceFinish(polished);
    }
    else if (newValue == "polishedfrontpainted") {
      fDetector->SetSurfaceFinish(polishedfrontpainted);
    }
    else if (newValue == "polishedbackpainted") {
      fDetector->SetSurfaceFinish(polishedbackpainted);
    }
    else if (newValue == "ground") {
      fDetector->SetSurfaceFinish(ground);
    }
    else if (newValue == "groundfrontpainted") {
      fDetector->SetSurfaceFinish(groundfrontpainted);
    }
    else if (newValue == "groundbackpainted") {
      fDetector->SetSurfaceFinish(groundbackpainted);
    }
    else if (newValue == "polishedlumirrorair") {
      fDetector->SetSurfaceFinish(polishedlumirrorair);
    }
    else if (newValue == "polishedlumirrorglue") {
      fDetector->SetSurfaceFinish(polishedlumirrorglue);
    }
    else if (newValue == "polishedair") {
      fDetector->SetSurfaceFinish(polishedair);
    }
    else if (newValue == "polishedteflonair") {
      fDetector->SetSurfaceFinish(polishedteflonair);
    }
    else if (newValue == "polishedtioair") {
      fDetector->SetSurfaceFinish(polishedtioair);
    }
    else if (newValue == "polishedtyvekair") {
      fDetector->SetSurfaceFinish(polishedtyvekair);
    }
    else if (newValue == "polishedvm2000air") {
      fDetector->SetSurfaceFinish(polishedvm2000air);
    }
    else if (newValue == "polishedvm2000glue") {
      fDetector->SetSurfaceFinish(polishedvm2000glue);
    }
    else if (newValue == "etchedlumirrorair") {
      fDetector->SetSurfaceFinish(etchedlumirrorair);
    }
    else if (newValue == "etchedlumirrorglue") {
      fDetector->SetSurfaceFinish(etchedlumirrorglue);
    }
    else if (newValue == "etchedair") {
      fDetector->SetSurfaceFinish(etchedair);
    }
    else if (newValue == "etchedteflonair") {
      fDetector->SetSurfaceFinish(etchedteflonair);
    }
    else if (newValue == "etchedtioair") {
      fDetector->SetSurfaceFinish(etchedtioair);
    }
    else if (newValue == "etchedtyvekair") {
      fDetector->SetSurfaceFinish(etchedtyvekair);
    }
    else if (newValue == "etchedvm2000air") {
      fDetector->SetSurfaceFinish(etchedvm2000air);
    }
    else if (newValue == "etchedvm2000glue") {
      fDetector->SetSurfaceFinish(etchedvm2000glue);
    }
    else if (newValue == "groundlumirrorair") {
      fDetector->SetSurfaceFinish(groundlumirrorair);
    }
    else if (newValue == "groundlumirrorglue") {
      fDetector->SetSurfaceFinish(groundlumirrorglue);
    }
    else if (newValue == "groundair") {
      fDetector->SetSurfaceFinish(groundair);
    }
    else if (newValue == "groundteflonair") {
      fDetector->SetSurfaceFinish(groundteflonair);
    }
    else if (newValue == "groundtioair") {
      fDetector->SetSurfaceFinish(groundtioair);
    }
    else if (newValue == "groundtyvekair") {
      fDetector->SetSurfaceFinish(groundtyvekair);
    }
    else if (newValue == "groundvm2000air") {
      fDetector->SetSurfaceFinish(groundvm2000air);
    }
    else if (newValue == "groundvm2000glue") {
      fDetector->SetSurfaceFinish(groundvm2000glue);
    }
    //         for Davis model
    else if (newValue == "Rough_LUT") {
      fDetector->SetSurfaceFinish(Rough_LUT);
    }
    else if (newValue == "RoughTeflon_LUT") {
      fDetector->SetSurfaceFinish(RoughTeflon_LUT);
    }
    else if (newValue == "RoughESR_LUT") {
      fDetector->SetSurfaceFinish(RoughESR_LUT);
    }
    else if (newValue == "RoughESRGrease_LUT") {
      fDetector->SetSurfaceFinish(RoughESRGrease_LUT);
    }
    else if (newValue == "Polished_LUT") {
      fDetector->SetSurfaceFinish(Polished_LUT);
    }
    else if (newValue == "PolishedTeflon_LUT") {
      fDetector->SetSurfaceFinish(PolishedTeflon_LUT);
    }
    else if (newValue == "PolishedESR_LUT") {
      fDetector->SetSurfaceFinish(PolishedESR_LUT);
    }
    else if (newValue == "PolishedESRGrease_LUT") {
      fDetector->SetSurfaceFinish(PolishedESRGrease_LUT);
    }
    else if (newValue == "Detector_LUT") {
      fDetector->SetSurfaceFinish(Detector_LUT);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Invalid surface finish: " << newValue;
      G4Exception("DetectorMessenger", "OpNovice2_003", FatalException,ed);
    }
  }

  //  MODEL
  else if (command == fSurfaceModelCmd) {
    if (newValue == "glisur") {
      fDetector->SetSurfaceModel(glisur);
    }
    else if (newValue == "unified") {
      fDetector->SetSurfaceModel(unified);
    }
    else if (newValue == "LUT") {
      fDetector->SetSurfaceModel(LUT);
    }
    else if (newValue == "DAVIS") {
      fDetector->SetSurfaceModel(DAVIS);
    }
    else if (newValue == "dichroic") {
      fDetector->SetSurfaceModel(dichroic);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Invalid surface model: " << newValue;
      G4Exception("DetectorMessenger", "ONovice2_001",
                  FatalException,ed);
    }
  }

  // TYPE
  else if (command == fSurfaceTypeCmd) {
    if (newValue == "dielectric_metal") {
      fDetector->SetSurfaceType(dielectric_metal);
    }
    else if (newValue == "dielectric_dielectric") {
      fDetector->SetSurfaceType(dielectric_dielectric);
    }
    else if (newValue == "dielectric_LUT") {
      fDetector->SetSurfaceType(dielectric_LUT);
    }
    else if (newValue == "dielectric_LUTDAVIS") {
      fDetector->SetSurfaceType(dielectric_LUTDAVIS);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Invalid surface type: " << newValue;
      G4Exception("DetectorMessenger", "OpNovice2_002", FatalException,ed);
    }
  }
  else if (command == fSurfaceSigmaAlphaCmd) {
    fDetector->SetSurfaceSigmaAlpha(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if (command == fBoxMatPropVectorCmd) {
    // got a string. need to convert it to physics vector.
    // string format is property name, then pairs of energy, value  
    // specify units for each value, eg 3.0*eV
    // space delimited
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    mpv->SetSpline(true);
    std::istringstream instring(newValue);
    G4String prop;
    instring >> prop;
    while (instring) {
      G4String tmp;
      instring >> tmp;
      if (tmp == "") { break; }
      G4double en = G4UIcommand::ConvertToDouble(tmp);
      instring >> tmp;
      G4double val;
      val = G4UIcommand::ConvertToDouble(tmp);
      mpv->InsertValues(en, val);
    } 
    const char* c = prop.c_str();

    fDetector->AddBoxMPV(c, mpv);
  }
  else if (command == fWorldMatPropVectorCmd) {
    // Convert string to physics vector
    // string format is property name, then pairs of energy, value  
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    std::istringstream instring(newValue);
    G4String prop;
    instring >> prop;
    while (instring) {
      G4String tmp;
      instring >> tmp;
      if (tmp == "") { break; }
      G4double en = G4UIcommand::ConvertToDouble(tmp);
      instring >> tmp;
      G4double val;
      val = G4UIcommand::ConvertToDouble(tmp);
      mpv->InsertValues(en, val);
    } 
    const char* c = prop.c_str();
    fDetector->AddWorldMPV(c, mpv);
  }
  else if (command == fSurfaceMatPropVectorCmd) {
    // Convert string to physics vector
    // string format is property name, then pairs of energy, value  
    // space delimited
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    G4cout << newValue << G4endl;
    std::istringstream instring(newValue);
    G4String prop;
    instring >> prop;
    while (instring) {
      G4String tmp;
      instring >> tmp;
      if (tmp == "") { break; }
      G4double en = G4UIcommand::ConvertToDouble(tmp);
      instring >> tmp;
      G4double val;
      val = G4UIcommand::ConvertToDouble(tmp);
      mpv->InsertValues(en, val);
    } 
    const char* c = prop.c_str();
    fDetector->AddSurfaceMPV(c, mpv);
  }

  else if (command == fBoxMatConstPropVectorCmd) {
    // Convert string to physics vector
    // string format is property name, then value
    // space delimited
    std::istringstream instring(newValue);
    G4String prop;
    G4String tmp;
    instring >> prop;
    instring >> tmp;
    G4double val = G4UIcommand::ConvertToDouble(tmp);
    const char* c = prop.c_str();
    fDetector->AddBoxMPCV(c, val);
  }
  else if (command == fWorldMatConstPropVectorCmd) {
    // Convert string to physics vector
    // string format is property name, then value
    // space delimited
    std::istringstream instring(newValue);
    G4String prop;
    G4String tmp;
    instring >> prop;
    instring >> tmp;
    G4double val = G4UIcommand::ConvertToDouble(tmp);
    const char* c = prop.c_str();
    fDetector->AddBoxMPCV(c, val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
