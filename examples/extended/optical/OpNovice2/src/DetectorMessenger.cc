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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4OpticalSurface.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include <sstream>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
  : G4UImessenger()
  , fDetector(Det)
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

  fSurfaceModelCmd = new G4UIcmdWithAString("/opnovice2/surfaceModel", this);
  fSurfaceModelCmd->SetGuidance("surface model.");
  fSurfaceModelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceModelCmd->SetToBeBroadcasted(false);

  fSurfaceSigmaAlphaCmd =
    new G4UIcmdWithADouble("/opnovice2/surfaceSigmaAlpha", this);
  fSurfaceSigmaAlphaCmd->SetGuidance("surface sigma alpha");
  fSurfaceSigmaAlphaCmd->SetGuidance(" parameter.");
  fSurfaceSigmaAlphaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceSigmaAlphaCmd->SetToBeBroadcasted(false);

  fSurfacePolishCmd = new G4UIcmdWithADouble("/opnovice2/surfacePolish", this);
  fSurfacePolishCmd->SetGuidance("surface polish");
  fSurfacePolishCmd->SetGuidance(" parameter (for Glisur model).");
  fSurfacePolishCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfacePolishCmd->SetToBeBroadcasted(false);

  fSurfaceMatPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/surfaceProperty", this);
  fSurfaceMatPropVectorCmd->SetGuidance("Set material property vector");
  fSurfaceMatPropVectorCmd->SetGuidance(" for the surface.");
  fSurfaceMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceMatPropVectorCmd->SetToBeBroadcasted(false);

  fSurfaceMatPropConstCmd =
    new G4UIcmdWithAString("/opnovice2/surfaceConstProperty", this);
  fSurfaceMatPropConstCmd->SetGuidance("Set material constant property");
  fSurfaceMatPropConstCmd->SetGuidance(" for the surface.");
  fSurfaceMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSurfaceMatPropConstCmd->SetToBeBroadcasted(false);

  fTankMatPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/boxProperty", this);
  fTankMatPropVectorCmd->SetGuidance("Set material property vector for ");
  fTankMatPropVectorCmd->SetGuidance("the box.");
  fTankMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fTankMatPropVectorCmd->SetToBeBroadcasted(false);

  fTankMatPropConstCmd =
    new G4UIcmdWithAString("/opnovice2/boxConstProperty", this);
  fTankMatPropConstCmd->SetGuidance("Set material constant property ");
  fTankMatPropConstCmd->SetGuidance("for the box.");
  fTankMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fTankMatPropConstCmd->SetToBeBroadcasted(false);

  fTankMaterialCmd = new G4UIcmdWithAString("/opnovice2/boxMaterial", this);
  fTankMaterialCmd->SetGuidance("Set material of box.");
  fTankMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fTankMaterialCmd->SetToBeBroadcasted(false);

  fWorldMatPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/worldProperty", this);
  fWorldMatPropVectorCmd->SetGuidance("Set material property vector ");
  fWorldMatPropVectorCmd->SetGuidance("for the world.");
  fWorldMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWorldMatPropVectorCmd->SetToBeBroadcasted(false);

  fWorldMatPropConstCmd =
    new G4UIcmdWithAString("/opnovice2/worldConstProperty", this);
  fWorldMatPropConstCmd->SetGuidance("Set material constant property");
  fWorldMatPropConstCmd->SetGuidance(" for the world.");
  fWorldMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWorldMatPropConstCmd->SetToBeBroadcasted(false);

  fWorldMaterialCmd = new G4UIcmdWithAString("/opnovice2/worldMaterial", this);
  fWorldMaterialCmd->SetGuidance("Set material of world.");
  fWorldMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWorldMaterialCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fOpticalDir;
  delete fSurfaceFinishCmd;
  delete fSurfaceTypeCmd;
  delete fSurfaceModelCmd;
  delete fSurfaceSigmaAlphaCmd;
  delete fSurfacePolishCmd;
  delete fSurfaceMatPropVectorCmd;
  delete fSurfaceMatPropConstCmd;
  delete fTankMatPropVectorCmd;
  delete fTankMatPropConstCmd;
  delete fTankMaterialCmd;
  delete fWorldMatPropVectorCmd;
  delete fWorldMatPropConstCmd;
  delete fWorldMaterialCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  //    FINISH
  if(command == fSurfaceFinishCmd)
  {
    if(newValue == "polished")
    {
      fDetector->SetSurfaceFinish(polished);
    }
    else if(newValue == "polishedfrontpainted")
    {
      fDetector->SetSurfaceFinish(polishedfrontpainted);
    }
    else if(newValue == "polishedbackpainted")
    {
      fDetector->SetSurfaceFinish(polishedbackpainted);
    }
    else if(newValue == "ground")
    {
      fDetector->SetSurfaceFinish(ground);
    }
    else if(newValue == "groundfrontpainted")
    {
      fDetector->SetSurfaceFinish(groundfrontpainted);
    }
    else if(newValue == "groundbackpainted")
    {
      fDetector->SetSurfaceFinish(groundbackpainted);
    }
    else if(newValue == "polishedlumirrorair")
    {
      fDetector->SetSurfaceFinish(polishedlumirrorair);
    }
    else if(newValue == "polishedlumirrorglue")
    {
      fDetector->SetSurfaceFinish(polishedlumirrorglue);
    }
    else if(newValue == "polishedair")
    {
      fDetector->SetSurfaceFinish(polishedair);
    }
    else if(newValue == "polishedteflonair")
    {
      fDetector->SetSurfaceFinish(polishedteflonair);
    }
    else if(newValue == "polishedtioair")
    {
      fDetector->SetSurfaceFinish(polishedtioair);
    }
    else if(newValue == "polishedtyvekair")
    {
      fDetector->SetSurfaceFinish(polishedtyvekair);
    }
    else if(newValue == "polishedvm2000air")
    {
      fDetector->SetSurfaceFinish(polishedvm2000air);
    }
    else if(newValue == "polishedvm2000glue")
    {
      fDetector->SetSurfaceFinish(polishedvm2000glue);
    }
    else if(newValue == "etchedlumirrorair")
    {
      fDetector->SetSurfaceFinish(etchedlumirrorair);
    }
    else if(newValue == "etchedlumirrorglue")
    {
      fDetector->SetSurfaceFinish(etchedlumirrorglue);
    }
    else if(newValue == "etchedair")
    {
      fDetector->SetSurfaceFinish(etchedair);
    }
    else if(newValue == "etchedteflonair")
    {
      fDetector->SetSurfaceFinish(etchedteflonair);
    }
    else if(newValue == "etchedtioair")
    {
      fDetector->SetSurfaceFinish(etchedtioair);
    }
    else if(newValue == "etchedtyvekair")
    {
      fDetector->SetSurfaceFinish(etchedtyvekair);
    }
    else if(newValue == "etchedvm2000air")
    {
      fDetector->SetSurfaceFinish(etchedvm2000air);
    }
    else if(newValue == "etchedvm2000glue")
    {
      fDetector->SetSurfaceFinish(etchedvm2000glue);
    }
    else if(newValue == "groundlumirrorair")
    {
      fDetector->SetSurfaceFinish(groundlumirrorair);
    }
    else if(newValue == "groundlumirrorglue")
    {
      fDetector->SetSurfaceFinish(groundlumirrorglue);
    }
    else if(newValue == "groundair")
    {
      fDetector->SetSurfaceFinish(groundair);
    }
    else if(newValue == "groundteflonair")
    {
      fDetector->SetSurfaceFinish(groundteflonair);
    }
    else if(newValue == "groundtioair")
    {
      fDetector->SetSurfaceFinish(groundtioair);
    }
    else if(newValue == "groundtyvekair")
    {
      fDetector->SetSurfaceFinish(groundtyvekair);
    }
    else if(newValue == "groundvm2000air")
    {
      fDetector->SetSurfaceFinish(groundvm2000air);
    }
    else if(newValue == "groundvm2000glue")
    {
      fDetector->SetSurfaceFinish(groundvm2000glue);
    }
    //         for Davis model
    else if(newValue == "Rough_LUT")
    {
      fDetector->SetSurfaceFinish(Rough_LUT);
    }
    else if(newValue == "RoughTeflon_LUT")
    {
      fDetector->SetSurfaceFinish(RoughTeflon_LUT);
    }
    else if(newValue == "RoughESR_LUT")
    {
      fDetector->SetSurfaceFinish(RoughESR_LUT);
    }
    else if(newValue == "RoughESRGrease_LUT")
    {
      fDetector->SetSurfaceFinish(RoughESRGrease_LUT);
    }
    else if(newValue == "Polished_LUT")
    {
      fDetector->SetSurfaceFinish(Polished_LUT);
    }
    else if(newValue == "PolishedTeflon_LUT")
    {
      fDetector->SetSurfaceFinish(PolishedTeflon_LUT);
    }
    else if(newValue == "PolishedESR_LUT")
    {
      fDetector->SetSurfaceFinish(PolishedESR_LUT);
    }
    else if(newValue == "PolishedESRGrease_LUT")
    {
      fDetector->SetSurfaceFinish(PolishedESRGrease_LUT);
    }
    else if(newValue == "Detector_LUT")
    {
      fDetector->SetSurfaceFinish(Detector_LUT);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Invalid surface finish: " << newValue;
      G4Exception("DetectorMessenger", "OpNovice2_003", FatalException, ed);
    }
  }

  //  MODEL
  else if(command == fSurfaceModelCmd)
  {
    if(newValue == "glisur")
    {
      fDetector->SetSurfaceModel(glisur);
    }
    else if(newValue == "unified")
    {
      fDetector->SetSurfaceModel(unified);
    }
    else if(newValue == "LUT")
    {
      fDetector->SetSurfaceModel(LUT);
    }
    else if(newValue == "DAVIS")
    {
      fDetector->SetSurfaceModel(DAVIS);
    }
    else if(newValue == "dichroic")
    {
      fDetector->SetSurfaceModel(dichroic);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Invalid surface model: " << newValue;
      G4Exception("DetectorMessenger", "ONovice2_001", FatalException, ed);
    }
  }

  // TYPE
  else if(command == fSurfaceTypeCmd)
  {
    if(newValue == "dielectric_metal")
    {
      fDetector->SetSurfaceType(dielectric_metal);
    }
    else if(newValue == "dielectric_dielectric")
    {
      fDetector->SetSurfaceType(dielectric_dielectric);
    }
    else if(newValue == "dielectric_LUT")
    {
      fDetector->SetSurfaceType(dielectric_LUT);
    }
    else if(newValue == "dielectric_LUTDAVIS")
    {
      fDetector->SetSurfaceType(dielectric_LUTDAVIS);
    }
    else if(newValue == "coated")
    {
      fDetector->SetSurfaceType(coated);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Invalid surface type: " << newValue;
      G4Exception("DetectorMessenger", "OpNovice2_002", FatalException, ed);
    }
  }
  else if(command == fSurfaceSigmaAlphaCmd)
  {
    fDetector->SetSurfaceSigmaAlpha(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fSurfacePolishCmd)
  {
    fDetector->SetSurfacePolish(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fTankMatPropVectorCmd)
  {
    // got a string. need to convert it to physics vector.
    // string format is property name, then pairs of energy, value
    // specify units for each value, eg 3.0*eV
    // space delimited
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    std::istringstream instring(newValue);
    G4String prop;
    instring >> prop;
    while(instring)
    {
      G4String tmp;
      instring >> tmp;
      if(tmp == "")
      {
        break;
      }
      G4double en = G4UIcommand::ConvertToDouble(tmp);
      instring >> tmp;
      G4double val = G4UIcommand::ConvertToDouble(tmp);
      mpv->InsertValues(en, val);
    }

    fDetector->AddTankMPV(prop, mpv);
  }
  else if(command == fWorldMatPropVectorCmd)
  {
    // Convert string to physics vector
    // string format is property name, then pairs of energy, value
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    std::istringstream instring(newValue);
    G4String prop;
    instring >> prop;
    while(instring)
    {
      G4String tmp;
      instring >> tmp;
      if(tmp == "")
      {
        break;
      }
      G4double en = G4UIcommand::ConvertToDouble(tmp);
      instring >> tmp;
      G4double val = G4UIcommand::ConvertToDouble(tmp);
      mpv->InsertValues(en, val);
    }
    fDetector->AddWorldMPV(prop, mpv);
  }
  else if(command == fSurfaceMatPropVectorCmd)
  {
    // Convert string to physics vector
    // string format is property name, then pairs of energy, value
    // space delimited
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    G4cout << newValue << G4endl;
    std::istringstream instring(newValue);
    G4String prop;
    instring >> prop;
    while(instring)
    {
      G4String tmp;
      instring >> tmp;
      if(tmp == "")
      {
        break;
      }
      G4double en = G4UIcommand::ConvertToDouble(tmp);
      instring >> tmp;
      G4double val = G4UIcommand::ConvertToDouble(tmp);
      mpv->InsertValues(en, val);
    }
    fDetector->AddSurfaceMPV(prop, mpv);
  }

  else if(command == fTankMatPropConstCmd)
  {
    // Convert string to physics vector
    // string format is property name, then value
    // space delimited
    std::istringstream instring(newValue);
    G4String prop;
    G4String tmp;
    instring >> prop;
    instring >> tmp;
    G4double val = G4UIcommand::ConvertToDouble(tmp);
    fDetector->AddTankMPC(prop, val);
  }
  else if(command == fWorldMatPropConstCmd)
  {
    // Convert string to physics vector
    // string format is property name, then value
    // space delimited
    std::istringstream instring(newValue);
    G4String prop;
    G4String tmp;
    instring >> prop;
    instring >> tmp;
    G4double val = G4UIcommand::ConvertToDouble(tmp);
    fDetector->AddWorldMPC(prop, val);
  }
  else if(command == fSurfaceMatPropConstCmd)
  {
    // Convert string to physics vector
    // string format is property name, then value
    // space delimited
    std::istringstream instring(newValue);
    G4String prop;
    G4String tmp;
    instring >> prop;
    instring >> tmp;
    G4double val = G4UIcommand::ConvertToDouble(tmp);
    fDetector->AddSurfaceMPC(prop, val);
  }
  else if(command == fWorldMaterialCmd)
  {
    fDetector->SetWorldMaterial(newValue);
  }
  else if(command == fTankMaterialCmd)
  {
    fDetector->SetTankMaterial(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
