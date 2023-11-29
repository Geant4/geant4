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
//----------------------------------------------------------------------------
//
// ClassName:   G4OpticalParametersMessenger
//
// Author:      P.Gumplinger 30.09.2009 //
//
// Modified:    P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
//----------------------------------------------------------------------------
//

#include "G4OpticalParametersMessenger.hh"
#include "G4OpticalParameters.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UImanager.hh"
#include "G4UIparameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalParametersMessenger::G4OpticalParametersMessenger(
  G4OpticalParameters* opticalParameters)
  : params(opticalParameters)

{
  G4bool toBeBroadcasted = false;
  fDir = new G4UIdirectory("/process/optical/", toBeBroadcasted);
  fDir->SetGuidance(
    "Commands related to the optical physics simulation engine.");

  fCerenkovDir =
    new G4UIdirectory("/process/optical/cerenkov/", toBeBroadcasted);
  fCerenkovDir->SetGuidance("Cerenkov process commands");
  fScintDir =
    new G4UIdirectory("/process/optical/scintillation/", toBeBroadcasted);
  fScintDir->SetGuidance("Scintillation process commands");
  fWlsDir = new G4UIdirectory("/process/optical/wls/", toBeBroadcasted);
  fWlsDir->SetGuidance("Wave length shifting process commands");
  fWls2Dir = new G4UIdirectory("/process/optical/wls2/", toBeBroadcasted);
  fWls2Dir->SetGuidance("Second Wave length shifting process commands");
  fBoundaryDir =
    new G4UIdirectory("/process/optical/boundary/", toBeBroadcasted);
  fBoundaryDir->SetGuidance("Boundary scattering commands");
  fMieDir = new G4UIdirectory("/process/optical/mie/", toBeBroadcasted);
  fMieDir->SetGuidance("Mie scattering process commands");
  fAbsDir = new G4UIdirectory("/process/optical/absorption/", toBeBroadcasted);
  fAbsDir->SetGuidance("absorption process commands");
  fRaylDir = new G4UIdirectory("/process/optical/rayleigh/", toBeBroadcasted);
  fRaylDir->SetGuidance("Rayleigh scattering commands");

  // general commands
  fActivateProcessCmd =
    new G4UIcommand("/process/optical/processActivation", this);
  fActivateProcessCmd->SetGuidance(
    "Activate/deactivate the specified optical process");
  auto par = new G4UIparameter("proc_name", 's', false);
  G4String candidates;
  for(G4int i = 0; i < kNoProcess; ++i)
  {
    candidates += G4OpticalProcessName(i);
    candidates += G4String(" ");
  }
  par->SetParameterCandidates(candidates);
  par->SetGuidance("the process name");
  fActivateProcessCmd->SetParameter(par);
  par = new G4UIparameter("flag", 'b', true);
  par->SetDefaultValue(true);
  par->SetGuidance("activation flag");
  fActivateProcessCmd->SetParameter(par);
  fActivateProcessCmd->AvailableForStates(G4State_PreInit);

  fVerboseCmd = new G4UIcmdWithAnInteger("/process/optical/verbose", this);
  fVerboseCmd->SetGuidance("Set default verbose level for optical processes");
  fVerboseCmd->SetParameterName("ver", true);
  fVerboseCmd->SetDefaultValue(1);
  fVerboseCmd->SetRange("ver>=0");
  fVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDumpCmd = new G4UIcommand("/process/optical/printParameters", this);
  fDumpCmd->SetGuidance("Print all optical parameters.");

  // Cerenkov ////////////////////
  fCerenkovMaxPhotonsCmd =
    new G4UIcmdWithAnInteger("/process/optical/cerenkov/setMaxPhotons", this);
  fCerenkovMaxPhotonsCmd->SetGuidance("Set maximum number of photons per step");
  fCerenkovMaxPhotonsCmd->SetParameterName("CerenkovMaxPhotons", false);
  fCerenkovMaxPhotonsCmd->SetRange("CerenkovMaxPhotons>=0");
  fCerenkovMaxPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCerenkovMaxBetaChangeCmd =
    new G4UIcmdWithADouble("/process/optical/cerenkov/setMaxBetaChange", this);
  fCerenkovMaxBetaChangeCmd->SetGuidance(
    "Set maximum change of beta of parent particle per step (in percent)");
  fCerenkovMaxBetaChangeCmd->SetParameterName("CerenkovMaxBetaChange", false);
  fCerenkovMaxBetaChangeCmd->SetRange("CerenkovMaxBetaChange>=0");
  fCerenkovMaxBetaChangeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCerenkovStackPhotonsCmd =
    new G4UIcmdWithABool("/process/optical/cerenkov/setStackPhotons", this);
  fCerenkovStackPhotonsCmd->SetGuidance(
    "Set whether or not to stack secondary Cerenkov photons");
  fCerenkovStackPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCerenkovTrackSecondariesFirstCmd = new G4UIcmdWithABool(
    "/process/optical/cerenkov/setTrackSecondariesFirst", this);
  fCerenkovTrackSecondariesFirstCmd->SetGuidance(
    "Whether to track secondary Cerenkov photons before the primary.");
  fCerenkovTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit,
                                                        G4State_Idle);

  fCerenkovVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/cerenkov/verbose", this);
  fCerenkovVerboseLevelCmd->SetGuidance("Verbose level for Cerenkov process.");
  fCerenkovVerboseLevelCmd->SetParameterName("verbose", true);
  fCerenkovVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fCerenkovVerboseLevelCmd->SetDefaultValue(2);
  fCerenkovVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Scintillation //////////////////////////
  fScintByParticleTypeCmd = new G4UIcmdWithABool(
    "/process/optical/scintillation/setByParticleType", this);
  fScintByParticleTypeCmd->SetGuidance(
    "Activate/Inactivate scintillation process by particle type");
  fScintByParticleTypeCmd->SetParameterName(
    "ScintillationByParticleTypeActivation", false);
  fScintByParticleTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fScintTrackInfoCmd =
    new G4UIcmdWithABool("/process/optical/scintillation/setTrackInfo", this);
  fScintTrackInfoCmd->SetGuidance(
    "Activate/Inactivate scintillation TrackInformation");
  fScintTrackInfoCmd->SetParameterName("ScintillationTrackInfo", false);
  fScintTrackInfoCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fScintFiniteRiseTimeCmd = new G4UIcmdWithABool(
    "/process/optical/scintillation/setFiniteRiseTime", this);
  fScintFiniteRiseTimeCmd->SetGuidance(
    "Set option of a finite rise-time for G4Scintillation");
  fScintFiniteRiseTimeCmd->SetGuidance(
    "If set, the G4Scintillation process expects the user to have set the");
  fScintFiniteRiseTimeCmd->SetGuidance(
    "constant material property SCINTILLATIONRISETIME{1,2,3}");
  fScintFiniteRiseTimeCmd->SetParameterName("FiniteRiseTime", false);
  fScintFiniteRiseTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fScintStackPhotonsCmd = new G4UIcmdWithABool(
    "/process/optical/scintillation/setStackPhotons", this);
  fScintStackPhotonsCmd->SetGuidance(
    "Set whether or not to stack secondary Scintillation photons");
  fScintStackPhotonsCmd->SetParameterName("ScintillationStackPhotons", true);
  fScintStackPhotonsCmd->SetDefaultValue(true);
  fScintStackPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fScintTrackSecondariesFirstCmd = new G4UIcmdWithABool(
    "/process/optical/scintillation/setTrackSecondariesFirst", this);
  fScintTrackSecondariesFirstCmd->SetGuidance(
    "Whether to track scintillation secondaries before primary.");
  fScintTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit,
                                                     G4State_Idle);

  fScintVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/scintillation/verbose", this);
  fScintVerboseLevelCmd->SetGuidance(
    "Verbose level for scintillation process.");
  fScintVerboseLevelCmd->SetParameterName("verbose", true);
  fScintVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fScintVerboseLevelCmd->AvailableForStates(G4State_Idle, G4State_PreInit);

  // WLS   //////////////////////////////////
  fWLSTimeProfileCmd =
    new G4UIcmdWithAString("/process/optical/wls/setTimeProfile", this);
  fWLSTimeProfileCmd->SetGuidance(
    "Set the WLS time profile (delta or exponential)");
  fWLSTimeProfileCmd->SetParameterName("WLSTimeProfile", false);
  fWLSTimeProfileCmd->SetCandidates("delta exponential");
  fWLSTimeProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fWLSVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/wls/verbose", this);
  fWLSVerboseLevelCmd->SetGuidance("Verbose level for WLS process.");
  fWLSVerboseLevelCmd->SetParameterName("verbose", true);
  fWLSVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fWLSVerboseLevelCmd->SetDefaultValue(1);
  fWLSVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // WLS2   //////////////////////////////////
  fWLS2TimeProfileCmd =
    new G4UIcmdWithAString("/process/optical/wls2/setTimeProfile", this);
  fWLS2TimeProfileCmd->SetGuidance(
    "Set the WLS2 time profile (delta or exponential)");
  fWLS2TimeProfileCmd->SetParameterName("WLS2TimeProfile", false);
  fWLS2TimeProfileCmd->SetCandidates("delta exponential");
  fWLS2TimeProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fWLS2VerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/wls2/verbose", this);
  fWLS2VerboseLevelCmd->SetGuidance("Verbose level for WLS2 process.");
  fWLS2VerboseLevelCmd->SetParameterName("verbose", true);
  fWLS2VerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fWLS2VerboseLevelCmd->SetDefaultValue(1);
  fWLS2VerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // boundary //////////////////////////////////////
  fBoundaryInvokeSDCmd =
    new G4UIcmdWithABool("/process/optical/boundary/setInvokeSD", this);
  fBoundaryInvokeSDCmd->SetGuidance(
    "Set option for calling InvokeSD in G4OpBoundaryProcess");
  fBoundaryInvokeSDCmd->SetParameterName("InvokeSD", false);
  fBoundaryInvokeSDCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBoundaryVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/boundary/verbose", this);
  fBoundaryVerboseLevelCmd->SetGuidance("Verbose level for boundary process.");
  fBoundaryVerboseLevelCmd->SetParameterName("verbose", true);
  fBoundaryVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fBoundaryVerboseLevelCmd->SetDefaultValue(1);
  fBoundaryVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // absorption //////////////////////////////////////
  fAbsorptionVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/absorption/verbose", this);
  fAbsorptionVerboseLevelCmd->SetGuidance(
    "Verbose level for absorption process.");
  fAbsorptionVerboseLevelCmd->SetParameterName("verbose", true);
  fAbsorptionVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fAbsorptionVerboseLevelCmd->SetDefaultValue(1);
  fAbsorptionVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // rayleigh //////////////////////////////////////
  fRayleighVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/rayleigh/verbose", this);
  fRayleighVerboseLevelCmd->SetGuidance("Verbose level for Rayleigh process.");
  fRayleighVerboseLevelCmd->SetParameterName("verbose", true);
  fRayleighVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fRayleighVerboseLevelCmd->SetDefaultValue(1);
  fRayleighVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // mie //////////////////////////////////////
  fMieVerboseLevelCmd =
    new G4UIcmdWithAnInteger("/process/optical/mie/verbose", this);
  fMieVerboseLevelCmd->SetGuidance("Verbose level for Mie process.");
  fMieVerboseLevelCmd->SetParameterName("verbose", true);
  fMieVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
  fMieVerboseLevelCmd->SetDefaultValue(1);
  fMieVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

G4OpticalParametersMessenger::~G4OpticalParametersMessenger()
{
  delete fDir;
  delete fCerenkovDir;
  delete fScintDir;
  delete fWlsDir;
  delete fBoundaryDir;
  delete fMieDir;
  delete fAbsDir;
  delete fRaylDir;
  delete fActivateProcessCmd;
  delete fVerboseCmd;
  delete fDumpCmd;
  delete fCerenkovMaxPhotonsCmd;
  delete fCerenkovMaxBetaChangeCmd;
  delete fCerenkovStackPhotonsCmd;
  delete fCerenkovTrackSecondariesFirstCmd;
  delete fCerenkovVerboseLevelCmd;
  delete fScintByParticleTypeCmd;
  delete fScintTrackInfoCmd;
  delete fScintStackPhotonsCmd;
  delete fScintVerboseLevelCmd;
  delete fScintFiniteRiseTimeCmd;
  delete fScintTrackSecondariesFirstCmd;
  delete fWLSTimeProfileCmd;
  delete fWLSVerboseLevelCmd;
  delete fWLS2TimeProfileCmd;
  delete fWLS2VerboseLevelCmd;
  delete fAbsorptionVerboseLevelCmd;
  delete fRayleighVerboseLevelCmd;
  delete fMieVerboseLevelCmd;
  delete fBoundaryVerboseLevelCmd;
  delete fBoundaryInvokeSDCmd;
}

void G4OpticalParametersMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{
  // physics needs to be rebuilt for all commands
  G4bool physicsModified = true;

  /// Apply command to the associated object.
  if(command == fActivateProcessCmd)
  {
    std::istringstream is(newValue.data());
    G4String pn;
    G4String flag;
    is >> pn >> flag;
    G4bool value = G4UIcommand::ConvertToBool(flag);
    params->SetProcessActivation(pn, value);
  }
  else if(command == fVerboseCmd)
  {
    params->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if(command == fDumpCmd)
  {
    params->Dump();
  }
  else if(command == fCerenkovMaxPhotonsCmd)
  {
    params->SetCerenkovMaxPhotonsPerStep(
      fCerenkovMaxPhotonsCmd->GetNewIntValue(newValue));
    G4cout << "Cerenkov max photons: " << params->GetCerenkovMaxPhotonsPerStep()
           << G4endl;
  }
  else if(command == fCerenkovMaxBetaChangeCmd)
  {
    params->SetCerenkovMaxBetaChange(
      fCerenkovMaxBetaChangeCmd->GetNewDoubleValue(newValue));
  }
  else if(command == fCerenkovStackPhotonsCmd)
  {
    params->SetCerenkovStackPhotons(
      fCerenkovStackPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if(command == fCerenkovTrackSecondariesFirstCmd)
  {
    params->SetCerenkovTrackSecondariesFirst(
      fCerenkovTrackSecondariesFirstCmd->GetNewBoolValue(newValue));
  }
  else if(command == fCerenkovVerboseLevelCmd)
  {
    params->SetCerenkovVerboseLevel(
      fCerenkovVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fScintByParticleTypeCmd)
  {
    params->SetScintByParticleType(
      fScintByParticleTypeCmd->GetNewBoolValue(newValue));
  }
  else if(command == fScintTrackInfoCmd)
  {
    params->SetScintTrackInfo(fScintTrackInfoCmd->GetNewBoolValue(newValue));
  }
  else if(command == fScintFiniteRiseTimeCmd)
  {
    params->SetScintFiniteRiseTime(
      fScintFiniteRiseTimeCmd->GetNewBoolValue(newValue));
  }
  else if(command == fScintStackPhotonsCmd)
  {
    params->SetScintStackPhotons(
      fScintStackPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if(command == fScintTrackSecondariesFirstCmd)
  {
    params->SetScintTrackSecondariesFirst(
      fScintTrackSecondariesFirstCmd->GetNewBoolValue(newValue));
  }
  else if(command == fScintVerboseLevelCmd)
  {
    params->SetScintVerboseLevel(
      fScintVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fWLSTimeProfileCmd)
  {
    params->SetWLSTimeProfile(newValue);
  }
  else if(command == fWLSVerboseLevelCmd)
  {
    params->SetWLSVerboseLevel(fWLSVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fWLS2TimeProfileCmd)
  {
    params->SetWLS2TimeProfile(newValue);
  }
  else if(command == fWLS2VerboseLevelCmd)
  {
    params->SetWLS2VerboseLevel(fWLS2VerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fAbsorptionVerboseLevelCmd)
  {
    params->SetAbsorptionVerboseLevel(
      fAbsorptionVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fRayleighVerboseLevelCmd)
  {
    params->SetRayleighVerboseLevel(
      fRayleighVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fMieVerboseLevelCmd)
  {
    params->SetMieVerboseLevel(fMieVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fBoundaryVerboseLevelCmd)
  {
    params->SetBoundaryVerboseLevel(
      fBoundaryVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if(command == fBoundaryInvokeSDCmd)
  {
    params->SetBoundaryInvokeSD(
      fBoundaryInvokeSDCmd->GetNewBoolValue(newValue));
  }
  if(physicsModified)
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }
}
