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


// Commands with '/defaults/' are duplicates and will be removed in
// the next major release of Geant4. Use commands with no /defaults/ instead


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalParametersMessenger::G4OpticalParametersMessenger(
                  G4OpticalParameters* opticalParameters)
  : params(opticalParameters)

{
    G4bool toBeBroadcasted = false;
    fDir = new G4UIdirectory("/process/optical/defaults/",toBeBroadcasted);
    fDir->SetGuidance("DEPRECATED Commands related to the optical physics simulation engine.");
    fDir2 = new G4UIdirectory("/process/optical/",toBeBroadcasted);
    fDir2->SetGuidance("Commands related to the optical physics simulation engine.");

    CreateDirectory("/process/optical/defaults/cerenkov/", "DEPRECATED Cerenkov process commands");
    CreateDirectory("/process/optical/defaults/scintillation/", "DEPRECATED Scintillation process commands");
    CreateDirectory("/process/optical/defaults/wls/", "DEPRECATED Wave length shifting process commands");
    CreateDirectory("/process/optical/defaults/boundary/", "DEPRECATED Boundary scattering commands");

    CreateDirectory("/process/optical/cerenkov/", "Cerenkov process commands");
    CreateDirectory("/process/optical/scintillation/", "Scintillation process commands");
    CreateDirectory("/process/optical/wls/", "Wave length shifting process commands");
    CreateDirectory("/process/optical/wls2/", "Second Wave length shifting process commands");
    CreateDirectory("/process/optical/boundary/", "Boundary scattering commands");
    CreateDirectory("/process/optical/mie/", "Mie scattering process commands");
    CreateDirectory("/process/optical/absorption/", "absorption process commands");
    CreateDirectory("/process/optical/rayleigh/", "Rayleigh scattering commands");

    // general commands
    fActivateProcessCmd= new G4UIcommand("/process/optical/processActivation", this);
    fActivateProcessCmd->SetGuidance("Activate/deactivate the specified optical process");
    G4UIparameter* par = new G4UIparameter("proc_name",'s',false);
    G4String candidates;
    for ( G4int i=0; i<kNoProcess; i++ ) {
        candidates += G4OpticalProcessName(i);
        candidates += G4String(" ");
    }
    par->SetParameterCandidates(candidates);
    par->SetGuidance("the process name");
    fActivateProcessCmd->SetParameter(par);
    par = new G4UIparameter("flag",'b',true);
    par->SetDefaultValue(true);
    par->SetGuidance("activation flag");
    fActivateProcessCmd->SetParameter(par);
    fActivateProcessCmd->AvailableForStates(G4State_PreInit);

    // DEPRECATED
    fTrackSecondariesFirstCmd = new G4UIcommand("/process/optical/setTrackSecondariesFirst", this);
    fTrackSecondariesFirstCmd->SetGuidance("Activate/deactivate tracking of secondaries before finishing their parent track");
    fTrackSecondariesFirstCmd->SetGuidance("DEPRECATED: Use /process/optical/cerenkov/setTrackSecondariesFirst and");
    fTrackSecondariesFirstCmd->SetGuidance("/process/optical/scintillation/setTrackSecondariesFirst and instead.");
    par = new G4UIparameter("proc_name",'s',false);
    par->SetParameterCandidates("Cerenkov Scintillation");
    fTrackSecondariesFirstCmd->SetParameter(par);
    par = new G4UIparameter("flag",'b',false);
    par->SetDefaultValue(true);
    fTrackSecondariesFirstCmd->SetParameter(par);
    fTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fVerboseCmd = new G4UIcmdWithAnInteger("/process/optical/verbose", this);
    fVerboseCmd->SetGuidance("Set default verbose level for optical processes");
    fVerboseCmd->SetParameterName("ver", true);
    fVerboseCmd->SetDefaultValue(1);
    fVerboseCmd->SetRange("ver>=0");
    fVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fDumpCmd = new G4UIcommand("/process/optical/printParameters", this);
    fDumpCmd->SetGuidance("Print all optical parameters.");

    // Cerenkov ////////////////////
    fCerenkovMaxPhotons1Cmd = new G4UIcmdWithAnInteger("/process/optical/defaults/cerenkov/setMaxPhotons", this);
    fCerenkovMaxPhotons1Cmd->SetGuidance("Set maximum number of photons per step");
    fCerenkovMaxPhotons1Cmd->SetGuidance("DEPRECATED: use /process/optical/cerenkov/setMaxPhotons instead.");
    fCerenkovMaxPhotons1Cmd->SetParameterName("CerenkovMaxPhotons", false);
    fCerenkovMaxPhotons1Cmd->SetRange("CerenkovMaxPhotons>=0");
    fCerenkovMaxPhotons1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovMaxPhotonsCmd = new G4UIcmdWithAnInteger("/process/optical/cerenkov/setMaxPhotons", this);
    fCerenkovMaxPhotonsCmd->SetGuidance("Set maximum number of photons per step");
    fCerenkovMaxPhotonsCmd->SetParameterName("CerenkovMaxPhotons", false);
    fCerenkovMaxPhotonsCmd->SetRange("CerenkovMaxPhotons>=0");
    fCerenkovMaxPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovMaxBetaChange1Cmd = new G4UIcmdWithADouble("/process/optical/defaults/cerenkov/setMaxBetaChange", this);
    fCerenkovMaxBetaChange1Cmd->SetGuidance("Set maximum change of beta of parent particle per step");
    fCerenkovMaxBetaChange1Cmd->SetGuidance("DEPRECATED: use /process/optical/cerenkov/setMaxBetaChange instead.");
    fCerenkovMaxBetaChange1Cmd->SetParameterName("CerenkovMaxBetaChange", false);
    fCerenkovMaxBetaChange1Cmd->SetRange("CerenkovMaxBetaChange>=0");
    fCerenkovMaxBetaChange1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovMaxBetaChangeCmd = new G4UIcmdWithADouble("/process/optical/cerenkov/setMaxBetaChange", this);
    fCerenkovMaxBetaChangeCmd->SetGuidance("Set maximum change of beta of parent particle per step");
    fCerenkovMaxBetaChangeCmd->SetParameterName("CerenkovMaxBetaChange", false);
    fCerenkovMaxBetaChangeCmd->SetRange("CerenkovMaxBetaChange>=0");
    fCerenkovMaxBetaChangeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovStackPhotons1Cmd = new G4UIcmdWithABool("/process/optical/defaults/cerenkov/setStackPhotons", this);
    fCerenkovStackPhotons1Cmd->SetGuidance("Set whether or not to stack secondary Cerenkov photons");
    fCerenkovStackPhotons1Cmd->SetGuidance("DEPRECATED: use /process/optical/cerenkov/setStackPhotons instead.");
    fCerenkovStackPhotons1Cmd->SetParameterName("CerenkovStackPhotons", true);
    fCerenkovStackPhotons1Cmd->SetDefaultValue(true);
    fCerenkovStackPhotons1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovStackPhotonsCmd = new G4UIcmdWithABool("/process/optical/cerenkov/setStackPhotons", this);
    fCerenkovStackPhotonsCmd->SetGuidance("Set whether or not to stack secondary Cerenkov photons");
    fCerenkovStackPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovTrackSecondariesFirstCmd = new G4UIcmdWithABool("/process/optical/cerenkov/setTrackSecondariesFirst", this);
    fCerenkovTrackSecondariesFirstCmd->SetGuidance("Whether to track secondary Cerenkov photons before the primary.");
    fCerenkovTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/cerenkov/verbose", this);
    fCerenkovVerboseLevelCmd->SetGuidance("Verbose level for Cerenkov process.");
    fCerenkovVerboseLevelCmd->SetParameterName("verbose", true);
    fCerenkovVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fCerenkovVerboseLevelCmd->SetDefaultValue(2);
    fCerenkovVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    // Scintillation //////////////////////////
    fScintYieldFactor1Cmd = new G4UIcmdWithADouble("/process/optical/defaults/scintillation/setYieldFactor", this);
    fScintYieldFactor1Cmd->SetGuidance("Set scintillation yield factor");
    fScintYieldFactor1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setYieldFactorinstead.");
    fScintYieldFactor1Cmd->SetParameterName("ScintillationYieldFactor", false);
    fScintYieldFactor1Cmd->SetRange("ScintillationYieldFactor>=0");
    fScintYieldFactor1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintYieldFactorCmd = new G4UIcmdWithADouble("/process/optical/scintillation/setYieldFactor", this);
    fScintYieldFactorCmd->SetGuidance("Set scintillation yield factor");
    fScintYieldFactorCmd->SetParameterName("ScintillationYieldFactor", false);
    fScintYieldFactorCmd->SetRange("ScintillationYieldFactor>=0");
    fScintYieldFactorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintExcitationRatioCmd = new G4UIcmdWithADouble("/process/optical/scintillation/setExcitationRatio", this);
    fScintExcitationRatioCmd->SetGuidance("Set scintillation excitation ratio");
    fScintExcitationRatioCmd->SetParameterName("ExcitationRatio", false);
    fScintExcitationRatioCmd->SetRange("ExcitationRatio >= 0 && ExcitationRatio <=1");
    fScintExcitationRatioCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintByParticleType1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setByParticleType", this);
    fScintByParticleType1Cmd->SetGuidance("Activate/Inactivate scintillation process by particle type");
    fScintByParticleType1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setByParticleType instead.");
    fScintByParticleType1Cmd->SetParameterName("ScintillationByParticleTypeActivation", false);
    fScintByParticleType1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintByParticleTypeCmd = new G4UIcmdWithABool("/process/optical/scintillation/setByParticleType", this);
    fScintByParticleTypeCmd->SetGuidance("Activate/Inactivate scintillation process by particle type");
    fScintByParticleTypeCmd->SetParameterName("ScintillationByParticleTypeActivation", false);
    fScintByParticleTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintEnhancedTimeConstantsCmd = new G4UIcmdWithABool("/process/optical/scintillation/setEnhancedTimeConstants", this);
    fScintEnhancedTimeConstantsCmd->SetGuidance("Activate/Inactivate enhanced time constants for scintillation.");
    fScintEnhancedTimeConstantsCmd->SetGuidance("This will be the default in the next major release.");
    fScintEnhancedTimeConstantsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintTrackInfo1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setTrackInfo", this);
    fScintTrackInfo1Cmd->SetGuidance("Activate/Inactivate scintillation TrackInformation");
    fScintTrackInfo1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setTrackInfo instead.");
    fScintTrackInfo1Cmd->SetParameterName("ScintillationTrackInfo", false);
    fScintTrackInfo1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintTrackInfoCmd = new G4UIcmdWithABool("/process/optical/scintillation/setTrackInfo", this);
    fScintTrackInfoCmd->SetGuidance("Activate/Inactivate scintillation TrackInformation");
    fScintTrackInfoCmd->SetParameterName("ScintillationTrackInfo", false);
    fScintTrackInfoCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintFiniteRiseTime1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setFiniteRiseTime", this);
    fScintFiniteRiseTime1Cmd->SetGuidance("Set option of a finite rise-time for G4Scintillation");
    fScintFiniteRiseTime1Cmd->SetGuidance("If set, the G4Scintillation process expects the user to have set the");
    fScintFiniteRiseTime1Cmd->SetGuidance("constant material property FAST/SLOWSCINTILLATIONRISETIME");
    fScintFiniteRiseTime1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setFiniteRiseTime instead.");
    fScintFiniteRiseTime1Cmd->SetParameterName("FiniteRiseTime", false);
    fScintFiniteRiseTime1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintFiniteRiseTimeCmd = new G4UIcmdWithABool("/process/optical/scintillation/setFiniteRiseTime", this);
    fScintFiniteRiseTimeCmd->SetGuidance("Set option of a finite rise-time for G4Scintillation");
    fScintFiniteRiseTimeCmd->SetGuidance("If set, the G4Scintillation process expects the user to have set the");
    fScintFiniteRiseTimeCmd->SetGuidance("constant material property FAST/SLOWSCINTILLATIONRISETIME");
    fScintFiniteRiseTimeCmd->SetParameterName("FiniteRiseTime", false);
    fScintFiniteRiseTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintStackPhotons1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setStackPhotons", this);
    fScintStackPhotons1Cmd->SetGuidance("Set whether or not to stack secondary Scintillation photons");
    fScintStackPhotons1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setStackPhotons instead.");
    fScintStackPhotons1Cmd->SetParameterName("ScintillationStackPhotons", true);
    fScintStackPhotons1Cmd->SetDefaultValue(true);
    fScintStackPhotons1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintStackPhotonsCmd = new G4UIcmdWithABool("/process/optical/scintillation/setStackPhotons", this);
    fScintStackPhotonsCmd->SetGuidance("Set whether or not to stack secondary Scintillation photons");
    fScintStackPhotonsCmd->SetParameterName("ScintillationStackPhotons", true);
    fScintStackPhotonsCmd->SetDefaultValue(true);
    fScintStackPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintTrackSecondariesFirstCmd = new G4UIcmdWithABool("/process/optical/scintillation/setTrackSecondariesFirst", this);
    fScintTrackSecondariesFirstCmd->SetGuidance("Whether to track scintillation secondaries before primary.");
    fScintTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/scintillation/verbose", this);
    fScintVerboseLevelCmd->SetGuidance("Verbose level for scintillation process.");
    fScintVerboseLevelCmd->SetParameterName("verbose", true);
    fScintVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fScintVerboseLevelCmd->AvailableForStates(G4State_Idle, G4State_PreInit);

    // WLS   //////////////////////////////////
    fWLSTimeProfile1Cmd = new G4UIcmdWithAString("/process/optical/defaults/wls/setTimeProfile", this);
    fWLSTimeProfile1Cmd->SetGuidance("Set the WLS time profile (delta or exponential)");
    fWLSTimeProfile1Cmd->SetGuidance("DEPRECATED: use /process/optical/wls/setTimeProfile instead.");
    fWLSTimeProfile1Cmd->SetParameterName("WLSTimeProfile", false);
    fWLSTimeProfile1Cmd->SetCandidates("delta exponential");
    fWLSTimeProfile1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fWLSTimeProfileCmd = new G4UIcmdWithAString("/process/optical/wls/setTimeProfile", this);
    fWLSTimeProfileCmd->SetGuidance("Set the WLS time profile (delta or exponential)");
    fWLSTimeProfileCmd->SetParameterName("WLSTimeProfile", false);
    fWLSTimeProfileCmd->SetCandidates("delta exponential");
    fWLSTimeProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fWLSVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/wls/verbose", this);
    fWLSVerboseLevelCmd->SetGuidance("Verbose level for WLS process.");
    fWLSVerboseLevelCmd->SetParameterName("verbose", true);
    fWLSVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fWLSVerboseLevelCmd->SetDefaultValue(1);
    fWLSVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    // WLS2   //////////////////////////////////
    fWLS2TimeProfileCmd = new G4UIcmdWithAString("/process/optical/wls2/setTimeProfile", this);
    fWLS2TimeProfileCmd->SetGuidance("Set the WLS2 time profile (delta or exponential)");
    fWLS2TimeProfileCmd->SetParameterName("WLS2TimeProfile", false);
    fWLS2TimeProfileCmd->SetCandidates("delta exponential");
    fWLS2TimeProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fWLS2VerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/wls2/verbose", this);
    fWLS2VerboseLevelCmd->SetGuidance("Verbose level for WLS2 process.");
    fWLS2VerboseLevelCmd->SetParameterName("verbose", true);
    fWLS2VerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fWLS2VerboseLevelCmd->SetDefaultValue(1);
    fWLS2VerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    // boundary //////////////////////////////////////
    fBoundaryInvokeSD1Cmd = new G4UIcmdWithABool("/process/optical/defaults/boundary/setInvokeSD", this);
    fBoundaryInvokeSD1Cmd->SetGuidance("Set option for calling InvokeSD in G4OpBoundaryProcess");
    fBoundaryInvokeSD1Cmd->SetGuidance("DEPRECATED: use /process/optical/boundary/setInvokeSD instead.");
    fBoundaryInvokeSD1Cmd->SetParameterName("InvokeSD", false);
    fBoundaryInvokeSD1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBoundaryInvokeSDCmd = new G4UIcmdWithABool("/process/optical/boundary/setInvokeSD", this);
    fBoundaryInvokeSDCmd->SetGuidance("Set option for calling InvokeSD in G4OpBoundaryProcess");
    fBoundaryInvokeSDCmd->SetParameterName("InvokeSD", false);
    fBoundaryInvokeSDCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBoundaryVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/boundary/verbose", this);
    fBoundaryVerboseLevelCmd->SetGuidance("Verbose level for boundary process.");
    fBoundaryVerboseLevelCmd->SetParameterName("verbose", true);
    fBoundaryVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fBoundaryVerboseLevelCmd->SetDefaultValue(1);
    fBoundaryVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    // absorption //////////////////////////////////////
    fBoundaryInvokeSD1Cmd = new G4UIcmdWithABool("/process/optical/defaults/boundary/setInvokeSD", this);
    fAbsorptionVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/absorption/verbose", this);
    fAbsorptionVerboseLevelCmd->SetGuidance("Verbose level for absorption process.");
    fAbsorptionVerboseLevelCmd->SetParameterName("verbose", true);
    fAbsorptionVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fAbsorptionVerboseLevelCmd->SetDefaultValue(1);
    fAbsorptionVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    // rayleigh //////////////////////////////////////
    fBoundaryInvokeSD1Cmd = new G4UIcmdWithABool("/process/optical/defaults/boundary/setInvokeSD", this);
    fRayleighVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/rayleigh/verbose", this);
    fRayleighVerboseLevelCmd->SetGuidance("Verbose level for Rayleigh process.");
    fRayleighVerboseLevelCmd->SetParameterName("verbose", true);
    fRayleighVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fRayleighVerboseLevelCmd->SetDefaultValue(1);
    fRayleighVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    // mie //////////////////////////////////////
    fMieVerboseLevelCmd = new G4UIcmdWithAnInteger("/process/optical/mie/verbose", this);
    fMieVerboseLevelCmd->SetGuidance("Verbose level for Mie process.");
    fMieVerboseLevelCmd->SetParameterName("verbose", true);
    fMieVerboseLevelCmd->SetRange("verbose >= 0 && verbose <= 2");
    fMieVerboseLevelCmd->SetDefaultValue(1);
    fMieVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

G4OpticalParametersMessenger::~G4OpticalParametersMessenger()
{
  delete fDir;
  delete fDir2;
  delete fActivateProcessCmd;
  delete fVerboseCmd;
  delete fDumpCmd;
  delete fCerenkovMaxPhotonsCmd;
  delete fCerenkovMaxPhotons1Cmd;
  delete fCerenkovMaxBetaChangeCmd;
  delete fCerenkovMaxBetaChange1Cmd;
  delete fCerenkovStackPhotonsCmd;
  delete fCerenkovStackPhotons1Cmd;
  delete fCerenkovTrackSecondariesFirstCmd;
  delete fCerenkovVerboseLevelCmd;
  delete fScintYieldFactorCmd;
  delete fScintYieldFactor1Cmd;
  delete fScintByParticleTypeCmd;
  delete fScintByParticleType1Cmd;
  delete fScintEnhancedTimeConstantsCmd;
  delete fScintTrackInfoCmd;
  delete fScintTrackInfo1Cmd;
  delete fScintStackPhotonsCmd;
  delete fScintStackPhotons1Cmd;
  delete fScintExcitationRatioCmd;
  delete fScintVerboseLevelCmd;
  delete fScintFiniteRiseTimeCmd;
  delete fScintFiniteRiseTime1Cmd;
  delete fScintTrackSecondariesFirstCmd;
  delete fWLSTimeProfileCmd;
  delete fWLSTimeProfile1Cmd;
  delete fWLSVerboseLevelCmd;
  delete fWLS2TimeProfileCmd;
  delete fWLS2VerboseLevelCmd;
  delete fAbsorptionVerboseLevelCmd;
  delete fRayleighVerboseLevelCmd;
  delete fMieVerboseLevelCmd;
  delete fBoundaryVerboseLevelCmd;
  delete fTrackSecondariesFirstCmd;
  delete fBoundaryInvokeSDCmd;
  delete fBoundaryInvokeSD1Cmd;
}

void G4OpticalParametersMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  // physics needs to be rebuilt for all commands
  G4bool physicsModified = true;

  /// Apply command to the associated object.
  if (command == fActivateProcessCmd) {
    std::istringstream is(newValue.data());
    G4String pn;
    G4String flag;
    is >> pn >> flag;
    G4bool value = G4UIcommand::ConvertToBool(flag);
    params->SetProcessActivation(pn, value);
  }
  else if (command == fTrackSecondariesFirstCmd) {
    std::istringstream is(newValue.data());
    G4String pn;
    G4String flag;
    is >> pn >> flag;
    G4bool value = G4UIcommand::ConvertToBool(flag);
    if (pn == "Cerenkov") params->SetCerenkovStackPhotons(value);
    else if (pn == "Scintillation") params->SetScintStackPhotons(value);
    else {
      G4ExceptionDescription msg;
      msg << "Process name not allowed:  "<<pn<<" (UI: "<<newValue<<")";
      G4Exception("G4OpticalParametersMessenger::SetNewValue(...)","Optical001",
                  FatalException,msg);
    }
  }
  else if (command == fVerboseCmd) {
    params->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if (command == fDumpCmd) {
    params->Dump();
  }
  else if (command == fCerenkovMaxPhotons1Cmd) {
    params->SetCerenkovMaxPhotonsPerStep(
          fCerenkovMaxPhotons1Cmd->GetNewIntValue(newValue));
    Deprecated();
  }
  else if (command == fCerenkovMaxPhotonsCmd) {
    params->SetCerenkovMaxPhotonsPerStep(
          fCerenkovMaxPhotonsCmd->GetNewIntValue(newValue));
    G4cout << "Cerenkov max photons: " << params->GetCerenkovMaxPhotonsPerStep() << G4endl;
  }
  else if (command == fCerenkovMaxBetaChange1Cmd) {
    params->SetCerenkovMaxBetaChange(
          fCerenkovMaxBetaChange1Cmd->GetNewDoubleValue(newValue));
    Deprecated();
  }
  else if (command == fCerenkovMaxBetaChangeCmd) {
    params->SetCerenkovMaxBetaChange(
          fCerenkovMaxBetaChangeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fCerenkovStackPhotons1Cmd) {
    params->SetCerenkovStackPhotons(
          fCerenkovStackPhotons1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fCerenkovStackPhotonsCmd) {
    params->SetCerenkovStackPhotons(
          fCerenkovStackPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fCerenkovTrackSecondariesFirstCmd) {
    params->SetCerenkovTrackSecondariesFirst(
          fCerenkovTrackSecondariesFirstCmd->GetNewBoolValue(newValue));
  }
  else if (command == fCerenkovVerboseLevelCmd) {
    params->SetCerenkovVerboseLevel(
          fCerenkovVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fScintYieldFactor1Cmd) {
    params->SetScintYieldFactor(
          fScintYieldFactor1Cmd->GetNewDoubleValue(newValue));
    Deprecated();
  }
  else if (command == fScintYieldFactorCmd) {
    params->SetScintYieldFactor(
          fScintYieldFactorCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fScintByParticleType1Cmd) {
    params->SetScintByParticleType(
         fScintByParticleType1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintByParticleTypeCmd) {
    params->SetScintByParticleType(
         fScintByParticleTypeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintEnhancedTimeConstantsCmd) {
    params->SetScintEnhancedTimeConstants(
         fScintEnhancedTimeConstantsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintTrackInfo1Cmd) {
    params->SetScintTrackInfo(
         fScintTrackInfo1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintTrackInfoCmd) {
    params->SetScintTrackInfo(
         fScintTrackInfoCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintFiniteRiseTime1Cmd) {
    params->SetScintFiniteRiseTime(
         fScintFiniteRiseTime1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintFiniteRiseTimeCmd) {
    params->SetScintFiniteRiseTime(
         fScintFiniteRiseTimeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintStackPhotons1Cmd) {
    params->SetScintStackPhotons(
          fScintStackPhotons1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintStackPhotonsCmd) {
    params->SetScintStackPhotons(
          fScintStackPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintExcitationRatioCmd) {
    params->SetScintExcitationRatio(
          fScintExcitationRatioCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fScintTrackSecondariesFirstCmd) {
    params->SetScintTrackSecondariesFirst(
          fScintTrackSecondariesFirstCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintVerboseLevelCmd) {
    params->SetScintVerboseLevel(
          fScintVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fWLSTimeProfile1Cmd) {
    params->SetWLSTimeProfile(newValue);
    Deprecated();
  }
  else if (command == fWLSTimeProfileCmd) {
    params->SetWLSTimeProfile(newValue);
  }
  else if (command == fWLSVerboseLevelCmd) {
    params->SetWLSVerboseLevel(fWLSVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fWLS2TimeProfileCmd) {
    params->SetWLS2TimeProfile(newValue);
  }
  else if (command == fWLS2VerboseLevelCmd) {
    params->SetWLS2VerboseLevel(fWLS2VerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fAbsorptionVerboseLevelCmd) {
    params->SetAbsorptionVerboseLevel(fAbsorptionVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fRayleighVerboseLevelCmd) {
    params->SetRayleighVerboseLevel(fRayleighVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fMieVerboseLevelCmd) {
    params->SetMieVerboseLevel(fMieVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fBoundaryVerboseLevelCmd) {
    params->SetBoundaryVerboseLevel(fBoundaryVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if (command == fBoundaryInvokeSD1Cmd) {
    params->SetBoundaryInvokeSD(fBoundaryInvokeSD1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fBoundaryInvokeSDCmd) {
    params->SetBoundaryInvokeSD(fBoundaryInvokeSDCmd->GetNewBoolValue(newValue));
  }
  if (physicsModified) {
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }
}

void G4OpticalParametersMessenger::Deprecated()
{
    G4ExceptionDescription ed;
    ed <<" This command has been deprecated and will be removed in the next" << G4endl
       << "major release. Use the same command without /defaults/ instead.";
    G4Exception("G4OpticalParametersMessenger", "optical001", JustWarning, ed);
}
