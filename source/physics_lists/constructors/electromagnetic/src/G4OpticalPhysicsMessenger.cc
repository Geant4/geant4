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
// ClassName:   G4OpticalPhysicsMessenger
//
// Author:      P.Gumplinger 30.09.2009 //
//
// Modified:    P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
//----------------------------------------------------------------------------
//

#include "G4OpticalPhysicsMessenger.hh"
#include "G4OpticalPhysics.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIparameter.hh"


// Commands with '/defaults/' are duplicates and will be removed in
// the next major release of Geant4. Use commands with no /defaults/ instead


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysicsMessenger::G4OpticalPhysicsMessenger(
                                            G4OpticalPhysics* opticalPhysics)
  : G4UImessenger(),
    fOpticalPhysics(opticalPhysics),
    fSelectedProcessIndex(kNoProcess),
    fActivateProcessCmd(nullptr),
    fVerboseCmd(nullptr),
    fTrackSecondariesFirstCmd(nullptr),

    fCerenkovMaxPhotonsCmd(nullptr),
    fCerenkovMaxPhotons1Cmd(nullptr),
    fCerenkovMaxBetaChangeCmd(nullptr),
    fCerenkovMaxBetaChange1Cmd(nullptr),
    fCerenkovStackPhotonsCmd(nullptr),
    fCerenkovStackPhotons1Cmd(nullptr),
    fCerenkovTrackSecondariesFirstCmd(nullptr),
    fCerenkovVerbosityCmd(nullptr),

    fScintYieldFactorCmd(nullptr),
    fScintYieldFactor1Cmd(nullptr),
    fScintByParticleTypeCmd(nullptr),
    fScintByParticleType1Cmd(nullptr),
    fScintTrackInfoCmd(nullptr),
    fScintTrackInfo1Cmd(nullptr),
    fScintStackPhotonsCmd(nullptr),
    fScintStackPhotons1Cmd(nullptr),
    fScintTrackSecondariesFirstCmd(nullptr),
    fScintFiniteRiseTimeCmd(nullptr),
    fScintFiniteRiseTime1Cmd(nullptr),
    fScintVerbosityCmd(nullptr),

    fWLSTimeProfileCmd(nullptr),
    fWLSTimeProfile1Cmd(nullptr),
    fWLSVerbosityCmd(nullptr),

    fBoundaryInvokeSDCmd(nullptr),
    fBoundaryInvokeSD1Cmd(nullptr),
    fBoundaryVerbosityCmd(nullptr),

    fAbsorptionVerbosityCmd(nullptr),
    fRayleighVerbosityCmd(nullptr),
    fMieVerbosityCmd(nullptr)

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

    fTrackSecondariesFirstCmd = new G4UIcommand("/process/optical/setTrackSecondariesFirst", this);
    fTrackSecondariesFirstCmd->SetGuidance("Activate/deactivate tracking of secondaries before finishing their parent track");
    par = new G4UIparameter("proc_name",'s',false);
    par->SetParameterCandidates(candidates);
    fTrackSecondariesFirstCmd->SetParameter(par);
    par = new G4UIparameter("flag",'b',false);
    par->SetDefaultValue(true);
    fTrackSecondariesFirstCmd->SetParameter(par);
    fTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit);

    fVerboseCmd = new G4UIcmdWithAnInteger("/process/optical/verbose", this);
    fVerboseCmd->SetGuidance("Set default verbosity level for optical processes");
    fVerboseCmd->SetParameterName("ver", true);
    fVerboseCmd->SetDefaultValue(1);
    fVerboseCmd->SetRange("ver>=0");
    fVerboseCmd->AvailableForStates(G4State_PreInit);

    //// Cerenkov ////////////////////
    fCerenkovMaxPhotons1Cmd = new G4UIcmdWithAnInteger("/process/optical/defaults/cerenkov/setMaxPhotons", this);
    fCerenkovMaxPhotons1Cmd->SetGuidance("Set default maximum number of photons per step");
    fCerenkovMaxPhotons1Cmd->SetGuidance("Note this command is used to set the default value,");
    fCerenkovMaxPhotons1Cmd->SetGuidance("if process is not active command will not have effect.");
    fCerenkovMaxPhotons1Cmd->SetGuidance("DEPRECATED: use /process/optical/cerenkov/setMaxPhotons instead.");
    fCerenkovMaxPhotons1Cmd->SetParameterName("CerenkovMaxPhotons", false);
    fCerenkovMaxPhotons1Cmd->SetRange("CerenkovMaxPhotons>=0");
    fCerenkovMaxPhotons1Cmd->AvailableForStates(G4State_PreInit);

    fCerenkovMaxPhotonsCmd = new G4UIcmdWithAnInteger("/process/optical/cerenkov/setMaxPhotons", this);
    fCerenkovMaxPhotonsCmd->SetGuidance("Set default maximum number of photons per step");
    fCerenkovMaxPhotonsCmd->SetGuidance("Note this command is used to set the default value,");
    fCerenkovMaxPhotonsCmd->SetGuidance("if process is not active command will not have effect.");
    fCerenkovMaxPhotonsCmd->SetParameterName("CerenkovMaxPhotons", false);
    fCerenkovMaxPhotonsCmd->SetRange("CerenkovMaxPhotons>=0");
    fCerenkovMaxPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovMaxBetaChange1Cmd = new G4UIcmdWithADouble("/process/optical/defaults/cerenkov/setMaxBetaChange", this);
    fCerenkovMaxBetaChange1Cmd->SetGuidance("Set default maximum change of beta of parent particle per step");
    fCerenkovMaxBetaChange1Cmd->SetGuidance("Note this command is used to set the default value,");
    fCerenkovMaxBetaChange1Cmd->SetGuidance("if process is not active command will not have effect.");
    fCerenkovMaxBetaChange1Cmd->SetGuidance("DEPRECATED: use /process/optical/cerenkov/setMaxBetaChange instead.");
    fCerenkovMaxBetaChange1Cmd->SetParameterName("CerenkovMaxBetaChange", false);
    fCerenkovMaxBetaChange1Cmd->SetRange("CerenkovMaxBetaChange>=0");
    fCerenkovMaxBetaChange1Cmd->AvailableForStates(G4State_PreInit);

    fCerenkovMaxBetaChangeCmd = new G4UIcmdWithADouble("/process/optical/cerenkov/setMaxBetaChange", this);
    fCerenkovMaxBetaChangeCmd->SetGuidance("Set default maximum change of beta of parent particle per step");
    fCerenkovMaxBetaChangeCmd->SetGuidance("Note this command is used to set the default value,");
    fCerenkovMaxBetaChangeCmd->SetGuidance("if process is not active command will not have effect.");
    fCerenkovMaxBetaChangeCmd->SetParameterName("CerenkovMaxBetaChange", false);
    fCerenkovMaxBetaChangeCmd->SetRange("CerenkovMaxBetaChange>=0");
    fCerenkovMaxBetaChangeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovStackPhotons1Cmd = new G4UIcmdWithABool("/process/optical/defaults/cerenkov/setStackPhotons", this);
    fCerenkovStackPhotons1Cmd->SetGuidance("Set default whether or not to stack secondary Cerenkov photons");
    fCerenkovStackPhotons1Cmd->SetGuidance("Note this command is used to set the default value,");
    fCerenkovStackPhotons1Cmd->SetGuidance("if process is not active command will not have effect.");
    fCerenkovStackPhotons1Cmd->SetGuidance("DEPRECATED: use /process/optical/cerenkov/setStackPhotons instead.");
    fCerenkovStackPhotons1Cmd->SetParameterName("CerenkovStackPhotons", true);
    fCerenkovStackPhotons1Cmd->AvailableForStates(G4State_PreInit);

    fCerenkovStackPhotonsCmd = new G4UIcmdWithABool("/process/optical/cerenkov/setStackPhotons", this);
    fCerenkovStackPhotonsCmd->SetGuidance("Set default whether or not to stack secondary Cerenkov photons");
    fCerenkovStackPhotonsCmd->SetGuidance("Note this command is used to set the default value,");
    fCerenkovStackPhotonsCmd->SetGuidance("if process is not active command will not have effect.");
    fCerenkovStackPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovTrackSecondariesFirstCmd = new G4UIcmdWithABool("/process/optical/cerenkov/setTrackSecondariesFirst", this);
    fCerenkovTrackSecondariesFirstCmd->SetGuidance("Whether to track secondary Cerenkov photons before the primary.");
    fCerenkovTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fCerenkovVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/cerenkov/verbose", this);
    fCerenkovVerbosityCmd->SetGuidance("Verbosity for Cerenkov process.");
    fCerenkovVerbosityCmd->SetParameterName("verbosity", true);
    fCerenkovVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fCerenkovVerbosityCmd->AvailableForStates(G4State_Idle);

    // Scintillation //////////////////////////
    fScintYieldFactor1Cmd = new G4UIcmdWithADouble("/process/optical/defaults/scintillation/setYieldFactor", this);
    fScintYieldFactor1Cmd->SetGuidance("Set scintillation yield factor");
    fScintYieldFactor1Cmd->SetGuidance("Note this command is used to set the default value,");
    fScintYieldFactor1Cmd->SetGuidance("if process is not active command will not have effect.");
    fScintYieldFactor1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setYieldFactorinstead.");
    fScintYieldFactor1Cmd->SetParameterName("ScintillationYieldFactor", false);
    fScintYieldFactor1Cmd->SetRange("ScintillationYieldFactor>=0");
    fScintYieldFactor1Cmd->AvailableForStates(G4State_PreInit);

    fScintYieldFactorCmd = new G4UIcmdWithADouble("/process/optical/scintillation/setYieldFactor", this);
    fScintYieldFactorCmd->SetGuidance("Set scintillation yield factor");
    fScintYieldFactorCmd->SetGuidance("Note this command is used to set the default value,");
    fScintYieldFactorCmd->SetGuidance("if process is not active command will not have effect.");
    fScintYieldFactorCmd->SetParameterName("ScintillationYieldFactor", false);
    fScintYieldFactorCmd->SetRange("ScintillationYieldFactor>=0");
    fScintYieldFactorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintExcitationRatioCmd = new G4UIcmdWithADouble("/process/optical/scintillation/setExcitationRatio", this);
    fScintExcitationRatioCmd->SetGuidance("Set scintillation excitation ratio");
    fScintExcitationRatioCmd->SetGuidance("Note this command is used to set the default value,");
    fScintExcitationRatioCmd->SetGuidance("if process is not active command will not have effect.");
    fScintExcitationRatioCmd->SetParameterName("ExcitationRatio", false);
    fScintExcitationRatioCmd->SetRange("ExcitationRatio >= 0 && ExcitationRatio <=1");
    fScintExcitationRatioCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintByParticleType1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setByParticleType", this);
    fScintByParticleType1Cmd->SetGuidance("Activate/Inactivate scintillation process by particle type");
    fScintByParticleType1Cmd->SetGuidance("Note this command is used to set the default value,");
    fScintByParticleType1Cmd->SetGuidance("if process is not active command will not have effect.");
    fScintByParticleType1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setByParticleType instead.");
    fScintByParticleType1Cmd->SetParameterName("ScintillationByParticleTypeActivation", false);
    fScintByParticleType1Cmd->AvailableForStates(G4State_PreInit);

    fScintByParticleTypeCmd = new G4UIcmdWithABool("/process/optical/scintillation/setByParticleType", this);
    fScintByParticleTypeCmd->SetGuidance("Activate/Inactivate scintillation process by particle type");
    fScintByParticleTypeCmd->SetGuidance("Note this command is used to set the default value,");
    fScintByParticleTypeCmd->SetGuidance("if process is not active command will not have effect.");
    fScintByParticleTypeCmd->SetParameterName("ScintillationByParticleTypeActivation", false);
    fScintByParticleTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintTrackInfo1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setTrackInfo", this);
    fScintTrackInfo1Cmd->SetGuidance("Activate/Inactivate scintillation TrackInformation");
    fScintTrackInfo1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setTrackInfo instead.");
    fScintTrackInfo1Cmd->SetParameterName("ScintillationTrackInfo", false);
    fScintTrackInfo1Cmd->AvailableForStates(G4State_PreInit);

    fScintTrackInfoCmd = new G4UIcmdWithABool("/process/optical/scintillation/setTrackInfo", this);
    fScintTrackInfoCmd->SetGuidance("Activate/Inactivate scintillation TrackInformation");
    fScintTrackInfoCmd->SetParameterName("ScintillationTrackInfo", false);
    fScintTrackInfoCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintFiniteRiseTime1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setFiniteRiseTime", this);
    fScintFiniteRiseTime1Cmd->SetGuidance("Set option of a finite rise-time for G4Scintillation");
    fScintFiniteRiseTime1Cmd->SetGuidance("If set, the G4Scintillation process expects the user to have set the");
    fScintFiniteRiseTime1Cmd->SetGuidance("constant material property FAST/SLOWSCINTILLATIONRISETIME");
    fScintFiniteRiseTime1Cmd->SetGuidance("Note this command is used to set the default value,");
    fScintFiniteRiseTime1Cmd->SetGuidance("if process is not active command will not have effect.");
    fScintFiniteRiseTime1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setFiniteRiseTime instead.");
    fScintFiniteRiseTime1Cmd->SetParameterName("FiniteRiseTime", false);
    fScintFiniteRiseTime1Cmd->AvailableForStates(G4State_PreInit);

    fScintFiniteRiseTimeCmd = new G4UIcmdWithABool("/process/optical/scintillation/setFiniteRiseTime", this);
    fScintFiniteRiseTimeCmd->SetGuidance("Set option of a finite rise-time for G4Scintillation");
    fScintFiniteRiseTimeCmd->SetGuidance("If set, the G4Scintillation process expects the user to have set the");
    fScintFiniteRiseTimeCmd->SetGuidance("constant material property FAST/SLOWSCINTILLATIONRISETIME");
    fScintFiniteRiseTimeCmd->SetGuidance("Note this command is used to set the default value,");
    fScintFiniteRiseTimeCmd->SetGuidance("if process is not active command will not have effect.");
    fScintFiniteRiseTimeCmd->SetParameterName("FiniteRiseTime", false);
    fScintFiniteRiseTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintStackPhotons1Cmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setStackPhotons", this);
    fScintStackPhotons1Cmd->SetGuidance("Set default whether or not to stack secondary Scintillation photons");
    fScintStackPhotons1Cmd->SetGuidance("Note this command is used to set the default value,");
    fScintStackPhotons1Cmd->SetGuidance("if process is not active command will not have effect.");
    fScintStackPhotons1Cmd->SetGuidance("DEPRECATED: use /process/optical/scintillation/setStackPhotons instead.");
    fScintStackPhotons1Cmd->SetParameterName("ScintillationStackPhotons", true);
    fScintStackPhotons1Cmd->AvailableForStates(G4State_PreInit);

    fScintStackPhotonsCmd = new G4UIcmdWithABool("/process/optical/scintillation/setStackPhotons", this);
    fScintStackPhotonsCmd->SetGuidance("Set default whether or not to stack secondary Scintillation photons");
    fScintStackPhotonsCmd->SetGuidance("Note this command is used to set the default value,");
    fScintStackPhotonsCmd->SetGuidance("if process is not active command will not have effect.");
    fScintStackPhotonsCmd->SetParameterName("ScintillationStackPhotons", true);
    fScintStackPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintTrackSecondariesFirstCmd = new G4UIcmdWithABool("/process/optical/scintillation/setTrackSecondariesFirst", this);
    fScintTrackSecondariesFirstCmd->SetGuidance("Whether to track scintillation secondaries before primary.");
    fScintTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fScintVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/scintillation/verbose", this);
    fScintVerbosityCmd->SetGuidance("Verbosity for scintillation process.");
    fScintVerbosityCmd->SetParameterName("verbosity", true);
    fScintVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fScintVerbosityCmd->AvailableForStates(G4State_Idle);

    // WLS   //////////////////////////////////
    fWLSTimeProfile1Cmd = new G4UIcmdWithAString("/process/optical/defaults/wls/setTimeProfile", this);
    fWLSTimeProfile1Cmd->SetGuidance("Set the WLS time profile (delta or exponential)");
    fWLSTimeProfile1Cmd->SetGuidance("DEPRECATED: use /process/optical/wls/setTimeProfile instead.");
    fWLSTimeProfile1Cmd->SetParameterName("WLSTimeProfile", false);
    fWLSTimeProfile1Cmd->SetCandidates("delta exponential");
    fWLSTimeProfile1Cmd->AvailableForStates(G4State_PreInit);

    fWLSTimeProfileCmd = new G4UIcmdWithAString("/process/optical/wls/setTimeProfile", this);
    fWLSTimeProfileCmd->SetGuidance("Set the WLS time profile (delta or exponential)");
    fWLSTimeProfileCmd->SetParameterName("WLSTimeProfile", false);
    fWLSTimeProfileCmd->SetCandidates("delta exponential");
    fWLSTimeProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fWLSVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/wls/verbose", this);
    fWLSVerbosityCmd->SetGuidance("Verbosity for WLS process.");
    fWLSVerbosityCmd->SetParameterName("verbosity", true);
    fWLSVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fWLSVerbosityCmd->AvailableForStates(G4State_Idle);

    // boundary //////////////////////////////////////
    fBoundaryInvokeSD1Cmd = new G4UIcmdWithABool("/process/optical/defaults/boundary/setInvokeSD", this);
    fBoundaryInvokeSD1Cmd->SetGuidance("Set option for calling InvokeSD in G4OpBoundaryProcess");
    fBoundaryInvokeSD1Cmd->SetGuidance("DEPRECATED: use /process/optical/boundary/setInvokeSD instead.");
    fBoundaryInvokeSD1Cmd->SetParameterName("InvokeSD", false);
    fBoundaryInvokeSD1Cmd->AvailableForStates(G4State_PreInit);

    fBoundaryInvokeSDCmd = new G4UIcmdWithABool("/process/optical/boundary/setInvokeSD", this);
    fBoundaryInvokeSDCmd->SetGuidance("Set option for calling InvokeSD in G4OpBoundaryProcess");
    fBoundaryInvokeSDCmd->SetParameterName("InvokeSD", false);
    fBoundaryInvokeSDCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBoundaryVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/boundary/verbose", this);
    fBoundaryVerbosityCmd->SetGuidance("Verbosity for boundary process.");
    fBoundaryVerbosityCmd->SetParameterName("verbosity", true);
    fBoundaryVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fBoundaryVerbosityCmd->AvailableForStates(G4State_Idle);

    // the others  ////////////////////////////////////
    fAbsorptionVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/absorption/verbose", this);
    fAbsorptionVerbosityCmd->SetGuidance("Verbosity for absorption process.");
    fAbsorptionVerbosityCmd->SetParameterName("verbosity", true);
    fAbsorptionVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fAbsorptionVerbosityCmd->AvailableForStates(G4State_Idle);

    fRayleighVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/rayleigh/verbose", this);
    fRayleighVerbosityCmd->SetGuidance("Verbosity for Rayleigh process.");
    fRayleighVerbosityCmd->SetParameterName("verbosity", true);
    fRayleighVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fRayleighVerbosityCmd->AvailableForStates(G4State_Idle);

    fMieVerbosityCmd = new G4UIcmdWithAnInteger("/process/optical/mie/verbose", this);
    fMieVerbosityCmd->SetGuidance("Verbosity for Mie process.");
    fMieVerbosityCmd->SetParameterName("verbosity", true);
    fMieVerbosityCmd->SetRange("verbosity >= 0 && verbosity <= 2");
    fMieVerbosityCmd->AvailableForStates(G4State_Idle);
}

G4OpticalPhysicsMessenger::~G4OpticalPhysicsMessenger()
{
  delete fDir;
  delete fDir2;
  delete fActivateProcessCmd;
  delete fVerboseCmd;
  delete fCerenkovMaxPhotonsCmd;
  delete fCerenkovMaxPhotons1Cmd;
  delete fCerenkovMaxBetaChangeCmd;
  delete fCerenkovMaxBetaChange1Cmd;
  delete fCerenkovStackPhotonsCmd;
  delete fCerenkovStackPhotons1Cmd;
  delete fCerenkovTrackSecondariesFirstCmd;
  delete fCerenkovVerbosityCmd;
  delete fScintYieldFactorCmd;
  delete fScintYieldFactor1Cmd;
  delete fScintByParticleTypeCmd;
  delete fScintByParticleType1Cmd;
  delete fScintTrackInfoCmd;
  delete fScintTrackInfo1Cmd;
  delete fScintStackPhotonsCmd;
  delete fScintStackPhotons1Cmd;
  delete fScintExcitationRatioCmd;
  delete fScintVerbosityCmd;
  delete fScintFiniteRiseTimeCmd;
  delete fScintFiniteRiseTime1Cmd;
  delete fScintTrackSecondariesFirstCmd;
  delete fWLSTimeProfileCmd;
  delete fWLSTimeProfile1Cmd;
  delete fWLSVerbosityCmd;
  delete fAbsorptionVerbosityCmd;
  delete fRayleighVerbosityCmd;
  delete fMieVerbosityCmd;
  delete fBoundaryVerbosityCmd;
  delete fTrackSecondariesFirstCmd;
  delete fBoundaryInvokeSDCmd;
  delete fBoundaryInvokeSD1Cmd;
}

void G4OpticalPhysicsMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
/// Apply command to the associated object.
  if (command == fActivateProcessCmd) {
    std::istringstream is(newValue.data());
    G4String pn;
    G4String flag;
    is >> pn >> flag;
    if  ( pn == "Cerenkov" )        {
        fSelectedProcessIndex = kCerenkov;
    } else if ( pn == "Scintillation" ) {
        fSelectedProcessIndex = kScintillation;
    } else if ( pn == "OpAbsorption" )  {
        fSelectedProcessIndex = kAbsorption;
    } else if ( pn == "OpRayleigh" )    {
        fSelectedProcessIndex = kRayleigh;
    } else if ( pn == "OpMieHG" )       {
        fSelectedProcessIndex = kMieHG;
    } else if ( pn == "OpBoundary" )    {
        fSelectedProcessIndex = kBoundary;
    } else if ( pn == "OpWLS" )         {
        fSelectedProcessIndex = kWLS;
    } else {
        G4ExceptionDescription msg;
        msg << "Not allowed process name: "<<pn<<" (UI: "<<newValue<<")";
        G4Exception("G4OpticalPhysicsMessenger::SetNewValue(...)","Optical001",FatalException,msg);
    }
    G4bool value = G4UIcommand::ConvertToBool(flag);
    fOpticalPhysics->Configure(fSelectedProcessIndex,value);
  }
  else if (command == fTrackSecondariesFirstCmd )
  {
      std::istringstream is(newValue.data());
      G4String pn;
      G4String flag;
      is >> pn >> flag;
      if ( pn == "Cerenkov" )        {
        fSelectedProcessIndex = kCerenkov;
      } else if ( pn == "Scintillation" ) {
        fSelectedProcessIndex = kScintillation;
      } else if ( pn == "OpAbsorption" )  {
        fSelectedProcessIndex = kAbsorption;
      } else if ( pn == "OpRayleigh" )    {
        fSelectedProcessIndex = kRayleigh;
      } else if ( pn == "OpMieHG" )       {
        fSelectedProcessIndex = kMieHG;
      } else if ( pn == "OpBoundary" )    {
        fSelectedProcessIndex = kBoundary;
      } else if ( pn == "OpWLS" )         {
        fSelectedProcessIndex = kWLS;
      } else {
          G4ExceptionDescription msg;
          msg << "Not allowed process name: "<<pn<<" (UI: "<<newValue<<")";
          G4Exception("G4OpticalPhysicsMessenger::SetNewValue(...)","Optical001",FatalException,msg);
      }
      G4bool value = G4UIcommand::ConvertToBool(flag);
      fOpticalPhysics->SetTrackSecondariesFirst(fSelectedProcessIndex,value);
  }
  else if (command == fVerboseCmd) {
        fOpticalPhysics->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if (command == fCerenkovMaxPhotons1Cmd) {
    fOpticalPhysics->SetMaxNumPhotonsPerStep(
          fCerenkovMaxPhotons1Cmd->GetNewIntValue(newValue));
    Deprecated();
  }
  else if (command == fCerenkovMaxPhotonsCmd) {
    fOpticalPhysics->SetMaxNumPhotonsPerStep(
          fCerenkovMaxPhotonsCmd->GetNewIntValue(newValue));
  }
  else if (command == fCerenkovMaxBetaChange1Cmd) {
    fOpticalPhysics->SetMaxBetaChangePerStep(
          fCerenkovMaxBetaChange1Cmd->GetNewDoubleValue(newValue));
    Deprecated();
  }
  else if (command == fCerenkovMaxBetaChangeCmd) {
    fOpticalPhysics->SetMaxBetaChangePerStep(
          fCerenkovMaxBetaChangeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fCerenkovStackPhotons1Cmd) {
    fOpticalPhysics->SetCerenkovStackPhotons(
          fCerenkovStackPhotons1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fCerenkovStackPhotonsCmd) {
    fOpticalPhysics->SetCerenkovStackPhotons(
          fCerenkovStackPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fCerenkovTrackSecondariesFirstCmd) {
    fOpticalPhysics->SetCerenkovTrackSecondariesFirst(
          fCerenkovTrackSecondariesFirstCmd->GetNewBoolValue(newValue));
  }
  else if (command == fCerenkovVerbosityCmd) {
    fOpticalPhysics->SetCerenkovVerbosity(
          fCerenkovVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fScintYieldFactor1Cmd) {
    fOpticalPhysics->SetScintillationYieldFactor(
          fScintYieldFactor1Cmd->GetNewDoubleValue(newValue));
    Deprecated();
  }
  else if (command == fScintYieldFactorCmd) {
    fOpticalPhysics->SetScintillationYieldFactor(
          fScintYieldFactorCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fScintByParticleType1Cmd) {
    fOpticalPhysics->SetScintillationByParticleType(
         fScintByParticleType1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintByParticleTypeCmd) {
    fOpticalPhysics->SetScintillationByParticleType(
         fScintByParticleTypeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintTrackInfo1Cmd) {
    fOpticalPhysics->SetScintillationTrackInfo(
         fScintTrackInfo1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintTrackInfoCmd) {
    fOpticalPhysics->SetScintillationTrackInfo(
         fScintTrackInfoCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintFiniteRiseTime1Cmd) {
    fOpticalPhysics->SetFiniteRiseTime(
         fScintFiniteRiseTime1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintFiniteRiseTimeCmd) {
    fOpticalPhysics->SetFiniteRiseTime(
         fScintFiniteRiseTimeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintStackPhotons1Cmd) {
    fOpticalPhysics->SetScintillationStackPhotons(
          fScintStackPhotons1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fScintStackPhotonsCmd) {
    fOpticalPhysics->SetScintillationStackPhotons(
          fScintStackPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintExcitationRatioCmd) {
    fOpticalPhysics->SetScintillationExcitationRatio(
          fScintExcitationRatioCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fScintTrackSecondariesFirstCmd) {
    fOpticalPhysics->SetScintillationTrackSecondariesFirst(
          fScintTrackSecondariesFirstCmd->GetNewBoolValue(newValue));
  }
  else if (command == fScintVerbosityCmd) {
    fOpticalPhysics->SetScintillationVerbosity(
          fScintVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fWLSTimeProfile1Cmd) {
    fOpticalPhysics->SetWLSTimeProfile(newValue);
    Deprecated();
  }
  else if (command == fWLSTimeProfileCmd) {
    fOpticalPhysics->SetWLSTimeProfile(newValue);
  }
  else if (command == fWLSVerbosityCmd) {
    fOpticalPhysics->SetWLSVerbosity(fWLSVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fAbsorptionVerbosityCmd) {
    fOpticalPhysics->SetAbsorptionVerbosity(fAbsorptionVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fRayleighVerbosityCmd) {
    fOpticalPhysics->SetRayleighVerbosity(fRayleighVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fMieVerbosityCmd) {
    fOpticalPhysics->SetMieVerbosity(fMieVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fBoundaryVerbosityCmd) {
    fOpticalPhysics->SetBoundaryVerbosity(fBoundaryVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fBoundaryInvokeSD1Cmd) {
    fOpticalPhysics->SetInvokeSD(fBoundaryInvokeSD1Cmd->GetNewBoolValue(newValue));
    Deprecated();
  }
  else if (command == fBoundaryInvokeSDCmd) {
    fOpticalPhysics
      ->SetInvokeSD(fBoundaryInvokeSDCmd->GetNewBoolValue(newValue));
  }
}

void G4OpticalPhysicsMessenger::Deprecated()
{
    G4ExceptionDescription ed;
    ed <<" This command has been deprecated and will be removed in the next" << G4endl
       << "major release. Use the same command without /defaults/ instead.";
    G4Exception("G4OpticalPhysicsMessenger", "optical001", JustWarning, ed);
}
