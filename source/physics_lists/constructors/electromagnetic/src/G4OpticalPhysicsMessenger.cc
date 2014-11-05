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

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysicsMessenger::G4OpticalPhysicsMessenger(
                                            G4OpticalPhysics* opticalPhysics)
  : G4UImessenger(),
    fOpticalPhysics(opticalPhysics),
    fSelectedProcessIndex(kNoProcess),
    fActivateProcessCmd(0),
    fSetOpProcessVerboseCmd(0),
    fSetCerenkovMaxPhotonsCmd(0),
    fSetCerenkovMaxBetaChangeCmd(0),
    fSetScintillationYieldFactorCmd(0),
    fSetScintillationByParticleTypeCmd(0),
    fSetWLSTimeProfileCmd(0),
    fSetTrackSecondariesFirstCmd(0),
    fSetFiniteRiseTimeCmd(0)
{
    G4bool toBeBroadcasted = false;
    fDir = new G4UIdirectory("/process/optical/defaults/",toBeBroadcasted);
    fDir->SetGuidance("Commands related to the optical physics simulation engine.");
    fDir2 = new G4UIdirectory("/process/optical/",toBeBroadcasted);
    fDir2->SetGuidance("Commands related to the optical physics simulation engine.");

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


    fSetOpProcessVerboseCmd = new G4UIcmdWithAnInteger("/process/optical/verbose", this);
    fSetOpProcessVerboseCmd->SetGuidance("Set default verbosity level for optical processes");
    fSetOpProcessVerboseCmd->SetParameterName("ver", true);
    fSetOpProcessVerboseCmd->SetDefaultValue(1);
    fSetOpProcessVerboseCmd->SetRange("ver>=0");
    fSetOpProcessVerboseCmd->AvailableForStates(G4State_PreInit);

    
    fSetTrackSecondariesFirstCmd = new G4UIcommand("/process/optical/setTrackSecondariesFirst", this);
    fSetTrackSecondariesFirstCmd->SetGuidance("Activate/deactivate tracking of secondaries before finishing their parent track");
    par = new G4UIparameter("proc_name",'s',false);
    par->SetParameterCandidates(candidates);
    fSetTrackSecondariesFirstCmd->SetParameter(par);
    par = new G4UIparameter("flag",'b',false);
    par->SetDefaultValue(true);
    fSetTrackSecondariesFirstCmd->SetParameter(par);
    fSetTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit);

    //This are repetition of process specific ui commands needed by ATLICE and VMC
    //(UI messages needed to exist in PreInit state before actual processes are instantiated)
    fSetCerenkovMaxPhotonsCmd = new G4UIcmdWithAnInteger("/process/optical/defaults/cerenkov/setMaxPhotons", this);
    fSetCerenkovMaxPhotonsCmd->SetGuidance("Set default maximum number of photons per step");
    fSetCerenkovMaxPhotonsCmd->SetGuidance("Note this command is used to set the default value,");
    fSetCerenkovMaxPhotonsCmd->SetGuidance("if process is not active command will not have effect.");
    fSetCerenkovMaxPhotonsCmd->SetParameterName("CerenkovMaxPhotons", false);
    fSetCerenkovMaxPhotonsCmd->SetRange("CerenkovMaxPhotons>=0");
    fSetCerenkovMaxPhotonsCmd->AvailableForStates(G4State_PreInit);

    fSetCerenkovMaxBetaChangeCmd = new G4UIcmdWithADouble("/process/optical/defaults/cerenkov/setMaxBetaChange", this);
    fSetCerenkovMaxBetaChangeCmd->SetGuidance("Set default maximum change of beta of parent particle per step");
    fSetCerenkovMaxBetaChangeCmd->SetGuidance("Note this command is used to set the default value,");
    fSetCerenkovMaxBetaChangeCmd->SetGuidance("if process is not active command will not have effect.");
    fSetCerenkovMaxBetaChangeCmd->SetParameterName("CerenkovMaxBetaChange", false);
    fSetCerenkovMaxBetaChangeCmd->SetRange("CerenkovMaxBetaChange>=0");
    fSetCerenkovMaxBetaChangeCmd->AvailableForStates(G4State_PreInit);

    fSetScintillationYieldFactorCmd = new G4UIcmdWithADouble("/process/optical/defaults/scintillation/setYieldFactor", this);
    fSetScintillationYieldFactorCmd->SetGuidance("Set scintillation yield factor");
    fSetScintillationYieldFactorCmd->SetGuidance("Note this command is used to set the default value,");
    fSetScintillationYieldFactorCmd->SetGuidance("if process is not active command will not have effect.");
    fSetScintillationYieldFactorCmd->SetParameterName("ScintillationYieldFactor", false);
    fSetScintillationYieldFactorCmd->SetRange("ScintillationYieldFactor>=0");
    fSetScintillationYieldFactorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetScintillationByParticleTypeCmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setByParticleType", this);
    fSetScintillationByParticleTypeCmd->SetGuidance("Activate/Inactivate scintillation process by particle type");
    fSetScintillationByParticleTypeCmd->SetGuidance("Note this command is used to set the default value,");
    fSetScintillationByParticleTypeCmd->SetGuidance("if process is not active command will not have effect.");
    fSetScintillationByParticleTypeCmd->SetParameterName("ScintillationByParticleTypeActivation", false);
    fSetScintillationByParticleTypeCmd->AvailableForStates(G4State_PreInit);

    fSetFiniteRiseTimeCmd = new G4UIcmdWithABool("/process/optical/defaults/scintillation/setFiniteRiseTime", this);
    fSetFiniteRiseTimeCmd->SetGuidance("Set option of a finite rise-time for G4Scintillation");
    fSetFiniteRiseTimeCmd->SetGuidance("If set, the G4Scintillation process expects the user to have set the constant material property FAST/SLOWSCINTILLATIONRISETIME");
    fSetFiniteRiseTimeCmd->SetGuidance("Note this command is used to set the default value,");
    fSetFiniteRiseTimeCmd->SetGuidance("if process is not active command will not have effect.");
    fSetFiniteRiseTimeCmd->SetParameterName("FiniteRiseTime", false);
    fSetFiniteRiseTimeCmd->AvailableForStates(G4State_PreInit);

    fSetWLSTimeProfileCmd = new G4UIcmdWithAString("/process/optical/defaults/wls/setTimeProfile", this);
    fSetWLSTimeProfileCmd->SetGuidance("Set the WLS time profile (delta or exponential)");
    fSetWLSTimeProfileCmd->SetParameterName("WLSTimeProfile", false);
    fSetWLSTimeProfileCmd->SetCandidates("delta exponential");
    fSetWLSTimeProfileCmd->AvailableForStates(G4State_PreInit);
}

G4OpticalPhysicsMessenger::~G4OpticalPhysicsMessenger()
{
// Destructor

  delete fDir;
  delete fDir2;
  delete fActivateProcessCmd;
  delete fSetOpProcessVerboseCmd;
  delete fSetCerenkovMaxPhotonsCmd;
  delete fSetCerenkovMaxBetaChangeCmd;
  delete fSetScintillationYieldFactorCmd;
  delete fSetScintillationByParticleTypeCmd;
  delete fSetWLSTimeProfileCmd;
  delete fSetFiniteRiseTimeCmd;
}

#include <iostream>
void G4OpticalPhysicsMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
/// Apply command to the associated object.
  if (command == fActivateProcessCmd) {
    std::istringstream is(newValue.data());
    G4String pn;
    G4int flag;
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
        G4Exception("G4OpticalPhysicsMessenger::SetNewValue(...)","Optical001",FatalException,pn);
    }
    fOpticalPhysics->Configure(fSelectedProcessIndex,flag);
  }
  else if (command == fSetTrackSecondariesFirstCmd )
  {
      std::istringstream is(newValue.data());
      G4String pn;
      G4int flag;
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
          G4Exception("G4OpticalPhysicsMessenger::SetNewValue(...)","Optical001",FatalException,pn);
      }
      fOpticalPhysics->SetTrackSecondariesFirst(fSelectedProcessIndex,flag);
  }
  else if (command == fSetOpProcessVerboseCmd) {
        fOpticalPhysics->SetVerboseLevel(fSetOpProcessVerboseCmd->GetNewIntValue(newValue));
  }
  else if (command == fSetCerenkovMaxPhotonsCmd) {
    fOpticalPhysics
      ->SetMaxNumPhotonsPerStep(
          fSetCerenkovMaxPhotonsCmd->GetNewIntValue(newValue));
  }
  else if (command == fSetCerenkovMaxBetaChangeCmd) {
    fOpticalPhysics
      ->SetMaxBetaChangePerStep(
          fSetCerenkovMaxBetaChangeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSetScintillationYieldFactorCmd) {
    fOpticalPhysics
      ->SetScintillationYieldFactor(
          fSetScintillationYieldFactorCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSetScintillationByParticleTypeCmd) {
    fOpticalPhysics
      ->SetScintillationByParticleType(
         fSetScintillationByParticleTypeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fSetFiniteRiseTimeCmd) {
    fOpticalPhysics
      ->SetFiniteRiseTime(
         fSetFiniteRiseTimeCmd->GetNewBoolValue(newValue));
  }
  else if (command == fSetWLSTimeProfileCmd) {
      fOpticalPhysics->SetWLSTimeProfile(newValue);
  }
}
