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

#include "G4ParticleTable.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4OpticalPhysicsMessenger::G4OpticalPhysicsMessenger(
                                            G4OpticalPhysics* opticalPhysics)
  : G4UImessenger(),
    fOpticalPhysics(opticalPhysics),
    fSelectedProcessIndex(kNoProcess),
    fSelectOpProcessCmd(0),
    fSetOpProcessUseCmd(0),
    fSetOpProcessVerboseCmd(0),
    fSetCerenkovMaxPhotonsCmd(0),
    fSetCerenkovMaxBetaChangeCmd(0),
    fSetScintillationYieldFactorCmd(0),
    fSetScintillationByParticleTypeCmd(0),
//    fSetOpticalSurfaceModelCmd(0),
    fSetWLSTimeProfileCmd(0),
    fSetTrackSecondariesFirstCmd(0),
    fSetFiniteRiseTimeCmd(0)
{
  fDir = new G4UIdirectory("/optics_engine/");
  fDir->SetGuidance("Commands related to the optical physics simulation engine.");
  
  fSelectOpProcessCmd
   = new G4UIcmdWithAString("/optics_engine/selectOpProcess", this);
  fSelectOpProcessCmd
   ->SetGuidance("Select optical process for applying use/verbose/trackfirst commands");
  fSelectOpProcessCmd->SetParameterName("OpProcess", false);
  G4String candidates;
  for ( G4int i=0; i<kNoProcess; i++ ) {
      candidates += G4OpticalProcessName(i);
      candidates += G4String(" ");
  }
  fSelectOpProcessCmd->SetCandidates(candidates);
  fSelectOpProcessCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  fSetOpProcessUseCmd
   = new G4UIcmdWithABool("/optics_engine/setOpProcessUse", this);
  fSetOpProcessUseCmd->SetGuidance("Use/Not use selected optical process");
  fSetOpProcessUseCmd->SetParameterName("OpProcessUse", false);
  fSetOpProcessUseCmd->AvailableForStates(G4State_PreInit);

  fSetOpProcessVerboseCmd 
    = new G4UIcmdWithAnInteger("/optics_engine/setOpProcessVerbose", this);  
  fSetOpProcessVerboseCmd->SetGuidance("Set verbosity level for selected optical process");
  fSetOpProcessVerboseCmd->SetParameterName("OpProcessVerbose", false);
  fSetOpProcessVerboseCmd->SetRange("OpProcessVerbose>=0");
  fSetOpProcessVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc); 

  fSetCerenkovMaxPhotonsCmd 
    = new G4UIcmdWithAnInteger("/optics_engine/setCerenkovMaxPhotons", this);  
  fSetCerenkovMaxPhotonsCmd->SetGuidance("Set maximum number of photons per step");
  fSetCerenkovMaxPhotonsCmd->SetParameterName("CerenkovMaxPhotons", false);
  fSetCerenkovMaxPhotonsCmd->SetRange("CerenkovMaxPhotons>=0");
  fSetCerenkovMaxPhotonsCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);  

  fSetCerenkovMaxBetaChangeCmd 
    = new G4UIcmdWithADouble("/optics_engine/setCerenkovMaxBetaChange", this);  
  fSetCerenkovMaxBetaChangeCmd
    ->SetGuidance("Set maximum change of beta of parent particle per step");
  fSetCerenkovMaxBetaChangeCmd->SetParameterName("CerenkovMaxBetaChange", false);
  fSetCerenkovMaxBetaChangeCmd->SetRange("CerenkovMaxBetaChange>=0");
  fSetCerenkovMaxBetaChangeCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  fSetScintillationYieldFactorCmd 
    = new G4UIcmdWithADouble("/optics_engine/setScintillationYieldFactor", this);  
  fSetScintillationYieldFactorCmd->SetGuidance("Set scintillation yield factor");
  fSetScintillationYieldFactorCmd->SetParameterName("ScintillationYieldFactor", false);
  fSetScintillationYieldFactorCmd->SetRange("ScintillationYieldFactor>=0");
  fSetScintillationYieldFactorCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  fSetScintillationByParticleTypeCmd
   = new G4UIcmdWithABool("/optics_engine/setScintillationByParticleType", this);
  fSetScintillationByParticleTypeCmd->SetGuidance("Activate/Inactivate scintillation process by particle type");
  fSetScintillationByParticleTypeCmd->SetParameterName("ScintillationByParticleTypeActivation", false);
  fSetScintillationByParticleTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

//  fSetOpticalSurfaceModelCmd 
//    = new G4UIcmdWithAString("/optics_engine/setOpticalSurfaceModel", this);  
//  fSetOpticalSurfaceModelCmd
//    ->SetGuidance("Set optical surface model (glisur or unified)");
//  fSetOpticalSurfaceModelCmd->SetParameterName("OpticalSurfaceModel", false);
//  fSetOpticalSurfaceModelCmd->SetCandidates("glisur unified");
//  fSetOpticalSurfaceModelCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  fSetWLSTimeProfileCmd
    = new G4UIcmdWithAString("/optics_engine/setWLSTimeProfile", this);
  fSetWLSTimeProfileCmd
    ->SetGuidance("Set the WLS time profile (delta or exponential)");
  fSetWLSTimeProfileCmd->SetParameterName("WLSTimeProfile", false);
  fSetWLSTimeProfileCmd->SetCandidates("delta exponential");
  fSetWLSTimeProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  fSetTrackSecondariesFirstCmd 
    = new G4UIcmdWithABool("/optics_engine/setTrackSecondariesFirst", this);  
  fSetTrackSecondariesFirstCmd
    ->SetGuidance("Set option to track secondaries before finishing their parent track");
  fSetTrackSecondariesFirstCmd->SetParameterName("TrackSecondariesFirst", false);
  fSetTrackSecondariesFirstCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  fSetFiniteRiseTimeCmd
    = new G4UIcmdWithABool("/optics_engine/setFiniteRiseTime", this);
  fSetFiniteRiseTimeCmd
     ->SetGuidance("Set option of a finite rise-time for G4Scintillation - If set, the G4Scintillation process expects the user to have set the constant material property FAST/SLOWSCINTILLATIONRISETIME");
  fSetFiniteRiseTimeCmd->SetParameterName("FiniteRiseTime", false);
  fSetFiniteRiseTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed, G4State_EventProc);
}

G4OpticalPhysicsMessenger::~G4OpticalPhysicsMessenger()
{
// Destructor

  delete fDir;
  delete fSelectOpProcessCmd; 
  delete fSetOpProcessUseCmd; 
  delete fSetOpProcessVerboseCmd;
  delete fSetCerenkovMaxPhotonsCmd;
  delete fSetCerenkovMaxBetaChangeCmd;
  delete fSetScintillationYieldFactorCmd;
  delete fSetScintillationByParticleTypeCmd;
//  delete fSetOpticalSurfaceModelCmd;
  delete fSetWLSTimeProfileCmd;
  delete fSetTrackSecondariesFirstCmd;
  delete fSetFiniteRiseTimeCmd;
}

void G4OpticalPhysicsMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
/// Apply command to the associated object.

  if (command == fSelectOpProcessCmd) {
    if      ( newValue == "Cerenkov" )        {
            fSelectedProcessIndex = kCerenkov;
    } else if ( newValue == "Scintillation" ) {
            fSelectedProcessIndex = kScintillation;
    } else if ( newValue == "OpAbsorption" )  {
            fSelectedProcessIndex = kAbsorption;
    } else if ( newValue == "OpRayleigh" )    {
            fSelectedProcessIndex = kRayleigh;
    } else if ( newValue == "OpMieHG" )       {
            fSelectedProcessIndex = kMieHG;
    } else if ( newValue == "OpBoundary" )    {
           fSelectedProcessIndex = kBoundary;
    } else if ( newValue == "OpWLS" )         {
           fSelectedProcessIndex = kWLS;
    }
  }
  else if (command == fSetOpProcessUseCmd) {
    fOpticalPhysics->
       Configure(fSelectedProcessIndex,
                 fSetOpProcessUseCmd->GetNewBoolValue(newValue));
  }  
  else if (command == fSetOpProcessVerboseCmd) {
    if ( fSelectedProcessIndex < kNoProcess ) {
       fOpticalPhysics->
          SetProcessVerbose(fSelectedProcessIndex,
                            fSetOpProcessVerboseCmd->GetNewIntValue(newValue));
    } else {
      for ( G4int i=0; i<kNoProcess; i++ ) {
          fOpticalPhysics->
            SetProcessVerbose(i,fSetOpProcessVerboseCmd->GetNewIntValue(newValue));
      }
    }
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
//  else if (command == fSetOpticalSurfaceModelCmd) {
//    if ( newValue == "glisur" ) {
//      fOpticalPhysics
//        ->SetOpticalSurfaceModel(glisur);
//    }    
//    if ( newValue == "unified" ) {
//      fOpticalPhysics
//        ->SetOpticalSurfaceModel(unified);
//    } 
//  }
  else if (command == fSetWLSTimeProfileCmd) {
    if ( newValue == "delta" ) {
      fOpticalPhysics
        ->SetWLSTimeProfile("delta");     }
    if ( newValue == "exponential" ) {
      fOpticalPhysics
        ->SetWLSTimeProfile("exponential");
    }
  } 
  else if (command == fSetTrackSecondariesFirstCmd) {
    fOpticalPhysics->SetTrackSecondariesFirst(fSelectedProcessIndex,
                                              fSetTrackSecondariesFirstCmd->
                                                    GetNewBoolValue(newValue));
  }
}
