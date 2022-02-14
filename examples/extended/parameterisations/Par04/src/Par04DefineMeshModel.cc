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
#include "Par04DefineMeshModel.hh"
#include <G4FastTrack.hh>              // for G4FastTrack
#include <G4Track.hh>                  // for G4Track
#include <G4VFastSimulationModel.hh>   // for G4VFastSimulationModel
#include <G4VUserEventInformation.hh>  // for G4VUserEventInformation
#include "G4Event.hh"                  // for G4Event
#include "G4EventManager.hh"           // for G4EventManager
#include "Par04EventInformation.hh"    // for Par04EventInformation
class G4FastStep;
class G4ParticleDefinition;
class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DefineMeshModel::Par04DefineMeshModel(G4String aModelName, G4Region* aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DefineMeshModel::Par04DefineMeshModel(G4String aModelName)
  : G4VFastSimulationModel(aModelName)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DefineMeshModel::~Par04DefineMeshModel() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04DefineMeshModel::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04DefineMeshModel::ModelTrigger(const G4FastTrack&)
{
  Par04EventInformation* info = dynamic_cast<Par04EventInformation*>(
    G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
  // check if particle direction and position were already set for this event
  if(info != nullptr)
    return !info->GetFlag();
  else
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DefineMeshModel::DoIt(const G4FastTrack& aFastTrack, G4FastStep&)
{
  Par04EventInformation* info = dynamic_cast<Par04EventInformation*>(
    G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
  if(info == nullptr)
  {
    info = new Par04EventInformation();
    G4EventManager::GetEventManager()->GetNonconstCurrentEvent()->SetUserInformation(info);
  }
  info->SetPosition(aFastTrack.GetPrimaryTrack()->GetPosition());
  info->SetDirection(aFastTrack.GetPrimaryTrack()->GetMomentumDirection());
  info->SetFlag(true);
}
