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
/// \file OpNoviceSteppingAction.cc
/// \brief Implementation of the OpNoviceSteppingAction class

#include "OpNoviceSteppingAction.hh"
#include "OpNoviceRun.hh"
#include "G4Event.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceSteppingAction::OpNoviceSteppingAction(OpNoviceEventAction* event)
  : G4UserSteppingAction()
  , fEventAction(event)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceSteppingAction::UserSteppingAction(const G4Step* step)
{
  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  const G4ParticleDefinition* particleDef =
    step->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  if(particleDef == opticalphoton)
  {
    G4StepPoint* endPoint = step->GetPostStepPoint();
    const G4VProcess* pds = endPoint->GetProcessDefinedStep();
    G4String procname     = pds->GetProcessName();
    if(procname == "OpRayleigh")
      fEventAction->AddRayleigh();
    else if(procname == "OpAbsorption")
      fEventAction->AddAbsorption();
    else if(procname == "OpMieHG")
      fEventAction->AddMie();

    // for boundary scattering, process name in 'transportation'.
    // Need to check differently:
    if(endPoint->GetStepStatus() == fGeomBoundary)
    {
      G4OpBoundaryProcessStatus theStatus = Undefined;
      G4ProcessManager* opManager         = opticalphoton->GetProcessManager();
      G4int n_proc = opManager->GetPostStepProcessVector(typeDoIt)->entries();
      G4ProcessVector* postStepDoItVector =
        opManager->GetPostStepProcessVector(typeDoIt);
      for(G4int i = 0; i < n_proc; ++i)
      {
        G4VProcess* currentProcess = (*postStepDoItVector)[i];

        auto opProc = dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
        if(opProc)
          theStatus = opProc->GetStatus();
      }
      if(theStatus != Undefined && theStatus != NotAtBoundary &&
         theStatus != StepTooSmall)
      {
        fEventAction->AddBoundary();
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
