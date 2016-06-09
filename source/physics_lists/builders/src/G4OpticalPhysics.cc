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
//---------------------------------------------------------------------------
//
// ClassName:   G4OpticalPhysics
//
// Author:      P.Gumplinger 30.09.2009
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4OpticalPhysics.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpticalPhysics::G4OpticalPhysics(G4int verbose) 
  : G4VPhysicsConstructor("Optical")
,
    wasActivated(false),

    fScintillationProcess(0),
    fCerenkovProcess(0),
    fOpWLSProcess(0),
    fOpAbsorptionProcess(0),
    fOpRayleighScatteringProcess(0),
    fOpMieHGScatteringProcess(0),
    fOpBoundaryProcess(0),
    fMaxNumPhotons(100),
    fMaxBetaChange(10.0),
    fYieldFactor(1.),
    fExcitationRatio(0.0),
    fSurfaceModel(unified),
    fProfile("delta"),
    fTrackSecondariesFirst(true),
    fScintillationByParticleType(false)
{
  verboseLevel = verbose;
  fMessenger = new G4OpticalPhysicsMessenger(this);
}

G4OpticalPhysics::G4OpticalPhysics(G4int verbose, const G4String& name)
  : G4VPhysicsConstructor(name),
    wasActivated(false),

    fScintillationProcess(0),
    fCerenkovProcess(0),
    fOpWLSProcess(0),
    fOpAbsorptionProcess(0),
    fOpRayleighScatteringProcess(0),
    fOpMieHGScatteringProcess(0),
    fOpBoundaryProcess(0),
    fMaxNumPhotons(100),
    fMaxBetaChange(10.0),
    fYieldFactor(1.),
    fExcitationRatio(0.0),
    fSurfaceModel(unified),
    fProfile("delta"),
    fTrackSecondariesFirst(true),
    fScintillationByParticleType(false)
{
  verboseLevel = verbose;
  fMessenger = new G4OpticalPhysicsMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysics::~G4OpticalPhysics()
{
  delete fMessenger;

  if (wasActivated) {

     delete fScintillationProcess;
     delete fCerenkovProcess;
     delete fOpWLSProcess;

     delete fOpAbsorptionProcess;
     delete fOpRayleighScatteringProcess;
     delete fOpMieHGScatteringProcess;
     delete fOpBoundaryProcess;

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4OpticalPhoton.hh"

void G4OpticalPhysics::ConstructParticle()
{
/// Instantiate particles.

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void G4OpticalPhysics::ConstructProcess()
{
// Construct optical processes.

  if (wasActivated) return;

  if(verboseLevel>0)
         G4cout <<"G4OpticalPhysics:: Add Optical Physics Processes"<< G4endl;

  // Add Optical Processes

  fOpAbsorptionProcess  = new G4OpAbsorption();
  fOpRayleighScatteringProcess = new G4OpRayleigh();
  fOpMieHGScatteringProcess = new G4OpMieHG();

  fOpBoundaryProcess    = new G4OpBoundaryProcess();
  fOpBoundaryProcess->SetModel(fSurfaceModel);

  fOpWLSProcess         = new G4OpWLS();
  fOpWLSProcess->UseTimeProfile(fProfile);

  G4ProcessManager * pManager = 0;
  pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (!pManager) {
     std::ostringstream o;
     o << "Optical Photon without a Process Manager";
     G4Exception("G4OpticalPhysics::ConstructProcess()","",
                  FatalException,o.str().c_str());
     return;
  }

  pManager->AddDiscreteProcess(fOpAbsorptionProcess);
  pManager->AddDiscreteProcess(fOpRayleighScatteringProcess);
  pManager->AddDiscreteProcess(fOpMieHGScatteringProcess);
  pManager->AddDiscreteProcess(fOpBoundaryProcess);
  pManager->AddDiscreteProcess(fOpWLSProcess);

  fScintillationProcess       = new G4Scintillation();
  fScintillationProcess->SetScintillationYieldFactor(fYieldFactor);
  fScintillationProcess->SetScintillationExcitationRatio(fExcitationRatio);
  fScintillationProcess->SetTrackSecondariesFirst(fTrackSecondariesFirst);
  fScintillationProcess->SetScintillationByParticleType(fScintillationByParticleType);

  // Use Birks Correction in the Scintillation process

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  fScintillationProcess->AddSaturation(emSaturation);

  fCerenkovProcess    = new G4Cerenkov();
  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotons);
  fCerenkovProcess->SetMaxBetaChangePerStep(fMaxBetaChange);
  fCerenkovProcess->SetTrackSecondariesFirst(fTrackSecondariesFirst);

  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
       std::ostringstream o;
       o << "Particle " << particleName << "without a Process Manager";
       G4Exception("G4OpticalPhysics::ConstructProcess()","",
                    FatalException,o.str().c_str());
    }

    if(fCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(fCerenkovProcess);
      pManager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
    if(fScintillationProcess->IsApplicable(*particle)){
      pManager->AddProcess(fScintillationProcess);
      pManager->SetProcessOrderingToLast(fScintillationProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(fScintillationProcess,idxPostStep);
    }

  }

  wasActivated = true;

}

void G4OpticalPhysics::SetScintillationYieldFactor(G4double yieldFactor)
{
/// Set the scintillation yield factor

  fYieldFactor = yieldFactor;

  if(fScintillationProcess)
    fScintillationProcess->SetScintillationYieldFactor(yieldFactor);
}

void G4OpticalPhysics::SetScintillationExcitationRatio(G4double excitationRatio)
{
/// Set the scintillation excitation ratio

  fExcitationRatio = excitationRatio;

  if(fScintillationProcess)
    fScintillationProcess->SetScintillationExcitationRatio(excitationRatio);
}

void G4OpticalPhysics::SetMaxNumPhotonsPerStep(G4int maxNumPhotons)
{
/// Limit step to the specified maximum number of Cherenkov photons

  fMaxNumPhotons = maxNumPhotons;

  if(fCerenkovProcess)
    fCerenkovProcess->SetMaxNumPhotonsPerStep(maxNumPhotons);
}

void G4OpticalPhysics::SetMaxBetaChangePerStep(G4double maxBetaChange)
{
/// Limit step to the specified maximum change of beta of the parent particle

  fMaxBetaChange = maxBetaChange;

  if(fCerenkovProcess)
    fCerenkovProcess->SetMaxBetaChangePerStep(maxBetaChange);
}

void G4OpticalPhysics::SetOpticalSurfaceModel(G4OpticalSurfaceModel model)
{
/// Set optical surface model (glisur or unified)

  fSurfaceModel = model;

  if(fOpBoundaryProcess)
    fOpBoundaryProcess->SetModel(model); 
}

void G4OpticalPhysics::SetWLSTimeProfile(G4String profile)
{
/// Set the WLS time profile (delta or exponential)

  fProfile = profile;

  if(fOpWLSProcess)
    fOpWLSProcess->UseTimeProfile(fProfile);
}

void G4OpticalPhysics::AddScintillationSaturation(G4EmSaturation* saturation)
{
/// Adds Birks Saturation to the G4Scintillation Process

  if(fScintillationProcess)
    fScintillationProcess->AddSaturation(saturation);
}

void G4OpticalPhysics::SetScintillationByParticleType(G4bool scintillationByParticleType)
{
  fScintillationByParticleType = scintillationByParticleType;

  if (fScintillationProcess)
     fScintillationProcess->SetScintillationByParticleType(scintillationByParticleType);
}

void G4OpticalPhysics::SetTrackSecondariesFirst(G4bool trackSecondariesFirst)
{
  fTrackSecondariesFirst = trackSecondariesFirst;

  if ( fCerenkovProcess )
     fCerenkovProcess->SetTrackSecondariesFirst(trackSecondariesFirst);

  if ( fScintillationProcess )
     fScintillationProcess->SetTrackSecondariesFirst(trackSecondariesFirst);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
