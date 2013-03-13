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
// Modified:    P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
//----------------------------------------------------------------------------
//

#include "G4OpticalPhysics.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"


// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4OpticalPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysics::G4OpticalPhysics(G4int verbose, const G4String& name)
  : G4VPhysicsConstructor(name),

    wasActivated(false),

    fScintillationProcess(NULL),
    fCerenkovProcess(NULL),
    fOpWLSProcess(NULL),
    fOpAbsorptionProcess(NULL),
    fOpRayleighScatteringProcess(NULL),
    fOpMieHGScatteringProcess(NULL),
    fOpBoundaryProcess(NULL),
    fMaxNumPhotons(100),
    fMaxBetaChange(10.0),
    fYieldFactor(1.),
    fExcitationRatio(0.0),
    fProfile("delta"),
    fFiniteRiseTime(false),
    fScintillationByParticleType(false)
{
  verboseLevel = verbose;
  fMessenger = new G4OpticalPhysicsMessenger(this);

  for ( G4int i=0; i<kNoProcess; i++ ) {
    fProcesses.push_back(NULL);
    fProcessUse.push_back(true);
    fProcessVerbose.push_back(verbose);
    fProcessTrackSecondariesFirst.push_back(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysics::~G4OpticalPhysics()
{
  delete fMessenger;

  if (wasActivated) {

     if (fScintillationProcess) delete fScintillationProcess;
     if (fCerenkovProcess) delete fCerenkovProcess;
     if (fOpWLSProcess) delete fOpWLSProcess;

     if (fOpAbsorptionProcess) delete fOpAbsorptionProcess;
     if (fOpRayleighScatteringProcess) delete fOpRayleighScatteringProcess;
     if (fOpMieHGScatteringProcess) delete fOpMieHGScatteringProcess;
     if (fOpBoundaryProcess) delete fOpBoundaryProcess;

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4OpticalPhysics::PrintStatistics() const
{
// Print all processes activation and their parameters

  for ( G4int i=0; i<kNoProcess; i++ ) {
    G4cout << "  " << G4OpticalProcessName(i) << " process:  ";
    if ( ! fProcessUse[i] ) {
      G4cout << "not used" << G4endl;
    }
    else {
      G4cout << "used" << G4endl;
      if ( i == kCerenkov ) {
        G4cout << "    Max number of photons per step: " << fMaxNumPhotons << G4endl;
        G4cout << "    Max beta change per step:       " << fMaxBetaChange << G4endl;
        if ( fProcessTrackSecondariesFirst[kCerenkov] ) G4cout << "  Track secondaries first:  activated" << G4endl;
      }
      if ( i == kScintillation ) {
        if (fScintillationByParticleType)
        G4cout << "    Scintillation by Particle Type:  activated " << G4endl;
        G4cout << "    Yield factor: "  << fYieldFactor << G4endl;
        G4cout << "    ExcitationRatio: " << fExcitationRatio << G4endl;
        if ( fProcessTrackSecondariesFirst[kScintillation] ) G4cout << "  Track secondaries first:  activated" << G4endl;
      }
      if ( i == kWLS ) {
        G4cout << "     WLS process time profile: " << fProfile << G4endl;
      }
    }
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

  fProcesses[kAbsorption] = fOpAbsorptionProcess  = new G4OpAbsorption();
  fProcesses[kRayleigh] = fOpRayleighScatteringProcess = new G4OpRayleigh();
  fProcesses[kMieHG] = fOpMieHGScatteringProcess = new G4OpMieHG();

  fProcesses[kBoundary] = fOpBoundaryProcess = new G4OpBoundaryProcess();

  fProcesses[kWLS] = fOpWLSProcess = new G4OpWLS();
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

  for ( G4int i=kAbsorption; i<=kWLS; i++ ) {
      if ( fProcessUse[i] ) {
         pManager->AddDiscreteProcess(fProcesses[i]);
      }
  }

  fProcesses[kScintillation] = fScintillationProcess = new G4Scintillation();
  fScintillationProcess->SetScintillationYieldFactor(fYieldFactor);
  fScintillationProcess->SetScintillationExcitationRatio(fExcitationRatio);
  fScintillationProcess->SetFiniteRiseTime(fFiniteRiseTime);
  fScintillationProcess->SetScintillationByParticleType(fScintillationByParticleType);
  fScintillationProcess->
       SetTrackSecondariesFirst(fProcessTrackSecondariesFirst[kScintillation]);

  // Use Birks Correction in the Scintillation process

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  fScintillationProcess->AddSaturation(emSaturation);

  fProcesses[kCerenkov] = fCerenkovProcess = new G4Cerenkov();
  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotons);
  fCerenkovProcess->SetMaxBetaChangePerStep(fMaxBetaChange);
  fCerenkovProcess->
       SetTrackSecondariesFirst(fProcessTrackSecondariesFirst[kCerenkov]);

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
       return;                 // else coverity complains for pManager use below	    
    }

    if( fCerenkovProcess->IsApplicable(*particle) &&
        fProcessUse[kCerenkov] ) {
          pManager->AddProcess(fCerenkovProcess);
          pManager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
    if( fScintillationProcess->IsApplicable(*particle) &&
        fProcessUse[kScintillation] ){
          pManager->AddProcess(fScintillationProcess);
          pManager->SetProcessOrderingToLast(fScintillationProcess,idxAtRest);
          pManager->SetProcessOrderingToLast(fScintillationProcess,idxPostStep);
    }

  }

  // Add verbose
  for ( G4int i=0; i<kNoProcess; i++ ) {
    fProcesses[i]->SetVerboseLevel(fProcessVerbose[i]);
  }

  if (verboseLevel > 1) PrintStatistics();
  if (verboseLevel > 0)
    G4cout << "### " << namePhysics << " physics constructed." << G4endl;

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

void G4OpticalPhysics::SetTrackSecondariesFirst(G4OpticalProcessIndex index,
                                                G4bool trackSecondariesFirst)
{
  if ( index >= kNoProcess ) return;
  if ( fProcessTrackSecondariesFirst[index] == trackSecondariesFirst ) return;
  fProcessTrackSecondariesFirst[index] = trackSecondariesFirst;

  if(fCerenkovProcess && index == kCerenkov )
    fCerenkovProcess->SetTrackSecondariesFirst(trackSecondariesFirst);

  if(fScintillationProcess && index == kScintillation)
    fScintillationProcess->SetTrackSecondariesFirst(trackSecondariesFirst);
}

void G4OpticalPhysics::SetFiniteRiseTime(G4bool finiteRiseTime)
{
  fFiniteRiseTime = finiteRiseTime;
  if(fScintillationProcess)
    fScintillationProcess->SetFiniteRiseTime(finiteRiseTime);
} 

void G4OpticalPhysics::Configure(G4OpticalProcessIndex index, G4bool isUse)
{
  // Configure the physics constructor to use/not use a selected process.
  // This method can only be called in PreInit> phase (before execution of
  // ConstructProcess). The process is not added to particle's process manager
  // and so it cannot be re-activated later in Idle> phase with the command
  // /process/activate.

  if ( index >= kNoProcess ) return;
  if ( fProcessUse[index] == isUse ) return;
  fProcessUse[index] = isUse;
}

void G4OpticalPhysics::SetProcessVerbose(G4int index,
                                         G4int inputVerboseLevel)
{
  // Set new verbose level to a selected process

  if ( index >= kNoProcess ) return;
  if ( fProcessVerbose[index] == inputVerboseLevel ) return;

  fProcessVerbose[index] = inputVerboseLevel;

  if ( fProcesses[index] ) fProcesses[index]->SetVerboseLevel(inputVerboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
