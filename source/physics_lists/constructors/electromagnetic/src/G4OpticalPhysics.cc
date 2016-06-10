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

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"

#include "G4OpBoundaryProcess.hh"

#include "G4OpWLS.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4OpticalPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysics::G4OpticalPhysics(G4int verbose, const G4String& name)
  : G4VPhysicsConstructor(name),

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
    fProcessUse.push_back(true);
    fProcessVerbose.push_back(verbose);
    fProcessTrackSecondariesFirst.push_back(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4OpticalPhysics::~G4OpticalPhysics()
{
  delete fMessenger;
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

#include "G4Threading.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void G4OpticalPhysics::ConstructProcess()
{
// Construct optical processes.

  if(verboseLevel>0)
         G4cout <<"G4OpticalPhysics:: Add Optical Physics Processes"<< G4endl;

  // A vector of optical processes
  std::vector<G4VProcess*> OpProcesses;

  for ( G4int i=0; i<kNoProcess; i++ ) OpProcesses.push_back(NULL);

  // Add Optical Processes

  G4OpAbsorption* OpAbsorptionProcess  = new G4OpAbsorption();
  G4OpRayleigh* OpRayleighScatteringProcess = new G4OpRayleigh();
  G4OpMieHG* OpMieHGScatteringProcess = new G4OpMieHG();

  OpProcesses[kAbsorption] = OpAbsorptionProcess;
  OpProcesses[kRayleigh] = OpRayleighScatteringProcess;
  OpProcesses[kMieHG] = OpMieHGScatteringProcess;

  G4OpBoundaryProcess* OpBoundaryProcess = new G4OpBoundaryProcess();

  OpProcesses[kBoundary] = OpBoundaryProcess;

  G4OpWLS* OpWLSProcess = new G4OpWLS();

  OpProcesses[kWLS] = OpWLSProcess;

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
         pManager->AddDiscreteProcess(OpProcesses[i]);
      }
  }

  G4Scintillation* ScintillationProcess = new G4Scintillation();

  OpProcesses[kScintillation] = ScintillationProcess;

  G4Cerenkov* CerenkovProcess = new G4Cerenkov();
  OpProcesses[kCerenkov] = CerenkovProcess;

  // static parameters - set only for masther thread
  if(!G4Threading::IsWorkerThread())
  {
    G4OpWLS::UseTimeProfile(fProfile);

    G4Scintillation::SetScintillationYieldFactor(fYieldFactor);
    G4Scintillation::SetScintillationExcitationRatio(fExcitationRatio);
    G4Scintillation::SetFiniteRiseTime(fFiniteRiseTime);
    G4Scintillation::SetScintillationByParticleType(fScintillationByParticleType);
    G4Scintillation::SetTrackSecondariesFirst(fProcessTrackSecondariesFirst[kScintillation]);

    // Use Birks Correction in the Scintillation process
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    G4Scintillation::AddSaturation(emSaturation);

    G4Cerenkov::SetMaxNumPhotonsPerStep(fMaxNumPhotons);
    G4Cerenkov::SetMaxBetaChangePerStep(fMaxBetaChange);
    G4Cerenkov::SetTrackSecondariesFirst(fProcessTrackSecondariesFirst[kCerenkov]);
  }

  aParticleIterator->reset();

  while( (*aParticleIterator)() ){

    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
       std::ostringstream o;
       o << "Particle " << particleName << "without a Process Manager";
       G4Exception("G4OpticalPhysics::ConstructProcess()","",
                    FatalException,o.str().c_str());
       return;                 // else coverity complains for pManager use below	    
    }

    if( CerenkovProcess->IsApplicable(*particle) &&
        fProcessUse[kCerenkov] ) {
          pManager->AddProcess(CerenkovProcess);
          pManager->SetProcessOrdering(CerenkovProcess,idxPostStep);
    }
    if( ScintillationProcess->IsApplicable(*particle) &&
        fProcessUse[kScintillation] ){
          pManager->AddProcess(ScintillationProcess);
          pManager->SetProcessOrderingToLast(ScintillationProcess,idxAtRest);
          pManager->SetProcessOrderingToLast(ScintillationProcess,idxPostStep);
    }

  }

  // Add verbose
  for ( G4int i=0; i<kNoProcess; i++ ) {
    if ( fProcessUse[i] ) OpProcesses[i]->SetVerboseLevel(fProcessVerbose[i]);
  }

  if (verboseLevel > 1) PrintStatistics();
  if (verboseLevel > 0)
    G4cout << "### " << namePhysics << " physics constructed." << G4endl;
}

void G4OpticalPhysics::SetScintillationYieldFactor(G4double yieldFactor)
{
/// Set the scintillation yield factor

  fYieldFactor = yieldFactor;
  G4Scintillation::SetScintillationYieldFactor(yieldFactor);
}

void G4OpticalPhysics::SetScintillationExcitationRatio(G4double excitationRatio)
{
/// Set the scintillation excitation ratio

  fExcitationRatio = excitationRatio;
  G4Scintillation::SetScintillationExcitationRatio(excitationRatio);
}

void G4OpticalPhysics::SetMaxNumPhotonsPerStep(G4int maxNumPhotons)
{
/// Limit step to the specified maximum number of Cherenkov photons

  fMaxNumPhotons = maxNumPhotons;
  G4Cerenkov::SetMaxNumPhotonsPerStep(maxNumPhotons);
}

void G4OpticalPhysics::SetMaxBetaChangePerStep(G4double maxBetaChange)
{
/// Limit step to the specified maximum change of beta of the parent particle

  fMaxBetaChange = maxBetaChange;
  G4Cerenkov::SetMaxBetaChangePerStep(maxBetaChange);
}

void G4OpticalPhysics::SetWLSTimeProfile(G4String profile)
{
/// Set the WLS time profile (delta or exponential)

  fProfile = profile;
  G4OpWLS::UseTimeProfile(profile);
}

void G4OpticalPhysics::AddScintillationSaturation(G4EmSaturation* saturation)
{
/// Adds Birks Saturation to the G4Scintillation Process
  G4Scintillation::AddSaturation(saturation);
}

void G4OpticalPhysics::
             SetScintillationByParticleType(G4bool scintillationByParticleType)
{
  fScintillationByParticleType = scintillationByParticleType;
  G4Scintillation::SetScintillationByParticleType(scintillationByParticleType);
}

void G4OpticalPhysics::SetTrackSecondariesFirst(G4OpticalProcessIndex index,
                                                G4bool trackSecondariesFirst)
{
  if ( index >= kNoProcess ) return;
  if ( fProcessTrackSecondariesFirst[index] == trackSecondariesFirst ) return;
  fProcessTrackSecondariesFirst[index] = trackSecondariesFirst;

  if ( index == kCerenkov )
     G4Cerenkov::SetTrackSecondariesFirst(trackSecondariesFirst);
  if ( index == kScintillation)
     G4Scintillation::SetTrackSecondariesFirst(trackSecondariesFirst);
}

void G4OpticalPhysics::SetFiniteRiseTime(G4bool finiteRiseTime)
{
  fFiniteRiseTime = finiteRiseTime;
  G4Scintillation::SetFiniteRiseTime(finiteRiseTime);
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
                                         G4int verbose)
{
  // Set new verbose level to a selected process

  if ( index >= kNoProcess ) return;
  if ( fProcessVerbose[index] == verbose ) return;
  fProcessVerbose[index] = verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
