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
/// \file field/field04/src/F04OpticalPhysics.cc
/// \brief Implementation of the F04OpticalPhysics class
//
//
#include "globals.hh"

#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "F04OpticalPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04OpticalPhysics::F04OpticalPhysics()
    : G4VPhysicsConstructor("Optical") {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04OpticalPhysics::~F04OpticalPhysics() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

void F04OpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ProcessManager.hh"

void F04OpticalPhysics::ConstructProcess()
{
    G4cout << "F04OpticalPhysics:: Add Optical Physics Processes"
           << G4endl;

  G4Scintillation* theScintProcess = new G4Scintillation();
  G4Cerenkov* theCerenkovProcess= new G4Cerenkov();

  G4OpAbsorption* theAbsorptionProcess= new G4OpAbsorption();
  G4OpRayleigh* theRayleighScattering = new G4OpRayleigh();
  G4OpMieHG* theMieHGScatteringProcess = new G4OpMieHG();
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  G4OpWLS* theWLSProcess=new G4OpWLS();

  G4ProcessManager* pManager =
                G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (!pManager) {
     std::ostringstream o;
     o << "Optical Photon without a Process Manager";
     G4Exception("F04OpticalPhysics::ConstructProcess()","",
                  FatalException,o.str().c_str());
  }

  pManager->AddDiscreteProcess(theAbsorptionProcess);
  pManager->AddDiscreteProcess(theRayleighScattering);
  pManager->AddDiscreteProcess(theMieHGScatteringProcess);

  pManager->AddDiscreteProcess(theBoundaryProcess);

  theWLSProcess->UseTimeProfile("delta");
//  theWLSProcess->UseTimeProfile("exponential");

  pManager->AddDiscreteProcess(theWLSProcess);

  theScintProcess->SetScintillationYieldFactor(1.);
  theScintProcess->SetScintillationExcitationRatio(0.0);
  theScintProcess->SetTrackSecondariesFirst(true);

  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
       std::ostringstream o;
       o << "Particle " << particleName << "without a Process Manager";
       G4Exception("F04OpticalPhysics::ConstructProcess()","",
                    FatalException,o.str().c_str());
    }

    if(theCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(theCerenkovProcess);
      pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    if(theScintProcess->IsApplicable(*particle)){
      pManager->AddProcess(theScintProcess);
      pManager->SetProcessOrderingToLast(theScintProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(theScintProcess,idxPostStep);
    }

  }
}
