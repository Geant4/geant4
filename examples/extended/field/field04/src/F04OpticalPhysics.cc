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
// $Id: F04OpticalPhysics.cc 79251 2014-02-20 16:16:23Z gcosmo $
//
/// \file field/field04/src/F04OpticalPhysics.cc
/// \brief Implementation of the F04OpticalPhysics class
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
#include "G4Threading.hh"

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

  pManager->AddDiscreteProcess(theWLSProcess);

  if(!G4Threading::IsWorkerThread())
  {
    G4OpWLS::UseTimeProfile("delta");
//  G4OpWLS::UseTimeProfile("exponential");

    G4Scintillation::SetScintillationYieldFactor(1.);
    G4Scintillation::SetScintillationExcitationRatio(0.0);
    G4Scintillation::SetTrackSecondariesFirst(true);
  }

  aParticleIterator->reset();

  while( (*aParticleIterator)() ){

    G4ParticleDefinition* particle = aParticleIterator->value();
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
