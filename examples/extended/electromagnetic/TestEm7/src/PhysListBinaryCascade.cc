//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// 

#include "PhysListBinaryCascade.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"

#include "G4BinaryCascade.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListBinaryCascade::PhysListBinaryCascade(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListBinaryCascade::~PhysListBinaryCascade()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListBinaryCascade::ConstructProcess()
{

  // Binary Cascade
  G4BinaryCascade * theBC = new G4BinaryCascade();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (particleName == "proton") {

      G4ProtonInelasticProcess* theInelasticProcess =
                                    new G4ProtonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBC);
      pmanager->AddDiscreteProcess(theInelasticProcess);
          G4cout << "Add hadronic physics for proton" << G4endl;

    } else if (particleName == "neutron") {

      G4NeutronInelasticProcess* theInelasticProcess =
                                    new G4NeutronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBC);
      pmanager->AddDiscreteProcess(theInelasticProcess);
          // fission
      G4HadronFissionProcess* theFissionProcess =
                                    new G4HadronFissionProcess;
      G4LFission* theFissionModel = new G4LFission;
      theFissionProcess->RegisterMe(theFissionModel);
      pmanager->AddDiscreteProcess(theFissionProcess);
         // capture
      G4HadronCaptureProcess* theCaptureProcess =
                                    new G4HadronCaptureProcess;
      G4LCapture* theCaptureModel = new G4LCapture;
      theCaptureProcess->RegisterMe(theCaptureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);
          G4cout << "Add hadronic physics for neutron" << G4endl;
      /*
    if (particle == G4IonC12::IonC12()) {
      pManager->AddDiscreteProcess(&theIonProcess);
      G4cout << "### IonSpecial process are registered for "
             << particle->GetParticleName()
             << G4endl;
    }
      */
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

