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

#include "test31PreCompound.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Proton.hh"   
#include "G4Neutron.hh"   
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31PreCompound::test31PreCompound(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31PreCompound::~test31PreCompound()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PreCompound::ConstructParticle()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PreCompound::ConstructProcess()
{
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (particle == G4Proton::Proton()) {

       G4ProtonInelasticProcess* theInelasticProcess = 
                 new G4ProtonInelasticProcess("inelastic");

       theInelasticProcess->GetCrossSectionDataStore()->AddDataSet(
                 new G4ProtonInelasticCrossSection());
    
       theInelasticProcess->RegisterMe(
                 new G4PreCompoundModel(new G4ExcitationHandler()));
       pmanager->AddDiscreteProcess(theInelasticProcess);
       G4cout << "PreCompaund model for protons are registered" << G4endl;

    } else if (particle == G4Neutron::Neutron()) {

       G4NeutronInelasticProcess* theInelasticProcess = 
                 new G4NeutronInelasticProcess("inelastic");

       theInelasticProcess->GetCrossSectionDataStore()->AddDataSet(
                 new G4NeutronInelasticCrossSection());
    
       theInelasticProcess->RegisterMe(
                 new G4PreCompoundModel(new G4ExcitationHandler()));

       pmanager->AddDiscreteProcess(theInelasticProcess);

       G4HadronFissionProcess* theFissionProcess =
                 new G4HadronFissionProcess();
       G4LFission* theFissionModel = new G4LFission;
       theFissionProcess->RegisterMe(theFissionModel);
       pmanager->AddDiscreteProcess(theFissionProcess);

         // capture
       G4HadronCaptureProcess* theCaptureProcess =
                 new G4HadronCaptureProcess;
       G4LCapture* theCaptureModel = new G4LCapture;
       theCaptureProcess->RegisterMe(theCaptureModel);
       pmanager->AddDiscreteProcess(theCaptureProcess);
       G4cout << "PreCompaund model for neutrons are registered" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

