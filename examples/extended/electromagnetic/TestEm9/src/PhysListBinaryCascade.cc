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
// $Id: PhysListBinaryCascade.cc,v 1.2 2003-10-13 16:31:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "PhysListBinaryCascade.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

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
  G4BinaryCascade * theBC = 0;

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (particleName == "proton") {

      theBC = new G4BinaryCascade();
      theIPproton.RegisterMe(theBC);
      G4CrossSectionDataStore * thePStore = theIPproton.GetCrossSectionDataStore();
      thePStore->AddDataSet(&theXSecProton);
      pmanager->AddDiscreteProcess(&theIPproton);

    } else if (particleName == "neutron") {

      theBC = new G4BinaryCascade();
      theIPneutron.RegisterMe(theBC);
      G4CrossSectionDataStore * thePStore = theIPneutron.GetCrossSectionDataStore();
      thePStore->AddDataSet(&theXSecNeutron);
      pmanager->AddDiscreteProcess(&theIPneutron);
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

