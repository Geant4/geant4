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
// $Id: G4EmDNAPhysics.cc,v 1.1 2009-03-25 11:28:07 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4EmDNAPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4DNAGenericIonsManager.hh"

// *** Processes and models

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"

#include "G4DNAIonisation.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNARuddIonisationModel.hh"

#include "G4DNAChargeDecrease.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"

#include "G4DNAChargeIncrease.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::G4EmDNAPhysics(
    G4int ver, const G4String& name)
  : G4VPhysicsConstructor(name), verbose(ver)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::~G4EmDNAPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics::ConstructParticle()
{

// leptons
  G4Electron::Electron();

// baryons
  G4Proton::Proton();

  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics::ConstructProcess()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") {

      // *** Elastic scattering (two alternative models available) ***
      
      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic();
      theDNAElasticProcess->SetModel(new G4DNAScreenedRutherfordElasticModel());
      
      // or alternative model
      // theDNAElasticProcess->SetModel(new G4DNAChampionElasticModel());     
      
      pmanager->AddDiscreteProcess(theDNAElasticProcess);

      // *** Excitation ***
      pmanager->AddDiscreteProcess(new G4DNAExcitation());

      // *** Ionisation ***
      pmanager->AddDiscreteProcess(new G4DNAIonisation());

    } else if ( particleName == "proton" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());

    } else if ( particleName == "hydrogen" ) {
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    } else if ( particleName == "alpha" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());

    } else if ( particleName == "alpha+" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    } else if ( particleName == "helium" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
