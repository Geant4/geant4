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
/// \file electromagnetic/TestEm13/src/PhysListEmStandard.cc
/// \brief Implementation of the PhysListEmStandard class
//
<<<<<<< HEAD
//
// $Id: PhysListEmStandard.cc 102356 2017-01-23 16:22:42Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4RayleighScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4KleinNishinaModel.hh"
#include "G4GammaConversion.hh"
<<<<<<< HEAD
=======
#include "G4GammaConversionToMuons.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
<<<<<<< HEAD
=======
#include "G4IonParametrisedLossModel.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#include "G4EmProcessOptions.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
<<<<<<< HEAD
  :  G4VPhysicsConstructor(name)
{ }
=======
   :  G4VPhysicsConstructor(name)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetMinEnergy(10*eV);
  param->SetMaxEnergy(10*TeV);
  param->SetNumberOfBinsPerDecade(10);
  
  param->SetVerbose(0);
  param->Dump();
}
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
<<<<<<< HEAD
{
=======
{ 
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  // Add standard EM Processes

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
<<<<<<< HEAD
      // gamma         
      ////pmanager->AddDiscreteProcess(new G4RayleighScattering);               
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      G4ComptonScattering* cs   = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel());
      pmanager->AddDiscreteProcess(cs);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
            
    } else if (particleName == "e-") {
      //electron
=======

      ////pmanager->AddDiscreteProcess(new G4RayleighScattering);               
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      G4ComptonScattering* compt   = new G4ComptonScattering;
      compt->SetEmModel(new G4KleinNishinaModel());
      pmanager->AddDiscreteProcess(compt);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      pmanager->AddDiscreteProcess(new G4GammaConversionToMuons);  
     
    } else if (particleName == "e-") {

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
      pmanager->AddProcess(new G4eIonisation,        -1,-1,1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,2);
            
    } else if (particleName == "e+") {
<<<<<<< HEAD
      //positron
      pmanager->AddProcess(new G4eIonisation,        -1,-1,1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,2);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,3);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MuIonisation,      -1,-1,1);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,2);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,3);       
=======

      pmanager->AddProcess(new G4eIonisation,        -1,-1,1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,2);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,3);
                  
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      pmanager->AddProcess(new G4MuIonisation,      -1,-1,1);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,2);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,3);       
                   
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {

      pmanager->AddProcess(new G4hIonisation,       -1,-1,1);
      pmanager->AddProcess(new G4hBremsstrahlung,   -1,-1,2);      
      pmanager->AddProcess(new G4hPairProduction,   -1,-1,3);        
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
     
    } else if( particleName == "alpha" || particleName == "GenericIon" ) { 
      pmanager->AddProcess(new G4ionIonisation,     -1,-1,1);

<<<<<<< HEAD
=======
      pmanager->AddProcess(new G4ionIonisation,     -1,-1,1);
            
    } else if( particleName == "GenericIon" ) {
 
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      pmanager->AddProcess(ionIoni,                 -1,-1,1);
      
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
<<<<<<< HEAD
      pmanager->AddProcess(new G4hIonisation,       -1,-1,1);
=======
      pmanager->AddProcess(new G4hIonisation,       -1,-1,1);      
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    }
  }
    
  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  de->SetFluo(true);
  de->SetAuger(false);  
  de->SetPIXE(false);  
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

