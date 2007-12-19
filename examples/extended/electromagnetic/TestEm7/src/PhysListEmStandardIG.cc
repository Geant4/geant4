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
// $Id: PhysListEmStandardIG.cc,v 1.6 2007-12-19 19:41:58 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmStandardIG.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4ScreenedNuclearRecoil.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4CoulombScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardIG::PhysListEmStandardIG(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardIG::~PhysListEmStandardIG()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandardIG::ConstructProcess()
{
  // Add standard EM Processes

  G4ScreenedNuclearRecoil* nucr = new G4ScreenedNuclearRecoil();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
            
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3,3);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 4,4);
             
    } else if (particleName == "alpha" || particleName == "He3") {
      G4ionIonisation* ion = new G4ionIonisation();
      ion->ActivateNuclearStopping(false);
      pmanager->AddProcess(ion, -1, 1, 1);
      pmanager->AddDiscreteProcess(nucr);      

    } else if (particleName == "GenericIon" ) { 
      G4ionIonisation* ion = new G4ionIonisation();
      ion->ActivateNuclearStopping(false);
      pmanager->AddProcess(ion, -1, 1, 1);
      pmanager->AddDiscreteProcess(nucr);      

    } else if (particleName == "proton" ||
	       particleName == "deuteron" ||
               particleName == "triton") { 
      pmanager->AddProcess(new G4hIonisation,      -1, 1,1);
      pmanager->AddDiscreteProcess(nucr);      
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4hIonisation,        -1, 2,2);      
    }
  }
  G4EmProcessOptions opt;
  opt.SetStepFunction(0.2, 100*um);
  opt.SetMinEnergy(0.1*keV);
  opt.SetMaxEnergy(100.*GeV);
  opt.SetDEDXBinning(360);
  opt.SetLambdaBinning(360);
  opt.SetLinearLossLimit(1.e-6);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

