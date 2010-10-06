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
// -------------------------------------------------------------------
// $Id: PhysicsList.cc,v 1.6 2010-10-06 12:16:59 sincerti Exp $
// -------------------------------------------------------------------

#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::PhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1*micrometer;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;
  
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructBosons()
{ 

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructBarions()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4StepLimiter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructEM()
{

// ****************************************************************
// Identical to G4EmStandardPhysics but added G4StepLimiter process
// ****************************************************************

  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

    } else if (particleName == "e-") {

      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1,-3, 3);

    } else if (particleName == "e+") {

      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      pmanager->AddProcess(new G4MuMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,    -1, 4, 4);       

    } else if (particleName == "alpha" ||
               particleName == "He3" ||
               particleName == "GenericIon") {
      // ions with charge >= +2
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2, 2);
     
    } else if (particleName == "proton") {
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4hPairProduction,    -1,-4, 4);

      pmanager->AddProcess(new G4StepLimiter(),-1,-1,5);
            
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  SetCutValue(cutForProton, "proton");
  SetCutValue(cutForProton, "anti_proton");
  
  //SetCutValueForOthers(defaultCutValue);
  
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  // G4Gamma::SetEnergyRange(lowcut,1e5); 
  SetGELowLimit(lowcut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  // G4Electron::SetEnergyRange(lowcut,1e5);
  SetGELowLimit(lowcut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetPositronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "Positron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4cerr << "PhysicsList::SetPositronLowLimit: Not currently able to set Positron LowLimit." << G4endl;
  G4Exception("Positron Low Limit: not implemented in PhysicsList"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetProtonLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "Proton cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;  
  }  

  G4cerr << "PhysicsList::SetProtonLowLimit: Not currently able to set Proton LowLimit." << G4endl;
  G4Exception("Proton Low Limit: not implemented in PhysicsList"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetGEPLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetGEPLowLimit:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  this->SetGELowLimit(lowcut); 

  G4cerr << " SetGEPLowLimit : Uncertain whether setting Positron low limit " << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void PhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetGELowLimit:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  
 
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowcut,1e5);
}
void PhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetElectronCut(G4double val)
{
  //  ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetPositronCut(G4double val)
{
  //  ResetCuts();
  cutForPositron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetProtonCut(G4double val)
{
  //ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
