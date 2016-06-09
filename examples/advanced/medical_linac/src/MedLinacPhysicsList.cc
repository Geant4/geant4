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
// $Id: MedLinacPhysicsList.cc,v 1.5 2006/06/29 16:04:39 gunter Exp $
//
//
// Code developed by: M. Piergentili
//
#include "globals.hh"
#include "MedLinacPhysicsList.hh"
#include "MedLinacPhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4UnitsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4StepLimiter.hh"
MedLinacPhysicsList::MedLinacPhysicsList():  G4VUserPhysicsList()
{   
  physicsListMessenger = new MedLinacPhysicsListMessenger(this);
  //defaultCutValue = defaultCut;
  //cutForGamma     = defaultCutValue;
  //cutForElectron  = defaultCutValue;
  //cutForPositron  = defaultCutValue;
   SetVerboseLevel(6);
}

MedLinacPhysicsList::~MedLinacPhysicsList()
{;}

void MedLinacPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
}
void MedLinacPhysicsList::ConstructBosons()
{
  // gamma
  G4Gamma::GammaDefinition();
}

void MedLinacPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void MedLinacPhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
  ConstructEM();
}

//---------------------------------------------------------------------

//#include "G4ComptonScattering.hh"
//#include "G4GammaConversion.hh"
//#include "G4PhotoElectricEffect.hh"

#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyRayleigh.hh" 

//#include "G4PenelopePhotoElectric.hh"
//#include "G4PenelopeCompton.hh"  
//#include "G4PenelopeGammaConversion.hh"
//#include "G4PenelopeRayleigh.hh" 

#include "G4MultipleScattering52.hh"

#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

//#include "G4PenelopeIonisation.hh" 
//#include "G4PenelopeBremsstrahlung.hh" 

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"


//---------------------------------------------------------------------

void MedLinacPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma  
 
      lowePhot = new  G4LowEnergyPhotoElectric("LowEnPhotoElec");
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager->AddDiscreteProcess(lowePhot);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
      pmanager -> AddProcess(new G4StepLimiter(), -1, -1, 3);

      //pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      //pmanager->AddDiscreteProcess(new G4ComptonScattering);
      //pmanager->AddDiscreteProcess(new G4GammaConversion);

      //pmanager->AddDiscreteProcess(new G4PenelopePhotoElectric);
      //pmanager->AddDiscreteProcess(new G4PenelopeCompton);
      //pmanager->AddDiscreteProcess(new G4PenelopeGammaConversion);      
      //pmanager->AddDiscreteProcess(new G4PenelopeRayleigh);

      
    } else if (particleName == "e-") {
      //electron
     
      loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
      loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");
    
      pmanager->AddProcess(new G4MultipleScattering52, -1, 1,1);
      pmanager->AddProcess(loweIon,     -1, 2,2);
      pmanager->AddProcess(loweBrem,    -1,-1,3);     
      pmanager -> AddProcess(new G4StepLimiter(), -1, -1, 3);
      
      //pmanager->AddProcess(new G4PenelopeIonisation,       -1, 2,2);
      //pmanager->AddProcess(new G4PenelopeBremsstrahlung,   -1,-1,3);
      
      //pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      //pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);    
  
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering52,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      pmanager -> AddProcess(new G4StepLimiter(), -1, -1, 3);
    }
  }
}


//---------------------------------------------------------------------

void MedLinacPhysicsList::SetCuts()
{
  // uppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                                
  SetVerboseLevel(6);               
                                            
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
 defaultCutValue = defaultCut;
 //G4cout <<"--------------------default cut  "<< defaultCutValue/mm << " mm " <<G4endl;
  SetCutsWithDefault();   

  //SetCutValue(cutForGamma, "gamma");
  //SetCutValue(cutForElectron, "e-");
  //SetCutValue(cutForPositron, "e+");

  // Retrieve verbose level
  SetVerboseLevel(temp);  

  // Production thresholds for detector regions

  //G4String regName[] = {"PrimaryCollimatorUp","PrimaryCollimatorLow"};
  //G4double cutValue[] = {1.*mm, 1.*mm};
  //for(G4int i=0;i<2;i++)
  // { 
  //G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName[i]);
  //G4ProductionCuts* cuts = new G4ProductionCuts;
  //cuts->SetProductionCut(cutValue[i]);
  //reg->SetProductionCuts(cuts);
  //}


  G4String regName = "PrimaryCollimatorLow";
  G4double cutValue = 8.*cm;  
  G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName);
  G4ProductionCuts* cuts = new G4ProductionCuts;
  cuts->SetProductionCut(cutValue);
  reg->SetProductionCuts(cuts);

  G4String regName1 = "PrimaryCollimatorUp";
  G4double cutValue1 = 8.*cm;  
  G4Region* reg1 = G4RegionStore::GetInstance()->GetRegion(regName1);
  G4ProductionCuts* cuts1 = new G4ProductionCuts;
  cuts1->SetProductionCut(cutValue1);
  reg1->SetProductionCuts(cuts1);

  

  if (verboseLevel >0){
    G4cout << "MedLinac::SetCuts: default cut length : "
         << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  

}



void MedLinacPhysicsList::SetCut (G4double val)
{ 
  defaultCut = val;
}
