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
// $Id: MedLinacPhysicsList.cc,v 1.2 2004-04-02 17:48:03 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
#include "globals.hh"
#include "MedLinacPhysicsList.hh"

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

MedLinacPhysicsList::MedLinacPhysicsList():  G4VUserPhysicsList()
{ 
  defaultCutValue = 0.1*mm;
   SetVerboseLevel(1);
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
  ConstructGeneral();
}

//---------------------------------------------------------------------

//#include "G4ComptonScattering.hh"
//#include "G4GammaConversion.hh"
//#include "G4PhotoElectricEffect.hh"

#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyRayleigh.hh" 
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 
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
      //pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      //pmanager->AddDiscreteProcess(new G4ComptonScattering);
      //pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
      //electron
     
      loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
      loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");
    
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(loweIon,     -1, 2,2);
      pmanager->AddProcess(loweBrem,    -1,-1,3);     


      //pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      //pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);    
  
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
    }
  }
}


//---------------------------------------------------------------------

#include "G4Decay.hh"
void MedLinacPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//-------------------------------------------------------------------

void MedLinacPhysicsList::SetCuts()
{
  // uppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                                SetVerboseLevel(0);                                                           
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   

  // Retrieve verbose level
  SetVerboseLevel(temp);  
}

