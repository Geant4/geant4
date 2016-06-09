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
//
// Code developed by: S.Guatelli
//
//    **********************************
//    *                                *
//    *     BrachyPhysicsList.cc       *
//    *                                *
//    **********************************
//
// $Id: BrachyPhysicsList.cc,v 1.3 2006/06/29 17:33:29 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
#include "BrachyPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"              


BrachyPhysicsList::BrachyPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 0.1*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  
  SetVerboseLevel(1);
}

BrachyPhysicsList::~BrachyPhysicsList()
{
}

void BrachyPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
}

void BrachyPhysicsList::ConstructBosons()
{ 
  // gamma
  G4Gamma::GammaDefinition();

}

void BrachyPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void BrachyPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

#include "G4MultipleScattering.hh"
// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 
// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 
// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

void BrachyPhysicsList::ConstructEM()
{
  theParticleIterator -> reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator -> value();
    G4ProcessManager* pmanager = particle -> GetProcessManager();
    G4String particleName = particle -> GetParticleName();
    
    //processes
    
    if (particleName == "gamma") {
      //gamma  
      lowePhot = new  G4LowEnergyPhotoElectric("LowEnPhotoElec");
      pmanager -> AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager -> AddDiscreteProcess(lowePhot);
      pmanager -> AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager -> AddDiscreteProcess(new G4LowEnergyGammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
      loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");
      loweBrem -> SetAngularGenerator("tsai");
    
      pmanager -> AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager -> AddProcess(loweIon,     -1, 2,2);
      pmanager -> AddProcess(loweBrem,    -1,-1,3);      
      
    } else if (particleName == "e+") {
      //positron      
      pmanager -> AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager -> AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager -> AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager -> AddProcess(new G4eplusAnnihilation,   0,-1,4);      
      
    }
  }  
}

void BrachyPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "BrachyPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  if (verboseLevel>0) DumpCutValuesTable();
}


