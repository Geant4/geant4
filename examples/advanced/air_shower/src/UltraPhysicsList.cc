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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraPhysicsList.cc
//    ****************************************************
//
//    Ultra Physics List class; Standard and Low Energy EM processes are defined for
//    the relevant particles. Optical processes are declared.
//
#include "G4ios.hh"
//#include "iomanip.h"
#include "globals.hh"

#include "UltraPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"

UltraPhysicsList::UltraPhysicsList() :  G4VModularPhysicsList(),
           fEmPhysicsList(0),
           fOpPhysicsList(0),
           fDecayPhysicsList(0),
           fVerboseLebel(1),
           fMaxNumPhotonStep(20)
{
  G4LossTableManager::Instance();
  SetDefaultCutValue(1*mm);

  // fMessenger = new UltraPhysicsListMessenger(this);
  // fStepMaxProcess = new StepMax();

  // Initilise flags

  SetVerboseLevel(1);


  // EM physics
  fEmName = G4String("emstandard");
  fEmPhysicsList = new G4EmStandardPhysics();
  fOpPhysicsList = new G4OpticalPhysics();

  // Decay Physics is always defined
  fDecayPhysicsList = new G4DecayPhysics();
}

UltraPhysicsList::~UltraPhysicsList() 
{
  delete fDecayPhysicsList;
  delete fEmPhysicsList;
  delete fOpPhysicsList;
//  delete fStepMaxProcess;
  for(size_t i=0; i<fHadronPhys.size(); i++) 
    delete fHadronPhys[i];
}


void UltraPhysicsList::ConstructParticle()
{
  fDecayPhysicsList->ConstructParticle();
}


void UltraPhysicsList::ConstructProcess()
{
  AddTransportation();
  if (fEmPhysicsList)
    fEmPhysicsList->ConstructProcess();

  if (fOpPhysicsList) 
    fOpPhysicsList->ConstructProcess();
  
  if (fDecayPhysicsList) 
    fDecayPhysicsList->ConstructProcess();

  for(size_t i=0; i<fHadronPhys.size(); ++i) {
    fHadronPhys[i]->ConstructProcess();
  }
}


void UltraPhysicsList::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "UltraPhysicsList::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
