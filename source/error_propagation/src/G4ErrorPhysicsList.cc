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
// $Id: G4ErrorPhysicsList.cc 99974 2016-10-13 07:22:33Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------

#include "globals.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ErrorPhysicsList.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
 
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "G4PhysicsTable.hh"
#include "G4Transportation.hh"

#include "G4ErrorEnergyLoss.hh"

//------------------------------------------------------------------------
G4ErrorPhysicsList::G4ErrorPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0E+9*cm; // set big step so that AlongStep computes all the energy
}


//------------------------------------------------------------------------
G4ErrorPhysicsList::~G4ErrorPhysicsList()
{
}


//------------------------------------------------------------------------
void G4ErrorPhysicsList::ConstructParticle()
{
// In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 
  //  gamma
  G4Gamma::GammaDefinition(); 
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  // pi+/-
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();

  // proton
  G4Proton::ProtonDefinition();

}


//------------------------------------------------------------------------
void G4ErrorPhysicsList::ConstructProcess()
{
  G4Transportation* theTransportationProcess= new G4Transportation();

#ifdef G4VERBOSE
    if (verboseLevel >= 4){
      G4cout << "G4VUserPhysicsList::ConstructProcess()  "<< G4endl;
    }
#endif

  // loop over all particles in G4ParticleTable
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ) {  // Loop checking, 06.08.2015, G.Cosmo
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!particle->IsShortLived()) {
      G4cout << particle << "G4ErrorPhysicsList:: particle process manager " << particle->GetParticleName() << " = " << particle->GetProcessManager() << G4endl;
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        G4String particleName = particle->GetParticleName();
        G4Exception("G4ErrorPhysicsList::ConstructProcess","No process manager",
                    RunMustBeAborted, particleName );
      } else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager ->AddProcess(theTransportationProcess);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
      }
    } else {
      // shortlived particle case
    }
  }

  ConstructEM();
}


//------------------------------------------------------------------------
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"

#include "G4eIonisation.hh"

#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuPairProduction.hh"

#include "G4PhysicsTable.hh"

#include "G4MuIonisation.hh"

#include "G4ErrorStepLengthLimitProcess.hh"
#include "G4ErrorMagFieldLimitProcess.hh"
#include "G4ErrorMessenger.hh"

void G4ErrorPhysicsList::ConstructEM()
{

  G4ErrorEnergyLoss* eLossProcess = new G4ErrorEnergyLoss;
  G4ErrorStepLengthLimitProcess* stepLengthLimitProcess = new G4ErrorStepLengthLimitProcess;
  G4ErrorMagFieldLimitProcess* magFieldLimitProcess = new G4ErrorMagFieldLimitProcess;
  new G4ErrorMessenger( stepLengthLimitProcess, magFieldLimitProcess, eLossProcess );

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ) {  // Loop checking, 06.08.2015, G.Cosmo
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

      //    } else if (particleName == "e-" || particleName == "e+"
      //               || particleName == "mu+" || particleName == "mu-" ) {
    }else if (!particle->IsShortLived() && particle->GetPDGCharge() != 0 ) {
  
      pmanager->AddContinuousProcess(eLossProcess,1);
      pmanager->AddDiscreteProcess( stepLengthLimitProcess, 2 ); 
      pmanager->AddDiscreteProcess( magFieldLimitProcess, 3 );     
      
      /*     } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
      //   G4VProcess* aMultipleScattering = new G4MultipleScattering();
     G4VProcess* anIonisation        = new G4hIonisation();     
     ////G4VProcess*  theUserCuts = new G4UserSpecialCuts();
     
     //
     // add processes
     pmanager->AddProcess(anIonisation);
     //   pmanager->AddProcess(aMultipleScattering);    
     ////pmanager->AddProcess(theUserCuts);
     
     //
     // set ordering for AlongStepDoIt
     //   pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
     pmanager->SetProcessOrdering(anIonisation, idxAlongStep,1);
     
     //
     // set ordering for PostStepDoIt
     //   pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
     pmanager->SetProcessOrdering(anIonisation, idxPostStep,1);
     ////pmanager->SetProcessOrdering(theUserCuts,     idxPostStep,2);
     */
    }
  }
}


//------------------------------------------------------------------------
void G4ErrorPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value or all particle types 
  SetCutsWithDefault(); 
  // if (verboseLevel>0) 
  //  DumpCutValuesTable();
}

