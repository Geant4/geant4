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
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//
//
//    ******      SUGGESTED PHYSICS FOR ACCURATE SIMULATIONS    *********
//
// At moment, if accurate simulations are necessary, we suggest the use of the
// Physics Lists 'HADRONTHERAPY_1';
// It can be activated inside any macro file using the command:
// /Physics/addPhysics HADRONTHERAPY_1

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyStepMax.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

// Local physic directly implemented in the Hadronthrapy directory
// Physic dedicated to the ion-ion inelastic processes
//
#include "LocalIonIonInelasticPhysic.hh"

#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4AutoDelete.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::HadrontherapyPhysicsList() : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();
    defaultCutValue = 1.*mm;
    cutForGamma     = defaultCutValue;
    cutForElectron  = defaultCutValue;
    cutForPositron  = defaultCutValue;
    
    pMessenger = new HadrontherapyPhysicsListMessenger(this);
    
    SetVerboseLevel(1);

    
    // ******     Definition of defaults for the physics processes *****
    // ******     in case no physics is called by the macro file   *****
    //
    // The default physics corresponds to the actual QGSP_BIC_HP list
    // but with the following differences:
    // --> G4EmStandardPhysics_option4 for the electromagnetic processes
    //     is used n place of the less accurate G4EmStandardPhysics
    // --> The G4RadioactiveDecayPhysics is add
    // --> G4HadronPhysicsQGSP_BIC is used in place of G4HadronPhysicsQGSP_BIC_HP
    // --> G4HadronElasticPhysics is used in place of G4HadronElasticPhysics_HP
    
    // Elecromagnetic physics
    //
    emPhysicsList = new G4EmStandardPhysics_option4();
    emName = G4String("emstandard_opt4");
    
    // Hadronic physics
    //
    hadronPhys.push_back( new G4DecayPhysics());
    hadronPhys.push_back( new G4RadioactiveDecayPhysics());
    hadronPhys.push_back( new G4IonBinaryCascadePhysics());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysics());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
    
    // Decay physics
    //
    decay_List = new G4DecayPhysics();
    radioactiveDecay_List = new G4RadioactiveDecayPhysics();
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{
    delete pMessenger;
    delete emPhysicsList;
    delete decay_List;  
    delete radioactiveDecay_List;
    hadronPhys.clear();
    for(size_t i=0; i<hadronPhys.size(); i++)
    {
        delete hadronPhys[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::ConstructParticle()
{
    decay_List -> ConstructParticle();
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::ConstructProcess()
{
    // Transportation
    //
    AddTransportation();
    
    // Electromagnetic physics
    //
    emPhysicsList -> ConstructProcess();
    em_config.AddModels();
    
    // Hadronic physics
    //
    for(size_t i=0; i < hadronPhys.size(); i++)
    {
        hadronPhys[i] -> ConstructProcess();
    }
    
    // step limitation (as a full process)
    //
    AddStepMax();
    
    //Parallel world sensitivity
    //
    G4ParallelWorldPhysics* pWorld = new G4ParallelWorldPhysics("DetectorROGeometry");
    pWorld->ConstructProcess();
    
    return;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddPhysicsList(const G4String& name)
{
    if (verboseLevel>1) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == emName) return;
    
    ///////////////////////////////////
    //   ELECTROMAGNETIC MODELS
    ///////////////////////////////////
    if (name == "standard_opt4") {
        emName = name;
        delete emPhysicsList;
        hadronPhys.clear();
        emPhysicsList = new G4EmStandardPhysics_option4();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
        
        // The following 'local_ion_ion_inelastic' is an example of implemenation of
        // inelastic hadronic models to be used when
        // mucleus-nulceus (ion-ion) interactions have to be
        // taken into account;
        // It must be used, of course, in connection with other lists: electromagnetic
        // plus hadronic elastic plus nucleon-nulcleon hadronic inelastic
        //
        // An example of coplete physics list using the 'local_ion_ion_inelastic' physics
        // is the one named HADRONTHERAPY_2, and defined below.
        //
    } else if (name == "standard_opt3") {
        emName = name;
        delete emPhysicsList;
        hadronPhys.clear();
        emPhysicsList = new G4EmStandardPhysics_option3();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;
    } else if (name == "local_ion_ion_inelastic") {
        hadronPhys.push_back(new LocalIonIonInelasticPhysic());
        locIonIonInelasticIsRegistered = true;
        
    } else if (name == "HADRONTHERAPY_1") {
        
        // The HADRONTHERAPY_1 physics list corresponds to the actual QGSP_BIC_HP list
        // but with the following differences:
        // --> G4EmStandardPhysics_option4 for the electromagnetic processes
        //     is used in place of the less accurate G4EmStandardPhysics
        // --> The G4RadioactiveDecayPhysics is added
        
        AddPhysicsList("standard_opt4");
        hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysicsHP());
        hadronPhys.push_back( new G4StoppingPhysics());
        hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
        
        G4cout << "HADRONTHERAPY_1 PHYSICS LIST has been activated" << G4endl;
    
        
    } else if (name == "HADRONTHERAPY_2") {
        
        // The HADRONTHERAPY_2 physics list corresponds to the actual QGSP_BIC_HP list
        // but with the following differences:
        // --> G4EmStandardPhysics_option4 for the electromagnetic processes
        //     is used in place of the less accurate G4EmStandardPhysics
        // --> The G4RadioactiveDecayPhysics is added
        // --> The 'local_ion_ion_inelastic' physics is used in place of the
        //     G4IonBinaryCascadePhysics(): it used the QMD model to treat
        //     the ion-ion inelastic interactions
        
        AddPhysicsList("standard_opt4");
        AddPhysicsList("local_ion_ion_inelastic");
        hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysicsHP());
        hadronPhys.push_back( new G4StoppingPhysics());
        hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
        
        G4cout << "HADRONTHERAPY_2 PHYSICS LIST has been acivated" << G4endl;
        
    } 
    else if (name == "QGSP_BIC_EMY") {
      AddPhysicsList("standard_opt3");
      //emPhysicsList = new G4EmLivermorePhysics();
      hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
      hadronPhys.push_back( new G4EmExtraPhysics());
      hadronPhys.push_back( new G4HadronElasticPhysics());
      hadronPhys.push_back( new G4StoppingPhysics());
      hadronPhys.push_back( new G4IonBinaryCascadePhysics());
      hadronPhys.push_back( new G4NeutronTrackingCut());
      hadronPhys.push_back( new G4DecayPhysics());
    }   
     else {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
        << " is not defined"
        << G4endl;
    }
    
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  // This process must exist in all threads.
  //
  HadrontherapyStepMax* stepMaxProcess  = new HadrontherapyStepMax();
  G4AutoDelete::Register( stepMaxProcess );
  
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (stepMaxProcess->IsApplicable(*particle) && pmanager)
      {
	pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::SetCuts()
{

    if (verboseLevel >0){
        G4cout << "PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }
    
    // Production thresholds for detector regions
    // The G4Regions, for which you want define a given cut via de macro command
    // '/run/setCutForRegion <G4Region name> <cut value>'
    // must be defined here
    //
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
    // At moment, only 'DetectorLog' is defined as G4Region
    //
    G4String regName[] = {"DetectorLog"};
    G4double fuc = 1.;
    for(G4int i=0;i<1;i++)
    {
        G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName[i]);
        G4ProductionCuts* cuts = new G4ProductionCuts;
        cuts->SetProductionCut(defaultCutValue*fuc);
        reg->SetProductionCuts(cuts);
        fuc *= 10.;
    }
}
