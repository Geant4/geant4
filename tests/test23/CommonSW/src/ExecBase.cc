
#include "ExecBase.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4GenericIon.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"

#include "G4ProcessManager.hh"

#include "G4Gamma.hh"

#include "G4StateManager.hh"
#include "G4ForceCondition.hh"

#include "TstReader.hh"



void ExecBase::Init( const TstReader* )
{

   if(!G4StateManager::GetStateManager()->SetNewState(G4State_PreInit))
       G4cout << "G4StateManager PROBLEM! " << G4endl;

   InitParticles();   
   
   //InitSetup( pset );
   
   //InitBeam( pset );

   if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

   return;

}

void ExecBase::InitParticles()
{

   // physics needs to be initialized before the 1st use of particle table,
   // because it constructs particles - otherwise the table is just empty
   //
   G4MesonConstructor pMesonConstructor;
   pMesonConstructor.ConstructParticle();

   G4BaryonConstructor pBaryonConstructor;
   pBaryonConstructor.ConstructParticle();  
  
   // This is needed because starting 9.6.ref05 G4IonTable::CtreateIon(...)
   // explicitly checks if generic ion has a process manager
   //
   G4GenericIon* gion = G4GenericIon::GenericIon();
   gion->SetProcessManager(new G4ProcessManager(gion));    
   //
   G4IonConstructor pIonConstructor;
   pIonConstructor.ConstructParticle();
  
   G4LeptonConstructor pLeptonConstructor;
   pLeptonConstructor.ConstructParticle();
  
   G4BosonConstructor pBosonConstructor;
   pBosonConstructor.ConstructParticle();

   G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
   partTable->SetReadiness();

   G4IonTable* ionTable = partTable->GetIonTable();
   ionTable->CreateAllIon();
   ionTable->CreateAllIsomer();

   return;

}
