#ifndef G4GPRBUILDERUTILS_HH
#define G4GPRBUILDERUTILS_HH

#include "G4ParticleDefinition.hh"
#include "G4GPRConverter.hh"
#include "G4GPRPhysicsListManagerSuperStore.hh"
#include "G4GPRPhysicsList.hh"
#include "G4GPRSingleProcessRelayT.hh"
#include "G4GPRKeySuperStore.hh"
//jane fixme - memory management for processes

namespace G4GPRBuilderUtils {

  template <typename Particle>
  void InitialiseGPR()
  {
    G4ParticleDefinition* def = Particle::Definition();
    if (0 == def->GetGPRManager()) def->SetGPRManager(new G4GPRManager(def));
  }

  template <typename Particle>
  void LoadObservers(G4VProcess* process)
  {
    G4GPRConverter::LoadObservers(&(*G4GPRObserverSuperStore::Instance())[Particle::Definition()], process);
  }

  template <typename Particle>
  void CreateDefaultPhysicsList(const G4String& id)
  {
    G4ParticleDefinition* def = Particle::Definition();
    
    // Load generalised processing from G4ProcessManager.
    G4GPRConverter::LoadGPR(def);

    G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];

    // Create new default - "true" means all processes are active by default
    physicsListManager->SetDefaultList(new G4GPRPhysicsList(id, true));
  }

  template <typename Particle, typename TriggerType, typename Trigger>
  void CreatePhysicsListWithTrigger(const G4String& id, const Trigger& trigger)
  {
    G4ParticleDefinition* def = Particle::Definition();
    
    // Load generalised processing from G4ProcessManager.
    G4GPRConverter::LoadGPR(def);

    G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];

    G4GPRPhysicsList* list = new G4GPRPhysicsList(id, false);
    physicsListManager->Register(list);

    // Register physics list triggers with physics list trigger manager
    G4GPRTriggerStore* physicsListTrigger = &(*G4GPRPhysicsListTriggerSuperStore::Instance())[def];
    
    physicsListTrigger->G4GPRTriggerManagerT<TriggerType>::Register(trigger, list, &G4GPRPhysicsList::FlipState, list->GetState());
    G4cout<<"jane done triggered physics list"<<G4endl;
  }

  template <typename Particle>
  G4GPRElementStore* GetStore(const G4String& physicsListName="")
  {
    G4ParticleDefinition* def = Particle::Definition();
    
    G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];

    G4GPRPhysicsList* physicsList = physicsListManager->GetPhysicsList(physicsListName);
    G4cout<<"jane adding to physics list "<<physicsList->GetName()<<G4endl;
    assert (0 != physicsList);

    G4GPRElementStore* elementStore = &(*G4GPRElementSuperStore::Instance())[def][physicsList];

    return elementStore;
  }

  template <typename Particle>
  void AddRestVProcess(G4VProcess* process, G4int idx, const G4String& physicsListName="")
  {
    typedef G4GPRSeedT<G4GPRProcessLists::AtRestGPIL> GPIL;
    typedef G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> DoIt;

    G4GPRElementStore* store = GetStore<Particle>(physicsListName);
    
    store->G4GPRManagerT<GPIL>::Register(new GPIL(process->GetProcessName(), process, &G4VProcess::AtRestGPIL, idx));
    store->G4GPRManagerT<DoIt>::Register(new DoIt(process->GetProcessName(), process, &G4VProcess::AtRestDoIt, idx));
  }

  template <typename Particle>
  void AddContinuousVProcess(G4VProcess* process, G4int idx, const G4String& physicsListName="")
  {
    typedef G4GPRSeedT<G4GPRProcessLists::ContinuousGPIL> GPIL;
    typedef G4GPRSeedT<G4GPRProcessLists::ContinuousDoIt> DoIt;

    G4GPRElementStore* store = GetStore<Particle>(physicsListName);
    
    store->G4GPRManagerT<GPIL>::Register(new GPIL(process->GetProcessName(), process, &G4VProcess::AlongStepGPIL, idx));
    store->G4GPRManagerT<DoIt>::Register(new DoIt(process->GetProcessName(), process, &G4VProcess::AlongStepDoIt, idx));
  }

  template <typename Particle>
  void AddDiscreteVProcess(G4VProcess* process, G4int idx, const G4String& physicsListName="")
  {
    typedef G4GPRSeedT<G4GPRProcessLists::DiscreteGPIL> GPIL;
    typedef G4GPRSeedT<G4GPRProcessLists::DiscreteDoIt> DoIt;

    G4GPRElementStore* store = GetStore<Particle>(physicsListName);
    
    store->G4GPRManagerT<GPIL>::Register(new GPIL(process->GetProcessName(), process, &G4VProcess::PostStepGPIL, idx));
    store->G4GPRManagerT<DoIt>::Register(new DoIt(process->GetProcessName(), process, &G4VProcess::PostStepDoIt, idx));
  }

  template <typename Particle, typename TriggerType, typename BiasingType, typename Mfn, typename Trigger>
  void RegisterTrigger(BiasingType*& biasing, Mfn mfn, const Trigger& trigger, const G4String& physicsListName="")
  {
    G4ParticleDefinition* def = Particle::Definition();

    G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
    G4GPRPhysicsList* physicsList = physicsListManager->GetPhysicsList(physicsListName);
    
    G4GPRTriggerStore* triggerStore = &(*G4GPRTriggerSuperStore::Instance())[def][physicsList];

    triggerStore->G4GPRTriggerManagerT<TriggerType>::Register(trigger, biasing, mfn);
    
    // Create and register key nodes with trigger manager so that know when an element has been activated or deactivated
    G4GPRNode* node1 = new G4GPRNode;

    triggerStore->G4GPRTriggerManagerT<TriggerType>::Register(trigger, node1, &G4GPRNode::FlipState);

    G4GPRKeyStore* keyStore = &(*G4GPRKeySuperStore::Instance())[def][physicsList];
    
    keyStore->G4GPRKeyManagerT<typename BiasingType::List>::AddNode(node1);
  }
  template <typename Particle, typename ProcessList, typename Biasing, typename TriggerType, typename Trigger>
  void AddRelayWithTrigger(const G4String& id, G4int idx, const Biasing& biasing, 
			   const Trigger& trigger, const G4String& physicsListName="")
  {
    G4ParticleDefinition* def = Particle::Definition();

    G4GPRConverter::LoadGPR(def);

    typedef G4GPRSingleProcessRelayT<ProcessList> Relay; 

    Relay* relay = new Relay(id, biasing, idx);
    G4GPRElementStore* store = GetStore<Particle>(physicsListName);

    store->G4GPRManagerT<Relay>::Register(relay);
    
    RegisterTrigger<Particle, TriggerType>(relay, &G4GPRSingleProcessRelayT<ProcessList>::ChangeState, trigger, physicsListName);

  }
  

}
#endif
