#ifndef G4GPRMANAGER_HH
#define G4GPRMANAGER_HH

#include "G4GPRProcessListGenerator.hh"
#include "G4GPRAssocT.hh"
#include "G4GPRTriggerSuperStore.hh"
#include "G4GPRPhysicsListTriggerSuperStore.hh"
#include "G4GPRPhysicsListManagerSuperStore.hh"

class G4GPRManager {

public:

  G4GPRManager(G4ParticleDefinition* def)
    :pDef(def)
  {
    pPhysicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
    pPhysicsListTrigger = &(*G4GPRPhysicsListTriggerSuperStore::Instance())[def];

    G4GPRPhysicsList* list = pPhysicsListManager->GetDefaultList();

    G4GPRProcessListGenerator* generator = new G4GPRProcessListGenerator(def, list);

    fGenerators.Register(list, generator);
    pCachedGenerator = generator;

    pCachedElementTrigger = &(*G4GPRTriggerSuperStore::Instance())[def][list];
  };

  G4GPRPhysicsList* GetActivePhysicsList() {return pPhysicsListManager->GetActiveList();}

  template <typename Scope, typename Arg1>
  void Fire(const Arg1& arg1) 
  {
    pPhysicsListTrigger->G4GPRTriggerManagerT<Scope>::Fire(arg1);

    UpdateCaches();
    
    pCachedElementTrigger->G4GPRTriggerManagerT<Scope>::Fire(arg1);
  }

  template <typename Scope, typename Arg1, typename Arg2>
  void Fire(const Arg1& arg1, const Arg2& arg2) 
  {
    pPhysicsListTrigger->G4GPRTriggerManagerT<Scope>::Fire(arg1, arg2);

    UpdateCaches();
    
    pCachedElementTrigger->G4GPRTriggerManagerT<Scope>::Fire(arg1, arg2);
  }
  
  template <typename List, typename Result>
  void GetList(Result*& result) 
  {
    UpdateCaches();
    pCachedGenerator->Generate<List>(result);
  }
  //jane make inline  
  void UpdateCaches()
  {

    if (pPhysicsListManager->ListChanged()) {
      G4GPRPhysicsList* newList = pPhysicsListManager->GetActiveList();
      
      if (!fGenerators.Retrieve(newList, pCachedGenerator)) {
	
	G4GPRProcessListGenerator* generator = new G4GPRProcessListGenerator(pDef, newList);
	fGenerators.Register(newList, generator);

	pCachedGenerator = generator;
      }
      pCachedElementTrigger = &(*G4GPRTriggerSuperStore::Instance())[pDef][newList];
      pPhysicsListManager->ResetFlags();
    }
  }

private:

  G4ParticleDefinition* pDef;
  G4GPRAssocT<G4GPRPhysicsList*, G4GPRProcessListGenerator*> fGenerators;

  G4GPRPhysicsListManager* pPhysicsListManager;
  G4GPRTriggerStore* pPhysicsListTrigger;

  G4GPRTriggerStore* pCachedElementTrigger;

  G4GPRProcessListGenerator* pCachedGenerator;

};
#endif
