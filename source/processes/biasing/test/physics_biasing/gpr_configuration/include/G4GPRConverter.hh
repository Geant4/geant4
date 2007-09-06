#ifndef G4GPRCONVERTER_HH
#define G4GPRCONVERTER_HH

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4GPRPhysicsListManagerSuperStore.hh"
#include "G4GPRPhysicsList.hh"
#include "G4GPRProcessLists.hh"
#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRManager.hh"
#include "G4GPRSeedT.hh"

namespace G4GPRConverter {

  template <typename List, typename Mfn>
  void LoadList(G4GPRElementStore* store, G4ProcessVector* vect, Mfn mfn) 
  {
    typedef G4GPRSeedT<List> Seed;

    G4int i(0);
    G4int size = vect->size();

    for (i=0; i<size; i++) {
      G4cout<<"jane adding seed for idx "<<i<<" "<<size<<G4endl;
      store->G4GPRManagerT<Seed>::Register(new Seed((*vect)[i]->GetProcessName(), (*vect)[i], mfn, i));
    }
  }

  void LoadGPR(G4ParticleDefinition* def)
  {
    G4cout<<"jane loading particle "<<def->GetParticleName()<<G4endl;
    if (0 != def->GetGPRManager()) {
      G4cout<<"jane already built"<<G4endl;
      return;
    }
    G4ProcessManager* processManager = def->GetProcessManager();

    G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
    G4GPRPhysicsList* physicsList = physicsListManager->GetDefaultList();
    
    G4GPRElementStore* store = &(*G4GPRElementSuperStore::Instance())[def][physicsList];
 
    LoadList<G4GPRProcessLists::AtRestGPIL>(store, processManager->GetAtRestProcessVector(typeGPIL), &G4VProcess::AtRestGPIL);
    LoadList<G4GPRProcessLists::AtRestDoIt>(store, processManager->GetAtRestProcessVector(typeDoIt), &G4VProcess::AtRestDoIt);

    LoadList<G4GPRProcessLists::ContinuousGPIL>(store, processManager->GetAlongStepProcessVector(typeGPIL), &G4VProcess::AlongStepGPIL);
    LoadList<G4GPRProcessLists::ContinuousDoIt>(store, processManager->GetAlongStepProcessVector(typeDoIt), &G4VProcess::AlongStepDoIt);   

    LoadList<G4GPRProcessLists::DiscreteGPIL>(store, processManager->GetPostStepProcessVector(typeGPIL), &G4VProcess::PostStepGPIL);
    LoadList<G4GPRProcessLists::DiscreteDoIt>(store, processManager->GetPostStepProcessVector(typeDoIt), &G4VProcess::PostStepDoIt);

    G4GPRManager* gprManager = new G4GPRManager(def);
    def->SetGPRManager(gprManager);
  }

  template <typename List>
  void CrossCheck(G4ProcessVector* pm, List* gm)
  {
    G4cout<<"jane list size: "<<pm->size()<<" "<<gm->size()<<G4endl;
    assert (pm->size() == static_cast<G4int>(gm->size()));

    for (G4int i=0; i<pm->size(); ++i) {
      G4cout<<"jane vproc: "<<(*pm)[i]->GetProcessName()<<" "<<", wrapper: "<<(*gm)[i].GetIdentifier()<<G4endl;
      assert((*pm)[i]->GetProcessName() == (*gm)[i].GetIdentifier());
    }
    
  }

  void TestLoading(G4ParticleDefinition* def)
  {
    G4cout<<"jane test loading "<<G4endl;
    G4ProcessManager* processManager = def->GetProcessManager();
    G4GPRManager* gprManager = def->GetGPRManager();

    assert (0 != processManager);
    assert (0 != gprManager);

    G4ProcessVector* atRestGPIL_pm = processManager->GetAtRestProcessVector(typeGPIL);
    G4ProcessVector* atRestDoIt_pm = processManager->GetAtRestProcessVector(typeDoIt);

    G4ProcessVector* continuousGPIL_pm = processManager->GetAlongStepProcessVector(typeGPIL);
    G4ProcessVector* continuousDoIt_pm = processManager->GetAlongStepProcessVector(typeDoIt);

    G4ProcessVector* discreteGPIL_pm = processManager->GetPostStepProcessVector(typeGPIL);
    G4ProcessVector* discreteDoIt_pm = processManager->GetPostStepProcessVector(typeDoIt);


    std::vector<G4GPRProcessWrappers::G4GPRAtRestGPIL>* atRestGPIL_gp(0);
    std::vector<G4GPRProcessWrappers::G4GPRAtRestDoIt>* atRestDoIt_gp(0);

    std::vector<G4GPRProcessWrappers::G4GPRContinuousGPIL>* continuousGPIL_gp(0);
    std::vector<G4GPRProcessWrappers::G4GPRContinuousDoIt>* continuousDoIt_gp(0);

    std::vector<G4GPRProcessWrappers::G4GPRDiscreteGPIL>* discreteGPIL_gp(0);
    std::vector<G4GPRProcessWrappers::G4GPRDiscreteDoIt>* discreteDoIt_gp(0);

    gprManager->GetList<G4GPRProcessLists::AtRestGPIL>(atRestGPIL_gp);
    gprManager->GetList<G4GPRProcessLists::AtRestDoIt>(atRestDoIt_gp);

    gprManager->GetList<G4GPRProcessLists::ContinuousGPIL>(continuousGPIL_gp);
    gprManager->GetList<G4GPRProcessLists::ContinuousDoIt>(continuousDoIt_gp);

    gprManager->GetList<G4GPRProcessLists::DiscreteGPIL>(discreteGPIL_gp);
    gprManager->GetList<G4GPRProcessLists::DiscreteDoIt>(discreteDoIt_gp);

    G4cout<<"jane checking atrestgpil:"<<G4endl;
    CrossCheck(atRestGPIL_pm, atRestGPIL_gp);
    G4cout<<"jane checking atrestdoit:"<<G4endl;
    CrossCheck(atRestDoIt_pm, atRestDoIt_gp);

    G4cout<<"jane checking continuousgpil:"<<G4endl;
    CrossCheck(continuousGPIL_pm, continuousGPIL_gp);
    G4cout<<"jane checking continuousdoit:"<<G4endl;
    CrossCheck(continuousDoIt_pm, continuousDoIt_gp);

    G4cout<<"jane checking discretegpil:"<<G4endl;
    CrossCheck(discreteGPIL_pm, discreteGPIL_gp);
    G4cout<<"jane checking discretedoit:"<<G4endl;
    CrossCheck(discreteDoIt_pm, discreteDoIt_gp);
  }

  void LoadGPR() {
    G4ParticleTable::G4PTblDicIterator* particles = G4ParticleTable::GetParticleTable()->GetIterator();
    
    // loop over all particles in G4ParticleTable
    particles->reset();
    while( (*particles)() ){
      G4ParticleDefinition* particle = particles->value();  
    G4GPRConverter::LoadGPR(particle);
    
    // Test loading was successful. Will also force default physics list
    // to be constructed for all particles
    G4cout<<"jane test loading"<<G4endl;
    G4GPRConverter::TestLoading(particle);
    }
  }
  

}
#endif
