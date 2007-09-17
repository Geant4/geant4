#ifndef G4GPRSELECTPROCESSES_HH
#define G4GPRSELECTPROCESSES_HH

#include "G4VProcess.hh"
#include "G4GPRProcessLists.hh"
#include "G4ProcessManager.hh"

// Private utilities
namespace {
  
  // jane fixme - yuk, yuk, yuk, rewrite all of this, just wanted
  // something quickly. Just want to store types, find a better way, maybe typelists

  template <typename List>
  struct GetVector {
    
    G4ProcessVector* operator()(G4ParticleDefinition*) {
      assert(0);
    }
  };

  template<>
  struct GetVector<G4GPRProcessLists::DiscreteDoIt> {

    G4ProcessVector* operator()(G4ParticleDefinition* def)
    {
      return def->GetProcessManager()->GetPostStepProcessVector(typeDoIt);
    }
  };

  struct VProcessTypeBase {

    template <typename List>
    std::vector<unsigned> GetIndices(G4ParticleDefinition* def) const
    {
      G4cout<<"jane in getindices"<<G4endl;
      std::vector<unsigned> result;
      GetVector<List> tmp;
      G4ProcessVector* processVector = tmp(def);
      G4cout<<"jane got process vector "<<processVector->size()<<G4endl;
      for (G4int i=0; i<processVector->size(); i++) {
	G4VProcess* process = (*processVector)[i];
	G4bool success = SameType(process);
	G4cout<<"jane getindices "<<def->GetParticleName()<<" "<<process->GetProcessName()<<" "<<success<<G4endl;
	if (success) result.push_back(i);
      }
      
      return result;
    }

    virtual G4bool SameType(G4VProcess*) const =0;


  };

  // just want to store type
  template <typename T>
  struct VProcessType : public VProcessTypeBase {

    typedef T Process;

    G4bool SameType(G4VProcess* incoming) const
    {
      G4cout<<"jane incoming type "<<incoming->GetProcessName()<<G4endl;
      return (0 != dynamic_cast<Process*>(incoming));
    }

  };
}


class G4GPRSelectProcesses {

public:

  typedef std::vector<unsigned> Result;

  G4GPRSelectProcesses()
    :fAllProcesses(false)
  {}

  void SelectAllProcesses() 
  {
    fAllProcesses = true;
  };

  void SelectProcess(unsigned idx)
  {
    fIndices.push_back(idx);
  }
  
  template <typename Type>
  void SelectVProcess()
  {
    fVProcesses.push_back(new VProcessType<Type>);
  }
  /*
  void SelectProcess(const G4String& id)
  {
    fIds.push_back(id);
  }
  */
  template <typename List>
  std::vector<unsigned> GetList(G4ParticleDefinition* def) const
  {
    G4cout<<"jane getting list for "<<def->GetParticleName()<<G4endl;
    //jane fixme - coded just to get it running - rewrite for coherence
    Result result;

    std::vector<unsigned> processTypes;

    result.insert(result.end(), fIndices.begin(), fIndices.end());
    Result vProcesses = GetVProcessIndices<List>(def);

    result.insert(result.end(), vProcesses.begin(), vProcesses.end());

    return result;
  }

  template <typename List>
  Result GetVProcessIndices(G4ParticleDefinition* def) const
  {
    G4cout<<"jane getting process for "<<def->GetParticleName()<<" "<<fVProcesses.size()<<G4endl;
    Result result;

    for (std::vector<VProcessTypeBase*>::const_iterator iter = fVProcesses.begin();
	 iter != fVProcesses.end(); ++iter) {
      G4cout<<"jane calling getindices"<<G4endl;
      Result indices = (*iter)->template GetIndices<List>(def);
      result.insert(result.end(), indices.begin(), indices.end());
    }

    return result;

  }
  
private:
  G4bool fAllProcesses;
  std::vector<unsigned> fIndices;
  std::vector<G4String> fIds;

  std::vector<VProcessTypeBase*> fVProcesses;
};

#endif
