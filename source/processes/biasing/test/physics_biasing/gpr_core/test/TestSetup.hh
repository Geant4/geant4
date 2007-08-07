#ifndef TESTSETUP_HH
#define TESTSETUP_HH

#include <deque>

namespace TestSetup {

  G4GPRTriggerSuperStore* triggerSuperStore = G4GPRTriggerSuperStore::Instance();

  typedef std::vector< G4GPRProcessWrappers::Wrappers<G4GPRProcessLists::DiscreteDoIt>::SeedWrapper > ProcessList;
  
  template <typename Generator>
  void GenerateList(Generator& generator, void (*config)(G4Track*), ProcessList*& result, G4Track* trk, G4Step* step) {

    config(trk);
    G4cout<<"jane generating"<<G4endl;
    triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Tracking::StartTracking>::Fire(trk);
    triggerSuperStore->G4GPRTriggerManagerT<G4GPRScopes::Stepping::StartStep>::Fire(*trk, *step);
    
    generator.template Generate<G4GPRProcessLists::DiscreteDoIt>(result);    
  }

  void PrintList(ProcessList*& list) 
  {
    // Iterate over process list
    for (ProcessList::iterator iter = list->begin(); iter != list->end(); iter++) {
      G4cout<<iter->GetIdentifier()<<" ";
    }
    G4cout<<G4endl;
  }

  // Regular G4VProcess
  struct VProcess : public G4VDiscreteProcess
  {
    VProcess():G4VDiscreteProcess("test"){}
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*){return 0;}
    
    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) {
    //    G4cout<<"Execute VProcess::PostStepDoIt"<<G4endl;
      return 0;
    }
  };
  
  // Alternative process
  struct OtherProcess
  {
    G4VParticleChange* Method(const G4Track&, const G4Step&)
    {
      //    G4cout<<"Execute OtherProcess::Method"<<G4endl;
    }
    
    G4VParticleChange* operator()(const G4Track&, const G4Step&)
    {
      //    G4cout<<"Execute OtherProcess::operator"<<G4endl;
      return 0;
    }
  };

  G4bool PrimaryTrackTrigger(G4Track* track) 
  {
    //  G4cout<<"jane executing MyTrigger "<<track->GetTrackID()<<G4endl;
    return (track->GetTrackID() == 0 ? true : false);
  }

  G4bool VolumeTrigger(const G4Track& track, const G4Step& step) 
  {
    //  G4cout<<"jane executing MyTrigger "<<track->GetTrackID()<<G4endl;
    return (track.GetTrackID() == 0 ? true : false);
  }
  
  class MaxEnergyTrigger {
    
  public:
    MaxEnergyTrigger()
      :fMaxEnergy(0)
    {}
    
    void SetMaxEnergy(G4double maxEnergy) {fMaxEnergy = maxEnergy;}
    G4double GetMaxEnergy() {return fMaxEnergy;}
    
    G4bool operator()(const G4Track& track, const G4Step& step) {
      //    G4cout<<"jane max energy trigger "<<track.GetKineticEnergy()<<G4endl;
      return (track.GetKineticEnergy() < fMaxEnergy ? true : false);
    }
    
  G4bool operator()(const G4Track& track) {
    //      G4cout<<"jane track energy "<<track.GetKineticEnergy()<<G4endl;
    //return (track.GetKineticEnergy() < fMaxEnergy ? true : false);
  }
    
  private:
    
    G4double fMaxEnergy;
  };
  
  
  void ConditionsA(G4Track* trk) {
    trk->SetTrackID(0);
    trk->SetKineticEnergy(99*MeV);
  }
  
  void ConditionsB(G4Track* trk) {
    trk->SetTrackID(0);
    trk->SetKineticEnergy(49*MeV);
  }
  
  void ConditionsC(G4Track* trk) {
    trk->SetTrackID(1);
    trk->SetKineticEnergy(49*MeV);
  }

  // Process function
  G4VParticleChange* MyFunction(const G4Track&, const G4Step&) {
    G4cout<<"Execute MyFunction"<<G4endl;
    return 0;
  }
  
  void Print(const std::deque<G4bool>& key) {
    G4cout<<"Printing key of length "<<key.size()<<": ";
    for (std::deque<G4bool>::const_iterator iter = key.begin(); iter != key.end(); ++iter) {
      G4cout<<*iter<<" ";
    }
    G4cout<<G4endl;
  }

}

#endif
