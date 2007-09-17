#ifndef G4GPRBUILDER_HH
#define G4GPRBUILDER_HH

// Collection of possible user configuration functions.
// jane fixme - memory management of added processes
#include "G4String.hh"
#include "G4GPRBuilderUtils.hh"

namespace G4GPRBuilder {

  template <typename Particle>
  void CreateDefaultPhysicsList(const G4String& id)
  {
    G4GPRBuilderUtils::CreateDefaultPhysicsList<Particle>(id);
  }

  // jane fixme default list idx
  template <typename Particle>
  void AddProcess(G4VProcess* process, G4int restIdx, G4int continuousIdx, G4int discreteIdx, 
		  const G4String& physicsList = "")
  {
    //jane fixme - memory management for vprocess
    //jane cleanup
    G4GPRBuilderUtils::LoadObservers<Particle>(process);
    G4GPRBuilderUtils::InitialiseGPR<Particle>();

    if (restIdx != -1) G4GPRBuilderUtils::AddRestVProcess<Particle>(process, restIdx, physicsList);
    if (continuousIdx != -1) G4GPRBuilderUtils::AddContinuousVProcess<Particle>(process, continuousIdx, physicsList);
    if (discreteIdx != -1) G4GPRBuilderUtils::AddDiscreteVProcess<Particle>(process, discreteIdx, physicsList);
  }

  template <typename Particle, typename TriggerType, typename Trigger>
  void CreatePhysicsListWithTrigger(const G4String& id, const Trigger& trigger)
  {
    G4GPRBuilderUtils::CreatePhysicsListWithTrigger<Particle, TriggerType>(id, trigger);
  }

  template <typename Particle, typename ProcessList, typename TriggerType, 
	    typename Biasing, typename Trigger>
  void AddTriggeredBiasing(const G4String id, G4int idx, const Biasing& biasing, 
			   const Trigger& trigger, const G4String& physicsList="")
  {
    G4GPRBuilderUtils::AddRelayWithTrigger<Particle, ProcessList, Biasing, TriggerType, Trigger>(id, idx, biasing, trigger, physicsList);
  }

  template <typename ProcessList, typename TriggerType, 
	    typename Biasing, typename Trigger>
  void AddTriggeredBiasing(const G4String id, const Biasing& biasing, 
			   const Trigger& trigger, G4GPRBiasingConfig& cfg, 
			   const G4String& physicsList="")
  {
    G4GPRBuilderUtils::AddRelayWithTrigger<ProcessList, TriggerType, Biasing, Trigger>
      (id, biasing, trigger, cfg, physicsList);
  }

  template <typename ProcessList, typename Biasing>
  void AddBiasing(const G4String id, const Biasing& biasing, G4GPRBiasingConfig& cfg, 
		  const G4String& physicsList="")
  {
    G4GPRBuilderUtils::AddRelay<ProcessList, Biasing>(id, biasing, cfg, physicsList);
  }

}

#endif
