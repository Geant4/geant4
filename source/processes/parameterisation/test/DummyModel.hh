#ifndef DummyModel_hh
#define DummyModel_hh

#include "G4VFastSimulationModel.hh"

class DummyModel : public G4VFastSimulationModel {
public:
  DummyModel(const G4String& aName, G4LogicalVolume* vol) : G4VFastSimulationModel(aName, vol) {}

  ~DummyModel() {}
  
  G4bool IsApplicable(const G4ParticleDefinition&) {return true;}
  G4bool ModelTrigger(const G4FastTrack &) {return true;}
  void DoIt(const G4FastTrack&, G4FastStep&) {}

};

#endif
