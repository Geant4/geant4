#ifndef A01_NEW_CALORIMETER_PHYSICS_LISTS_HH
#define A01_NEW_CALORIMETER_PHYSICS_LISTS_HH

#include "A01Triggers.hh"

#include "G4VUserPhysicsBiasing.hh"
#include "G4GPRBuilder.hh"
#include "G4GPRTriggerTypes.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Transportation.hh"

using namespace G4GPRBuilder;
using namespace A01Triggers;

class A01Biasing_New_Calorimeter_PhysicsLists : public G4VUserPhysicsBiasing {

public:
  
  void ConstructBiasing() 
  {
    // CreatePhysicsList<Particle>(Trigger type)(List name, TriggeR)
    // AddProcess<Particle>(process, atRestIdx, alongStepIdx, postStepIdx)
    G4String caloListName("Calorimeter_PhysicsList");
    
    CreatePhysicsListWithTrigger<G4Gamma, G4GPRTriggerTypes::Geometry::NewVolume>
      (caloListName, &CalorimeterTrigger);
    AddProcess<G4Gamma>(new G4Transportation,  -1, 0, 0, caloListName);
    
    CreatePhysicsListWithTrigger<G4Electron, G4GPRTriggerTypes::Geometry::NewVolume>
      (caloListName, &CalorimeterTrigger);
    AddProcess<G4Electron>(new G4Transportation,  -1, 0, 0, caloListName);
    
    CreatePhysicsListWithTrigger<G4Positron, G4GPRTriggerTypes::Geometry::NewVolume>
      (caloListName, &CalorimeterTrigger);
    AddProcess<G4Positron>(new G4Transportation,  -1, 0, 0, caloListName);
  }
};
#endif

