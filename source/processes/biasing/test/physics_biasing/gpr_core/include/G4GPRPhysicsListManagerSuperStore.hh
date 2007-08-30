#ifndef G4GPRPHYSICSLISTMANAGERSUPERSTORE_HH
#define G4GPRPHYSICSLISTMANAGERSUPERSTORE_HH

#include "G4GPRPhysicsListManager.hh"
#include "G4GPRAssocT.hh"
#include "G4ParticleDefinition.hh"
#include "G4GPRSingletonHierarchyT.hh"

typedef G4GPRAssocT<G4ParticleDefinition*, G4GPRPhysicsListManager>  Tmp;

typedef G4GPRSingletonHierarchyT< G4GPRTypeList_1( Tmp ) > G4GPRPhysicsListManagerSuperStore;

#endif
