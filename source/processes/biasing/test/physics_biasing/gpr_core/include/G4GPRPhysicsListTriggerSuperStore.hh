#ifndef G4GPRPHYSICSLISTTRIGGERSUPERSTORE_HH
#define G4GPRPHYSICSLISTTRIGGERSUPERSTORE_HH

#include "G4GPRTriggerSuperStore.hh"

typedef G4GPRAssocT<G4ParticleDefinition*, G4GPRTriggerStore> TmpAssoc;

typedef G4GPRSingletonHierarchyT< G4GPRTypeList_1(TmpAssoc) > G4GPRPhysicsListTriggerSuperStore;

#endif
