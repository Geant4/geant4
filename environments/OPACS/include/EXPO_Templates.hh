//
// I need that for SunOS-CC !!!
//    G.Barrand
//
// RW :
#include "G4ApplicationState.hh"
#include "G4Element.hh"
#include "G4Event.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4MPVEntry.hh"
#include "G4Material.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4NavigationLevel.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include "G4ProcTblElement.hh"
#include "G4ProcessAttribute.hh"
#include "G4ProcessManager.hh"
#include "G4SDStructure.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelProxy.hh"
#include "G4Track.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UImessenger.hh"
#include "G4UIparameter.hh"
#include "G4VDecayChannel.hh"
#include "G4VDigiCollection.hh"
#include "G4VHitsCollection.hh"
#include "G4VModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"
#include "G4VStateDependent.hh"
#include "G4VTrajectory.hh"
#include "G4UnitsTable.hh"
#include "G4UItokenNum.hh"
#include "G4Trajectory.hh"

template class RWTPtrSlistDictionary<RWCString, G4MaterialPropertyVector>;
template class RWTPtrSlistDictionary<RWCString, G4ParticleDefinition>;
template class RWTPtrSlistDictionary<int, G4ParticleDefinition>;
template class RWTPtrSlistDictionaryIterator<RWCString, G4MaterialPropertyVector>;
template class RWTPtrSlistDictionaryIterator<RWCString, G4ParticleDefinition>;
template class RWTPtrSlistDictionaryIterator<int, G4ParticleDefinition>;
template class RWTValSlistDictionary<RWCString, void*>;
template class RWTValSlistDictionary<const G4ParticleDefinition*, G4EnergyLossTablesHelper>;
template class RWTValSlistDictionary<unsigned long, int>;
template class RWTValVector<yystype>;

template class RWTPtrHashDictionary<RWCString, G4MaterialPropertyVector>;
template class RWTPtrHashDictionary<RWCString, G4ParticleDefinition>;
template class RWTPtrHashDictionary<int, G4ParticleDefinition>;
template class RWTPtrHashDictionaryIterator<RWCString, G4MaterialPropertyVector>;
template class RWTPtrHashDictionaryIterator<RWCString, G4ParticleDefinition>;
template class RWTPtrHashDictionaryIterator<int, G4ParticleDefinition>;
template class RWTPtrOrderedVector<G4Element>;
template class RWTPtrOrderedVector<G4Event>;
template class RWTPtrOrderedVector<G4Isotope>;
template class RWTPtrOrderedVector<G4LogicalVolume>;
template class RWTPtrOrderedVector<G4MPVEntry>;
template class RWTPtrOrderedVector<G4Material>;
template class RWTPtrOrderedVector<G4ParticleDefinition>;
template class RWTPtrOrderedVector<G4PhysicsVector>;
template class RWTPtrOrderedVector<G4ProcTblElement>;
template class RWTPtrOrderedVector<G4ProcessAttribute>;
template class RWTPtrOrderedVector<G4ProcessManager>;
template class RWTPtrOrderedVector<G4SDStructure>;
template class RWTPtrOrderedVector<G4SmartVoxelNode>;
template class RWTPtrOrderedVector<G4SmartVoxelProxy>;
template class RWTPtrOrderedVector<G4Track>;
template class RWTPtrOrderedVector<G4UIcommand>;
template class RWTPtrOrderedVector<G4UIcommandTree>;
template class RWTPtrOrderedVector<G4UImessenger>;
template class RWTPtrOrderedVector<G4UIparameter>;
template class RWTPtrOrderedVector<G4UnitDefinition>;
template class RWTPtrOrderedVector<G4UnitsCategory>;
template class RWTPtrOrderedVector<G4VDecayChannel>;
template class RWTPtrOrderedVector<G4VDigiCollection>;
template class RWTPtrOrderedVector<G4VHitsCollection>;
template class RWTPtrOrderedVector<G4VModel>;
template class RWTPtrOrderedVector<G4VPhysicalVolume>;
template class RWTPtrOrderedVector<G4VProcess>;
template class RWTPtrOrderedVector<G4VSensitiveDetector>;
template class RWTPtrOrderedVector<G4VSolid>;
template class RWTPtrOrderedVector<G4VStateDependent>;
template class RWTPtrOrderedVector<G4VTrajectory>;
template class RWTPtrOrderedVector<G4ValVector>;
template class RWTPtrSortedVector<G4MPVEntry>;
template class RWTPtrSortedVector<G4VDecayChannel>;
template class RWTValHashDictionary<RWCString, void*>;
template class RWTValHashDictionary<const G4ParticleDefinition*, G4EnergyLossTablesHelper>;
template class RWTValHashDictionary<unsigned long, int>;
template class RWTValOrderedVector<G4ApplicationState>;
template class RWTValOrderedVector<Hep3Vector>;
//template class RWTValOrderedVector<HepPlane3D>;
template class RWTValOrderedVector<HepPoint3D>;
template class RWTValOrderedVector<RWCString>;
template class RWTValOrderedVector<double>;
template class RWTValOrderedVector<int (*)(void*)>;
template class RWTValOrderedVector<int>;
template class RWTValOrderedVector<void (*)(void)>;
template class RWTValOrderedVector<void*>;
template class RWTValOrderedVector<yystype>;
template class RWTValVector<G4ApplicationState>;
template class RWTValVector<G4NavigationLevel>;
template class RWTValVector<Hep3Vector>;
//template class RWTValVector<HepPlane3D>;
template class RWTValVector<HepPoint3D>;
template class RWTValVector<RWCString>;
template class RWTValVector<double>;
template class RWTValVector<int (*)(void*)>;
template class RWTValVector<int>;
template class RWTValVector<void (*)(void)>;
template class RWTValVector<void*>;
template class RWTPtrOrderedVector<G4VTrajectoryPoint>;

//
// G4 :
#include "G4DecayProducts.hh"
#include "G4Run.hh"

template class G4Allocator<G4DCofThisEvent>;
template class G4Allocator<G4DecayProducts>;
template class G4Allocator<G4DynamicParticle>;
template class G4Allocator<G4Event>;
template class G4Allocator<G4HCofThisEvent>;
template class G4Allocator<G4NavigationLevel>;
template class G4Allocator<G4NavigationLevelRep>;
template class G4Allocator<G4PrimaryParticle>;
template class G4Allocator<G4PrimaryVertex>;
template class G4Allocator<G4Run>;
template class G4Allocator<G4StackedTrack>;
template class G4Allocator<G4Track>;
template class G4Allocator<G4Trajectory>;
template class G4Allocator<G4TrajectoryPoint>;




