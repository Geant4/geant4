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
#include "G4Scene.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelProxy.hh"
#include "G4Track.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UImessenger.hh"
#include "G4UIparameter.hh"
#include "G4VDecayChannel.hh"
#include "G4VDigiCollection.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VHitsCollection.hh"
#include "G4VModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VSceneHandler.hh"
#include "G4VSolid.hh"
#include "G4VStateDependent.hh"
#include "G4VTrajectory.hh"
#include "G4VViewer.hh"
#include "G4UnitsTable.hh"
#include "G4UItokenNum.hh"
#include "G4Trajectory.hh"

class SoSeparator;

template class RWTPtrSlistDictionary<RWCString, G4MaterialPropertyVector>;
template class RWTPtrSlistDictionary<RWCString, G4ParticleDefinition>;
template class RWTPtrSlistDictionary<int, G4ParticleDefinition>;
template class RWTPtrSlistDictionaryIterator<RWCString, G4MaterialPropertyVector>;
template class RWTPtrSlistDictionaryIterator<RWCString, G4ParticleDefinition>;
template class RWTPtrSlistDictionaryIterator<int, G4ParticleDefinition>;
template class RWTValSlistDictionary<RWCString, void*>;
template class RWTValSlistDictionary<_WidgetRec*, RWCString>;
template class RWTValSlistDictionary<const G4ParticleDefinition*, G4EnergyLossTablesHelper>;
template class RWTValSlistDictionary<const G4VPhysicalVolume*, SoSeparator*>;
template class RWTValSlistDictionary<unsigned long, int>;
template class G4RWTValVector<yystype>;

template class G4RWTPtrHashDictionary<RWCString, G4MaterialPropertyVector>;
template class G4RWTPtrHashDictionary<RWCString, G4ParticleDefinition>;
template class G4RWTPtrHashDictionary<int, G4ParticleDefinition>;
template class G4RWTPtrHashDictionaryIterator<RWCString, G4MaterialPropertyVector>;
template class G4RWTPtrHashDictionaryIterator<RWCString, G4ParticleDefinition>;
template class G4RWTPtrHashDictionaryIterator<int, G4ParticleDefinition>;
template class G4RWTPtrOrderedVector<G4Element>;
template class G4RWTPtrOrderedVector<G4Event>;
template class G4RWTPtrOrderedVector<G4Isotope>;
template class G4RWTPtrOrderedVector<G4LogicalVolume>;
template class G4RWTPtrOrderedVector<G4MPVEntry>;
template class G4RWTPtrOrderedVector<G4Material>;
template class G4RWTPtrOrderedVector<G4ParticleDefinition>;
template class G4RWTPtrOrderedVector<G4PhysicsVector>;
template class G4RWTPtrOrderedVector<G4ProcTblElement>;
template class G4RWTPtrOrderedVector<G4ProcessAttribute>;
template class G4RWTPtrOrderedVector<G4ProcessManager>;
template class G4RWTPtrOrderedVector<G4SDStructure>;
template class G4RWTPtrOrderedVector<G4Scene>;
template class G4RWTPtrOrderedVector<G4SmartVoxelNode>;
template class G4RWTPtrOrderedVector<G4SmartVoxelProxy>;
template class G4RWTPtrOrderedVector<G4Track>;
template class G4RWTPtrOrderedVector<G4UIcommand>;
template class G4RWTPtrOrderedVector<G4UIcommandTree>;
template class G4RWTPtrOrderedVector<G4UImessenger>;
template class G4RWTPtrOrderedVector<G4UIparameter>;
template class G4RWTPtrOrderedVector<G4UnitDefinition>;
template class G4RWTPtrOrderedVector<G4UnitsCategory>;
template class G4RWTPtrOrderedVector<G4VDecayChannel>;
template class G4RWTPtrOrderedVector<G4VDigiCollection>;
template class G4RWTPtrOrderedVector<G4VGraphicsSystem>;
template class G4RWTPtrOrderedVector<G4VHitsCollection>;
template class G4RWTPtrOrderedVector<G4VModel>;
template class G4RWTPtrOrderedVector<G4VPhysicalVolume>;
template class G4RWTPtrOrderedVector<G4VProcess>;
template class G4RWTPtrOrderedVector<G4VSceneHandler>;
template class G4RWTPtrOrderedVector<G4VSensitiveDetector>;
template class G4RWTPtrOrderedVector<G4VSolid>;
template class G4RWTPtrOrderedVector<G4VStateDependent>;
template class G4RWTPtrOrderedVector<G4VTrajectory>;
template class G4RWTPtrOrderedVector<G4VViewer>;
template class G4RWTPtrOrderedVector<G4ValVector>;
template class G4RWTPtrSortedVector<G4MPVEntry>;
template class G4RWTPtrSortedVector<G4VDecayChannel>;
template class G4RWTValHashDictionary<RWCString, void*>;
template class G4RWTValHashDictionary<_WidgetRec*, RWCString>;
template class G4RWTValHashDictionary<const G4ParticleDefinition*, G4EnergyLossTablesHelper>;
template class G4RWTValHashDictionary<const G4VPhysicalVolume*, SoSeparator*>;
template class G4RWTValHashDictionary<unsigned long, int>;
template class G4RWTValOrderedVector<G4ApplicationState>;
template class G4RWTValOrderedVector<Hep3Vector>;
template class G4RWTValOrderedVector<HepPlane3D>;
template class G4RWTValOrderedVector<HepPoint3D>;
template class G4RWTValOrderedVector<RWCString>;
template class G4RWTValOrderedVector<double>;
template class G4RWTValOrderedVector<int (*)(void*)>;
template class G4RWTValOrderedVector<int>;
template class G4RWTValOrderedVector<void (*)(void)>;
template class G4RWTValOrderedVector<void*>;
template class G4RWTValOrderedVector<yystype>;
template class G4RWTValVector<G4ApplicationState>;
template class G4RWTValVector<G4NavigationLevel>;
template class G4RWTValVector<Hep3Vector>;
template class G4RWTValVector<HepPlane3D>;
template class G4RWTValVector<HepPoint3D>;
template class G4RWTValVector<RWCString>;
template class G4RWTValVector<double>;
template class G4RWTValVector<int (*)(void*)>;
template class G4RWTValVector<int>;
template class G4RWTValVector<void (*)(void)>;
template class G4RWTValVector<void*>;
template class G4RWTPtrOrderedVector<G4VTrajectoryPoint>;

//
// G4 :
#include "G4DecayProducts.hh"
#include "G4Run.hh"
#include "G4VisCommandTemplates.hh"
#include "G4VisCommandsCamera.hh"
#include "G4VisCommandsViewer.hh"
#include "G4VisCommandsScene.hh"
#include "G4VisCommandsClear.hh"
#include "G4VisCommandsCreateScene.hh"
#include "G4VisCommandsCreateView.hh"
#include "G4VisCommandsCopy.hh"
#include "G4VisCommandsDraw.hh"
#include "G4VisCommandsDelete.hh"
#include "G4VisCommandsSet.hh"
#include "G4VisCommandsLights.hh"
#include "G4VisCommandsPrint.hh"
#include "G4VisCommandsRefresh.hh"
#include "G4VisCommandsShow.hh"

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
template class G4VisButtonCommandMessenger<G4VisCommandLightsMoveWithCamera>;
template class G4VisButtonCommandMessenger<G4VisCommandSetCullCoveredDaughters>;
template class G4VisButtonCommandMessenger<G4VisCommandSetCullInvisible>;
template class G4VisButtonCommandMessenger<G4VisCommandSetCulling>;
template class G4VisCommandDirectoryMessenger<G4VisCommandCamera>;
template class G4VisCommandDirectoryMessenger<G4VisCommandClear>;
template class G4VisCommandDirectoryMessenger<G4VisCommandCopy>;
template class G4VisCommandDirectoryMessenger<G4VisCommandCreateScene>;
template class G4VisCommandDirectoryMessenger<G4VisCommandCreateView>;
template class G4VisCommandDirectoryMessenger<G4VisCommandDelete>;
template class G4VisCommandDirectoryMessenger<G4VisCommandDraw>;
template class G4VisCommandDirectoryMessenger<G4VisCommandLights>;
template class G4VisCommandDirectoryMessenger<G4VisCommandPrint>;
template class G4VisCommandDirectoryMessenger<G4VisCommandRefresh>;
template class G4VisCommandDirectoryMessenger<G4VisCommandSet>;
template class G4VisCommandDirectoryMessenger<G4VisCommandShow>;
template class G4VisSimpleCommandMessenger<G4VisCommandCameraReset>;
template class G4VisSimpleCommandMessenger<G4VisCommandClearView>;
template class G4VisSimpleCommandMessenger<G4VisCommandCopyView>;
template class G4VisSimpleCommandMessenger<G4VisCommandCreateSceneClear>;
template class G4VisSimpleCommandMessenger<G4VisCommandCreateViewNewScene>;
template class G4VisSimpleCommandMessenger<G4VisCommandCreateViewNewView>;
template class G4VisSimpleCommandMessenger<G4VisCommandDeleteScene>;
template class G4VisSimpleCommandMessenger<G4VisCommandDeleteView>;
template class G4VisSimpleCommandMessenger<G4VisCommandDrawCurrent>;
template class G4VisSimpleCommandMessenger<G4VisCommandRefreshView>;
template class G4VisSimpleCommandMessenger<G4VisCommandShowView>;


