//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Author: Christian Velten (2025)

#include "G4MoleculeCounterManager.hh"

#include "G4AutoLock.hh"
#include "G4MoleculeCounterManagerMessenger.hh"
#include "G4MoleculeCounterTemplates.hh"
#include "G4StateManager.hh"
#include "G4Step.hh"
#include "G4Threading.hh"

//------------------------------------------------------------------------------

static G4Mutex managerInstance;
static G4Mutex masterInstanceMutex;
static G4Mutex workerInstancesMutex;

const G4MoleculeCounterManager* G4MoleculeCounterManager::fpMasterInstance = nullptr;
std::vector<const G4MoleculeCounterManager*> G4MoleculeCounterManager::fWorkerInstances = {};
G4ThreadLocal std::unique_ptr<G4MoleculeCounterManager> G4MoleculeCounterManager::fpInstance = nullptr;

//------------------------------------------------------------------------------

G4ThreadLocal std::atomic<G4bool> G4MoleculeCounterManager::fBeginOfEventTriggered(false);
std::atomic<G4bool> G4MoleculeCounterManager::fBeginOfRunTriggered(false);

//------------------------------------------------------------------------------

//
// CTOR & INSTANCE
//

G4MoleculeCounterManager::G4MoleculeCounterManager(G4MoleculeCounterManager::Private)
  : fVerbosity(0), fIsInitialized(false), fIsActive(true)
{
  fpMessenger = std::make_unique<G4MoleculeCounterManagerMessenger>(this);
}

G4MoleculeCounterManager::~G4MoleculeCounterManager()
{
  // The manager owns all the counters as raw pointers, clean up:
  DeregisterAllCounters();

  if (fVerbosity > 0) {
    if (GetResetCountersBeforeRun() && !fBeginOfRunTriggered.load()) {
      G4Exception("G4MoleculeCounterManager::~G4MoleculeCounterManager", "MOLMAN000", JustWarning,
                  "The molecule counter manager was configured to reset counters before each run"
                  " but the BeginOfRunAction was never triggered!\n"
                  "Ensure that the user code calls either the G4DNAChemistryManager's or "
                  "G4MoleculeCounterManager's Run methods!");
    }
    if (GetResetCountersBeforeEvent() && !fBeginOfEventTriggered.load()
        && (!G4Threading::IsMultithreadedApplication() || G4Threading::IsWorkerThread()))
    {
      // ignore, if this is the master of a MT application
      G4Exception("G4MoleculeCounterManager::~G4MoleculeCounterManager", "MOLMAN000", JustWarning,
                  "The molecule counter manager was configured to reset counters before each event"
                  " but the BeginOfEventAction was never triggered!\n"
                  "This can occurr if this thread is never processing an event but:\n"
                  "Ensure that the user code calls either the G4DNAChemistryManager's or "
                  "G4MoleculeCounterManager's Event methods!");
    }
  }

  // Handle the recorded master and worker instances
  if (G4Threading::IsMasterThread()) {
    // This is the master's dtor: set the master const* to nullptr (as it will point to nothing of use)
    G4AutoLock lockMaster(&masterInstanceMutex);
    fpMasterInstance = nullptr;
  }
  else {
    // This is a worker's dtor: make sure that this worker is still in the list of instances.
    // If not, throw a warning (something fishy going on?)
    // Either way, set this worker's const* to nullptr in the list of instances.
    G4AutoLock lock(&workerInstancesMutex);
    auto it = std::find(fWorkerInstances.begin(), fWorkerInstances.end(), this);
    if (it == fWorkerInstances.end()) {
      G4Exception(
        "G4MoleculeCounterManager::~G4MoleculeCounterManager", "MOLMAN_DTOR", JustWarning,
        "The destroyed instance of G4MoleculeCounterManager has not been registered as a worker!");
    }
    else {
      (*it) = nullptr;
    }
  }
}

void G4MoleculeCounterManager::RegisterInstance()
{
  if (fInstancesRegistered) {
    G4Exception("G4MoleculeCounterManager::RegisterInstance", "MOLMAN000", FatalException,
                "Instances were already registered once!");
  }
  else {
    if (G4Threading::IsMasterThread()) {
      G4AutoLock lock(&masterInstanceMutex);
      if (fpMasterInstance != nullptr) {
        G4Exception("G4MoleculeCounterManager::RegisterInstance", "MOLMAN000", FatalException,
                    "Master instance was set already!");
      }
      fpMasterInstance = Instance();
    }
    else {
      G4AutoLock lock(&workerInstancesMutex);
      fWorkerInstances.push_back(Instance());
    }
    fInstancesRegistered = true;
  }
}

G4MoleculeCounterManager* G4MoleculeCounterManager::Instance()
{
  if (fpInstance == nullptr) {
    G4AutoLock lock(&managerInstance);
    fpInstance = std::make_unique<G4MoleculeCounterManager>(G4MoleculeCounterManager::Private());
    fpInstance->RegisterInstance();
  }
  return fpInstance.get();
}

G4MoleculeCounterManager* G4MoleculeCounterManager::GetInstanceIfExists()
{
  return fpInstance.get();
}

void G4MoleculeCounterManager::DeleteInstance()
{
  G4AutoLock lock(&managerInstance);

  if (fpInstance != nullptr) {
    fpInstance.reset();
    // this should (test!) trigger the dtor, which will delete/deregister the counters
  }
}

// ---------------------------------------------------------------------

//
// Initialization
//

void G4MoleculeCounterManager::Initialize()
{
  if (fVerbosity > 0) {
    G4cout << "G4MoleculeCounterManager::Initialize ("
           << (G4Threading::IsMasterThread() ? "master" : "worker") << ")" << G4endl;
  }

  if (G4Threading::IsMultithreadedApplication()) {
    if (G4Threading::IsWorkerThread())
      InitializeWorker();
    else
      InitializeMaster();
  }
  else {
    InitializeMaster();
  }
}

void G4MoleculeCounterManager::InitializeMaster()
{
  if (fIsInitialized) return;
  for (auto& counter : fCounters)
    counter.second->Initialize();
}

void G4MoleculeCounterManager::InitializeWorker()
{
  if (fIsInitialized) return;
  for (auto& counter : fCounters)
    counter.second->Initialize();
}

// ---------------------------------------------------------------------

//
// Add & Remove Molecules
//

void G4MoleculeCounterManager::AddMoleculeWithoutTrack(const G4MolecularConfiguration* molecule,
                                                       G4double time, G4int n)
{
  for (auto& [id, counter] : fCounters) {
    counter->AddMolecule(counter->BuildSimpleIndex(molecule), time, n);
  }
}

void G4MoleculeCounterManager::RemoveMoleculeWithoutTrack(const G4MolecularConfiguration* molecule,
                                                          G4double time, G4int n)
{
  for (auto& [id, counter] : fCounters) {
    counter->RemoveMolecule(counter->BuildSimpleIndex(molecule), time, n);
  }
}

void G4MoleculeCounterManager::AddMolecule(const G4Track* aTrack, G4double time, G4int n)
{
  for (auto& [id, counter] : fCounters) {
    counter->AddMolecule(counter->BuildIndex(aTrack), time, n);
  }
}

void G4MoleculeCounterManager::RemoveMolecule(const G4Track* aTrack, G4double time, G4int n)
{
  for (auto& [id, counter] : fCounters) {
    counter->RemoveMolecule(counter->BuildIndex(aTrack), time, n);
  }
}

void G4MoleculeCounterManager::AddMolecule(const G4Track* aTrack, const G4StepPoint* aStepPoint,
                                           G4double time, G4int n)
{
  for (auto& [id, counter] : fCounters) {
    if (counter->GetSensitiveToStepping()) {
      counter->AddMolecule(counter->BuildIndex(aTrack, aStepPoint), time, n);
    }
  }
}

void G4MoleculeCounterManager::RemoveMolecule(const G4Track* aTrack, const G4StepPoint* aStepPoint,
                                              G4double time, G4int n)
{
  for (auto& [id, counter] : fCounters) {
    if (counter->GetSensitiveToStepping()) {
      counter->RemoveMolecule(counter->BuildIndex(aTrack, aStepPoint), time, n);
    }
  }
}

// ---------------------------------------------------------------------

void G4MoleculeCounterManager::RecordReaction(const G4DNAMolecularReactionData* reactionData,
                                              G4double time, G4int n)
{
  for (auto& [id, counter] : fReactionCounters) {
    counter->RecordReaction(counter->BuildSimpleIndex(reactionData), time, n);
  }
}

// ---------------------------------------------------------------------

//
// Notification of Step or Finalization (must be implemented in a custom G4ITTrackingInteractivity!
//

/* Example:
 void ChemistrySteppingAction::UserSteppingAction(const G4Step* aStep)
 {
   if (G4MoleculeCounterManager::Instance()->GetIsActive())
     G4MoleculeCounterManager::Instance()->NotifyOfStep(aStep);
 }
 */

void G4MoleculeCounterManager::NotifyOfStep(const G4Step* aStep)
{
  const G4Track* aTrack = aStep->GetTrack();

  for (auto& [id, counter] : fCounters) {
    auto preStepIndex = counter->BuildIndex(aTrack, aStep->GetPreStepPoint());
    auto postStepIndex = counter->BuildIndex(aTrack, aStep->GetPostStepPoint());

    if (!(*preStepIndex == *postStepIndex)) {
#ifdef G4VERBOSE
      if (GetVerbosity() > 1) {
        G4cout << "G4MoleculeCounterManager::NotifyOfStep for counter " << counter->GetName()
               << ":\n-- Pre = " << preStepIndex->GetInfo() << "\n"
               << "-- Post = " << postStepIndex->GetInfo() << G4endl;
      }
#endif
      counter->RemoveMolecule(std::move(preStepIndex), aTrack->GetGlobalTime(), 1);
      counter->AddMolecule(std::move(postStepIndex), aTrack->GetGlobalTime(), 1);
    }
  }
}

void G4MoleculeCounterManager::NotifyOfFinalize()
{
  for (auto& [id, counter] : fCounters) {
    counter->SchedulerFinalizedTracking();
  }
}

//
// Broadcast Functions
//

void G4MoleculeCounterManager::BroadcastIgnoreMolecule(const G4MoleculeDefinition* molecule)
{
  for (auto& [id, counter] : fCounters) {
    counter->IgnoreMolecule(molecule);
  }
}

void G4MoleculeCounterManager::BroadcastIgnoreReactant(const G4MolecularConfiguration* molecule)
{
  for (auto& [id, counter] : fCounters) {
    counter->IgnoreReactant(molecule);
  }
}

void G4MoleculeCounterManager::BroadcastRegisterAllMoleculesAndReactants()
{
  for (auto& [id, counter] : fCounters) {
    counter->RegisterAll();
  }
}

//------------------------------------------------------------------------------

//
// Manager Functions
//

G4int G4MoleculeCounterManager::RegisterCounter(std::unique_ptr<G4VMoleculeCounter> counter)
{
  auto idProvider = [&]() {
    if (fCounters.size() == 0) return 0;
    auto indices = G4::MoleculeCounter::GetMapIndices(fCounters);
    auto lastIndex = *indices.rbegin();
    return ++lastIndex;
  };

  return RegisterCounter(fCounters, std::move(counter), idProvider);
}

//------------------------------------------------------------------------------

G4int G4MoleculeCounterManager::RegisterCounter(std::unique_ptr<G4VMoleculeReactionCounter> counter)
{
  auto idProvider = [&] {
    if (fReactionCounters.size() == 0) return 0;
    auto indices = G4::MoleculeCounter::GetMapIndices(fReactionCounters);
    auto lastIndex = *indices.rbegin();
    return ++lastIndex;
  };

  return RegisterCounter(fReactionCounters, std::move(counter), idProvider);
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::DeregisterAllCounters()
{
  for (auto& [id, ctr] : fCounters) {
    delete ctr;
  }
  fCounters.clear();

  for (auto& [id, ctr] : fReactionCounters) {
    delete ctr;
  }
  fReactionCounters.clear();
}

//------------------------------------------------------------------------------

//
// Manipulate Counters
//

void G4MoleculeCounterManager::ResetCounters()
{
  if (fVerbosity > 0) {
    G4cout << "G4MoleculeCounterManager::ResetCounters ("
           << (G4Threading::IsMasterThread() ? "master" : "worker") << ")" << G4endl;
  }

  for (auto& [id, counter] : fCounters)
    counter->ResetCounter();
  for (auto& [id, counter] : fReactionCounters)
    counter->ResetCounter();
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::ActivateCounterAtTimes(G4int id, G4double aboveTime,
                                                      G4double belowTime, G4bool aboveTimeInclusive,
                                                      G4bool belowTimeInclusive)
{
  if (fVerbosity > 0) {
    G4cout << "G4MoleculeCounterManager::ActivateCounterAtTimes ("
           << (G4Threading::IsMasterThread() ? "master" : "worker") << ")" << G4endl;
  }

  auto counter = GetEditableMoleculeCounter(id);
  counter->SetActiveLowerBound(aboveTime, aboveTimeInclusive);
  counter->SetActiveUpperBound(belowTime, belowTimeInclusive);
}

void G4MoleculeCounterManager::ActivateReactionCounterAtTimes(G4int id, G4double aboveTime,
                                                              G4double belowTime,
                                                              G4bool aboveTimeInclusive,
                                                              G4bool belowTimeInclusive)
{
  if (fVerbosity > 0) {
    G4cout << "G4MoleculeCounterManager::ActivateReactionCounterAtTimes ("
           << (G4Threading::IsMasterThread() ? "master" : "worker") << ")" << G4endl;
  }

  auto counter = GetEditableMoleculeReactionCounter(id);
  counter->SetActiveLowerBound(aboveTime, aboveTimeInclusive);
  counter->SetActiveUpperBound(belowTime, belowTimeInclusive);
}

//------------------------------------------------------------------------------

G4VMoleculeCounter* G4MoleculeCounterManager::GetEditableMoleculeCounter(G4int id) const
{
  auto it = fCounters.find(id);
  if (it == fCounters.end()) {
    G4ExceptionDescription description;
    description << "No molecule counter with Id = " << id << " was found!\n";
    G4Exception("G4MoleculeCounterManager::GetMoleculeCounter", "MOLMAN001", FatalErrorInArgument,
                description);
    return nullptr;
  }else{
    return it->second;
  }
}

std::vector<const G4VMoleculeCounter*> G4MoleculeCounterManager::GetMoleculeCounters() const
{
  std::vector<const G4VMoleculeCounter*> output;
  output.reserve(fCounters.size());
  for (auto& [id, counter] : fCounters)
    output.push_back(counter);
  return output;
}

std::vector<const G4VMoleculeCounter*> G4MoleculeCounterManager::GetMoleculeCounters(G4String name) const
{
  std::vector<const G4VMoleculeCounter*> output;
  for (auto& [id, counter] : fCounters) {
    if (name == counter->GetName()) output.push_back(counter);
  }
  return output;
}

//------------------------------------------------------------------------------

G4VMoleculeReactionCounter*
G4MoleculeCounterManager::GetEditableMoleculeReactionCounter(G4int id) const
{
  auto it = fReactionCounters.find(id);
  if (it == fReactionCounters.end()) {
    G4ExceptionDescription description;
    description << "No molecule reaction counter with Id = " << id << " was found!\n";
    G4Exception("G4MoleculeCounterManager::GetMoleculeReactionCounter", "MOLMAN001",
                FatalErrorInArgument, description);
    return nullptr;
  }else
  {
    return it->second;
  }
}

std::vector<const G4VMoleculeReactionCounter*>
G4MoleculeCounterManager::GetMoleculeReactionCounters() const
{
  std::vector<const G4VMoleculeReactionCounter*> output;
  output.reserve(fReactionCounters.size());
  for (auto& [id, counter] : fReactionCounters)
    output.push_back(counter);
  return output;
}

std::vector<const G4VMoleculeReactionCounter*>
G4MoleculeCounterManager::GetMoleculeReactionCounters(G4String name) const
{
  std::vector<const G4VMoleculeReactionCounter*> output;
  for (auto& [id, counter] : fReactionCounters) {
    if (name == counter->GetName()) output.push_back(counter);
  }
  return output;
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::BeginOfEventAction(const G4Event*)
{
  fBeginOfEventTriggered = true;

  if (GetResetCountersBeforeEvent()) {
    // trigger reset if:
    // * is not an MT app (master = worker)
    // * is an MT app, master, and we want to reset the master
    // * is an MT app and worker
    if (!G4Threading::IsMultithreadedApplication()
        || (G4Threading::IsMultithreadedApplication() && G4Threading::IsMasterThread()
            && GetResetMasterCounterWithWorkers())
        || (G4Threading::IsMultithreadedApplication() && G4Threading::IsWorkerThread()))
      ResetCounters();
  }
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::BeginOfRunAction(const G4Run*)
{
  fBeginOfRunTriggered = true;

  if (GetResetCountersBeforeRun()) {
    // trigger reset if:
    // * is not an MT app (master = worker)
    // * is an MT app, master, and we want to reset the master
    // * is an MT app and worker
    if (!G4Threading::IsMultithreadedApplication()
        || (G4Threading::IsMultithreadedApplication() && G4Threading::IsMasterThread()
            && GetResetMasterCounterWithWorkers())
        || (G4Threading::IsMultithreadedApplication() && G4Threading::IsWorkerThread()))
      ResetCounters();
  }
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::EndOfEventAction(const G4Event*)
{
  // EndOfEvent is never triggered on the master of a G4MT
  if (GetAccumulateCounterIntoMaster() && GetResetCountersBeforeEvent()
      && G4Threading::IsMultithreadedApplication() && G4Threading::IsWorkerThread())
  {
    AbsorbWorkerManagerCounters(this);
  }
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::EndOfRunAction(const G4Run*)
{
  if (GetAccumulateCounterIntoMaster() && GetResetCountersBeforeRun()
      && !GetResetCountersBeforeEvent() &&
      // if ResetBeforeEvent, AbsorbWorkerManagerCounters will have been triggered already
      G4Threading::IsMultithreadedApplication() && G4Threading::IsMasterThread())
    AbsorbWorkerManagerCounters();
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::AbsorbWorkerManagerCounters(
  const G4MoleculeCounterManager* selectedWorker)
{
  if (selectedWorker == nullptr && !G4Threading::IsMasterThread()) {
    // Can only call without worker from master thread!
    G4ExceptionDescription description;
    description << "This method may only be called from the master thread!";
    G4Exception("G4MoleculeCounterManager::AbsorbWorkerManagerCounters", "MOLMAN999",
                FatalException, description);
  }

  // prevent changes to any of the instances
  G4AutoLock lockMaster(&masterInstanceMutex);
  G4AutoLock lockWorker(&workerInstancesMutex);

  for (auto const& worker : fWorkerInstances) {
    if (selectedWorker == nullptr || worker != selectedWorker) continue;

    for (auto& [id, masterCounter] : fpMasterInstance->fCounters) {
      // acquire worker counter
      auto workerCounter = worker->GetMoleculeCounter(id);
      masterCounter->AbsorbCounter(workerCounter);
    }

    for (auto& [id, masterCounter] : fpMasterInstance->fReactionCounters) {
      // acquire worker counter
      auto workerCounter = worker->GetMoleculeReactionCounter(id);
      masterCounter->AbsorbCounter(workerCounter);
    }
  }
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::DumpMasterCounters() const
{
  G4AutoLock lock(&masterInstanceMutex);

  for (auto const& pCounter : fpMasterInstance->GetMoleculeCounters()) {
    G4cout << "===========================================================================  \n"
           << " >> [MASTER] Dumping Molecule Counter `" << pCounter->GetName() << "`\n"
           << G4endl;
    pCounter->Dump();
    G4cout << "\n===========================================================================  "
           << G4endl;
  }

  for (auto const& pCounter : fpMasterInstance->GetMoleculeReactionCounters()) {
    G4cout << "===========================================================================  \n"
           << " >> [MASTER] Dumping Molecule Reaction Counter `" << pCounter->GetName() << "`\n"
           << G4endl;
    pCounter->Dump();
    G4cout << "\n===========================================================================  "
           << G4endl;
  }
}

//------------------------------------------------------------------------------

void G4MoleculeCounterManager::DumpWorkerCounters() const
{
  G4AutoLock lock(&workerInstancesMutex);

  for (auto const& worker : G4MoleculeCounterManager::fWorkerInstances) {
    for (auto const& pCounter : worker->GetMoleculeCounters()) {
      G4cout << "===========================================================================  \n"
             << " >> [WORKER<" << worker << ">] Dumping Molecule Counter `" << pCounter->GetName()
             << "`\n"
             << G4endl;
      pCounter->Dump();
      G4cout << "\n===========================================================================  "
             << G4endl;
    }

    for (auto const& pCounter : worker->GetMoleculeReactionCounters()) {
      G4cout << "===========================================================================  \n"
             << " >> [WORKER<" << worker << ">] Dumping Molecule Reaction Counter `"
             << pCounter->GetName() << "`\n"
             << G4endl;
      pCounter->Dump();
      G4cout << "\n===========================================================================  "
             << G4endl;
    }
  }
}

//------------------------------------------------------------------------------

std::atomic<G4bool> G4MoleculeCounterManager::fResetCountersBeforeEvent(false);
std::atomic<G4bool> G4MoleculeCounterManager::fResetCountersBeforeRun(false);
std::atomic<G4bool> G4MoleculeCounterManager::fAccumulateCounterIntoMaster(true);
std::atomic<G4bool> G4MoleculeCounterManager::fResetMasterCounterWithWorkers(false);

G4bool G4MoleculeCounterManager::GetResetCountersBeforeEvent() const
{
  return fResetCountersBeforeEvent.load();
}
void G4MoleculeCounterManager::SetResetCountersBeforeEvent(G4bool flag)
{
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit) {
    fResetCountersBeforeEvent = flag;
  }
  else
    G4Exception("G4DNAChemistryManager::SetResetCountersBeforeEvent", "WRONG_STATE", FatalException,
                "This flag may only be set during the PreInit state!");
}

G4bool G4MoleculeCounterManager::GetResetCountersBeforeRun() const
{
  return fResetCountersBeforeRun.load();
}
void G4MoleculeCounterManager::SetResetCountersBeforeRun(G4bool flag)
{
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit)
    fResetCountersBeforeRun = flag;
  else
    G4Exception("G4DNAChemistryManager::SetResetCountersBeforeRun", "WRONG_STATE", FatalException,
                "This flag may only be set during the PreInit state!");
}

G4bool G4MoleculeCounterManager::GetAccumulateCounterIntoMaster() const
{
  return fAccumulateCounterIntoMaster.load();
}
void G4MoleculeCounterManager::SetAccumulateCounterIntoMaster(G4bool flag)
{
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit)
    fAccumulateCounterIntoMaster = flag;
  else
    G4Exception("G4DNAChemistryManager::SetAccumulateCounterIntoMaster", "WRONG_STATE",
                FatalException, "This flag may only be set during the PreInit state!");
}

G4bool G4MoleculeCounterManager::GetResetMasterCounterWithWorkers() const
{
  return fResetMasterCounterWithWorkers.load();
}
void G4MoleculeCounterManager::SetResetMasterCounterWithWorkers(G4bool flag)
{
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit)
    fResetMasterCounterWithWorkers = flag;
  else
    G4Exception("G4DNAChemistryManager::SetResetMasterCounterWithWorkers", "WRONG_STATE",
                FatalException, "This flag may only be set during the PreInit state!");
}

//------------------------------------------------------------------------------
