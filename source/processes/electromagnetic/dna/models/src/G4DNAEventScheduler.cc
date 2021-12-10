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
//
#include <memory>
#include "G4DNAEventScheduler.hh"
#include "G4DNAGillespieDirectMethod.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4DNAUpdateSystemModel.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4Timer.hh"
#include "G4Scheduler.hh"
#include "G4UserMeshAction.hh"
#include "G4MoleculeCounter.hh"
#include "G4DNAScavengerMaterial.hh"

G4DNAEventScheduler::G4DNAEventScheduler(const G4DNABoundingBox& boundingBox,
                                         G4int pixel)
  : IEventScheduler()
  , fVerbose(0)
  , fInitialized(false)
  , fStartTime(1 * ps)
  , fEndTime(10000 * s)
  , fStepNumber(0)
  , fMaxStep(INT_MAX)
  , fRunning(true)
  , fTimeStep(DBL_MAX)
  , fGlobalTime(fStartTime)
  , fJumpingNumber(0)
  , fReactionNumber(0)
  , fPixel(pixel)
  , fIsChangeMesh(false)
  , fSetChangeMesh(true)
  , fStepNumberInMesh(0)
  , fInitialPixels(fPixel)
  , fpMesh(new G4DNAMesh(boundingBox, fPixel))
  , fpGillespieReaction(new G4DNAGillespieDirectMethod())
  , fpEventSet(new G4DNAEventSet())
  , fpUpdateSystem(new G4DNAUpdateSystemModel())
  , fpUserMeshAction(nullptr)
{
  if(!CheckingReactionRadius(fpMesh->GetResolution()))
  {
    G4String WarMessage = "resolution is not good : " +
                          std::to_string(fpMesh->GetResolution() / nm);
    G4Exception("G4DNAEventScheduler::InitializeInMesh()", "WrongResolution",
                JustWarning, WarMessage);
  }
}

void G4DNAEventScheduler::ClearAndReChargeCounter()
{
  fCounterMap.clear();
  if(fTimeToRecord.empty())
  {
    G4cout << "fTimeToRecord is empty " << G4endl;
  }
  fLastRecoredTime = fTimeToRecord.begin();

  if(G4VMoleculeCounter::Instance()->InUse())  // copy from MoleculeCounter
  {
    G4MoleculeCounter::RecordedMolecules species;
    species = G4MoleculeCounter::Instance()->GetRecordedMolecules();
    if(species.get() == nullptr)
    {
      return;
    }
    else if(species->empty())
    {
      G4MoleculeCounter::Instance()->ResetCounter();
      return;
    }
    for(auto time_mol : fTimeToRecord)
    {
      if(time_mol > fStartTime)
      {
        continue;
      }

      for(auto molecule : *species)
      {
        G4int n_mol = G4MoleculeCounter::Instance()->GetNMoleculesAtTime(
          molecule, time_mol);

        if(n_mol < 0)
        {
          G4cerr << "G4DNAEventScheduler::ClearAndReChargeCounter() ::N "
                    "molecules not valid < 0 "
                 << G4endl;
          G4Exception("", "N<0", FatalException, "");
        }
        fCounterMap[time_mol][molecule] = n_mol;
      }

      fLastRecoredTime++;
    }
    G4MoleculeCounter::Instance()->ResetCounter();  // reset
    G4MoleculeCounter::Instance()->Use(false);      // no more used
  }
}

[[maybe_unused]] void G4DNAEventScheduler::AddTimeToRecord(const G4double& time)
{
  if(fTimeToRecord.find(time) == fTimeToRecord.end())
  {
    fTimeToRecord.insert(time);
  }
}

G4DNAEventScheduler::~G4DNAEventScheduler() = default;

void G4DNAEventScheduler::Voxelizing()
{
  auto pMainList = G4ITTrackHolder::Instance()->GetMainList();
  std::map<G4DNAMesh::Key, MapList> TrackKeyMap;
  for(auto track : *pMainList)
  {
    auto molType = GetMolecule(track)->GetMolecularConfiguration();

    auto pScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial*>(
      G4Scheduler::Instance()->GetScavengerMaterial());
    if(pScavengerMaterial != nullptr &&
       pScavengerMaterial->find(molType))  // avoid voxelize the scavenger
    {
      continue;
    }

    auto key = fpMesh->GetKey(track->GetPosition());
    if(TrackKeyMap.find(key) != TrackKeyMap.end())
    {
      std::map<MolType, size_t>& TrackTypeMap = TrackKeyMap[key];
      if(TrackTypeMap.find(molType) != TrackTypeMap.end())
      {
        TrackTypeMap[molType]++;
      }
      else
      {
        TrackTypeMap[molType] = 1;
      }
    }
    else
    {
      TrackKeyMap[key][molType] = 1;
    }
  }

  for(auto& it : TrackKeyMap)
  {
    fpMesh->SetVoxelMapList(it.first, std::move(it.second));
  }
}

void G4DNAEventScheduler::ReVoxelizing(G4int pixel)
{
  fPixel       = pixel;
  auto newMesh = new G4DNAMesh(fpMesh->GetBoundingBox(), fPixel);

  auto begin = fpMesh->begin();
  auto end   = fpMesh->end();
  std::map<G4DNAMesh::Key, MapList> TrackKeyMap;
  for(; begin != end; begin++)
  {
    auto index  = fpMesh->GetIndex(begin->first);
    auto newKey = newMesh->GetKey(fpMesh->GetIndex(index, fPixel));
    auto node   = begin->second;
    // if (node == nullptr) continue;
    if(TrackKeyMap.find(newKey) == TrackKeyMap.end())
    {
      TrackKeyMap[newKey] = node->GetMapList();
    }
    else
    {
      for(const auto& it : node->GetMapList())
      {
        TrackKeyMap[newKey][it.first] += it.second;
      }
      if(fVerbose > 1)
      {
        G4cout << "key : " << begin->first << " index : " << index
               << " new index : " << fpMesh->GetIndex(index, fPixel)
               << " new key : " << newKey
               << " number: " << node->GetMapList().begin()->second << G4endl;
      }
    }
  }
  fpMesh.reset(newMesh);

  for(auto& it : TrackKeyMap)
  {
    fpMesh->SetVoxelMapList(it.first, std::move(it.second));
  }
}
void G4DNAEventScheduler::Reset()
{
  // find another solultion
  fGlobalTime = fEndTime;

  //
  // RecordTime();//Last register for counter

  if(fVerbose > 0)
  {
    G4cout << "End Processing and reset Gird, ScavengerTable, EventSet for new "
              "simulation!!!!"
           << G4endl;
  }
  fInitialized    = false;
  fTimeStep       = 0;
  fStepNumber     = 0;
  fGlobalTime     = fStartTime;
  fRunning        = true;
  fReactionNumber = 0;
  fJumpingNumber  = 0;

  fpEventSet->RemoveEventSet();
  fpMesh->Reset();
}

void G4DNAEventScheduler::Initialize()
{
  if(!fInitialized)
  {
    fPixel = fInitialPixels;
    fpMesh = std::make_unique<G4DNAMesh>(fpMesh->GetBoundingBox(), fPixel);

    // Scavenger();

    auto pScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial*>(
      G4Scheduler::Instance()->GetScavengerMaterial());
    if(pScavengerMaterial == nullptr)
    {
      G4cout << "pScavengerMaterial == nullptr" << G4endl;
    }
    else
    {
      if(fVerbose > 1)
      {
        pScavengerMaterial->PrintInfo();
      }
    }

    Voxelizing();
    fpGillespieReaction->SetVoxelMesh(*fpMesh);
    fpGillespieReaction->SetEventSet(fpEventSet.get());
    fpGillespieReaction->SetTimeStep(
      0);  // reset fTimeStep = 0 in fpGillespieReaction
    fpGillespieReaction->Initialize();
    fpUpdateSystem->SetMesh(fpMesh.get());
    ClearAndReChargeCounter();
    fInitialized = true;
  }

  if(fVerbose > 0)
  {
    fpUpdateSystem->SetVerbose(1);
  }

  if(fVerbose > 2)
  {
    fpMesh->PrintMesh();
  }
}
void G4DNAEventScheduler::InitializeInMesh()
{
  if(fPixel <= 1)
  {
    fRunning = false;
    return;
  }
  // TEST /3
  ReVoxelizing(fPixel / 2);  //
  // ReVoxelizing(fPixel/3);//

  fpGillespieReaction->SetVoxelMesh(*fpMesh);
  fpUpdateSystem->SetMesh(fpMesh.get());
  fpGillespieReaction->Initialize();
}

void G4DNAEventScheduler::ResetInMesh()
{
  if(fVerbose > 0)
  {
    G4cout
      << "*** End Processing In Mesh and reset Mesh, EventSet for new Mesh!!!!"
      << G4endl;
  }
  fpEventSet->RemoveEventSet();
  fInitialized      = false;
  fIsChangeMesh     = false;
  fReactionNumber   = 0;
  fJumpingNumber    = 0;
  fStepNumberInMesh = 0;
}

G4double G4DNAEventScheduler::GetStartTime() const { return fStartTime; }

G4double G4DNAEventScheduler::GetEndTime() const { return fEndTime; }

[[maybe_unused]] G4double G4DNAEventScheduler::GetTimeStep() const { return fTimeStep; }

G4int G4DNAEventScheduler::GetVerbose() const { return fVerbose; }

[[maybe_unused]] void G4DNAEventScheduler::SetMaxNbSteps(G4int max) { fMaxStep = max; }

[[maybe_unused]] void G4DNAEventScheduler::SetStartTime(G4double time)
{
  fStartTime = time;
}

void G4DNAEventScheduler::Stop() { fRunning = false; }
void G4DNAEventScheduler::Run()
{
  G4Timer localtimer;
  if(fVerbose > 0)
  {
    localtimer.Start();
    G4cout << "***G4DNAEventScheduler::Run*** for Pixel : " << fPixel << G4endl;
  }
  while(fEndTime > fGlobalTime && fRunning)
  {
    RunInMesh();
  }
  if(fVerbose > 0)
  {
    if(!fRunning)
    {
      G4cout << " StepNumber(" << fStepNumber << ") = MaxStep(" << fMaxStep
             << ")" << G4endl;
    }
    else if(fEndTime <= fGlobalTime)
    {
      G4cout << " GlobalTime(" << fGlobalTime << ") > EndTime(" << fEndTime
             << ")"
             << " StepNumber : " << fStepNumber << G4endl;
    }
    localtimer.Stop();
    G4cout << "***G4DNAEventScheduler::Ending::"
           << G4BestUnit(fGlobalTime, "Time")
           << " Events left : " << fpEventSet->size() << G4endl;
    if(fVerbose > 1) {
      fpMesh->PrintMesh();
    }
    G4cout << " Computing Time : " << localtimer << G4endl;
  }
  Reset();
}

void G4DNAEventScheduler::RunInMesh()
{
  if(!fInitialized)
  {
    InitializeInMesh();
  }
  G4Timer localtimerInMesh;
  // if (fVerbose > 0)
  {
    localtimerInMesh.Start();
    G4double C = 20;
    G4double D = G4MoleculeTable::Instance()
                 ->GetConfiguration("H2O2")
                 ->GetDiffusionCoefficient();
    G4double transferTime = std::pow(fpMesh->GetResolution(), 2) * C / (6 * D);
    G4cout << "***G4DNAEventScheduler::RunInMesh*** for Pixel : " << fPixel
           << "  transferTime : " << G4BestUnit(transferTime, "Time") << G4endl;
    G4cout << "  resolution : " << G4BestUnit(fpMesh->GetResolution(), "Length")
           << G4endl;
  }

  if(fVerbose > 2)
  {
    fpMesh->PrintMesh();
  }

  if(fpUserMeshAction != nullptr)
  {
    fpUserMeshAction->BeginOfMesh(fpMesh.get(), fGlobalTime);
  }

  // if diffusive jumping is avaiable, EventSet is never empty
  while(!fpEventSet->Empty() && !fIsChangeMesh && fEndTime > fGlobalTime)
  {
    Stepping();
    fGlobalTime = fTimeStep + fStartTime;

    if(fpUserMeshAction != nullptr)
    {
      fpUserMeshAction->InMesh(fpMesh.get(), fGlobalTime);
    }

    if(fVerbose > 2)
    {
      G4cout << "fGlobalTime : " << G4BestUnit(fGlobalTime, "Time")
             << " fTimeStep : " << G4BestUnit(fTimeStep, "Time") << G4endl;
    }

    G4double C = 20;
    G4double D = G4MoleculeTable::Instance()
                 ->GetConfiguration("H2O2")
                 ->GetDiffusionCoefficient();
    if(D == 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "D == 0";
      G4Exception("G4DNAEventScheduler::RunInMesh", "G4DNAEventScheduler001",
                  FatalErrorInArgument, exceptionDescription);
    }
    G4double transferTime = std::pow(fpMesh->GetResolution(), 2) * C / (6 * D);

    // if(fStepNumberInMesh > 40000 && fPixel != 1)
    if(transferTime < fTimeStep &&
       fPixel != 1)  // dont change Mesj if fPixel == 1
    {
      if(fVerbose > 1)
      {
        G4cout << " Pixels : " << fPixel << "  resolution : "
               << G4BestUnit(fpMesh->GetResolution(), "Length")
               << "  fStepNumberInMesh : " << fStepNumberInMesh
               << " at fGlobalTime : " << G4BestUnit(fGlobalTime, "Time")
               << " at fTimeStep : " << G4BestUnit(fTimeStep, "Time")
               << "  fReactionNumber : " << fReactionNumber
               << " transferTime : " << G4BestUnit(transferTime, "Time")
               << G4endl;
      }
      if(fSetChangeMesh)
      {
        fIsChangeMesh = true;
      }
    }
  }

  if(fVerbose > 1)
  {
    localtimerInMesh.Stop();
    G4cout << "***G4DNAEventScheduler::Ending::"
           << G4BestUnit(fGlobalTime, "Time")
           << " Event left : " << fpEventSet->size() << G4endl;
    G4cout << " Computing Time : " << localtimerInMesh << " Due to : ";
    if(fpEventSet->Empty())
    {
      G4cout << "EventSet is Empty" << G4endl;
    }
    else if(fIsChangeMesh)
    {
      G4cout << "Changing Mesh from : " << fPixel
             << " pixels to : " << fPixel / 2 << " pixels" << G4endl;
      G4cout << "Info : ReactionNumber : " << fReactionNumber
             << "   JumpingNumber : " << fJumpingNumber << G4endl;
    }
    else if(fEndTime > fGlobalTime)
    {
      G4cout << " GlobalTime(" << fGlobalTime << ") > EndTime(" << fEndTime
             << ")"
             << " StepNumber : " << fStepNumber << G4endl;
    }
    if(fVerbose > 2)
    {
      fpMesh->PrintMesh();
    }
    G4cout << G4endl;
  }

  if(fpUserMeshAction != nullptr)
  {
    fpUserMeshAction->EndOfMesh(fpMesh.get(), fGlobalTime);
  }
  ResetInMesh();
}

void G4DNAEventScheduler::Stepping()  // this event loop
{
  fStepNumber < fMaxStep ? fStepNumber++ : fRunning = false;

  if(fpEventSet->size() > fpMesh->size())
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "fpEventSet->size() > fpMesh->size()";
    G4Exception("G4DNAEventScheduler::Stepping", "G4DNAEventScheduler002",
                FatalErrorInArgument, exceptionDescription);
  };

  auto selected   = fpEventSet->begin();
  const auto& key = (*selected)->GetKey();
  auto index      = fpMesh->GetIndex(key);

  if(fVerbose > 1)
  {
    G4cout << "G4DNAEventScheduler::Stepping()*********************************"
              "*******"
           << G4endl;
    (*selected)->PrintEvent();
  }

  // get selected time step
  fTimeStep = (*selected)->GetTime();

  // selected data
  auto pJumping  = (*selected)->GetJumpingData();
  auto pReaction = (*selected)->GetReactionData();

  fpUpdateSystem->SetGlobalTime(fTimeStep +
                                fStartTime);  // this is just for printing

  if(pJumping == nullptr && pReaction == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "pJumping == nullptr && pReaction == nullptr";
    G4Exception("G4DNAEventScheduler::Stepping", "G4DNAEventScheduler003",
                FatalErrorInArgument, exceptionDescription);
  }

  fpGillespieReaction->SetTimeStep(fTimeStep);

  if(pJumping == nullptr)
  {
    fpUpdateSystem->UpdateSystem(index, *pReaction);

    fpEventSet->RemoveEvent(selected);
    // create new event
    fpGillespieReaction->CreateEvent(key);
    fReactionNumber++;

    // recordTime in reaction
    RecordTime();
  }
  else if(pReaction == nullptr)
  {
    // dont change this
    fpUpdateSystem->UpdateSystem(index, *pJumping);
    auto jumpingKey = fpMesh->GetKey(pJumping->second);
    fpEventSet->RemoveEvent(selected);

    // create new event
    // should create Jumping before key
    fpGillespieReaction->CreateEvent(jumpingKey);
    fpGillespieReaction->CreateEvent(key);

    fJumpingNumber++;
  }
  if(fVerbose > 1)
  {
    G4cout << "G4DNAEventScheduler::Stepping::end "
              "Print***********************************"
           << G4endl;
    G4cout << G4endl;
  }
  fStepNumberInMesh++;
}

void G4DNAEventScheduler::SetEndTime(const G4double& endTime)
{
  fEndTime = endTime;
}

void G4DNAEventScheduler::RecordTime()
{
  auto recordTime = *fLastRecoredTime;
  if(fGlobalTime >= recordTime && fCounterMap[recordTime].empty())
  {
    auto begin = fpMesh->begin();
    auto end   = fpMesh->end();
    for(; begin != end; begin++)
    {
      auto node = begin->second;
      if(node == nullptr) {
        continue;
      }
      for(const auto& it : node->GetMapList())
      {
        fCounterMap[recordTime][it.first] += it.second;
      }
    }
    fLastRecoredTime++;

#ifdef DEBUG
    PrintRecordTime();
    G4MoleculeTable* pMoleculeTable = G4MoleculeTable::Instance();
    auto iter = pMoleculeTable->GetConfigurationIterator();
    iter.reset();
    while(iter())
    {
      auto conf = iter.value();

      G4cout << "GlobalTime : " << G4BestUnit(fGlobalTime, "Time")
             << "  recordTime : " << G4BestUnit(recordTime, "Time") << "  "
             << conf->GetName()
             << "  number : " << fCounterMap[recordTime][conf]
             << "  MoleculeCounter : "
             << G4MoleculeCounter::Instance()->GetCurrentNumberOf(conf)
             << G4endl;

      assert(G4MoleculeCounter::Instance()->GetCurrentNumberOf(conf) ==
             fCounterMap[recordTime][conf]);
    }
#endif
  }
}

void G4DNAEventScheduler::PrintRecordTime()
{
  G4cout << "fCounterMap.size : " << fCounterMap.size() << G4endl;

  for(const auto& i : fCounterMap)
  {
    auto map   = i.second;
    auto begin = map.begin();  //
    auto end   = map.end();    //
    for(; begin != end; begin++)
    {
      auto molecule = begin->first;
      auto number   = begin->second;
      if(number == 0)
      {
        continue;
      }
      G4cout << "molecule : " << molecule->GetName() << " number : " << number
             << G4endl;
    }
  }
}

std::map<G4double /*time*/, G4DNAEventScheduler::MapCounter>
G4DNAEventScheduler::GetCounterMap() const
{
  return fCounterMap;
}

void G4DNAEventScheduler::SetUserMeshAction(
  std::unique_ptr<G4UserMeshAction> pUserMeshAction)
{
  fpUserMeshAction = std::move(pUserMeshAction);
}

G4DNAMesh* G4DNAEventScheduler::GetMesh() const { return fpMesh.get(); }

G4int G4DNAEventScheduler::GetPixels() const { return fPixel; }

G4bool G4DNAEventScheduler::CheckingReactionRadius(G4double resolution)
{
  auto pMolecularReactionTable = G4DNAMolecularReactionTable::Instance();
  auto reactionDataList = pMolecularReactionTable->GetVectorOfReactionData();
  if(reactionDataList.empty())
  {
    G4cout << "reactionDataList.empty()" << G4endl;
    return true;
  }
  else
  {
    for(auto it : reactionDataList)
    {
      if(it->GetEffectiveReactionRadius() >= resolution / CLHEP::pi)
      {
        G4cout << it->GetReactant1()->GetName() << " + "
               << it->GetReactant2()->GetName() << G4endl;
        G4cout << "G4DNAEventScheduler::ReactionRadius : "
               << G4BestUnit(it->GetEffectiveReactionRadius(), "Length")
               << G4endl;
        G4cout << "resolution : " << G4BestUnit(resolution, "Length") << G4endl;
        return false;
      }
    }
    return true;
  }
}