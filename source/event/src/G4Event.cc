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
// G4Event class implementation
//
// Author: M.Asai, SLAC/JLAB
// --------------------------------------------------------------------

#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4VHitsCollection.hh"
#include "G4VDigiCollection.hh"
#include "G4ios.hh"
#include "G4SubEvent.hh"

G4Allocator<G4Event>*& anEventAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4Event>* _instance = nullptr;
  return _instance;
}

G4Event::G4Event(G4int evID)
  : eventID(evID)
{
}

G4Event::~G4Event()
{
  G4PrimaryVertex* nextVertex = thePrimaryVertex;
  while(nextVertex != nullptr)
  {
    G4PrimaryVertex* thisVertex = nextVertex;
    nextVertex = thisVertex->GetNext();
    thisVertex->ClearNext();
    delete thisVertex;
  }
  thePrimaryVertex = nullptr;
  delete HC;
  delete DC;
  if(trajectoryContainer != nullptr)
  {
    trajectoryContainer->clearAndDestroy();
    delete trajectoryContainer;
  }
  delete userInfo;
  delete randomNumberStatus;
  delete randomNumberStatusForProcessing;

  // Following G4Exception are temporally issuing JustWarning to delete 
  // unprocessed tracks in sub-events. Once sub-event mechanism is completely 
  // implemented, G4Exception should cause FatalException.

  G4int remainingSE = 0;
  for(auto& sem : fSubEvtStackMap)
  {
    if((sem.second!=nullptr)&&!(sem.second->empty()))
    {
      remainingSE += sem.second->size();
      for(auto& se : *(sem.second))
      {
        se->clearAndDestroy();
      }
      sem.second->clear();
    }
  }
  if(remainingSE>0)
  {
    G4ExceptionDescription ed;
    ed << "Deleting G4Event (id:" << eventID << ") that still has "
       << remainingSE << " sub-events un-processed.";
    G4Exception("G4Event::~G4Event()","SubEvt0001",JustWarning,ed);
  }

  if(!(fSubEvtVector.empty()))
  {
    G4ExceptionDescription ed;
    ed << "Deleting G4Event (id:" << eventID << ") that has "
       << fSubEvtVector.size() << " sub-events still processing.";
    G4Exception("G4Event::~G4Event()","SubEvt0001",JustWarning,ed);
    for(auto& se : fSubEvtVector)
    {
      se->clearAndDestroy();
      delete se;
    }
  }
}

G4bool G4Event::operator==(const G4Event& right) const
{
  return ( eventID == right.eventID );
}

G4bool G4Event::operator!=(const G4Event& right) const
{
  return ( eventID != right.eventID );
}

void G4Event::Print() const
{
  G4cout << "G4Event " << eventID << G4endl;
}

void G4Event::Draw() const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager == nullptr) return;

  if(trajectoryContainer != nullptr)
  {
    std::size_t n_traj = trajectoryContainer->entries();
    for(std::size_t i=0; i<n_traj; ++i)
    { (*trajectoryContainer)[i]->DrawTrajectory(); }
  }

  if(HC != nullptr)
  {
    std::size_t n_HC = HC->GetCapacity();
    for(std::size_t j=0; j<n_HC; ++j)
    {
      G4VHitsCollection* VHC = HC->GetHC((G4int)j);
      if(VHC != nullptr) VHC->DrawAllHits();
    }
  }

  if(DC != nullptr)
  {
    std::size_t n_DC = DC->GetCapacity();
    for(std::size_t j=0; j<n_DC; ++j)
    {
      G4VDigiCollection* VDC = DC->GetDC((G4int)j);
      if(VDC != nullptr) VDC->DrawAllDigi();
    }
  }
}

G4int G4Event::StoreSubEvent(G4int ty,G4SubEvent* se)
{
  std::set<G4SubEvent*>* sev = nullptr;
  auto ses = fSubEvtStackMap.find(ty);
  if(ses==fSubEvtStackMap.end())
  { 
    sev = new std::set<G4SubEvent*>;
    fSubEvtStackMap[ty] = sev;
  }
  else
  { sev = ses->second; }
  sev->insert(se);
  return (G4int)sev->size();
}

G4SubEvent* G4Event::PopSubEvent(G4int ty)
{
  G4SubEvent* se = nullptr;
  auto ses = fSubEvtStackMap.find(ty);
  if(ses!=fSubEvtStackMap.end())
  {
    auto sev = ses->second;
    if(!(sev->empty()))
    {
      se = sev->extract(sev->begin()).value();
      SpawnSubEvent(se);
    }
  }
  return se;
}

G4int G4Event::SpawnSubEvent(G4SubEvent* se)
{
  auto ss = fSubEvtVector.find(se);
  if(ss!=fSubEvtVector.end())
  {
    G4ExceptionDescription ed;
    ed << "Sub-event " << se << " of type " << se->GetSubEventType()
       << " with " << se->GetNTrack() << " tracks has already spawned.";
    G4Exception("G4Event::SpawnSubEvent","SubEvent9001",
                FatalException,ed);
  }
  fSubEvtVector.insert(se);
  return (G4int)fSubEvtVector.size();
}

void G4Event::MergeSubEventResults(const G4Event* se)
{
#ifdef G4_STORE_TRAJECTORY
  if(se->trajectoryContainer!=nullptr && se->trajectoryContainer->size()>0)
  {
    if(trajectoryContainer==nullptr) trajectoryContainer = new G4TrajectoryContainer;
    for(auto& trj : *(se->trajectoryContainer->GetVector()))
    { trajectoryContainer->push_back(trj); }
  }
#endif
  // Note:
  //     - scores are merged directly to the scoring manager
  //     - hits collections should be merged by the user event action
}

G4int G4Event::TerminateSubEvent(G4SubEvent* se)
{
  auto ss = fSubEvtVector.find(se);
  if(ss==fSubEvtVector.end())
  {
    G4ExceptionDescription ed;
    ed << "Sub-event " << se << " of type " << se->GetSubEventType()
       << " with " << se->GetNTrack() << " tracks has never been spawned.";
    G4Exception("G4Event::TerminateSubEvent","SubEvent9002",
                FatalException,ed);
  }

  fSubEvtVector.erase(ss);

  ss = fSubEvtVector.find(se);
  if(ss!=fSubEvtVector.end())
  {
    G4ExceptionDescription ed;
    ed << "Sub-event " << se << " of type " << se->GetSubEventType()
       << " with " << se->GetNTrack() << " appears more than once. PANIC!";
    G4Exception("G4Event::TerminateSubEvent","SubEvent9003",
                FatalException,ed);
  }

  se->clearAndDestroy();
  delete se;
  return (G4int)fSubEvtVector.size();
}

