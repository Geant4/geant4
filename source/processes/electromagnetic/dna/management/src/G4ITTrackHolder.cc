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
/*
 * G4MIMolecularTracks.cc
 *
 *  Created on: 17 sept. 2014
 *      Author: kara
 */

#include <G4Scheduler.hh>
#include <G4VScheduler.hh>
#include "G4ITTrackHolder.hh"
#include "G4IT.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"

using namespace std;

PriorityList::PriorityList() :
     fpMainList(nullptr), fpWaitingList(nullptr)
{
}

PriorityList::PriorityList(G4TrackManyList& allMainList):
     fpMainList(nullptr), fpWaitingList(nullptr)
{
  NewMainList(allMainList);
}

PriorityList::PriorityList(const PriorityList& right) :
    
    fpMainList(right.fpMainList),
    fpWaitingList(right.fpWaitingList)
{
}

PriorityList::~PriorityList()
{
  if (fpMainList != nullptr)
  {
    delete fpMainList;
    fpMainList = nullptr;
  }
  if (fpWaitingList != nullptr)
  {
    delete fpWaitingList;
    fpWaitingList = nullptr;
  }
}

void PriorityList::NotifyDeletingList(G4TrackList* __list)
{
  if (__list == fpMainList)
  {
//    StopWatching(fpMainList);
    fpMainList = nullptr;
  }
  else if (__list == fpWaitingList)
  {
//    StopWatching(fpWaitingList);
    fpWaitingList = nullptr;
  }
}

void PriorityList::NewMainList(G4TrackList* __list,
                               G4TrackManyList& allMainList)
{
  fpMainList = __list;
  allMainList.Add(__list);
  Watch(fpMainList);
}

G4TrackList* PriorityList::NewMainList(G4TrackManyList& allMainList)
{
  auto  trackList = new G4TrackList();
  NewMainList(trackList, allMainList);
  return fpMainList;
}

void PriorityList::PushToMainList(G4Track* __track,
                                  G4TrackManyList& allMainList)
{
  if (fpMainList == nullptr)
  {
    NewMainList(allMainList);
  }
  fpMainList->push_back(__track);
}

void PriorityList::TransferToMainList(G4TrackList*& __list,
                                      G4TrackManyList& allMainList)
{
  if (fpMainList != nullptr)
  {
    __list->transferTo(fpMainList);
    delete __list;
    __list = nullptr;
  }
  else
  {
    NewMainList(__list, allMainList);
  }
}

void PriorityList::PushToListOfSecondaries(G4Track* __track,
                                           G4TrackManyList& listOfAllSecondaries)
{
  //      if (priorityList->fSecondaries.empty())
  if (fSecondaries.GetListNode() != nullptr)
  {
    listOfAllSecondaries.Add(&fSecondaries);
  }
  fSecondaries.push_back(__track);
}

void PriorityList::PushToWaitingList(G4Track* __track)
{
  if (fpWaitingList == nullptr)
  {
    fpWaitingList = new G4TrackList();
  }
  fpWaitingList->push_back(__track);
}

void PriorityList::TransferSecondariesToMainList()
{
  fSecondaries.transferTo(fpMainList);
}

void PriorityList::PushToMainList(G4Track* track)
{
  if (fpMainList == nullptr) fpMainList = new G4TrackList();
  fpMainList->push_back(track);
}

void PriorityList::MergeWithMainList(G4TrackList* trackList)
{
  if (fpMainList == nullptr) fpMainList = new G4TrackList();
  trackList->transferTo(trackList);
}

int PriorityList::GetNTracks()
{
  int nTracks = 0;

  if (fpMainList != nullptr)
  {
    nTracks += fpMainList->size();
  }

  if (fpWaitingList != nullptr)
  {
    nTracks += fpWaitingList->size();
  }

  nTracks += fSecondaries.size();

  return nTracks;
}

//=============================================================================
// G4ITTrackHolder
//=============================================================================

G4ThreadLocal G4ITTrackHolder* G4ITTrackHolder::fgInstance(nullptr);
G4ITTrackHolder* G4ITTrackHolder::fgMasterInstance(nullptr);

G4Mutex creationOfTheMasterInstance;
G4Mutex pushToTheMasterInstance;

G4ITTrackHolder* G4ITTrackHolder::Instance()
{
  if (fgInstance == nullptr)
  {
    fgInstance = new G4ITTrackHolder();
    if(G4Threading::IsMasterThread() ||
       !G4Threading::IsMultithreadedApplication() 
    ) 
    {
      fgMasterInstance = fgInstance;
    }
    
  }
  return fgInstance;
}

G4ITTrackHolder* G4ITTrackHolder::MasterInstance()
{
  G4AutoLock lock(&creationOfTheMasterInstance);
  if (fgMasterInstance == nullptr)
  {
    fgMasterInstance = new G4ITTrackHolder();
  }
  lock.unlock();
  return fgMasterInstance;
}

G4ITTrackHolder::G4ITTrackHolder()
{
  fNbTracks = -1;
  fMainListHaveBeenSet = false;
  fVerbose = 0;

  fPostActivityGlobalTime = -1;
//  fPreActivityGlobalTime = -1;
}

G4ITTrackHolder::~G4ITTrackHolder()
{
  auto end = fLists.end();

  for (auto it = fLists.begin(); it != end;
      it++)
  {
    delete it->second;
    it->second = 0;
  }

  if (!fDelayedList.empty())
  {
    auto fDelayedList_i = fDelayedList.begin();
    auto fDelayedList_end = fDelayedList.end();

    for (; fDelayedList_i != fDelayedList_end; fDelayedList_i++)
    {
      auto it = fDelayedList_i->second.begin();
      auto __end =
          fDelayedList_i->second.end();

      for (; it != __end; it++)
      {
        delete (it->second);
        it->second = 0;
      }
    }
    fDelayedList.clear();
  }

  fAllMainList.RemoveLists();
  fAllSecondariesList.RemoveLists();
  fNbTracks = -1;
}

/*
 void G4MIMolecularTracks::Decide()
 {
 cout << "G4MIMolecularTracks::Decide" << endl;

 if (fDelayedList.empty())
 {
 cout << "fDelayedList.empty()" << endl;
 return;
 }
 fPostActivityGlobalTime = GetNextTime();
 //  PushActivity(workspace->GetScheduler(), this);
 }
 */

/*
 * param time = time of the merged list
 * returned = was there actually merged data ?
 */
bool G4ITTrackHolder::MergeNextTimeToMainList(double& time)
{
//  G4cout << "G4ITTrackHolder::MergeNextTimeToMainList" << G4endl;
  if (fDelayedList.empty())
  {
    return false;
  }

//  G4cout << "fDelayedList.size = " << fDelayedList.size() <<G4endl;

  auto it =
      fDelayedList.begin()->second.begin();
  auto end =
      fDelayedList.begin()->second.end();
  if (it == end) return false;

  bool output = false;
  for (; it != end; it++)
  {
    PriorityList* right_listUnion(nullptr);

    auto it_listUnion = fLists.find(
        it->first);
    if (it_listUnion == fLists.end())
    {
      right_listUnion = (fLists[it->first] = new PriorityList());
    }
    else
    {
      if (it_listUnion->second == 0)
      {
        it_listUnion->second = new PriorityList();
      }
      right_listUnion = it_listUnion->second;
    }

    if (it->second == 0) continue;

    /*
     if (right_listUnion->GetMainList() == 0)
     {
     //      right_listUnion->fpMainList = new G4TrackList();
     //      if(it->second)
     //      {
     right_listUnion->NewMainList(it->second, fAllMainList);
     //      }
     }
     else
     {
     right_listUnion->TransferToMainList(it->second);
     delete it->second;
     }*/

    right_listUnion->TransferToMainList(it->second, fAllMainList);

    if (!output)
    {
      if (static_cast<int>(!right_listUnion->GetMainList()->empty()) != 0)
      {
        output = true;
      }
    }
    it->second = 0;
  }

  if (output) time = fDelayedList.begin()->first;
  fDelayedList.erase(fDelayedList.begin());
  return output;
}

void G4ITTrackHolder::MergeSecondariesWithMainList()
{
  auto it = fLists.begin();
  auto end = fLists.end();

  for (; it != end; it++)
  {
    if (it->second->GetMainList() == nullptr)
    {
      it->second->NewMainList(fAllMainList);
    }

    it->second->TransferSecondariesToMainList();
  }
}

//_________________________________________________________________________

void G4ITTrackHolder::AddTrackID(G4Track* track)
{
  //if(fNbTracks == 0) fNbTracks = -1;
  track->SetTrackID(fNbTracks);
  fNbTracks--;
}

//_________________________________________________________________________

void G4ITTrackHolder::Push(G4Track* track)
{
//  if (G4VScheduler::Instance()->IsRunning())
//  {
//    G4ExceptionDescription exceptionDescription;
//    exceptionDescription
//        << "G4ITTrackHolder::PushTrack : You are trying to push tracks while the "
//        "ITStepManager is running";
//    G4Exception("G4ITTrackHolder::PushTrack", "ITStepManager012",
//                FatalErrorInArgument, exceptionDescription);
//  }
  _PushTrack(track);

//  G4MIConstituent::NotifyEntityAdded(track);
}
//_________________________________________________________________________
void G4ITTrackHolder::PushTo(G4Track* track, PriorityList::Type type)
{
  int moleculeID = GetIT(track)->GetITSubType();
  auto it = fLists.find(moleculeID);

  PriorityList* priorityList(nullptr);

  if (it == fLists.end())
  {
    priorityList = new PriorityList(fAllMainList);
    fLists[moleculeID] = priorityList;
  }
  else
  {
    priorityList = it->second;
  }

  switch (type)
  {
    case PriorityList::MainList:
    {
      priorityList->PushToMainList(track, fAllMainList);
      break;
    }
    case PriorityList::SecondariesList:
    {
      priorityList->PushToListOfSecondaries(track, fAllSecondariesList);
      break;
    }
    case PriorityList::WaitingList:
    {
      priorityList->PushToWaitingList(track);
      return;
      break;
    }

    default:
    {
      return;
      break;
    }
  }
}
//_________________________________________________________________________

void G4ITTrackHolder::_PushTrack(G4Track* track)
{
  if (track == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "You are trying to push a non-existing track (track pointer is null)"
        << G4endl;

    G4Exception("G4ITTrackHolder::_PushTrack", "ITStepManager014",
                FatalErrorInArgument, exceptionDescription);
  }

  G4double globalTime = track->GetGlobalTime();

  if (track->GetTrackID() == 0)
  {
    // Set track ID
    AddTrackID(track);
  }

  double currentGlobalTime = G4Scheduler::Instance()->GetGlobalTime();

#ifdef G4VERBOSE
  if (fVerbose != 0)
  {
    G4cout << G4endl;
    G4cout << "\t"<< ">> Pushing a track -->  ";
    G4cout << GetIT(track)->GetName() << " (" << track->GetTrackID() <<")"
    << " -- ";
    G4cout << "Global current time: " << G4BestUnit(currentGlobalTime,"Time")
    << "\t";
    G4cout << "Track's time: " << G4BestUnit(track->GetGlobalTime(),"Time")
    << G4endl;
  }
#endif

  if (!G4Scheduler::Instance()->IsRunning())
  {
    if (globalTime < currentGlobalTime)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "You are trying to push a track with a global time"
          << " inferior to the current simulation time." << G4endl<< "The time is going back : " << G4endl
      << "The time in the step manager : "
      << G4BestUnit(currentGlobalTime,"Time")
      << G4endl
      << "The time of the track : "
      << G4BestUnit(globalTime,"Time")
      << G4endl
      << "(ITStepManager is not yet running)"
      << G4endl;

      G4Exception("G4ITTrackHolder::_PushTrack", "ITStepManager014",
                  FatalErrorInArgument, exceptionDescription);
    }

    // Push the track to the rigth track list :
    // If the track time is the same as the main track list,
    // it will be push to the main track list
    // otherwise, it will be pushed to the delayed track list.
    if (!fMainListHaveBeenSet)
    {
      PushDelayed(track);
    }
    else
    {
      if (globalTime == currentGlobalTime)
      {
        #ifdef G4VERBOSE
          if (fVerbose != 0)
          {
            G4cout << G4endl;
            G4cout << "\t"<< ">> Pushing to *main* list -->  ";
            G4cout << GetIT(track)->GetName() << " (" << track->GetTrackID() <<")"
            << " -- ";
            G4cout << "Global current time: " << G4BestUnit(currentGlobalTime,"Time")
            << "\t";
            G4cout << "Track's time: " << G4BestUnit(track->GetGlobalTime(),"Time")
            << G4endl;
          }
        #endif
        PushTo(track, PriorityList::MainList);
      }
      else
      {
        // if(currentGlobalTime > 1*CLHEP::picosecond) abort();
        #ifdef G4VERBOSE
          if (fVerbose != 0)
          {
            G4cout << G4endl;
            G4cout << "\t"<< ">> Pushing to *delayed* list -->  ";
            G4cout << GetIT(track)->GetName() << " (" << track->GetTrackID() <<")"
            << " -- ";
            G4cout << "Global current time: " << G4BestUnit(currentGlobalTime,"Time")
            << "\t";
            G4cout << "Track's time: " << G4BestUnit(track->GetGlobalTime(),"Time")
            << G4endl;
          }
        #endif
        PushDelayed(track);
      }
    }
  }
  else // Is running
  {
    double timeDifference = globalTime - currentGlobalTime;
    double timeTolerance = G4Scheduler::Instance()->GetTimeTolerance();

    if (timeDifference < -1 * timeTolerance)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "You are trying to push a track with a global time"
          << " inferior to the current simulation time." << G4endl<< "The time is going back : "
      << G4endl
      << "The time in the step manager : "
      << G4BestUnit(timeDifference,"Time")
      << G4endl
      << "The time of the track : " << G4BestUnit(globalTime,"Time")
      << G4endl
      << "(ITStepManager is running)"
      << G4endl;

      G4Exception("G4ITTrackHolder::_PushTrack", "ITStepManager015",
                  FatalErrorInArgument, exceptionDescription);
    }

    // Push the track to the rigth track list :
    // If the track time is the same as the main track list,
    // it will be push to the secondary list
    // otherwise, it will be pushed to the delayed track list.
    if (fabs(timeDifference) < timeTolerance)
    {
//      G4cout << "Is pushing " << GetIT(track)->GetName() << G4endl;

      #ifdef G4VERBOSE
        if (fVerbose != 0)
        {
          G4cout << G4endl;
          G4cout << "\t"<< ">> Pushing to *secondary* list -->  ";
          G4cout << GetIT(track)->GetName() << " (" << track->GetTrackID() <<")"
          << " -- ";
          G4cout << "Global current time: " << G4BestUnit(currentGlobalTime,"Time")
          << "\t";
          G4cout << "Track's time: " << G4BestUnit(track->GetGlobalTime(),"Time")
          << G4endl;
        }
      #endif
      PushTo(track, PriorityList::SecondariesList);
    }
    else // globalTime < fGlobalTime already taken into account above
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "While running you cannot push a track"
          << " with a bigger global time than the current global time" << G4endl<< "The time in the step manager : "
      << G4BestUnit(currentGlobalTime,"Time")
      << G4endl
      << "The time of the track : " << G4BestUnit(globalTime,"Time")
      << G4endl
      << "(ITStepManager is running)"
      << G4endl;

      G4Exception("G4ITTrackHolder::_PushTrack", "ITStepManager016",
                  FatalErrorInArgument, exceptionDescription);
      // PushDelayed(track, globalTime);
    }
  }
}

//_________________________________________________________________________

void G4ITTrackHolder::PushDelayed(G4Track* track)
{
#ifdef G4VERBOSE
  if (fVerbose != 0)
  {
    G4cout << "\t" << ">> Pushing a delayed track" << G4endl;
  }
#endif

  int moleculeID = GetIT(track)->GetITSubType();
  //  std::map<int, PriorityList>::iterator it = fLists.find(moleculeID);

  G4double globalTime = track->GetGlobalTime();

  auto it_delayed =
  fDelayedList.find(globalTime);

  if (it_delayed == fDelayedList.end())
  {
    (fDelayedList[globalTime][moleculeID] = new G4TrackList())->push_back(
        track);
  }
  else
  {
    auto it_trackList =
    it_delayed->second.find(moleculeID);

    if (it_trackList == it_delayed->second.end())
    {
      (it_delayed->second[moleculeID] = new G4TrackList())->push_back(track);
    }
    else
    {
      if (it_trackList->second != 0)
      {
        it_trackList->second->push_back(track);
      }
    }
  }

  //  fDelayedList[globalTime][moleculeID]

  /*
   std::map<double,std::map<int, G4TrackList* > >::iterator it_delayed =
   fDelayedList.begin();

   std::map<double,std::map<int, G4TrackList* > >::iterator end_delayed =
   fDelayedList.end();

   for(it_delayed != end_delayed ; it_delayed++)
   {
   std::map<int, G4TrackList*> & trackListMap = it->second;


   }
   */
  /*
   std::map<double,G4TrackList* > :: iterator
   fDelayedList_i = fDelayedList.find(globalTime) ;

   if(fDelayedList_i == fDelayedList.end())
   {

   G4TrackList* newList = new G4TrackList ;
   newList -> push_back(track);
   fDelayedList[globalTime] = newList ;
   }
   else
   {
   fDelayedList_i->second-> push_back(track);
   }
   */
}
//______________________________________________________________________________

void G4ITTrackHolder::KillTracks()
{
  if (fToBeKilledList.empty()) return;
#ifdef G4VERBOSE
  if (fVerbose > 1)
  {
    G4cout << "*** G4ITTrackHolder::KillTracks , step #"
           << G4VScheduler::Instance()->GetNbSteps()
           << " ***" << G4endl;
    G4cout << "Nb of tracks to kill "<< fToBeKilledList.size() << G4endl;
    G4cout << setw(25) << left << "#Name"
           << setw(25) << "track ID"<< G4endl;

    G4TrackList::iterator it = fToBeKilledList.begin();
    for(; it != fToBeKilledList.end();)
    {
      G4Track* toBeErased = *it;

      G4cout << setw(25) << GetIT(toBeErased)->GetName()
      << setw(25) << toBeErased->GetTrackID()
      << G4endl;

      it = fToBeKilledList.erase(toBeErased);
    }
  }
  else
#endif
  fToBeKilledList.erase(fToBeKilledList.begin(), fToBeKilledList.end());
}

void G4ITTrackHolder::Clear()
{
  fAllMainList.ClearLists();
  fAllSecondariesList.ClearLists();
//  fAllMainList.RemoveLists();
//  fAllSecondariesList.RemoveLists();

  auto it = fLists.begin();

  for (; it != fLists.end(); it++)
  {
    delete it->second;
    it->second = 0;
  }
  fLists.clear();

  auto it1 = fDelayedList.begin();

  for (; it1 != fDelayedList.end(); it1++)
  {
    auto it2 = it1->second.begin();

    for (; it2 != it1->second.end(); it2++)
    {
      delete it2->second;
      it2->second = 0;
    }
  }

  fDelayedList.clear();

//  fAllMainList.ClearLists();
//  fAllSecondariesList.ClearLists();
  fAllMainList.RemoveLists();
  fAllSecondariesList.RemoveLists();
  KillTracks();

  fNbTracks = -1;
}

PriorityList* G4ITTrackHolder::GetPriorityList(Key i)
{
  auto it = fLists.find(i);
  if (it == fLists.end()) return nullptr;
  return it->second;
}

G4TrackList* G4ITTrackHolder::GetMainList(G4int i)
{
  PriorityList* priorityList = GetPriorityList(i);
  if (priorityList != nullptr)
  {
    return priorityList->GetMainList();
  }
  return nullptr;
}

bool G4ITTrackHolder::AddWatcher(Key id,
                                 G4TrackList::Watcher* watcher,
                                 PriorityList::Type type)
{
  auto it = fLists.find(id);
  if (it == fLists.end()) return false;

  G4TrackList* trackList = it->second->Get(type);
  if (trackList == nullptr) return false;
  trackList->AddWatcher(watcher);
  return true;
}

void G4ITTrackHolder::AddWatcherForMainList(G4TrackList::Watcher* watcher)
{
  fAllMainList.AddGlobalWatcher(watcher);
}

void G4ITTrackHolder::AddWatcherForKillList(G4TrackList::Watcher* watcher)
{
  watcher->Watch(&fToBeKilledList);
}

void G4ITTrackHolder::PushToMaster(G4Track* track)
{
  G4ITTrackHolder* trackHolder = MasterInstance();

  G4AutoLock lock(&pushToTheMasterInstance);
  trackHolder->PushDelayed(track);
  lock.unlock();
}

size_t G4ITTrackHolder::GetNTracks()
{
  size_t nTracks(0);
  nTracks += fAllMainList.size();
  nTracks += fAllSecondariesList.size();

  //    G4cout << "nTracks = " << nTracks << G4endl;

  auto delayedmap_it = fDelayedList.begin();
  auto delayedmap_end = fDelayedList.end();

  for (; delayedmap_it != delayedmap_end; delayedmap_it++)
  {
    auto it = delayedmap_it->second.begin();
    auto end = delayedmap_it->second.end();

    for (; it != end; it++)
    {
      if (it->second != nullptr) nTracks += it->second->size();
    }
  }

  //    G4cout << "nTracks = " << nTracks << G4endl;

  return nTracks;
}

void G4ITTrackHolder::MoveMainToWaitingList()
{
  auto it = fLists.begin();
  auto end = fLists.end();
  for (; it != end; it++)
  {
    if (PriorityList* lists = it->second)
    {
      lists->SetWaitingList(lists->GetMainList());
      //TODO
    }
  }
  fAllMainList.RemoveLists();
}

bool G4ITTrackHolder::DelayListsNOTEmpty()
{
  auto __it = fDelayedList.begin();
  auto __end = fDelayedList.end();
  for (; __it != __end; __it++)
  {
    std::map<Key, G4TrackList*>& mapOfLists = __it->second;
    if (!mapOfLists.empty())
    {
      auto it = mapOfLists.begin();
      auto end = mapOfLists.end();
      for (; it != end; it++)
      {
        if (G4TrackList* mainList = it->second)
        {
          if (!(mainList->empty())) return true;
        }
      }
    }
  }
  return false;
}

bool G4ITTrackHolder::CheckMapIsNOTEmpty(MapOfPriorityLists& mapOfLists,
                                         PriorityList::Type type)
{
  auto it = mapOfLists.begin();
  auto end = mapOfLists.end();
  for (; it != end; it++)
  {
    if (PriorityList* lists = it->second)
    {
      if (G4TrackList* trackList = lists->Get(type))
      {
        if (!(trackList->empty())) return true;
      }
    }
  }
  return false;
}
