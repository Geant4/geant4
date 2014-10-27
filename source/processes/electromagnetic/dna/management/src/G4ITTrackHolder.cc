/*
 * G4MIMolecularTracks.cc
 *
 *  Created on: 17 sept. 2014
 *      Author: kara
 */

#include <G4ITScheduler.hh>
#include <G4VScheduler.hh>
#include "G4ITTrackHolder.hh"
#include "G4IT.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"

G4ThreadLocal G4ITTrackHolder* G4ITTrackHolder::fgInstance(0);

G4ITTrackHolder* G4ITTrackHolder::Instance()
{
  if (fgInstance == 0) new G4ITTrackHolder();
  return fgInstance;
}

G4ITTrackHolder::G4ITTrackHolder() :
    G4VITTrackHolder()
{
  fNbTracks = -1;
  fMainListHaveBeenSet = false;
  fVerbose = 0;

  fPostActivityGlobalTime = -1;
//  fPreActivityGlobalTime = -1;7
  fgInstance = this;
}

G4ITTrackHolder::~G4ITTrackHolder()
{
  std::map<Key, PriorityList*>::iterator end = fLists.end();

  for (std::map<Key, PriorityList*>::iterator it = fLists.begin(); it != end;
      it++)
  {
    delete it->second;
    it->second = 0;
  }

  if (!fDelayedList.empty())
  {
    MapOfDelayedLists::iterator fDelayedList_i = fDelayedList.begin();
    MapOfDelayedLists::iterator fDelayedList_end = fDelayedList.end();

    for (; fDelayedList_i != fDelayedList_end; fDelayedList_i++)
    {
      std::map<Key, G4TrackList*>::iterator it = fDelayedList_i->second.begin();
      std::map<Key, G4TrackList*>::iterator __end =
          fDelayedList_i->second.end();

      for (; it != __end; it++)
      {
        if (it->second) delete (it->second);
        it->second = 0;
      }
    }
    fDelayedList.clear();
  }

  fAllMainList.clear();
  fAllSecondariesList.clear();
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
//  G4cout << "G4MIMolecularTracks::MergeNextTimeToMainList" << G4endl;
  if (fDelayedList.empty()) return false;

//  G4cout << "fDelayedList.size = " << fDelayedList.size() <<G4endl;

  std::map<Key, G4TrackList*>::iterator it =
      fDelayedList.begin()->second.begin();
  std::map<Key, G4TrackList*>::iterator end =
      fDelayedList.begin()->second.end();
  if (it == end) return false;

  bool output = false;
  for (; it != end; it++)
  {
    PriorityList* right_listUnion(0);

    std::map<Key, PriorityList*>::iterator it_listUnion = fLists.find(
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
    if (right_listUnion->fpMainList == 0)
    {
      right_listUnion->fpMainList = new G4TrackList();
      fAllMainList.Add(right_listUnion->fpMainList);
    }
//    cout << it->second << endl;
//    cout << right_listUnion->fpMainList << endl;

//    G4cout << "Transfering : "
//        << GetIT(*(it->second->begin()))->GetName()
//        << G4endl;
//    G4cout << "At time : "
//            << fDelayedList.begin()->first
//            << G4endl;

    it->second->transferTo(right_listUnion->fpMainList);

    if (output == false && right_listUnion->fpMainList->size()) output = true;
    delete it->second;
  }

  if (output) time = fDelayedList.begin()->first;
  fDelayedList.erase(fDelayedList.begin());
  return output;
}

void G4ITTrackHolder::MergeSecondariesWithMainList()
{
  std::map<Key, PriorityList*>::iterator it = fLists.begin();
  std::map<Key, PriorityList*>::iterator end = fLists.end();

  for (; it != end; it++)
  {
    if (it->second->fpMainList == 0)
    {
      it->second->NewMainList(fAllMainList);
    }

    it->second->fSecondaries.transferTo(it->second->fpMainList);
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
  if (G4VScheduler::Instance()->IsRunning())
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "G4ITStepManager::PushTrack : You are trying to push tracks while the "
        "ITStepManager is running";
    G4Exception("G4ITStepManager::PushTrack", "ITStepManager012",
                FatalErrorInArgument, exceptionDescription);
  }
  _PushTrack(track);

//  G4MIConstituent::NotifyEntityAdded(track);
}
//_________________________________________________________________________
void G4ITTrackHolder::PushTo(G4Track* track, PriorityList::Type type)
{
  int moleculeID = GetIT(track)->GetITSubType();
  std::map<Key, PriorityList*>::iterator it = fLists.find(moleculeID);

  PriorityList* priorityList(0);

  if (it == fLists.end())
  {
    priorityList = new PriorityList();
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
      G4TrackList* trackList = 0;

      if (priorityList->fpMainList == 0)
      {
        trackList = priorityList->NewMainList(fAllMainList);
      }
      else
      {
        trackList = priorityList->fpMainList;
      }

      trackList->push_back(track);
      break;
    }
    case PriorityList::SecondariesList:
    {
      if (priorityList->fSecondaries.empty())
      {
        fAllSecondariesList.Add(&priorityList->fSecondaries);
      }

      priorityList->fSecondaries.push_back(track);
      break;
    }
    case PriorityList::WaitingList:
    {
      G4TrackList* trackList = 0;

      if (priorityList->fpWaitingList == 0)
      {
        trackList = new G4TrackList();
        priorityList->fpWaitingList = trackList;
      }
      else
      {
        trackList = priorityList->fpWaitingList;
      }

      trackList->push_back(track);

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
void G4ITTrackHolder::PushFront(G4Track* track, PriorityList::Type type)
{
  int moleculeID = GetIT(track)->GetITSubType();
  std::map<Key, PriorityList*>::iterator it = fLists.find(moleculeID);

  PriorityList* priorityList(0);

  if (it == fLists.end())
  {
    priorityList = new PriorityList();
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
      G4TrackList* trackList = 0;

      if (priorityList->fpMainList == 0)
      {
        trackList = priorityList->NewMainList(fAllMainList);
      }
      else
      {
        trackList = priorityList->fpMainList;
      }

      trackList->push_front(track);
      break;
    }
    case PriorityList::SecondariesList:
    {
      if (priorityList->fSecondaries.empty())
      {
        fAllSecondariesList.Add(&priorityList->fSecondaries);
      }
      priorityList->fSecondaries.push_front(track);
      break;
    }
    case PriorityList::WaitingList:
    {
      G4TrackList* trackList = 0;

      if (priorityList->fpWaitingList == 0)
      {
        trackList = new G4TrackList();
        priorityList->fpWaitingList = trackList;
      }
      else
      {
        trackList = priorityList->fpWaitingList;
      }

      trackList->push_front(track);

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
  if (track == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "You are trying to push a non-existing track (track pointer is null)"
        << G4endl;

    G4Exception("G4ITStepManager::_PushTrack", "ITStepManager014",
                FatalErrorInArgument, exceptionDescription);
  }

  G4double globalTime = track->GetGlobalTime();

  if (track->GetTrackID() == 0)
  {
    // Set track ID
    AddTrackID(track);
  }

  double currentGlobalTime = G4ITScheduler::Instance()->GetGlobalTime();

#ifdef G4VERBOSE
  if (fVerbose > 5)
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

  if (G4ITScheduler::Instance()->IsRunning() == false)
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

      G4Exception("G4ITStepManager::_PushTrack", "ITStepManager014",
                  FatalErrorInArgument, exceptionDescription);
    }

    // Push the track to the rigth track list :
    // If the track time is the same as the main track list,
    // it will be push to the main track list
    // otherwise, it will be pushed to the delayed track list.
    if (fMainListHaveBeenSet == false)
    {
      PushDelayed(track, globalTime);
    }
    else
    {
      if (globalTime == currentGlobalTime)
      {
        PushTo(track, PriorityList::MainList);
      }
      else
      {
        PushDelayed(track, globalTime);
      }
    }
  }
  else // Is running
  {
    double timeDifference = globalTime - currentGlobalTime;
    double timeTolerance = G4ITScheduler::Instance()->GetTimeTolerance();

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

      G4Exception("G4ITStepManager::_PushTrack", "ITStepManager015",
                  FatalErrorInArgument, exceptionDescription);
    }

    // Push the track to the rigth track list :
    // If the track time is the same as the main track list,
    // it will be push to the secondary list
    // otherwise, it will be pushed to the delayed track list.
    if (fabs(timeDifference) < timeTolerance)
    {
//      G4cout << "Is pushing " << GetIT(track)->GetName() << G4endl;

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

      G4Exception("G4ITStepManager::_PushTrack", "ITStepManager016",
                  FatalErrorInArgument, exceptionDescription);
      // PushDelayed(track, globalTime);
    }
  }
}

//_________________________________________________________________________

void G4ITTrackHolder::PushDelayed(G4Track* track, const G4double& globalTime)
{
#ifdef G4VERBOSE
  if (fVerbose > 5)
  {
    G4cout << "\t" << ">> Pushing a delayed track" << G4endl;
  }
#endif

  int moleculeID = GetIT(track)->GetITSubType();
  //  std::map<int, PriorityList>::iterator it = fLists.find(moleculeID);

  std::map<double, std::map<Key, G4TrackList*> >::iterator it_delayed =
  fDelayedList.find(globalTime);

  if (it_delayed == fDelayedList.end())
  {
    (fDelayedList[globalTime][moleculeID] = new G4TrackList())->push_back(
        track);
  }
  else
  {
    std::map<Key, G4TrackList*>::iterator it_trackList =
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
  if (fToBeKilledList.size() == 0) return;
#ifdef G4VERBOSE
  if (fVerbose > 5)
  {
    G4cout << "*** G4ITStepManager::KillTracks , step #"
//           << G4MIWorldEngine::Instance()->GetScheduler()->GetNbSteps()
           << G4VScheduler::Instance()->GetNbSteps() << " ***" << G4endl;
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
  std::map<Key, PriorityList*>::iterator it = fLists.begin();

  for (; it != fLists.end(); it++)
  {
    delete it->second;
    it->second = 0;
  }
  fLists.clear();

  MapOfDelayedLists::iterator it1 = fDelayedList.begin();

  for (; it1 != fDelayedList.end(); it1++)
  {
    std::map<Key, G4TrackList*>::iterator it2 = it1->second.begin();

    for (; it2 != it1->second.end(); it2++)
    {
      delete it2->second;
      it2->second = 0;
    }
  }

  fDelayedList.clear();

  fAllMainList.clear();
  fAllSecondariesList.clear();

  fNbTracks = -1;
}
