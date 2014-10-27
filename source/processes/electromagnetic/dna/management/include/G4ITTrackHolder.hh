/*
 * G4ITTrackHolder.hh
 *
 *  Created on: 17 sept. 2014
 *      Author: kara
 */

#ifndef G4ITTRACKHOLDER_HH
#define G4ITTRACKHOLDER_HH

#include "G4IT.hh"
#include "G4TrackList.hh"
#include "G4VITTrackHolder.hh"
#include <iostream>

using namespace std;

typedef int Key;

class G4ITScheduler;

class PriorityList
{
  friend class G4ITTrackHolder;
public:
  enum Type
  {
    MainList = 0,
    SecondariesList = 1,
    WaitingList = 2,
    Undefined = -1
  };

  PriorityList() :
      fpMainList(0), fpWaitingList(0)
  {
  }

  PriorityList(const PriorityList& right) :
      fpMainList(right.fpMainList), fpWaitingList(right.fpWaitingList)
  {
  }

  ~PriorityList()
  {
    if (fpMainList) delete fpMainList;
    if (fpWaitingList) delete fpWaitingList;
  }

  G4TrackList* GetMainList()
  {
    return fpMainList;
  }

  G4TrackList* NewMainList(G4TrackManyList& allMainList)
  {
    G4TrackList* trackList = new G4TrackList();
    fpMainList = trackList;
    allMainList.Add(trackList);
    return fpMainList;
  }

  void PushToMainList(G4Track* track)
  {
    if (fpMainList == 0) fpMainList = new G4TrackList();
    fpMainList->push_back(track);
  }

  void MergeWithMainList(G4TrackList* trackList)
  {
    if (fpMainList == 0) fpMainList = new G4TrackList();
    trackList->transferTo(trackList);
  }

  G4TrackList* Get(Type type)
  {
    switch (type)
    {
      case MainList:
        return fpMainList;
        break;
      case SecondariesList:
        return &fSecondaries;
        break;
      case WaitingList:
        return fpWaitingList;
        break;
      case Undefined:
        return 0;
    }
    return 0;
  }

private:
  G4TrackList* fpMainList;
  G4TrackList fSecondaries;
  // to merge with fpMainList
  G4TrackList* fpWaitingList;
  // Waiting queue of currentList
};

class G4ITTrackHolder : public G4VITTrackHolder
{
  /* UR:
   * Push on time
   * Push delayed
   * Exception when going back
   * Get all tracks
   */

  static G4ThreadLocal G4ITTrackHolder* fgInstance;

  friend class G4ITScheduler;

public:

  static G4ITTrackHolder* Instance();

  G4ITTrackHolder();
  virtual
  ~G4ITTrackHolder();

  //----- Time of the next set of tracks -----
  inline double GetNextTime();

  virtual void Push(G4Track*);

  void Clear();

  inline void PuskToKill(G4Track* track)
  {
    fToBeKilledList.push_back(track);
  }
  void PushFront(G4Track*, PriorityList::Type type);
  bool MergeNextTimeToMainList(double& time);
  void MergeSecondariesWithMainList();

  // ----- To call at the end of the step
  void KillTracks();

  typedef std::map<Key, PriorityList*> MapOfLists;
  typedef std::map<double, std::map<Key, G4TrackList*> > MapOfDelayedLists;

  std::map<Key, PriorityList*>& GetLists()
  {
    return fLists;
  }

  MapOfDelayedLists& GetDelayedLists()
  {
    return fDelayedList;
  }

  bool MainListsNOTEmpty()
  {
    return CheckMapIsNOTEmpty(fLists, PriorityList::MainList);
  }

  bool SecondaryListsNOTEmpty()
  {
    return CheckMapIsNOTEmpty(fLists, PriorityList::SecondariesList);
  }

  void MoveMainToWaitingList()
  {
    MapOfLists::iterator it = fLists.begin();
    MapOfLists::iterator end = fLists.end();
    for (; it != end; it++)
    {
      if (PriorityList* lists = it->second)
      {
        lists->fpWaitingList = lists->fpMainList;
      }
    }
    fAllMainList.clear();
  }

  bool DelayListsNOTEmpty()
  {
    MapOfDelayedLists::iterator __it = fDelayedList.begin();
    MapOfDelayedLists::iterator __end = fDelayedList.end();
    for (; __it != __end; __it++)
    {
      std::map<Key, G4TrackList*>& mapOfLists = __it->second;
      if (mapOfLists.empty() == false)
      {
        std::map<Key, G4TrackList*>::iterator it = mapOfLists.begin();
        std::map<Key, G4TrackList*>::iterator end = mapOfLists.end();
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

  bool CheckMapIsNOTEmpty(MapOfLists& mapOfLists, PriorityList::Type type)
  {
    MapOfLists::iterator it = mapOfLists.begin();
    MapOfLists::iterator end = mapOfLists.end();
    for (; it != end; it++)
    {
      if (PriorityList* lists = it->second)
      {
        if (G4TrackList* mainList = lists->Get(type))
        {
          if (!(mainList->empty())) return true;
        }
      }
    }
    return false;
  }

  G4TrackManyList* GetMainList()
  {
    return &fAllMainList;
  }

  G4TrackManyList* GetSecondariesList()
  {
    return &fAllSecondariesList;
  }

  virtual size_t GetNTracks()
  {
    size_t nTracks(0);
    nTracks += fAllMainList.size();
    nTracks += fAllSecondariesList.size();

    MapOfDelayedLists::iterator delayedmap_it = fDelayedList.begin();
    MapOfDelayedLists::iterator delayedmap_end = fDelayedList.end();

    for (; delayedmap_it != delayedmap_end; delayedmap_it++)
    {
      std::map<Key, G4TrackList*>::iterator it = delayedmap_it->second.begin();
      std::map<Key, G4TrackList*>::iterator end = delayedmap_it->second.end();

      for (; it != end; it++)
      {
        if (it->second) nTracks += it->second->size();
      }
    }

    return nTracks;
  }

protected:
  void AddTrackID(G4Track* track);
  void _PushTrack(G4Track* track);
  void PushTo(G4Track*, PriorityList::Type);
  void PushDelayed(G4Track* track, const G4double& globalTime);

protected:
  std::map<Key, PriorityList*> fLists;
  MapOfDelayedLists fDelayedList;
  G4TrackList fToBeKilledList;
  bool fMainListHaveBeenSet;
  int fVerbose;
  int fNbTracks;

  double fPostActivityGlobalTime;
  //  double fPreActivityGlobalTime ;

  G4TrackManyList fAllMainList;
  G4TrackManyList fAllSecondariesList;
};

inline double G4ITTrackHolder::GetNextTime()
{
  if (fDelayedList.empty()) return DBL_MAX;
  return fDelayedList.begin()->first;
}

#endif /* G4MIMOLECULARTRACKS_HH_ */
