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
 * G4ITReactionInfo.hh
 *
 *  Created on: 1 févr. 2015
 *      Author: matkara
 *      updated: HoangTran
 */

#ifndef G4ITREACTIONINFO_HH_
#define G4ITREACTIONINFO_HH_

#include "tls.hh"
#include <list>
#include <map>
#include <G4memory.hh>
#include "G4Track.hh"
#include <set>

using G4TrackVectorHandle = G4shared_ptr< std::vector<G4Track*> >;

#ifndef compTrackPerID__
#define compTrackPerID__
struct compTrackPerID
{
  G4bool operator()(G4Track* rhs, G4Track* lhs) const
  {
    return rhs->GetTrackID() < lhs->GetTrackID();
  }
};
#endif

class G4Track;
class G4ITReactionSet;
class G4ITReactionPerTrack;
class G4ITReaction;
using G4ITReactionPtr = G4shared_ptr<G4ITReaction>;
using G4ITReactionPerTrackPtr = G4shared_ptr<G4ITReactionPerTrack>;

using G4ITReactionList = std::list<G4ITReactionPtr>;
using G4ITReactionPerTrackMap = std::map<G4Track*,
                 G4ITReactionPerTrackPtr,
                 compTrackPerID>;
using G4ReactionPerTrackIt = std::list<std::pair<G4ITReactionPerTrackPtr,
                            G4ITReactionList::iterator> >;

struct compReactionPerTime
{
     G4bool operator()(G4ITReactionPtr rhs,
                     G4ITReactionPtr lhs) const;
};

using G4ITReactionPerTime = std::multiset<G4ITReactionPtr, compReactionPerTime>;
using G4ITReactionPerTimeIt = std::multiset<G4ITReactionPtr, compReactionPerTime>::iterator;

class G4ITReaction : public G4enable_shared_from_this<G4ITReaction>
{
  G4ITReaction(G4double time, G4Track*, G4Track*);
public:
  static G4ITReactionPtr New(G4double time, G4Track* trackA, G4Track* trackB)
  {
    return G4ITReactionPtr(new G4ITReaction(time, trackA, trackB));
  }
  virtual ~G4ITReaction();

  G4Track* GetReactant(G4Track* trackA) const
  {
    if(fReactants.first != trackA) return fReactants.first;
    return fReactants.second;
  }

  std::pair<G4Track*, G4Track*> GetReactants() const{return fReactants;}
  std::size_t GetHash() const;
  G4double GetTime() const
  {
    return fTime;
  }

  void RemoveMe();

  void AddIterator(G4ITReactionPerTrackPtr reactionPerTrack,
                   G4ITReactionList::iterator it)
  {
    fReactionPerTrack.emplace_back(reactionPerTrack, it);
  }

  void AddIterator(G4ITReactionPerTimeIt it)
  {
    fReactionPerTimeIt = new G4ITReactionPerTimeIt(it);
  }

  G4double fTime;
  std::pair<G4Track*, G4Track*> fReactants;
  G4ReactionPerTrackIt fReactionPerTrack;
  G4ITReactionPerTimeIt* fReactionPerTimeIt;

  //static G4ThreadLocal std::set<G4ITReaction*>* gAll;
};

class G4ITReactionPerTrack  : public G4enable_shared_from_this<G4ITReactionPerTrack>
{
  G4ITReactionPerTrack() = default;
public:
  static G4ITReactionPerTrackPtr New()
  {
    return G4ITReactionPerTrackPtr(new G4ITReactionPerTrack());
  }

  virtual ~G4ITReactionPerTrack()
  {
    fReactions.clear();
  }

  void AddReaction(G4ITReactionPtr reaction)
  {
    auto it =
      fReactions.insert(fReactions.end(), reaction);
    reaction->AddIterator(this->shared_from_this(), it);
  }

  void AddIterator(G4ITReactionPerTrackMap::iterator it)
  {
    fReactionSetIt.push_back(it);
  }

  G4bool RemoveThisReaction(G4ITReactionList::iterator it);
  void RemoveMe()
  {
    G4ITReactionPerTrackPtr backMeUp = this->shared_from_this();

    G4ITReactionList::iterator next;
    for(auto it = fReactions.begin() ;
        it !=  fReactions.end() ; it = next)
    {
      next = it;
      ++next;
       (*it)->RemoveMe();
    }
    fReactions.clear();
    fReactionSetIt.clear();
  }

  G4ITReactionList& GetReactionList()
  {
    return fReactions;
  }

  std::list<G4ITReactionPerTrackMap::iterator>& GetListOfIterators()
  {
    return fReactionSetIt;
  }

protected:
  G4ITReactionList fReactions;
  std::list<G4ITReactionPerTrackMap::iterator> fReactionSetIt;
};

class G4ITReactionSet
{
  G4ITReactionSet() : fReactionPerTime(compReactionPerTime())
  {
    fpInstance = this;
    fSortByTime = false;
  }
public:
  virtual ~G4ITReactionSet()
  {
    fReactionPerTrack.clear();
    fReactionPerTime.clear();
  }
  
  static G4ITReactionSet* Instance()
  {
    if(fpInstance == nullptr) new G4ITReactionSet();

    return fpInstance;
  }

  //------------------------------------------------------------------------------------

  void AddReaction(G4double time, G4Track* trackA, G4Track* trackB)
  {
    if(CanAddThisReaction(trackA, trackB))
    {
      G4ITReactionPtr reaction(G4ITReaction::New(time, trackA, trackB));
      AddReaction(trackA, reaction);
      AddReaction(trackB, reaction);

      if(fSortByTime)
      {
        auto it = fReactionPerTime.insert(reaction);
        reaction->AddIterator(it);
      }
    }
  }

//Hoang: this function checks if this reaction is added
  G4bool CanAddThisReaction(G4Track* trackA, G4Track* trackB)
  {
      auto it_track = fReactionPerTrack.find(trackA);
      G4ITReactionPerTrackPtr reactionPerTrack;
      if(it_track == fReactionPerTrack.end())
      {
          return true;
      }
      
      reactionPerTrack = it_track->second;
      auto list = reactionPerTrack->GetReactionList();
      //for(auto it_list = list.begin(); it_list != list.end(); ++it_list)
      for(const auto& it_list:list)
      {
          if ((*it_list).GetReactant(trackA)->GetTrackID() == trackB->GetTrackID())
          {
              return false;
          }
      }
      return true;
  }

  void AddReactions(G4double time, G4Track* trackA, G4TrackVectorHandle reactants)
  {
    auto it = reactants->begin();
    for(;it != reactants->end() ; ++it)
    {
      AddReaction(time, trackA, *it);
    }
  }

  void RemoveReactionSet(G4Track* track)
  {
    auto it = fReactionPerTrack.find(track);
    if(it != fReactionPerTrack.end())
    {
      G4ITReactionPerTrackPtr backItUp = it->second->shared_from_this();
      backItUp->RemoveMe();
      //fReactionPerTrack.erase(it); // not needed : once empty ==> auto-erase
      it = fReactionPerTrack.find(track);
      if(it != fReactionPerTrack.end())
      {
        fReactionPerTrack.erase(it);
      }
    }
  }

  void SelectThisReaction(G4ITReactionPtr reaction)
  {
    reaction->RemoveMe();
    RemoveReactionSet(reaction->GetReactants().first);
    RemoveReactionSet(reaction->GetReactants().second);
  }

  G4ITReactionPerTrackMap& GetReactionMap()
  {
    return fReactionPerTrack;
  }

  void RemoveReactionPerTrack(G4ITReactionPerTrackPtr reactionPerTrack)
  {
    for(auto & it : reactionPerTrack->GetListOfIterators())
    {
      fReactionPerTrack.erase(it);
    }
    reactionPerTrack->GetListOfIterators().clear();
    reactionPerTrack->GetReactionList().clear();
  }

  void CleanAllReaction()
  {
    for(auto it = fReactionPerTrack.begin();
        it != fReactionPerTrack.end() ;
        it = fReactionPerTrack.begin())
    {
      it->second->RemoveMe();
    }
    fReactionPerTrack.clear();
    fReactionPerTime.clear();
  }

  G4bool Empty()
  {
    return fReactionPerTrack.empty();
  }

  G4ITReactionPerTime& GetReactionsPerTime()
  {
    return fReactionPerTime;
  }

  void SortByTime(){
	  fSortByTime = true;
  }

protected:
  void AddReaction(G4Track* track, G4ITReactionPtr reaction)
  {
    auto it = fReactionPerTrack.find(track);

    G4ITReactionPerTrackPtr reactionPerTrack;

    if(it == fReactionPerTrack.end())
    {
     reactionPerTrack = G4ITReactionPerTrack::New();
     std::pair< G4ITReactionPerTrackMap::iterator,bool> pos =
          fReactionPerTrack.insert(std::make_pair(track, reactionPerTrack));
     reactionPerTrack->AddIterator(pos.first);
    }
    else
    {
     reactionPerTrack = it->second;
    }

    reactionPerTrack->AddReaction(reaction);
  }
  G4ITReactionPerTrackMap fReactionPerTrack;
  G4ITReactionPerTime fReactionPerTime;

  G4bool fSortByTime;
  static G4ThreadLocal G4ITReactionSet* fpInstance;
};

#endif /* G4ITREACTIONINFO_HH_ */
