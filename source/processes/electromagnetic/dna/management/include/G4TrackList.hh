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
// $Id$
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4TRACKLIST_H
#define G4TRACKLIST_H

#include "globals.hh"
#include "G4ReferenceCountedHandle.hh"
#include <CLHEP/Utility/memory.h>
#include <vector>
#include "G4IT.hh"

class G4Track;
class G4TrackList;
class G4TrackManyList;
class G4TrackList_Boundary;

/** Comments :
 * - A track cannot belong to two different track lists
 * - Erase a given track is constant complexity
 * - This development was thought to be used together with G4IT
 */

struct _ListRef
{
  friend class G4TrackList;
  friend class G4TrackManyList;
protected:
  inline _ListRef(G4TrackList* __list) :
      fpTrackList(__list)
  {
    ;
  }

  G4TrackList* fpTrackList;
};

/**
 * G4TrackListNode is the entity actually stored
 * by the G4TrackList. A G4TrackListNode should
 * belong only to one list. Also, a track
 * should belong only to one list.
 */

class G4TrackListNode
{
  friend class G4TrackList;
  friend class G4TrackManyList;

public:
  G4Track* GetTrack()
  {
    return fpTrack;
  }
  G4TrackListNode* GetNext()
  {
    return fpNext;
  }
  G4TrackListNode* GetPrevious()
  {
    return fpPrevious;
  }
  bool IsAttached()
  {
    return fAttachedToList;
  }

protected:
  /** Default constructor */
  G4TrackListNode(G4Track* track = 0);
  /** Default destructor */
  ~G4TrackListNode();
  void SetNext(G4TrackListNode* node)
  {
    fpNext = node;
  }
  void SetPrevious(G4TrackListNode* node)
  {
    fpPrevious = node;
  }
  void SetAttachedToList(bool flag)
  {
    fAttachedToList = flag;
  }

  bool fAttachedToList;
  CLHEP::shared_ptr<_ListRef> fListRef;
  G4Track* fpTrack;
  G4TrackListNode* fpPrevious;
  G4TrackListNode* fpNext;
};

struct G4TrackList_iterator;
struct G4TrackManyList_iterator;

/**
 * G4TrackList is used by G4ITStepManager to save
 * G4IT tracks only. Its advantage lies to a fast
 * search of a track in this list.
 */

class G4TrackList
{
protected:
  G4int fNbTracks;
  G4TrackListNode * fpStart;
  G4TrackListNode * fpFinish;
  CLHEP::shared_ptr<_ListRef> fListRef;

  G4TrackListNode fBoundary;
  // Must be empty and link to the last non-empty node of the list
  // and to the first non-empty node of the list (begin())
  // The iterator returned by end() is linked to this empty node

public:
  typedef G4TrackList_iterator iterator;

  G4TrackList();
  ~G4TrackList();

  inline G4Track* back()
  {
    if (fNbTracks != 0) return fpFinish->GetTrack();
    else return 0;
  }
  inline G4int size() const
  {
    return fNbTracks;
  }
  inline bool empty() const;
  iterator insert(iterator /*position*/, G4Track*);
  inline iterator begin();
  inline iterator end();
  /**
   * return an iterator that contains an empty node
   * use for boundary checking only
   */

  bool Holds(const G4Track*) const;

  inline void push_front(G4Track* __track);
  inline void push_back(G4Track* __track);
  G4Track* pop_back();

  void remove(G4Track*);

  iterator pop(G4Track*);
  iterator pop(G4TrackListNode*);
  iterator pop(iterator __first, iterator __last);
  iterator erase(G4Track*);
  /**
   * Complexity = constant
   * By storing the node inside the object, we avoid
   * searching through all the container.
   */

  iterator erase(iterator __first, iterator __last);
  /**
   * Complexity = linear in size between __first and __last
   */
  void transferTo(G4TrackList*);
  /**
   * Complexity = constant
   */

  static G4TrackListNode* GetNode(G4Track*);
  static G4TrackList* GetTrackList(G4Track*);
  static G4TrackList* GetTrackList(G4TrackListNode* __trackListNode);
  static void Pop(G4Track*);

protected:
  G4TrackListNode* CreateNode(G4Track*);
  static G4TrackListNode* __GetNode(G4Track*);
  G4TrackListNode* Flag(G4Track*);
  G4TrackListNode* Unflag(G4Track*);
  void Unflag(G4TrackListNode* __trackListNode);
  void CheckFlag(G4TrackListNode*);
  void DeleteTrack(G4Track*);

  void Hook(G4TrackListNode* /*position*/, G4TrackListNode* /*toHook*/);
  void Unhook(G4TrackListNode*);
  G4TrackListNode* EraseTrackListNode(G4Track*);

private:
  G4TrackList(const G4TrackList& other);
  G4TrackList & operator=(const G4TrackList &right);
  G4int operator==(const G4TrackList &right) const;
  G4int operator!=(const G4TrackList &right) const;
};


/*
 * Roll over many list as if it was one
 */
class G4TrackManyList
{
protected:
  std::vector<G4TrackList*> fAssociatedLists;
  // TODO use "marked list" insted of vector

public:
  typedef G4TrackManyList_iterator iterator;

  G4TrackManyList() :
      fAssociatedLists()
  {
  }

  inline void Add(G4TrackList* __list)
  {
    if (__list == 0) return;
    fAssociatedLists.push_back(__list); // TODO use the table doubling tech
  }

  inline bool Holds(G4Track* __track) const
  {
    std::vector<G4TrackList*>::const_iterator __it = fAssociatedLists.begin();
    std::vector<G4TrackList*>::const_iterator __end = fAssociatedLists.end();
    for(;__it != __end ; __it++) if((*__it)->Holds(__track)) return true;
    return false;
  }

  inline size_t size() const
  {
    size_t __size(0);
    std::vector<G4TrackList*>::const_iterator __it = fAssociatedLists.begin();
    std::vector<G4TrackList*>::const_iterator __end = fAssociatedLists.end();
    for(;__it != __end ; __it++) __size+= (*__it)->size();
    return __size;
  }

  inline void clear()
  {
    fAssociatedLists.clear();
  }

  inline iterator begin();
  inline iterator end();

  void pop(G4Track*);

};

/**
 * G4TrackList_iterator enables to go through
 * the tracks contained by a list.
 */
struct G4TrackList_iterator
{
  friend class G4TrackList;
  typedef G4TrackList_iterator _Self;
  typedef G4TrackListNode _Node;

  G4TrackList_iterator() :
      fpNode(0)
  {
  }

  explicit G4TrackList_iterator(_Node* __x) :
      fpNode(__x)
  {
  }

  _Node* GetNode()
  {
    return fpNode;
  }

  G4Track*
  operator*();

  const G4Track*
  operator*() const;

  G4Track*
  operator->();

  const G4Track*
  operator->() const;

  _Self&
  operator++()
  {
    fpNode = fpNode->GetNext();
    return *this;
  }

  _Self operator++(int)
  {
    _Self __tmp = *this;
    fpNode = fpNode->GetNext();
    return __tmp;
  }

  _Self&
  operator--()
  {
    fpNode = fpNode->GetPrevious();
    return *this;
  }

  _Self operator--(int)
  {
    _Self __tmp = *this;
    fpNode = fpNode->GetPrevious();
    return __tmp;
  }

  bool operator==(const _Self& __x) const
  {
    return (fpNode == __x.fpNode);
  }

  bool operator!=(const _Self& __x) const
  {
    return (fpNode != __x.fpNode);
  }

private:
  // The only member points to the G4TrackList_iterator element.
  _Node* fpNode;
};

struct G4TrackManyList_iterator
{
  friend class G4TrackManyList;
  typedef G4TrackManyList_iterator _Self;
  typedef G4TrackListNode _Node;

  G4TrackList_iterator fIterator;
  std::vector<G4TrackList*>::iterator fCurrentListIt;
  std::vector<G4TrackList*>* fLists;

  G4TrackManyList_iterator() :
      fIterator(0), fLists(0)
  {
  }

  explicit G4TrackManyList_iterator(G4TrackList_iterator __x,
                                    std::vector<G4TrackList*>::iterator __it,
                                    std::vector<G4TrackList*>* __lists) :
      fIterator(__x), fCurrentListIt(__it), fLists(__lists)
  {
  }

  G4TrackManyList_iterator(const G4TrackManyList_iterator& __x) :
      fIterator(__x.fIterator),
      fCurrentListIt(__x.fCurrentListIt),
      fLists(__x.fLists)
  {
  }

  _Node* GetNode()
  {
    return fIterator.GetNode();
  }

  G4TrackList* GetTrackList()
  {
    return *fCurrentListIt;
  }

  G4Track* operator*()
  {
    return *fIterator;
  }
  const G4Track* operator*() const
  {
      return *fIterator;
  }
  G4Track* operator->()
  {
      return *fIterator;
  }
  const G4Track* operator->() const
  {
      return *fIterator;
  }

  _Self UpdateToNextValidList();
  _Self& operator++();

  _Self operator++(int)
  {
    return operator++();
  }

  _Self&
  operator--()
  {
    if(fLists->empty())
    {
//      G4cout << "Vide" << G4endl;
      fIterator = G4TrackList_iterator();
      return *this;
    }
    if(fCurrentListIt == fLists->begin())
    {
//      G4cout << "Sur le premier paquet" << G4endl;
      if(fIterator == (*fCurrentListIt)->begin())
      {
//        G4cout << "retourne NULL" << G4endl;
        fIterator = G4TrackList_iterator();
        return *this;
      }
    }

    if(fCurrentListIt == fLists->end())
    {
//      G4cout << "Sur le dernier paquet" << G4endl;
      fCurrentListIt--;
      fIterator = (*fCurrentListIt)->end();
    }
    else if(fIterator == (*fCurrentListIt)->begin())
    {
      fCurrentListIt--;
      fIterator = (*fCurrentListIt)->end();
    }

    fIterator--;

    while(( (*fCurrentListIt)->empty()|| fIterator.GetNode() == 0
        || fIterator.GetNode()->GetTrack() == 0)
        && fCurrentListIt != fLists->begin())
    {
//      G4cout << "!:::!"<< G4endl;
      fIterator = (*fCurrentListIt)->begin();
      fCurrentListIt--;
      fIterator = (*fCurrentListIt)->end();
      fIterator--;
    }

    if(fIterator.GetNode() == 0
            && fCurrentListIt == fLists->begin())
    {
      fIterator = G4TrackList_iterator();
      return *this;
    }

     return *this;
  }

  _Self operator--(int)
  {
    return operator--();
  }

  bool operator==(const _Self& __x) const
  {
    return (fIterator == __x.fIterator && fCurrentListIt == __x.fCurrentListIt);
  } // Fast check

  bool operator!=(const _Self& __x) const
  {
    return !(this->operator ==(__x));
  }

protected:
  void HasReachedEnd()
  {
    if(fLists->empty() == false)
    {
      fIterator = (*fLists->rbegin())->end();
    }
    else
    {
      fIterator = G4TrackList_iterator();
    }
  }
};

inline bool G4TrackList::empty() const
{
  return (fNbTracks == 0);
}

inline G4TrackList::iterator G4TrackList::begin()
{
  return iterator(fpStart);
}

inline G4TrackList::iterator G4TrackList::end()
{
  return iterator(&(fBoundary));
}
// return an iterator that contains an empty node
// use for boundary checking only

inline void G4TrackList::push_front(G4Track* track)
{
  insert(begin(), track);
}

inline void G4TrackList::push_back(G4Track* track)
{
  insert(end(), track);
}

inline G4TrackManyList::iterator G4TrackManyList::begin()
{
  if (fAssociatedLists.empty())
  {
    return G4TrackManyList_iterator(G4TrackList_iterator(),
                                    fAssociatedLists.end(), &fAssociatedLists);
  }

  return G4TrackManyList_iterator(fAssociatedLists[0]->begin(),
                                  fAssociatedLists.begin(), &fAssociatedLists);
}

inline G4TrackManyList::iterator G4TrackManyList::end()
{
  if (fAssociatedLists.empty())
  {
    return G4TrackManyList_iterator(G4TrackList_iterator(),
                                    fAssociatedLists.end(), &fAssociatedLists);
  }

  return G4TrackManyList_iterator((*fAssociatedLists.rbegin())->end(),
                                  fAssociatedLists.end(), &fAssociatedLists);
}

#endif // G4TRACKLIST_H
