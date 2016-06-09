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
// $Id: G4TrackList.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4TRACKLIST_H
#define G4TRACKLIST_H

#include "globals.hh"
#include "G4ReferenceCountedHandle.hh"

class G4Track;
class G4TrackList;
class G4TrackList_Boundary;

/** Comments :
* - A track cannot belong to two different track lists
* - Erase a given track is constant complexity
* - This development was thought to be used together with G4IT
*/

struct _ListRef
{
    friend class G4TrackList ;
protected :
    inline _ListRef(G4TrackList* __list) :fpTrackList(__list)
    {;}

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

public :
    G4Track* GetTrack() { return fpTrack; }
    G4TrackListNode* GetNext() { return fpNext;}
    G4TrackListNode* GetPrevious() { return fpPrevious;}
    bool IsAttached() { return fAttachedToList;}

protected:
    /** Default constructor */
    G4TrackListNode(G4Track* track= 0);
    /** Default destructor */
    ~G4TrackListNode();
    void SetNext(G4TrackListNode* node) { fpNext = node;}
    void SetPrevious(G4TrackListNode* node) { fpPrevious = node;}
    void SetAttachedToList(bool flag) { fAttachedToList = flag;}

    bool fAttachedToList;
    G4ReferenceCountedHandle<_ListRef> fListRef;
    G4Track*         fpTrack;
    G4TrackListNode* fpPrevious;
    G4TrackListNode* fpNext;
};

struct G4TrackList_iterator ;

/**
  * G4TrackList is used by G4ITStepManager to save
  * G4IT tracks only. Its advantage lies to a fast
  * search of a track in this list.
  */

class G4TrackList
{
private :
    G4int                               fNbTracks;
    G4TrackListNode *                   fpStart;
    G4TrackListNode *                   fpFinish;
    G4ReferenceCountedHandle<_ListRef>  fListRef;

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
        if(fNbTracks != 0)
            return fpFinish->GetTrack();
        else return 0 ;
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

protected:
    G4TrackListNode* CreateNode(G4Track*);
    G4TrackListNode* Flag(G4Track*);
    G4TrackListNode* Unflag(G4Track*);
    void CheckFlag(G4TrackListNode*);
    void DeleteTrack(G4Track*);

    void Hook(G4TrackListNode* /*position*/, G4TrackListNode* /*toHook*/);
    void Unhook(G4TrackListNode*);
    G4TrackListNode* EraseTrackListNode(G4Track*);

private:
    G4TrackList(const G4TrackList& other);
    G4TrackList & operator=
    (const G4TrackList &right);
    G4int operator==(const G4TrackList &right) const;
    G4int operator!=(const G4TrackList &right) const;
};

/**
  * G4TrackList_iterator enables to go through
  * the tracks contained by a list.
  */

struct G4TrackList_iterator
{
    friend class G4TrackList;
    typedef G4TrackList_iterator                    _Self;
    typedef G4TrackListNode                          _Node;

    G4TrackList_iterator()
        : fpNode() { }

    explicit
    G4TrackList_iterator(_Node* __x)
        : fpNode(__x) { }

    _Node* GetNode()
    { return fpNode; }

    G4Track*
    operator*();

    const G4Track*
    operator*() const;

    G4Track*
    operator->() ;

    const G4Track*
    operator->() const;

    _Self&
    operator++()
    {
        fpNode = fpNode->GetNext();
        return *this;
    }

    _Self
    operator++(int)
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

    _Self
    operator--(int)
    {
        _Self __tmp = *this;
        fpNode = fpNode->GetPrevious();
        return __tmp;
    }

    bool
    operator==(const _Self& __x) const
    { return (fpNode == __x.fpNode); }

    bool
    operator!=(const _Self& __x) const
    {
        return (fpNode != __x.fpNode);
    }

private:
    // The only member points to the G4TrackList_iterator element.
    _Node* fpNode;
};

inline bool G4TrackList::empty() const
{ return (fNbTracks == 0); }


inline G4TrackList::iterator G4TrackList::begin()
{ return iterator(fpStart); }

inline G4TrackList::iterator G4TrackList::end()
{ return iterator( &(fBoundary) ); }
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

#endif // G4TRACKLIST_H
