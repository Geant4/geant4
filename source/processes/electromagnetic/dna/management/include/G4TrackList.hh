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
    inline _ListRef(G4TrackList* __list) :fTrackList(__list)
    {;}

    G4TrackList* fTrackList;
};

class G4TrackListNode
{
     friend class G4TrackList;
     friend class G4TrackList_Boundary;

public :
     G4Track* GetTrack() { return fTrack; }
     G4TrackListNode* GetNext() { return fNext;}
     G4TrackListNode* GetPrevious() { return fPrevious;}

protected:
        /** Default constructor */
        G4TrackListNode(G4Track* track= 0);
        /** Default destructor */
        ~G4TrackListNode();

        bool fAttachedToList;
        G4ReferenceCountedHandle<_ListRef> fListRef;
        void SetNext(G4TrackListNode* node) { fNext = node;}
        void SetPrevious(G4TrackListNode* node) { fPrevious = node;}

        G4Track*         fTrack;
        G4TrackListNode* fPrevious;
        G4TrackListNode* fNext;
};

struct G4TrackList_iterator ;
class G4TrackList
{
    private :
        G4int                               fNbTracks;
        G4TrackListNode *                   fStart;
        G4TrackListNode *                   fFinish;
        G4ReferenceCountedHandle<_ListRef>  fListRef;

        struct G4TrackList_Boundary
        {
            G4TrackListNode fEmptyNode;
            // Must be empty and link to the last non-empty node of the list
            // and to the first non-empty node of the list (begin())
            inline G4TrackList_Boundary():fEmptyNode(){;}
            ~G4TrackList_Boundary(){;}
        };

        G4TrackList_Boundary fBoundary;

    public:
        typedef G4TrackList_iterator iterator;

        G4TrackList();
        ~G4TrackList();

        inline G4Track* back()
        {
            if(fNbTracks != 0)
             return fFinish->GetTrack();
            else return 0 ;
        }
        inline G4int size()
        {
            return fNbTracks;
        }
        inline bool empty();
        iterator insert(iterator /*position*/, G4Track*);
        inline iterator begin();
        inline iterator end();
        // return an iterator that contains an empty node
        // use for boundary checking only

        // TODO
        /*inline reverse_iterator rbegin()
        inline reverse_iterator rend();
        // return an iterator that contains an empty node
        // use for boundary checking only*/

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

struct G4TrackList_iterator
{
    typedef G4TrackList_iterator                    _Self;
    typedef G4TrackListNode                          _Node;

    G4TrackList_iterator()
    : fNode() { }

    explicit
    G4TrackList_iterator(_Node* __x)
    : fNode(__x) { }

    _Node* GetStackedTrack()
    { return fNode; }

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
        fNode = fNode->GetNext();
        return *this;
    }

    _Self
    operator++(int)
    {
        _Self __tmp = *this;
        fNode = fNode->GetNext();
        return __tmp;
    }

    _Self&
    operator--()
    {
        fNode = fNode->GetPrevious();
        return *this;
    }

    _Self
    operator--(int)
    {
        _Self __tmp = *this;
        fNode = fNode->GetPrevious();
        return __tmp;
    }

    bool
    operator==(const _Self& __x) const
    { return (fNode == __x.fNode); }

    bool
    operator!=(const _Self& __x) const
    {
         return (fNode != __x.fNode);
    }

    // The only member points to the G4TrackList_iterator element.
    _Node* fNode;
};

inline bool G4TrackList::empty()
{ return (fNbTracks == 0); }


inline G4TrackList::iterator G4TrackList::begin()
{ return iterator(fStart); }

G4TrackList::iterator G4TrackList::end()
{ return iterator( &(fBoundary.fEmptyNode) ); }
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
