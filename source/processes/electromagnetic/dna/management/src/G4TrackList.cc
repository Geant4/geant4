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
// $Id: G4TrackList.cc 65022 2012-11-12 16:43:12Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4TrackList.hh"
#include "G4IT.hh"
#include "G4Track.hh"

using namespace std;

//***********************************************************
// TrackList_iterator
G4Track*
G4TrackList_iterator::operator*()
{ return fpNode->GetTrack(); }

G4Track*
G4TrackList_iterator::operator->()
{ return fpNode->GetTrack(); }

const G4Track*
G4TrackList_iterator::operator*() const
{ return fpNode->GetTrack(); }

const G4Track*
G4TrackList_iterator::operator->() const
{ return fpNode->GetTrack(); }


//***********************************************************
// TrackNodeList

G4TrackListNode::G4TrackListNode(G4Track* track) :
    fpTrack(track),
    fpPrevious(0),
    fpNext(0)
{
    fAttachedToList = false;
}

G4TrackListNode::~G4TrackListNode()
{;}

//***********************************************************

G4TrackList::G4TrackList() : fBoundary()
{
    fListRef    = new _ListRef(this);
    fpStart      = 0;
    fpFinish     = 0;
    fNbTracks   = 0 ;
    fBoundary.SetPrevious(&fBoundary);
    fBoundary.SetNext(&fBoundary);
    fBoundary.fAttachedToList = true;
}

// should not be used
G4TrackList::G4TrackList(const G4TrackList& /*other*/) : fBoundary()
{
    // One track should not belong to two different trackLists

    fpFinish = 0;
    fpStart = 0;
    fNbTracks = 0;
    fListRef = 0;
}

G4TrackList& G4TrackList::operator=(const G4TrackList& other)
{
    // One track should not belong to two different trackList
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    return *this;
}

G4TrackList::~G4TrackList()
{
    if( fNbTracks != 0 )
    {
        G4TrackListNode * __stackedTrack = fpStart;
        G4TrackListNode * __nextStackedTrack;

        // delete tracks in the stack
        while(  __stackedTrack && __stackedTrack != &(fBoundary) )
        {
            __nextStackedTrack = __stackedTrack->GetNext();
            G4Track* __track = __stackedTrack->GetTrack();

            delete __stackedTrack;
            __stackedTrack = 0;

            if(__track)
            {
                //////////////
                DeleteTrack(__track);
                __track = 0;
                //////////////
            }

            __stackedTrack = __nextStackedTrack;
        }
    }
    fNbTracks = 0;
}

bool G4TrackList::Holds(const G4Track* track) const
{
    return (GetIT(track)->GetTrackListNode()->fListRef->fpTrackList == this)  ;
}

G4TrackListNode* G4TrackList::Flag(G4Track* __track)
{
    G4IT* __iTrack = GetIT(__track);
    G4TrackListNode* __trackListNode = __iTrack->GetTrackListNode();

    if(__trackListNode != 0)
    {
        // Suggestion move the node to this list
        if(__trackListNode->fAttachedToList)
        {
            G4ExceptionDescription exceptionDescription ;
            exceptionDescription << "This track "<< __iTrack->GetName() ;
            exceptionDescription << " is already attached to a TrackList ";
            G4Exception("G4TrackList::Flag","G4TrackList001",
                        FatalErrorInArgument,exceptionDescription);
         }
    }
    else
    {
        __trackListNode = new G4TrackListNode(__track);
        __iTrack->SetTrackListNode(__trackListNode);
    }

    __trackListNode->fAttachedToList = true;
    __trackListNode->fListRef = fListRef;
    return __trackListNode;
}

G4TrackListNode* G4TrackList::CreateNode(G4Track* __track)
{
    G4TrackListNode* __trackListNode = Flag(__track);
    return __trackListNode;
}

void G4TrackList::Hook(G4TrackListNode* __position, G4TrackListNode* __toHook)
{
    if(fNbTracks == 0)
    {
        // DEBUG
        //        G4cout << "fNbTracks == 0" << G4endl;
        fpStart = __toHook;
        fpFinish = __toHook;
        __toHook->SetNext(&fBoundary);
        __toHook->SetPrevious(&fBoundary);
        fBoundary.SetNext(__toHook);
        fBoundary.SetPrevious(__toHook);
    }
    else if( __position == &fBoundary)
    {
        // DEBUG
        //        G4cout << "__position == &fBoundary" << G4endl;
        fpFinish->SetNext( __toHook );
        __toHook->SetPrevious( fpFinish );

        __toHook->SetNext(&fBoundary);
        fBoundary.SetPrevious( __toHook );

        fpFinish = __toHook;
    }
    else if( __position == fpStart )
    {
        // DEBUG
        //        G4cout << "__position == fStart" << G4endl;
        __toHook->SetPrevious( &fBoundary );
        fBoundary.SetNext(__toHook);
        __toHook->SetNext(fpStart);
        fpStart->SetPrevious(__toHook);
        fpStart = __toHook;
    }
    else
    {
        // DEBUG
        //        G4cout << "else" << G4endl;
        G4TrackListNode* __previous = __position->GetPrevious();
        __toHook->SetPrevious(__previous);
        __toHook->SetNext(__position);
        __position->SetPrevious(__toHook);
        __previous->SetNext(__toHook);
    }

    fNbTracks++;
}

void G4TrackList::Unhook(G4TrackListNode* __toUnHook)
{
    G4TrackListNode* __previous = __toUnHook->GetPrevious();
    G4TrackListNode* __next = __toUnHook->GetNext();

    __toUnHook->SetPrevious(0);
    __toUnHook->SetNext(0);

    if( fNbTracks == 1 )
    {
        fpStart = 0;
        fpFinish = 0;
    }
    else
    {
        if(__toUnHook == fpFinish)
        {
            fpFinish = __previous;
        }
        if(__toUnHook == fpStart)
        {
            fpStart = __next;
        }
    }

    // There should be always a __next and a __previous
    // because it is a circle link
    __next->SetPrevious(__previous);
    __previous->SetNext(__next);

    fNbTracks--;
}

G4TrackList::iterator G4TrackList::insert(G4TrackList::iterator __position, G4Track* __track)
{
    G4TrackListNode* __node = CreateNode(__track);
    Hook(__position.fpNode, __node);
    return iterator(__node);
}

//____________________________________________________________________
//
//                      WITHDRAW FROM LIST
//____________________________________________________________________
void G4TrackList::CheckFlag(G4TrackListNode* __trackListNode)
{
    if(__trackListNode -> fListRef->fpTrackList != this)
    {
        G4Track* track = __trackListNode->GetTrack();
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription
                << "The track "<< GetIT(track)->GetName()
                << " with trackID " << track->GetTrackID()
                << " is not correctly linked to a TrackList."
                << G4endl
                << "You are probably trying to withdraw this track "
                << "from the list but it probably does not belong to "
                << "this track list." << G4endl;
        G4Exception("G4TrackList::CheckFlag","G4TrackList002",
                    FatalErrorInArgument,exceptionDescription);
    }
}

G4TrackListNode* G4TrackList::Unflag(G4Track* __track)
{
    G4IT* __IT = GetIT(__track);
    G4TrackListNode* __trackListNode = __IT->GetTrackListNode();
    // TODO : complete the exception
    if(__trackListNode == 0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "This track "<< GetIT(__track)->GetName() ;
        exceptionDescription << " was not connected to any trackList ";
        G4Exception("G4TrackList::Unflag","G4TrackList003",
                    FatalErrorInArgument,exceptionDescription);
        return 0;
    }
    CheckFlag(__trackListNode);
    __trackListNode->fAttachedToList = false;
    __trackListNode->fListRef = 0;
    return __trackListNode;
}

G4Track* G4TrackList::pop_back()
{
    if( fNbTracks == 0 ) return 0;
    G4TrackListNode * __aStackedTrack = fpFinish;
    Unhook( __aStackedTrack );
    Unflag( __aStackedTrack->GetTrack() );
    return __aStackedTrack->GetTrack();
}

G4TrackList::iterator G4TrackList::pop(G4Track* __track)
{
    G4TrackListNode* __node = Unflag(__track);
    iterator __next(__node->GetNext());
    Unhook(__node);
    return __next;
}

G4TrackListNode* G4TrackList::EraseTrackListNode(G4Track* __track)
{
    G4TrackListNode* __node = Unflag(__track);
    GetIT(__track)->SetTrackListNode(0);
    G4TrackListNode* __next = __node->GetNext();
    Unhook(__node);
    delete __node;
    return __next;
}

void G4TrackList::DeleteTrack(G4Track* __track)
{
    G4Step* __step = const_cast<G4Step*>(__track->GetStep());
    if(__step)
    {
        if(__step->GetfSecondary()) __step->DeleteSecondaryVector();
        delete __step;
    }
    delete __track;
}

G4TrackList::iterator G4TrackList::erase(G4Track* __track)
{
    G4TrackListNode* __next_node = EraseTrackListNode(__track);
    //////////////////
    DeleteTrack(__track);
    __track = 0;
    //////////////////
    iterator __next(__next_node);
    return __next;
}

void G4TrackList::remove(G4Track* __track)
{
    this->erase(__track);
}

G4TrackList::iterator
G4TrackList::pop(iterator __first, iterator __last)
{
    if(fNbTracks == 0) return iterator(&fBoundary);

    while (__first != __last)
    {
        if(__first . fpNode)
            __first = pop(*__first);
    }
    return __last;
}


G4TrackList::iterator
G4TrackList::erase(iterator __first, iterator __last)
{
    if(fNbTracks == 0) return iterator(&fBoundary);

    while (__first != __last)
    {
        if(__first . fpNode)
            __first = erase(*__first);
    }
    return __last;
}

void G4TrackList::transferTo(G4TrackList* __destination)
{
    if(fNbTracks==0) return;

    if(__destination->fNbTracks == 0)
    {
        __destination->fpStart       =    this->fpStart ;
        __destination->fpFinish      =    this->fpFinish ;
        __destination->fNbTracks    =    this->fNbTracks;

        __destination->fBoundary.SetNext(fpStart);
        __destination->fBoundary.SetPrevious(fpFinish);

        __destination->fpFinish->SetNext(&__destination->fBoundary);
        __destination->fpStart->SetPrevious(&__destination->fBoundary);
    }
    else
    {
        this->fpStart->SetPrevious(__destination->fpFinish);
        __destination->fpFinish->SetNext(this->fpStart);
        __destination->fBoundary.SetPrevious(this->fpFinish);
        this->fpFinish->SetNext(&__destination->fBoundary);

        __destination->fpFinish = this->fpFinish;
        __destination->fNbTracks += this->fNbTracks;
    }

    fNbTracks = 0;
    fpStart = 0;
    fpFinish = 0;
    this->fBoundary.SetPrevious(&this->fBoundary);
    this->fBoundary.SetNext(&this->fBoundary);

    fListRef->fpTrackList = __destination;
}
