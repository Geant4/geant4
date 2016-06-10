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
// $Id: G4TrackList.cc 93616 2015-10-27 08:59:17Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

/*
 * G4FastList.cc
 *
 *  Created on: 18 nov. 2014
 *      Author: kara
 */

#include "G4TrackList.hh"

//! SPECIFIC TO TRACKS
template<>
G4FastListNode<G4Track>* G4FastList<G4Track>::__GetNode(G4Track* __track)
{
  G4IT* __IT = GetIT(__track);
  G4FastListNode<G4Track>* __trackListNode = __IT->GetListNode();
  // TODO : complete the exception
  if (__trackListNode == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "This track " << GetIT(__track)->GetName();
    exceptionDescription << " was not connected to any trackList ";
    G4Exception("G4FastList<OBJECT>::Unflag", "G4TrackList003", FatalErrorInArgument,
        exceptionDescription);
    return 0;
  }
  return __trackListNode;
}

//! SPECIFIC TO TRACKS
template<>
void G4FastList<G4Track>::DeleteObject(G4Track* __track)
{
  if (!G4AllocatorList::GetAllocatorListIfExist()) return;

  G4Step* __step = const_cast<G4Step*>(__track->GetStep());
  if (__step)
  {
    if (__step->GetfSecondary())
      __step->DeleteSecondaryVector();
    delete __step;
  }
  delete __track;
}


template<>
  void G4FastListNode<G4Track>::DetachYourSelf()
  {
    if(fpObject)
    {
      GetIT(fpObject)->SetListNode(0);
    }
  }


//! SPECIFIC TO TRACKS
template<>
bool G4FastList<G4Track>::Holds(const G4Track* track) const
{
  return (GetIT(track)->GetListNode()->fListRef->fpList == this);
}

//! SPECIFIC TO TRACKS
template<>
G4FastListNode<G4Track>* G4FastList<G4Track>::Flag(G4Track* __track)
{
  G4IT* __iTrack = GetIT(__track);
  G4FastListNode<G4Track>* __trackListNode = __iTrack->GetListNode();

  if (__trackListNode != 0)
  {
    // Suggestion move the node to this list
    if (__trackListNode->fAttachedToList)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "This track " << __iTrack->GetName();
      exceptionDescription << " is already attached to a TrackList ";
      G4Exception("G4FastList<OBJECT>::Flag", "G4TrackList001",
          FatalErrorInArgument,
          exceptionDescription);
    }
  } else
  {
    __trackListNode = new G4FastListNode<G4Track>(__track);
    __iTrack->SetListNode(__trackListNode);
  }

  __trackListNode->fAttachedToList = true;
  __trackListNode->fListRef = fListRef;
  return __trackListNode;
}

//! SPECIFIC TO TRACKS
template<>
void G4FastList<G4Track>::CheckFlag(G4FastListNode<G4Track>* __trackListNode)
{
  if (__trackListNode->fListRef->fpList != this)
  {
    G4Track* track = __trackListNode->GetObject();
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "The track " << GetIT(track)->GetName()
        << " with trackID " << track->GetTrackID()
        << " is not correctly linked to a TrackList."
        << G4endl
        << "You are probably trying to withdraw this track "
        << "from the list but it probably does not belong to "
        << "this track list." << G4endl;
    G4Exception("G4FastList<OBJECT>::CheckFlag", "G4FastList002",
        FatalErrorInArgument, exceptionDescription);
  }
}

//! SPECIFIC TO TRACKS
template<>
G4FastListNode<G4Track>* G4FastList<G4Track>::EraseListNode(G4Track* __track)
{
  G4FastListNode<G4Track>* __node = Unflag(__track);
  G4FastListNode<G4Track>* __next = __node->GetNext();
  Unhook(__node);
  __node->DetachYourSelf();
  delete __node;
  return __next;
}


//! SPECIFIC TO TRACKS
template<>
G4FastListNode<G4Track>* G4FastList<G4Track>::GetNode(G4Track* __track)
{
  G4IT* __IT = GetIT(__track);
  G4FastListNode<G4Track>* __trackListNode = __IT->GetListNode();
  // TODO : complete the exception
  if (__trackListNode == 0)
  {
    return 0;
  }
  return __trackListNode;
}

//! SPECIFIC TO TRACKS
template<>
G4FastList<G4Track>* G4FastList<G4Track>::GetList(G4Track* __track)
{
  G4IT* __IT = GetIT(__track);
  G4FastListNode<G4Track>* __trackListNode = __IT->GetListNode();

  if(__trackListNode == 0) return 0;
  if(__trackListNode->fListRef == nullptr) return 0;

  return __trackListNode->fListRef->fpList;
}

