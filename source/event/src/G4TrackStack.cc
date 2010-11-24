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
// $Id: G4TrackStack.cc,v 1.9 2010-11-24 22:56:57 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4TrackStack.hh"
#include "G4SmartTrackStack.hh"
#include "G4VTrajectory.hh"

G4TrackStack::G4TrackStack()
:n_stackedTrack(0),firstStackedTrack(0),lastStackedTrack(0)
{
  maxNTracks = 0;
}

G4TrackStack::~G4TrackStack()
{
  if( n_stackedTrack != 0 )
  {
    G4StackedTrack * aStackedTrack = firstStackedTrack;
    G4StackedTrack * nextStackedTrack;

    // delete tracks in the stack
    while( aStackedTrack != 0 )
    {
      nextStackedTrack = aStackedTrack->GetNext();
      delete aStackedTrack->GetTrack();
      delete aStackedTrack->GetTrajectory();
      delete aStackedTrack;
      aStackedTrack = nextStackedTrack;
    }
  }
}

const G4TrackStack & G4TrackStack::operator=(const G4TrackStack &right) 
{
  n_stackedTrack = right.n_stackedTrack;
  firstStackedTrack = right.firstStackedTrack;
  lastStackedTrack = right.lastStackedTrack;
  return *this; 
}

int G4TrackStack::operator==(const G4TrackStack &right) const
{ return (firstStackedTrack==right.firstStackedTrack); }
int G4TrackStack::operator!=(const G4TrackStack &right) const
{ return (firstStackedTrack!=right.firstStackedTrack); }

void G4TrackStack::TransferTo(G4TrackStack * aStack)
{
  if(n_stackedTrack==0) return;

  if(aStack->n_stackedTrack == 0)
  {
    *aStack = *this;
  }
  else
  {
    aStack->lastStackedTrack->SetNext( firstStackedTrack );
    firstStackedTrack->SetPrevious( aStack->lastStackedTrack );
    aStack->lastStackedTrack = lastStackedTrack;
    aStack->n_stackedTrack += n_stackedTrack;
  }

  n_stackedTrack = 0;
  firstStackedTrack = 0;
  lastStackedTrack = 0;
}

void G4TrackStack::TransferTo(G4SmartTrackStack * aStack)
{
  while(n_stackedTrack)
  { aStack->PushToStack(PopFromStack()); }
}

G4StackedTrack * G4TrackStack::PopFromStack()
{
  if( n_stackedTrack == 0 ) return 0;
  G4StackedTrack * aStackedTrack = lastStackedTrack;
  GrabFromStack( aStackedTrack );
  return aStackedTrack;
}

void G4TrackStack::PushToStack( G4StackedTrack * aStackedTrack )
{
  if(aStackedTrack)
  {
    if( n_stackedTrack == 0 )
    {
      aStackedTrack->SetPrevious( 0 );
      firstStackedTrack = aStackedTrack;
    }
    else
    {
      lastStackedTrack->SetNext( aStackedTrack );
      aStackedTrack->SetPrevious( lastStackedTrack );
    }
    lastStackedTrack = aStackedTrack;
    n_stackedTrack++;
    if(n_stackedTrack>maxNTracks) maxNTracks = n_stackedTrack;
  }
}

void G4TrackStack::GrabFromStack( G4StackedTrack * aStackedTrack )
{
  if( n_stackedTrack == 1 )
  {
    firstStackedTrack = 0;
    lastStackedTrack = 0;
  }
  else
  {
    if( aStackedTrack == firstStackedTrack )
    {
      firstStackedTrack = aStackedTrack->GetNext();
      firstStackedTrack->SetPrevious( 0 );
    }
    else
    {
      if( aStackedTrack == lastStackedTrack )
      {
        lastStackedTrack = aStackedTrack->GetPrevious();
        lastStackedTrack->SetNext( 0 );
      }
      else
      {
        aStackedTrack->GetPrevious()
          ->SetNext( aStackedTrack->GetNext() );
        aStackedTrack->GetNext()
          ->SetPrevious( aStackedTrack->GetPrevious() );
      }
    }
  }
  n_stackedTrack--;
}

void G4TrackStack::clear()
{
  G4StackedTrack * aStackedTrack = firstStackedTrack;
  G4StackedTrack * nextStackedTrack;

  if ( n_stackedTrack == 0 ) return;

  // delete tracks in the stack
  while( aStackedTrack != 0 )
  {
    nextStackedTrack = aStackedTrack->GetNext();
    delete aStackedTrack->GetTrack();
    delete aStackedTrack->GetTrajectory();
    delete aStackedTrack;
    aStackedTrack = nextStackedTrack;
  }
  n_stackedTrack = 0;
  firstStackedTrack = 0;
  lastStackedTrack = 0;
}


