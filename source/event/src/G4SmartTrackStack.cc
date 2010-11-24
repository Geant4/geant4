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
// $Id: G4SmartTrackStack.cc,v 1.5 2010-11-24 22:56:57 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4SmartTrackStack.hh"
#include "G4VTrajectory.hh"

G4SmartTrackStack::G4SmartTrackStack()
:fTurn(0),nTurn(5)
{
  for(int i=0;i<nTurn;i++)
  { stacks[i] = new G4TrackStack(); }
  // If entry of one sub-stack exceeds safetyValve1, we will stick
  // to that sub-stack until entry of that sub-stack goes down
  // to safetyValve2.
  nStick = 100;
  safetyValve1 = 3000; 
  safetyValve2 = safetyValve1 - nStick;
  maxNTracks = 0;
}

G4SmartTrackStack::~G4SmartTrackStack()
{
  for(int i=0;i<nTurn;i++)
  { delete stacks[i]; }
}

const G4SmartTrackStack & G4SmartTrackStack::operator=(const G4SmartTrackStack &) 
{ return *this; }
int G4SmartTrackStack::operator==(const G4SmartTrackStack &right) const
{ return (this==&right); }
int G4SmartTrackStack::operator!=(const G4SmartTrackStack &right) const
{ return (this!=&right); }

void G4SmartTrackStack::TransferTo(G4TrackStack * aStack)
{
  for(int i=0;i<nTurn;i++)
  { stacks[i]->TransferTo(aStack); }
}

G4StackedTrack * G4SmartTrackStack::PopFromStack()
{
  if( n_stackedTrack() == 0 ) return 0;
  G4StackedTrack * aStackedTrack = 0;
  while(!aStackedTrack)
  {
    if(stacks[fTurn]->GetNTrack()==0)
    {
      fTurn = (fTurn+1)%nTurn;
      //G4cout<<"++++++++ Shift to Stack ["<<fTurn<<"] with "<<stacks[fTurn]->GetNTrack()<<" stacked tracks."<<G4endl;
    }
    else
    { aStackedTrack = stacks[fTurn]->PopFromStack(); }
  }
  return aStackedTrack;
}

#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

void G4SmartTrackStack::PushToStack( G4StackedTrack * aStackedTrack )
{
  static G4ParticleDefinition* neutDef = G4Neutron::Definition();
  static G4ParticleDefinition* elecDef = G4Electron::Definition();
  static G4ParticleDefinition* gammDef = G4Gamma::Definition();
  static G4ParticleDefinition* posiDef = G4Positron::Definition();

  if(!aStackedTrack) return;

  G4int iDest = 0;
  if( aStackedTrack->GetTrack()->GetParentID() == 0 )
  {
    // We have a primary track, which should go first.
    fTurn = 0; // reseting the turn
  }
  else
  {
    G4ParticleDefinition* partDef = aStackedTrack->GetTrack()->GetDefinition();
    if(partDef==neutDef)
    { iDest = 1; }
    else if(partDef==elecDef)
    { iDest = 2; }
    else if(partDef==gammDef)
    { iDest = 3; }
    else if(partDef==posiDef)
    { iDest = 4; }
  }

  stacks[iDest]->PushToStack(aStackedTrack);
  if(stacks[iDest]->GetNTrack()>safetyValve1)
  {
    // Too many tracks in the stack. Process tracks in this stack first
    // unless the current stack also have too many tracks.
    if(stacks[fTurn]->GetNTrack()<safetyValve2)
    {
      fTurn = iDest;
      safetyValve2 = stacks[iDest]->GetNTrack() - nStick;
      //G4cout<<"++++++++ Shift to Stack ["<<fTurn<<"] with "<<stacks[fTurn]->GetNTrack()<<" stacked tracks."<<G4endl;
    }
  }

  if(n_stackedTrack()>maxNTracks) maxNTracks = n_stackedTrack();
}

void G4SmartTrackStack::clear()
{
  for(int i=0;i<nTurn;i++)
  { stacks[i]->clear(); }
}


