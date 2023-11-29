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
// G4SmartTrackStack class implementation
//
// Author: S.Kamperis - 4 October 2012
// --------------------------------------------------------------------

#include "G4SmartTrackStack.hh"
#include "G4VTrajectory.hh"
#include "G4Track.hh"

void G4SmartTrackStack::dumpStatistics()
{
  // Print to stderr so that we can split stats output from normal
  // output of Geant4 which is typically being printed to stdout
  for (G4int i=0; i<nTurn; ++i)
  {
    G4cerr << stacks[i]->GetNTrack() << " ";
    G4cerr << stacks[i]->getTotalEnergy() << " ";
  }
  G4cerr << G4endl;
}

G4SmartTrackStack::G4SmartTrackStack() 
{
  for(G4int i=0; i<nTurn; ++i)
  {
    stacks[i] = new G4TrackStack(5000);
    energies[i] = 0.;
  }
}

G4SmartTrackStack::~G4SmartTrackStack()
{
  for (auto* sp : stacks)
  {
    delete sp;
  }
}

void G4SmartTrackStack::TransferTo(G4TrackStack* aStack)
{
  for (auto* sp : stacks)
  {
    sp->TransferTo(aStack);
  }
  nTracks = 0;
}

G4StackedTrack G4SmartTrackStack::PopFromStack()
{
  G4StackedTrack aStackedTrack;

  if (nTracks != 0)
  {
    while (true)
    {
      if (stacks[fTurn]->GetNTrack() != 0u)
      {
        aStackedTrack = stacks[fTurn]->PopFromStack();
        energies[fTurn] -= aStackedTrack.GetTrack()->GetDynamicParticle()->GetTotalEnergy();
        --nTracks;
        break;
      }
      
      fTurn = (fTurn+1) % nTurn;
    }
  }

  return aStackedTrack;
}

enum
{
  electronCode = 11, positronCode = -11, gammaCode = 22, neutronCode = 2112
};

void G4SmartTrackStack::PushToStack( const G4StackedTrack& aStackedTrack )
{

  G4int iDest = 0;
  if (aStackedTrack.GetTrack()->GetParentID() != 0)
  {
    G4int code = aStackedTrack.GetTrack()->GetDynamicParticle()->GetPDGcode();
    if (code == electronCode)
      iDest = 2;
    else if (code == gammaCode)
      iDest = 3;
    else if (code == positronCode)
      iDest = 4;
    else if (code == neutronCode)
      iDest = 1;
  }
  else
  {
    // We have a primary track, which should go first.
    fTurn = 0; // reseting the turn
  }
  stacks[iDest]->PushToStack(aStackedTrack);
  energies[iDest] += aStackedTrack.GetTrack()->GetDynamicParticle()->GetTotalEnergy();
  ++nTracks;
  
  G4long dy1 = stacks[iDest]->GetNTrack() - stacks[iDest]->GetSafetyValue1();
  G4long dy2 = stacks[fTurn]->GetNTrack() - stacks[fTurn]->GetSafetyValue2();
  
  if (dy1 > 0 || dy1 > dy2 ||
      (iDest == 2 &&
       stacks[iDest]->GetNTrack() < 50 && energies[iDest] < energies[fTurn]))
  {
    fTurn = iDest;
  }
  
  if (nTracks > maxNTracks) maxNTracks = nTracks;
}

void G4SmartTrackStack::clear()
{
  for (G4int i = 0; i < nTurn; ++i)
  {
    stacks[i]->clear();
    energies[i] = 0.0;
    fTurn = 0;
  }
  nTracks = 0;
}

void G4SmartTrackStack::clearAndDestroy()
{
  for (G4int i = 0; i < nTurn; ++i)
  {
    stacks[i]->clearAndDestroy();
    energies[i] = 0.0;
    fTurn = 0;
  }
  nTracks = 0;
}
