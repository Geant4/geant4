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
// Author: Makoto Asai (SLAC)
// --------------------------------------------------------------------

#include "G4TrackStack.hh"
#include "G4SmartTrackStack.hh"
#include "G4VTrajectory.hh"
#include "G4Track.hh"

G4TrackStack::~G4TrackStack()
{
  clearAndDestroy();
}

void G4TrackStack::clearAndDestroy()
{
  for( auto i = begin(); i != end(); ++i )
  {
    delete (*i).GetTrack();
    delete (*i).GetTrajectory();
  }
  clear();
}

void G4TrackStack::TransferTo(G4TrackStack* aStack)
{
  for(auto i = begin(); i != end(); ++i)
  {
    aStack->push_back(*i);
  }
  clear();
}

void G4TrackStack::TransferTo(G4SmartTrackStack * aStack)
{
  while (size())
  {
    aStack->PushToStack(PopFromStack());
  }
}

G4double G4TrackStack::getTotalEnergy(void) const
{
  G4double totalEnergy = 0.0;
  for ( auto i = cbegin(); i != cend(); ++i )
  {
    totalEnergy += (*i).GetTrack()->GetDynamicParticle()->GetTotalEnergy();
  }
  return totalEnergy;
}
