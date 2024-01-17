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
// G4SubEvent class implementation
//
// Author: Makoto Asai (JLAB) - 23/Aug/23
// --------------------------------------------------------------------

#include "G4SubEvent.hh"
#include "G4VTrajectory.hh"
#include "G4Track.hh"

G4SubEvent::~G4SubEvent()
{
  clearAndDestroy();
}

void G4SubEvent::clearAndDestroy()
{
  for(auto & i : *this)
  {
    delete i.GetTrack();
    delete i.GetTrajectory();
  }
  clear();
}

G4StackedTrack G4SubEvent::PopFromStack()
{
  G4Exception("G4SubEvent::PopFromStack","EventStack0001",
       FatalException,"This method must not be invoked.");
  return G4StackedTrack(nullptr);
}

G4double G4SubEvent::getTotalEnergy() const
{
  G4double totalEnergy = 0.0;
  for (const auto & i : *this)
  {
    totalEnergy += i.GetTrack()->GetDynamicParticle()->GetTotalEnergy();
  }
  return totalEnergy;
}
