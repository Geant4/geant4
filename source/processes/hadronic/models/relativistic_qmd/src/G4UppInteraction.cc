//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//


#include "G4UppInteraction.hh"
#include "G4KineticTrack.hh"


G4UppInteraction::G4UppInteraction(const G4double Time, const G4int i, const G4int j)
{
  IncomingParticleIndex.push_back(i);
  IncomingParticleIndex.push_back(j);
  InteractionTime = Time;
}


G4UppTrackChange G4UppInteraction::Perform()
{
  G4UppTrackChange aChange;
  return aChange;
}
