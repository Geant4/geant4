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

#ifndef G4UPPINTERACTION_H
#define G4UPPINTERACTION_H


#include "G4UppTrackVector.hh"
#include "G4UppTrackChange.hh"
#include "g4std/vector"


class G4UppInteraction
{
public:

  G4UppInteraction() : InteractionTime(-1) {}
  G4UppInteraction(const G4double Time, const G4int i, const G4int j);

  G4UppTrackChange Perform();

  G4double getInteractionTime() 
    { return InteractionTime; }
  G4UppTrackVector* getOutgoingParticles() 
    { return &OutgoingParticles; }

private:

  G4std::vector<G4int> IncomingParticleIndex;
  G4UppTrackVector OutgoingParticles;
  G4double InteractionTime;

};


#endif // G4UPPINTERACTION_H



