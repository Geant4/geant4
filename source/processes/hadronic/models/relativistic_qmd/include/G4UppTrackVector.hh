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

#ifndef G4UPPTRACKVECTOR_H
#define G4UPPTRACKVECTOR_H


#include "g4std/vector"
#include "G4UppTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Nucleon.hh"


class G4UppTrackChange;


class G4UppTrackVector : public G4std::vector<G4UppTrack*>
{
public:

  void add(const G4KineticTrackVector& aTrackVector);
  void add(const G4RWTPtrOrderedVector<G4Nucleon>& aNucleonVector, 
	   const G4int nonInteractionGroup=0); 

  G4int getIndex(const G4UppTrack* trackPtr) const;

  G4UppTrackChange update(const G4UppTrackChange& aTrackChange);

  void setGlobalTime(const G4double newGlobalTime) 
    { globalTime = newGlobalTime; }
  G4double getGlobalTime() const 
    { return globalTime; }

  G4bool isPauliBlocked(const G4UppTrackVector& particlesToCheck) const;
  G4bool isPauliBlocked(const G4UppTrack* particleToCheck) const;

  void resetChangedFlags();

  void dump() const;

private:

  G4double globalTime;

};


#endif // G4UPPTRACKVECTOR_H
