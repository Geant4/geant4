
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
