
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

  void add(const G4KineticTrackVector& v);
  void add(const G4RWTPtrOrderedVector<G4Nucleon>& v, const G4int g=0); 
  G4int getIndex(const G4UppTrack* tPtr) const;
  G4UppTrackChange Update(const G4UppTrackChange& aChange);
  G4double getGlobalTime() const 
    { return globalTime; }
  void setGlobalTime(const G4double t) { globalTime=t; }
  G4bool isPauliBlocked(const G4UppTrackVector& t) const;
  G4bool isPauliBlocked(const G4UppTrack* t) const;
  void resetChanged();
  void dump() const;

private:

  G4double globalTime;

};


#endif // G4UPPTRACKVECTOR_H
