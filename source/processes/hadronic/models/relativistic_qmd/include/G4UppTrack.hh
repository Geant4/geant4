
#ifndef G4UPPTRACK_H
#define G4UPPTRACK_H


#include "G4KineticTrack.hh"
#include "G4Nucleon.hh"
#include "G4VKineticNucleon.hh"
#include "G4LorentzVector.hh"

//class G4UppTrackVector;

class G4UppTrack : public G4KineticTrack
{
public:

  G4UppTrack() : 
    NumberOfCollisions(0), changed(true), nonInteractionGroup(0) {};
  G4UppTrack(const G4KineticTrack& aTrack) :
    G4KineticTrack(aTrack), NumberOfCollisions(0), 
    changed(true), nonInteractionGroup(0) {};
  G4UppTrack(const G4Nucleon& aNucleon, const G4int g=0) :
    NumberOfCollisions(0), changed(true), nonInteractionGroup(g)
    {
      SetDefinition(aNucleon.GetDefinition());
      SetPosition(aNucleon.GetPosition());
      Set4Momentum(aNucleon.Get4Momentum());
    };
  
  G4bool hasChanged() const { return changed; }
  void setChanged(G4bool c) { changed=c; } 
  G4int getNonInteractionGroup() const { return nonInteractionGroup; }
  void setNonInteractionGroup(G4int g) { nonInteractionGroup=g; }
  G4int getNumberOfCollisions() const { return NumberOfCollisions; }
  void setNumberOfCollisions(G4int n) { NumberOfCollisions=n; }
  G4bool isLastPartner(const G4UppTrack* PartPtr) const;
  void addLastPartner(const G4UppTrack* PartPtr) { lastPartners.push_back(PartPtr); }
  void clearLastPartner() { lastPartners.clear(); }
  G4LorentzVector Get4Position() const
     { return G4LorentzVector(GetPosition(),LocalTime); }
  void Set4Position(G4LorentzVector& a4Pos) 
     { SetPosition(a4Pos.vect()); LocalTime = a4Pos.t(); }
  void dump() const;

private:

  G4double LocalTime;
  G4int NumberOfCollisions;
  G4bool changed;
  vector<const G4UppTrack*> lastPartners;
  G4int nonInteractionGroup;

};


#endif // G4UPPTRACK_H


