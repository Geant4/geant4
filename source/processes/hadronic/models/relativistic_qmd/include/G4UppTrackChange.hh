
#ifndef G4UPPTRACKCHANGE_H
#define G4UPPTRACKCHANGE_H


#include "g4std/vector"
#include "G4UppTrackVector.hh"


struct G4UppTrackChange
{
public:

  G4UppTrackVector oldParticles;
  G4UppTrackVector newParticles;

  void dump() const;

};


#endif // G4UPPTRACKCHANGE_H

