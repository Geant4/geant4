
#ifndef G4UPPTRACKCHANGE_H
#define G4UPPTRACKCHANGE_H


#include "g4std/vector"
#include "G4UppTrackVector.hh"


struct G4UppTrackChange
{
public:

  G4std::vector<G4UppTrackVector::iterator> oldParticles;
  G4UppTrackVector newParticles;

};


#endif // G4UPPTRACKCHANGE_H

