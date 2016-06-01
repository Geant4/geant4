#ifndef G4VAnnihilator_h
#define G4VAnnihilator_h 1

#include "G4KineticTrackVector.hh"


class G4VAnnihilator 
{
public:
   G4VAnnihilator();
  virtual ~G4VAnnihilator();
public:
    virtual G4KineticTrackVector* Scatter(const G4KineticTrack &aProjectile, const G4KineticTrack &aTarget) = 0;
};

#endif 


