#ifndef G4VAnnihilator_h
#define G4VAnnihilator_h 1

class G4KineticTrackVector;
class G4KineticTrack;


class G4VAnnihilator 
{
public:
   G4VAnnihilator();
  virtual ~G4VAnnihilator();
public:
    virtual G4KineticTrackVector* Scatter(const G4KineticTrack &aProjectile, const G4KineticTrack &aTarget) = 0;
};

#endif 


