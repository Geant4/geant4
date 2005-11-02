#ifndef G4TRAJECTORYDRAWBYPARTICLEID
#define G4TRAJECTORYDRAWBYPARTICLEID

#include "G4VTrajectoryDrawer.hh"
#include "G4Colour.hh"
#include "G4String.hh"
#include <map>

using std::map;

class G4TrajectoryDrawByParticleID : public G4VTrajectoryDrawer {

public:
 
  G4TrajectoryDrawByParticleID(const G4String& name = "NULL");
  
  virtual ~G4TrajectoryDrawByParticleID();

  virtual void Draw(const G4VTrajectory&, G4int);

  virtual void Print() const;

  void SetDefault(const G4Colour&);
  void Set(const G4String& particle, const G4String& colour);

private:

  map<G4String, G4Colour> fMap;
  G4Colour fDefault;
};

#endif
