#ifndef G4MesonAbsorption_hh
#define G4MesonAbsorption_hh

#include <vector>
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include "G4BCAction.hh"

class G4MesonAbsorption : public G4BCAction
{
  public:
  G4MesonAbsorption()
  {
    std::cout << "Please enter the microscopic pion absorption cross-section"
              << std::endl;
    std::cin >> theCross;
  }
  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & someCandidates,
		       G4double aCurrentTime);

  virtual G4KineticTrackVector * 
         GetFinalState(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & theTargets);

  G4CollisionInitialState * GetCollision(G4KineticTrack * projectile, 
                                         std::vector<G4KineticTrack *> targets);

  private:
  G4double GetTimeToAbsorption(const G4KineticTrack& trk1, const G4KineticTrack& trk2);
  
  void FindAndFillCluster(G4KineticTrackVector & result, 
                         G4KineticTrack * aProjectile,
			 std::vector<G4KineticTrack *> & someCandidates);
  
  G4double AbsorptionCrossSection(const G4KineticTrack & trk1, const G4KineticTrack & trk2);  
  private:
  std::vector<G4CollisionInitialState *> theCollisions;
  G4double theCross;
};

#endif
