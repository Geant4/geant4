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
#ifndef G4MesonAbsorption_hh
#define G4MesonAbsorption_hh

#include <vector>
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include "G4BCAction.hh"

class G4MesonAbsorption : public G4BCAction
{
  public:
  G4MesonAbsorption(){}
  virtual ~G4MesonAbsorption(){}

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
