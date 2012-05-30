//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
    G4MesonAbsorption(const G4MesonAbsorption &);
    G4MesonAbsorption & operator= (const G4MesonAbsorption &);

  private:
  G4double GetTimeToAbsorption(const G4KineticTrack& trk1, const G4KineticTrack& trk2);
  
  void FindAndFillCluster(G4KineticTrackVector & result, 
                         G4KineticTrack * aProjectile,
			 std::vector<G4KineticTrack *> & someCandidates);
  
  G4double AbsorptionCrossSection(const G4KineticTrack & trk1, const G4KineticTrack & trk2);  
  private:
  std::vector<G4CollisionInitialState *> theCollisions;
    
};

#endif
