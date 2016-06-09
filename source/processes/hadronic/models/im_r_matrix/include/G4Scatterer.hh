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
// $Id: G4Scatterer.hh,v 1.2 2003/10/10 14:08:22 hpw Exp $ //
//
//
// removing a auther spec that was part of a template.

#ifndef G4SCATTERER_HH
#define G4SCATTERER_HH

#include "globals.hh"
#include <vector>
#include "G4VScatterer.hh"
#include "G4VCollision.hh"
#include "G4KineticTrackVector.hh"
#include "G4CollisionVector.hh"
#include "G4BCAction.hh"

class G4KineticTrack;

class G4Scatterer : public G4VScatterer, public G4BCAction
{
public:

  G4Scatterer();
  
  virtual ~G4Scatterer();
  
  virtual G4double GetTimeToInteraction(const G4KineticTrack& trk1, 
					const G4KineticTrack& trk2);
  
  G4double GetCrossSection(const G4KineticTrack& trk1, 	
			   const G4KineticTrack& trk2);
			   
  virtual G4KineticTrackVector* Scatter(const G4KineticTrack& trk1, 	
					   const G4KineticTrack& trk2);

  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & someCandidates,
		       G4double aCurrentTime);

  virtual G4KineticTrackVector * 
         GetFinalState(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & theTargets);

private:

  G4VCollision* FindCollision(const G4KineticTrack& trk1, 
			      const G4KineticTrack& trk2);
  
  G4CollisionVector collisions;
  std::vector<G4CollisionInitialState *> theCollisions;

}; 
#endif 


