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
// $Id: G4Scatterer.hh,v 1.5 2006-06-29 20:35:58 gunter Exp $ //
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
					const G4KineticTrack& trk2) const;
  
  G4double GetCrossSection(const G4KineticTrack& trk1, 	
			   const G4KineticTrack& trk2) const;
			   
  virtual G4KineticTrackVector* Scatter(const G4KineticTrack& trk1, 	
					   const G4KineticTrack& trk2) const;

  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & someCandidates,
		       G4double aCurrentTime);

  virtual G4KineticTrackVector * 
         GetFinalState(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & theTargets);

private:

  const G4VCollision* FindCollision(const G4KineticTrack& trk1,
			      const G4KineticTrack& trk2) const;
  
  struct Register
  {
    template<class T> void operator()(T*, G4CollisionVector * aC)
    {
      G4VCollision* aT = new T; 
      aC->push_back(aT);
    }
  };

  static G4CollisionVector collisions;
  std::vector<G4CollisionInitialState *> theCollisions;

}; 
#endif 


