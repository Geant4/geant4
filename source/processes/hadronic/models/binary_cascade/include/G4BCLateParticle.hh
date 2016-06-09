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
#ifndef G4BCLateParticle_h
#define G4BCLateParticle_h 1

#include <vector>

#include "G4BCAction.hh"
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

class G4BCLateParticle : public G4BCAction
{
  public:
  
  G4BCLateParticle(){}
  virtual ~G4BCLateParticle(){}
  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & ,
		       G4double theCurrentTime)
  {
    theColl.clear();
    G4double ctime=std::max(0.,aProjectile->GetFormationTime()) + theCurrentTime;

    G4KineticTrackVector noTarget;
    G4CollisionInitialState * aLate = 
            new G4CollisionInitialState(ctime,
	    				 aProjectile, noTarget, this);
    theColl.push_back(aLate);    
    return theColl;
  }

  virtual G4KineticTrackVector * GetFinalState(G4KineticTrack * aProjectile, 
	                           std::vector<G4KineticTrack *> & )
  {
    G4KineticTrackVector * result = new G4KineticTrackVector;
    G4KineticTrack * lateParticle=new G4KineticTrack(*aProjectile);
    result->push_back(lateParticle);
    return result;
  }
  
  private:
  std::vector<G4CollisionInitialState *> theColl;
};

#endif
