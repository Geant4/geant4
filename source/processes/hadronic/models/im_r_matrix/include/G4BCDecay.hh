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
#ifndef G4BCDecay_h
#define G4BCDecay_h 1
#include "G4BCAction.hh"
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include <vector>

class G4BCDecay : public G4BCAction
{
  public:
  
  G4BCDecay(){}
  virtual ~G4BCDecay(){}
  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & ,
		       G4double theCurrentTime)
  {
    theColl.clear();
    if(aProjectile->GetDefinition()->IsShortLived())
    {
      G4double aTime = theCurrentTime+aProjectile->SampleResidualLifetime();
      G4KineticTrackVector noTarget;
      G4CollisionInitialState * aDecay = 
            new G4CollisionInitialState(aTime, aProjectile, noTarget, this);
      theColl.push_back(aDecay);
    }
    return theColl;
  }

  virtual G4KineticTrackVector * GetFinalState(G4KineticTrack * aProjectile, 
	                           std::vector<G4KineticTrack *> & )
  {
    return aProjectile->Decay();
  }
  
  private:
  
  std::vector<G4CollisionInitialState *> theColl;
};

#endif
