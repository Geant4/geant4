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
#ifndef G4BCLateParticle_h
#define G4BCLateParticle_h 1
#include "G4BCAction.hh"
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include <vector>

//#define debug_BCLateParticle 1

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
    if(aProjectile->GetFormationTime() > 0 )
    {
      G4KineticTrackVector noTarget;
      G4CollisionInitialState * aLate = 
            new G4CollisionInitialState(aProjectile->GetFormationTime(),
	    				 aProjectile, noTarget, this);
      theColl.push_back(aLate);
#ifdef debug_BCLateParticle
  G4cout << "Particle starting late " << aProjectile << " "  
         <<  aProjectile->GetDefinition()->GetParticleName() << " "
	 <<  1/MeV*aProjectile->Get4Momentum() << " " 
	 <<  1/fermi*aProjectile->GetPosition() << " " 
	 <<  aProjectile->GetFormationTime()
	 << G4endl;
#endif	 
    }
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
