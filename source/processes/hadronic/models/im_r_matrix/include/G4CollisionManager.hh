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
#ifndef G4CollisionManager_hh
#define G4CollisionManager_hh

#include "globals.hh"
#include "G4ListOfCollisions.hh"
#include "G4KineticTrackVector.hh"

class G4KineticTrack;
class G4CollisionInitialState;

class G4CollisionManager
{
public:
  G4CollisionManager();
  ~G4CollisionManager();

  G4int Entries();
  void AddCollision(G4double time, G4KineticTrack * proj,
		    G4KineticTrack * target = NULL);
  void AddCollision(G4CollisionInitialState * collision)
  {
    theCollisionList->push_back(collision);
  }
  void RemoveCollision(G4CollisionInitialState * collision);
  void RemoveTracksCollisions(G4KineticTrackVector * ktv);
  void ClearAndDestroy();
  G4CollisionInitialState * GetNextCollision();
  void Print();

private:
  G4ListOfCollisions * theCollisionList;   //was sorted (by time) vector...
};



inline G4int G4CollisionManager::Entries()
{
  return theCollisionList->size();
}


#endif














