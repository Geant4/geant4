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
		    G4KineticTrack * target = nullptr);
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
  G4CollisionManager(const G4CollisionManager &);
  G4CollisionManager & operator= (const G4CollisionManager &);

private:
  G4ListOfCollisions * theCollisionList;   //was sorted (by time) vector...
};

inline G4int G4CollisionManager::Entries()
{
  return (G4int)theCollisionList->size();
}


#endif
