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
//
// $Id: G4CollisionInitialState.hh 89688 2015-04-27 10:06:44Z gcosmo $
//
// $Id: G4CollisionInitialState.hh,v 1.0 1998/06/30
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 30th June 1998

// -----------------------------------------------------------------------------

#ifndef G4CollisionInitialState_hh
#define G4CollisionInitialState_hh

#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4HadTmpUtil.hh"

class G4BCAction;

class G4CollisionInitialState
{

public:
  G4CollisionInitialState();
  G4CollisionInitialState(G4double time, G4KineticTrack * aPrimary,
			  G4KineticTrack * aTarget);
  G4CollisionInitialState(G4double time, G4KineticTrack * aPrimary,
			  const G4KineticTrackVector & aTarget,
			  G4BCAction * aFSGenerator);

  ~G4CollisionInitialState() { }

  G4CollisionInitialState(const G4CollisionInitialState & right);
  G4CollisionInitialState & operator=(const G4CollisionInitialState & right);

  int operator<(const G4CollisionInitialState & right) const
    {return (theCollisionTime < right.theCollisionTime);}

  int operator==(const G4CollisionInitialState& right) const
    {return (theCollisionTime == right.theCollisionTime);}


  G4KineticTrack * GetPrimary(void)
    {return thePrimary;}
  void SetPrimary(G4KineticTrack * aPrimary)
    {thePrimary = aPrimary;}

  G4KineticTrack * GetTarget(void)
    {return theTarget;}
  void SetTarget(G4KineticTrack * aTarget)
    {theTarget = aTarget;}

  void AddTarget(G4KineticTrack * aTarget)
    {theTs.push_back(aTarget);}
  G4KineticTrackVector  & GetTargetCollection(void)
    {return theTs;}
  G4KineticTrackVector * GetFinalState();
  G4int GetTargetBaryonNumber()
  {
    G4double result=0;
    for(size_t i=0; i<theTs.size(); i++)
    {
      result += theTs[i]->GetDefinition()->GetBaryonNumber();
    }
    return G4lrint(result);
  }
  G4int GetTargetCharge()
  {
    G4double result=0;
    for(size_t i=0; i<theTs.size(); i++)
    {
      result += theTs[i]->GetDefinition()->GetPDGCharge();
    }
    return G4lrint(result);
  }

// -new interface post pion:

  G4double GetCollisionTime(void)
    {return theCollisionTime;}
  void SetCollisionTime(G4double value)
    {theCollisionTime = value;}

// for debugging only
const G4BCAction * GetGenerator()
  {
    return theFSGenerator;
  }


  void Print() const;

//  friend std::ostream& operator<<(std::ostream & out, const G4CollisionInitialState & collision){
//  out=1;
//  return out;
//  }


private:

  G4double theCollisionTime;
  G4KineticTrack * thePrimary;
  G4KineticTrack * theTarget;
  G4KineticTrackVector theTs;
  G4BCAction * theFSGenerator;
};

#endif


