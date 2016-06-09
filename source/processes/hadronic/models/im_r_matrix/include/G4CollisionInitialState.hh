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
//
// $Id: G4CollisionInitialState.hh,v 1.4.4.1 2004/03/24 13:18:28 hpw Exp $
// GEANT4 tag $Name: geant4-06-01 $
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
// +new interface post pion:
  G4CollisionInitialState(G4double time, G4KineticTrack * aPrimary,
			  const G4KineticTrackVector & aTarget,
			  G4BCAction * aFSGenerator);
// -new interface post pion:

  G4CollisionInitialState(G4CollisionInitialState & right);

  ~G4CollisionInitialState() { }

  const G4CollisionInitialState & operator=(const G4CollisionInitialState & right);
      
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

// +new interface post pion:
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
  G4BCAction * GetGenerator()
  {
    return theFSGenerator;
  }
private:

  G4double theCollisionTime;
  G4KineticTrack * thePrimary;
  G4KineticTrack * theTarget;
  G4KineticTrackVector theTs;
  G4BCAction * theFSGenerator;
};

#endif


