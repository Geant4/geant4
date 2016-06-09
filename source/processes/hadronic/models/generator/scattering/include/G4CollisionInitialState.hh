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
// $Id: G4CollisionInitialState.hh,v 1.2 2002/12/12 19:17:37 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
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

class G4CollisionInitialState 
{

public:
  G4CollisionInitialState();
  G4CollisionInitialState(G4double time, G4KineticTrack * aPrimary,
			  G4KineticTrack * aTarget);

  G4CollisionInitialState(G4CollisionInitialState & right);

  ~G4CollisionInitialState() { }

  const G4CollisionInitialState & operator=(const G4CollisionInitialState & right);
      
  int operator<(const G4CollisionInitialState & right) const
    {return (theCollisionTime < right.theCollisionTime);}
      
  int operator==(const G4CollisionInitialState& right) const
    {return (theCollisionTime == right.theCollisionTime);}

  G4KineticTrack * GetPrimary(void)            {return thePrimary;}
  void SetPrimary(G4KineticTrack * aPrimary)   {thePrimary = aPrimary;}
 
  G4KineticTrack * GetTarget(void)             {return theTarget;}
  void SetTarget(G4KineticTrack * aTarget)     {theTarget = aTarget;}

  G4double GetCollisionTime(void)             {return theCollisionTime;}
  void SetCollisionTime(G4double value)       {theCollisionTime = value;}

private:

  G4double theCollisionTime;
  G4KineticTrack * thePrimary;
  G4KineticTrack * theTarget;
};

#endif


