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
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, A. Feliciello, 30th June 1998
// -----------------------------------------------------------------------------

#include "G4CollisionInitialState.hh"

// Class G4CollisionInitialState 

G4CollisionInitialState::G4CollisionInitialState() :
   theCollisionTime(DBL_MAX), thePrimary(NULL), theTarget(NULL)
{
}



G4CollisionInitialState::G4CollisionInitialState(G4double time, 
   G4KineticTrack * aPrimary, G4KineticTrack * aTarget)
{
   theCollisionTime = time;
   thePrimary       = aPrimary;
   theTarget        = aTarget;
}



G4CollisionInitialState::G4CollisionInitialState(G4CollisionInitialState & right)
{
   theCollisionTime = right.theCollisionTime;
   thePrimary       = right.thePrimary;
   theTarget        = right.theTarget;
}

const G4CollisionInitialState & G4CollisionInitialState::
      operator=(const G4CollisionInitialState& right)
{
   if (this != &right)
   {
      theCollisionTime = right.theCollisionTime;
      thePrimary       = right.thePrimary;
      theTarget        = right.theTarget;
   }

   return *this;
}






