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
// $Id: ParRandomState.hh,v 1.1 2002-03-05 15:22:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------
// Each of the three members below must be modified for the application.
// This default version simply sets the random seed on the slave using
// the eventID.
// --------------------------------------------------------------------

#ifndef ParRandomState_h
#define ParRandomState_h 1

#include "ParMarshaledObj.hh"

class ParMarshaledRandomState : public MarshaledObj
{
  public:
  // If 0, GetNextRandomStateForSlave not called.
  static const int RANDOM_STATE_SIZE = sizeof(long);

  // Called on master
  inline void *GetNextRandomStateForSlave()
  { static long * nextSeed = new long;
    *nextSeed = (long) (100000000L * HepRandom::getTheGenerator()->flat());
#ifdef TOPC_DEBUG
    G4cout << "GET SEED: " << *nextSeed << G4endl;
#endif
    return nextSeed;
  }

  // Called on slave at beginning of new event.
  inline void SetNextRandomStateForSlave( G4int eventID, void *randomState )
  { HepRandom::setTheSeed( *(long *)randomState );
#ifdef TOPC_DEBUG
    G4cout << "SET SEED: " << *(long *)randomState << G4endl;
#endif
  }

  public:
  // The constructors do not need to be modified.
  inline ParMarshaledRandomState( G4int eventID )
  : MarshaledObj(RANDOM_STATE_SIZE+sizeof(int)+sizeof(G4int))
  { Marshal(eventID);
    Marshal( GetNextRandomStateForSlave(), RANDOM_STATE_SIZE );
  }
  inline ParMarshaledRandomState( void *buf )
  : MarshaledObj(buf)
  // MarshaledObj() called with buf and 0 size, creates object for unmarshaling
  { }

  inline void unmarshalEventIDandSetState( G4int & eventID )
  { Unmarshal(eventID);
    // SharedObjPtr is pointer to random state.
    SetNextRandomStateForSlave( eventID, UnmarshalAsSharedObjPtr() );
  }
};

#endif
