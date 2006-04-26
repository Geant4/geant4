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
// $Id: G4VUserParallelWorld.hh,v 1.1 2006-04-26 15:24:24 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VUserParallelWorld_h
#define G4VUserParallelWorld_h 1

class G4VPhysicalVolume;

#include "globals.hh"

// class description:
//
//  This is the abstract base class of the user's parallel world.
// The user MUST NOT create the worls volume, but should get it from GetWorld()
// protected method, and fill inside as the Construct() method.
// The constructor must take a unique name of the parallel world, which is used
// for the name of the world physical volume of this parallel world.
//

class G4VUserParallelWorld
{
  public:
    G4VUserParallelWorld(G4String worldName);
    virtual ~G4VUserParallelWorld();

  public:
    virtual void Construct() = 0;

  protected:
    G4String fWorldName;
    
  protected:
    G4VPhysicalVolume* GetWorld();

  public:
    inline G4String GetName() { return fWorldName; }

};

#endif

