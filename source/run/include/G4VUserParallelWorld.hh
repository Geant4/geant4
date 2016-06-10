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
// $Id: G4VUserParallelWorld.hh 69917 2013-05-17 13:28:02Z gcosmo $
//

#ifndef G4VUserParallelWorld_h
#define G4VUserParallelWorld_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;

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
    virtual void ConstructSD();

  protected:
    G4String fWorldName;
    
  protected:
    G4VPhysicalVolume* GetWorld();

  public:
    inline G4String GetName() { return fWorldName; }

  protected:
    void SetSensitiveDetector(const G4String& logVolName,
                G4VSensitiveDetector* aSD,G4bool multi=false);
    void SetSensitiveDetector(G4LogicalVolume* logVol,
                G4VSensitiveDetector* aSD);
};

#endif

