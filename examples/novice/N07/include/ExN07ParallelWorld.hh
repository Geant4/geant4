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
// $Id$
//
// 

#ifndef ExN07ParallelWorld_h
#define ExN07ParallelWorld_h 1

#include "G4VUserParallelWorld.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class ExN07ParallelWorld : public G4VUserParallelWorld
{
  public:
    ExN07ParallelWorld(G4String worldName);
    virtual ~ExN07ParallelWorld();

  public:
    virtual void Construct();

  private:
    void SetupGeometry();
    void SetupDetectors();
     
  public:
    void SetSerialGeometry(G4bool ser);
    inline G4bool IsSerial() const
    { return serial; }

  private:
    G4LogicalVolume*   calorLogical[3];
    G4LogicalVolume*   layerLogical[3];
    G4VPhysicalVolume* calorPhysical[3];
    G4VPhysicalVolume* layerPhysical[3];
    G4String           calName[3];
    G4bool             constructed;
    G4bool             serial;
    G4double           totalThickness;
    G4int              numberOfLayers;
};


#endif

