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

#include "Test2ParallelWorld.hh"
#include "Test2GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"


Test2ParallelWorld::Test2ParallelWorld(G4String worldName)
  :G4VUserParallelWorld(worldName), fbConstructed(false),
   SDC(0),PSC(0),fSensitivityType(0)
{;}

Test2ParallelWorld::~Test2ParallelWorld()
{
  if(SDC) delete SDC;
  if(PSC) delete PSC;
}

void Test2ParallelWorld::Construct() {

  if(!fbConstructed)
  { 
    fbConstructed = true;
    SetupGeometry();
    switch ( fSensitivityType ){
    case 0:
      break;
    case 1:
      SetupPSDetectors();
      break;
    case 2:
      SetupSDDetectors();
      break;
    default:
      break;
    }
  }
}

void Test2ParallelWorld::SetupGeometry()
{
  //     
  // World
  //
  G4VPhysicalVolume * ghostWorld = GetWorld();

  G4cout << "Test2ParallelWorld::SetupGeometry"<<G4endl;
  //                               
  // 3D nested phantom
  //
  // parameters
  phantomSize[0] = 1.*m;
  phantomSize[1] = 1.*m;
  phantomSize[2] = 1.*m;
  nSegment[0] = 10;
  nSegment[1] = 10;
  nSegment[2] = 10;

  GEOM = 
    new Test2GeometryConstruction(phantomSize,NULL,nSegment);
  fWorldPhys = GEOM->ConstructGeometry(ghostWorld);

}


#include "G4SDManager.hh"
#include "Test2SDConstruction.hh"

void Test2ParallelWorld::SetupSDDetectors() {
  G4LogicalVolume* phantomLogical = GEOM->GetSensitiveLogical();

  SDC = new Test2SDConstruction("ParallelWorldSD",nSegment);
  SDC->SetupSensitivity(phantomLogical);
}

#include "Test2PSConstruction.hh"

void Test2ParallelWorld::SetupPSDetectors() {
  G4LogicalVolume* phantomLogical = GEOM->GetSensitiveLogical();

  PSC = new Test2PSConstruction("ParallelWorldPS",nSegment);
  PSC->SetupSensitivity(phantomLogical);
}

