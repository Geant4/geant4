// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorConstruction.hh,v 1.3 2000-05-02 13:50:14 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "g4std/vector"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

  private:
     G4double expHall_x;
     G4double expHall_y;
     G4double expHall_z;

     G4double calBox_x;
     G4double calBox_y;
     G4double calBox_z;
     G4double rotAngle;
     G4double calPos;
     G4String calMaterialName;

     G4double trackerRadius;
     G4double trackerHight;
     G4double trackerPos;

     G4std::vector<G4Material*> materialPointerStore;

  public:
     inline void SetCalMaterial(G4String name)
     { calMaterialName = name; };
};

#endif

