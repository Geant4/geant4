// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BuildGeom.cc,v 1.5 1999-05-22 06:31:34 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//    cltog4
//
//    Convert Geant3 call List to Geant4 init file
//
//    T. Wenaus  LLNL  6/95
//
// 1-27-97 Akbar Mokhtarani
// change the name to BuildGeom for use with visualization
//code. This routine is called by test18.cc in prototype/visualization/test
// directory. It reads the call List 'g3calls.dat' produced by rztog4
// code from geant 3 geometry and builds the g4 geometry
// It returns a pointer to the logical volume of the mother of all worlds.

#include <iomanip.h>
#include <fstream.h>
#include "G4ios.hh"
#include "G4GeometryManager.hh"
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"
#include "G3RotTable.hh"
#include "G3VolTable.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "VolTableEntry.hh"

extern ofstream ofile;

void G3CLRead(G4String &, char *);
void checkVol(G4LogicalVolume*, G4int);
void checkVol();


G4LogicalVolume* G4BuildGeom(G4String& inFile){

  // Read the call List and interpret to Generate Geant4 geometry

  G4cout << "Reading the call List file " << inFile << "..." << endl;

  G3CLRead(inFile, NULL);

  G4cout << "Call List file read completed. Build geometry" << endl;

  // Build the geometry

  G4cout << "G3toG4 top level volume is " << G3Vol.GetFirstVTE()->GetName()
	 << endl;

  /*
        // Retrieve the top-level G3toG4 logical mother volume pointer

  G4LogicalVolume* lG3toG4 = G3Vol.GetG3toG4Mother();

  // mark as invisible

  lG3toG4->SetVisAttributes(G4VisAttributes::Invisible);
    
  G4cout << "Top-level G3toG4 logical volume " << lG3toG4->GetName() << " "
         << *(lG3toG4 -> GetVisAttributes()) << endl;
        
        // check the geometry here

  G4int debug=0;
        
  if (debug){
    G4cout << "scan through G4LogicalVolumeStore:" << endl;
    checkVol();
  }
  return lG3toG4; */
  return 0;
}

void checkVol()
{
  G4LogicalVolumeStore* theStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* ll = (*theStore)[0];
  G4int level=0;
  checkVol(ll, level);
}

void checkVol(G4LogicalVolume* _lvol, G4int level)
{
  G4LogicalVolume* _ldvol;
  G4VPhysicalVolume* _pdvol;
  level++;
  
  G4int ndau = _lvol -> GetNoDaughters();
  
  G4cout << "G44LogicalVolume " << _lvol->GetName() << " at level " << level
	 << " contains " << ndau << " daughters." << endl;
  for (int idau=0; idau<ndau; idau++){
    _pdvol = _lvol-> GetDaughter(idau);
    _ldvol = _pdvol -> GetLogicalVolume();
    G4cout << "G4VPhysical volume " << setw(5) << _pdvol -> GetName() 
	 << " (G4LogicalVolume " << setw(5) << _ldvol->GetName() << ")" 
	 << endl;
    checkVol(_ldvol, level);
  }
  return;
}





