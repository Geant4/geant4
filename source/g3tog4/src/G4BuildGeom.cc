// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BuildGeom.cc,v 1.8 1999-07-21 08:40:13 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $


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
#include "G3toG4BuildTree.hh"

extern ofstream ofile;

void G3CLRead(G4String &, char *);
void checkVol(G4LogicalVolume*, G4int);
void checkVol();


G4LogicalVolume* G4BuildGeom(G4String& inFile){

  G4int irot=0;
  G4gsrotm(0, 90, 0, 90, 90, 0, 0);

  G4cout << "Instantiated unit rotation matrix irot=" << irot << endl;
 
  // Read the call List and interpret to Generate Geant4 geometry

  G4cout << "Reading the call List file " << inFile << "..." << endl;

  G3CLRead(inFile, NULL);

  G3Vol.VTEStat();

  G4cout << "Call List file read completed. Build geometry" << endl;

  // Print the materials

  // G3Mat.print();

  // Build the geometry

  G4cout << "G3toG4 top level volume is " << G3Vol.GetFirstVTE()->GetName()
	 << endl;

  G3toG4BuildTree(G3Vol.GetFirstVTE());

  // Retrieve the top-level G3toG4 logical mother volume pointer

  G4LogicalVolume* lG3toG4 = G3Vol.GetFirstVTE()->GetLV();

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
  return lG3toG4;
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





