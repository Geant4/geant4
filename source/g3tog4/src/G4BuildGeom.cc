// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BuildGeom.cc,v 1.2 1999-05-01 21:42:01 lockman Exp $
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

#include "G4ios.hh"
#include <fstream.h>
#include "G4GeometryManager.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
extern ofstream ofile;

void G3CLRead(G4String &, char *);
void checkVol(G4LogicalVolume*, G4int);


G4LogicalVolume* G4BuildGeom(G4String& inFile)
{
        // Read the call List and interpret to Generate Geant4 geometry

    G4cout << "Reading the call List file " << inFile << "..." << endl;

    G3CLRead(inFile, NULL);

    G4cout << "Call List file read completed." << endl;

        // Retrieve the top-level G3toG4 logical mother volume pointer

    G4LogicalVolume* lG3toG4 = G3Vol.GetLV();

        // make the top-level volume invisible
    
    lG3toG4 -> SetVisAttributes(G4VisAttributes::Invisible);
    G4cout << "Top-level G3toG4 logical volume " << lG3toG4->GetName() << " "
         << *(lG3toG4 -> GetVisAttributes()) << endl;
        
        // check the geometry here

    G4int debug=1;
        
    if (debug){
        G4int level=0;
        checkVol(lG3toG4, level);
    }
        
    return lG3toG4;
}

void checkVol(G4LogicalVolume* _lvol, G4int level)
{
    G4LogicalVolume* _ldvol;
    G4VPhysicalVolume* _pdvol;
    level++;
    
    G4int ndau = _lvol -> GetNoDaughters();
    
    for (int idau=0; idau<ndau; idau++){
        _pdvol = _lvol-> GetDaughter(idau);
        
        G4cout << "logical volume " << _lvol->GetName() << " at level " << level
             << " contains Physical Daughter volume "
             << _pdvol -> GetName() << endl;
	_ldvol = _pdvol -> GetLogicalVolume();
        checkVol(_ldvol, level);
    }
    return;
}





