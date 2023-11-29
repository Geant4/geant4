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
//
// modified by I. Hrivnacova, 13.10.99 

#include <iomanip>
#include <fstream>
#include "globals.hh"
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"
#include "G3RotTable.hh"
#include "G3VolTable.hh"
#include "G3PartTable.hh"
#include "G3DetTable.hh"
#include "G3toG4BuildTree.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

extern std::ofstream ofile;

void G3CLRead(G4String &, char *);
void checkVol(G4LogicalVolume*, G4int);
void checkVol();


G4LogicalVolume* G4BuildGeom(G4String& inFile){

  G4int irot=0;
  G4gsrotm(0, 90, 0, 90, 90, 0, 0);

  G4cout << "Instantiated unit rotation matrix irot=" << irot << G4endl;
 
  // Read the call List and interpret to Generate Geant4 geometry

  G4cout << "Reading the call List file " << inFile << "..." << G4endl;

  G3CLRead(inFile, 0);

  G3Part.PrintAll();

  G3Det.PrintAll();

  G3Vol.PrintAll();

  G4cout << "Call List file read completed. Build geometry" << G4endl;

  // Build the geometry

  G3VolTableEntry* topVTE = G3Vol.GetFirstVTE();
  G4cout << "G3toG4 top level volume is " << topVTE->GetName() << G4endl;

  // modified
  G3toG4BuildTree(topVTE, 0);

  // Retrieve the top-level G3toG4 logical mother volume pointer

  G4LogicalVolume* topLV = topVTE->GetLV();

  // position the top logical volume
  // (in Geant3 the top volume is not positioned)
  // 
  new G4PVPlacement(0, G4ThreeVector(), topLV->GetName(), topLV, 0, false, 0);

  // mark as invisible

  topLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    
  G4cout << "Top-level G3toG4 logical volume " << topLV->GetName() << " "
         << *(topLV->GetVisAttributes()) << G4endl;
        
        // check the geometry here

  #ifdef G3G4DEBUG
    G4cout << "scan through G4LogicalVolumeStore:" << G4endl;
    checkVol();
  #endif

  return topLV;
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
  ++level;
  
  std::size_t ndau = _lvol -> GetNoDaughters();
  
  G4cout << "G44LogicalVolume " << _lvol->GetName() << " at level " << level
	 << " contains " << ndau << " daughters." << G4endl;
  for (std::size_t idau=0; idau<ndau; ++idau)
  {
    _pdvol = _lvol-> GetDaughter(idau);
    _ldvol = _pdvol -> GetLogicalVolume();
    G4cout << "G4VPhysical volume " << std::setw(5) << _pdvol -> GetName() 
	 << " (G4LogicalVolume " << std::setw(5) << _ldvol->GetName() << ")" 
	 << G4endl;
    checkVol(_ldvol, level);
  }
  return;
}











