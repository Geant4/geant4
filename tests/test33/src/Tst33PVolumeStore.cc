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
// $Id: Tst33PVolumeStore.cc,v 1.4 2002-11-20 09:38:26 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33PVolumeStore.cc
//
// ----------------------------------------------------------------------

#include "Tst33PVolumeStore.hh"
#include "g4std/strstream"


#include "G4VPhysicalVolume.hh"

Tst33PVolumeStore::Tst33PVolumeStore(){}
Tst33PVolumeStore::~Tst33PVolumeStore(){}
  
void Tst33PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  Tst33SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  if (it != fSetGeometryCell.end()) {
    G4cout << "Tst33PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

G4GeometryCell Tst33PVolumeStore::
GetGeometryCell(G4int i) const {
  Tst33SetGeometryCell::const_iterator  pCell(fSetGeometryCell.end());
  G4String cellName(GetCellName(i));

  for (Tst33SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    if (vol.GetName() == cellName) {
      pCell = it;
    } 
  }
  if (pCell == fSetGeometryCell.end()) {
    G4cout << "Tst33PVolumeStore::GetGeometryCell: no G4GeometryCell named: " 
	   << cellName << ", found" << G4endl;
    G4Exception("G4GeometryCell not found");
  }
  return *pCell;
}

G4String Tst33PVolumeStore::GetCellName(G4int i) const {
  char st[200];
  G4std::ostrstream os(st,200);
  os << "cell_";
  if (i<10) {
    os << "0";
  }
  os << i 
     << '\0';
  G4String name(st);
  return name;
}
