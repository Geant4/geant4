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
// $Id: Tst33PVolumeStore.cc,v 1.9 2006-06-29 22:01:09 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33PVolumeStore.cc
//
// ----------------------------------------------------------------------

#include "Tst33PVolumeStore.hh"
#include <sstream>


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
GetGeometryCell(G4int i, const G4String &nameExt) const {
  Tst33SetGeometryCell::const_iterator  pCell(fSetGeometryCell.end());
  G4String cellName(GetCellName(i));
  cellName += nameExt;

  for (Tst33SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();

    if (vol.GetName() == cellName) {
      pCell = it;
      break;
    } 
  }
  if (pCell == fSetGeometryCell.end()) {
    G4cout << "Tst33PVolumeStore::GetGeometryCell: no G4GeometryCell named: " 
	   << cellName << ", found" << G4endl;
    G4Exception("Tst33PVolumeStore::GetGeometryCell()",
        "TST33-07", FatalException, "G4GeometryCell not found");
  }
  return *pCell;
}

G4String Tst33PVolumeStore::GetCellName(G4int i) const {
  std::ostringstream os;
  os << "cell_";
  if (i<10) {
    os << "0";
  }
  os << i ;
  G4String name = os.str();
  return name;
}
