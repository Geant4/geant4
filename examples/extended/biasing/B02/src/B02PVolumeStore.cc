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
// $Id: B02PVolumeStore.cc,v 1.1 2002-11-08 14:52:17 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// B02PVolumeStore.cc
//
// ----------------------------------------------------------------------

#include "B02PVolumeStore.hh"
#include "g4std/strstream"


#include "G4VPhysicalVolume.hh"

B02PVolumeStore::B02PVolumeStore(){}
B02PVolumeStore::~B02PVolumeStore(){}
  
void B02PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  B02SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  if (it != fSetGeometryCell.end()) {
    G4cout << "B02PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

const G4VPhysicalVolume *B02PVolumeStore::
GetPVolume(const G4String &name) const {
  const G4VPhysicalVolume *pvol = 0;
  for (B02SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    if (vol.GetName() == name) {
      pvol =  &vol;
    } 
  }
  if (!pvol) {
    G4cout << "B02PVolumeStore::GetPVolume: no physical volume named: " 
	   << name << ", found" << G4endl;
  }
  return pvol;
}

G4String B02PVolumeStore::GetPNames() const {
  G4String NameString;
  for (B02SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    char st[200];
    G4std::ostrstream os(st,200);
    os << vol.GetName() << "_" << it->GetReplicaNumber() 
       << "\n" << '\0';
    G4String cellname(st);
    
    //    G4String cellname(vol.GetName());
    //    cellname += G4String("_");
    //    cellname += G4std::str(it->GetReplicaNumber());

    NameString += cellname;
  }
  return NameString;
}
