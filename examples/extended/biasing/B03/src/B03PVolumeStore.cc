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
// $Id: B03PVolumeStore.cc,v 1.1 2002-11-08 17:35:18 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// B03PVolumeStore.cc
//
// ----------------------------------------------------------------------

#include "B03PVolumeStore.hh"
#include "g4std/strstream"


#include "G4VPhysicalVolume.hh"

B03PVolumeStore::B03PVolumeStore(){}
B03PVolumeStore::~B03PVolumeStore(){}
  
void B03PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  B03SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  if (it != fSetGeometryCell.end()) {
    G4cout << "B03PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

const G4VPhysicalVolume *B03PVolumeStore::
GetPVolume(const G4String &name) const {
  const G4VPhysicalVolume *pvol = 0;
  for (B03SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    if (vol.GetName() == name) {
      pvol =  &vol;
    } 
  }
  if (!pvol) {
    G4cout << "B03PVolumeStore::GetPVolume: no physical volume named: " 
	   << name << ", found" << G4endl;
  }
  return pvol;
}

G4String B03PVolumeStore::GetPNames() const {
  G4String NameString;
  for (B03SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
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
