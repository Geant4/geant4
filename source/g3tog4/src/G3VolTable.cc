// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.3 1999-05-12 08:09:55 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3VolTable.hh"

// This "class" remembers the top-level G3toG4 mother volume. It also returns
// logical and physical volume pointers, given their names and PV copy numbers.
// If no LV/PV name is given, the top-level LV/PV pointer is returned.

#include "globals.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

G3VolTable::G3VolTable() 
  : G3toG4LogicalMother(0), G3toG4PhysicalMother(0) {;}

G3VolTable::~G3VolTable(){;}

void
G3VolTable::SetMother(G4LogicalVolume* LV){
  if (G3toG4LogicalMother == 0) {
    G3toG4LogicalMother = LV;
  }
};

G4LogicalVolume* 
G3VolTable::GetLV(const G4String& name){
  if (name == "-") {
    return G3toG4LogicalMother;
  } else { // loop over the logical volume store, look for name
    G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
    G4int nlv = lvs->entries();
    for (int i=0; i<nlv; i++){
      G4String vn = (*lvs)[i]->GetName();
      if (name == vn) return (*lvs)[i];
    }
    // G4cerr << "G3VolTable: no such logical volume '" << name << "'" << endl;
    return 0;
  }
}    

G4VPhysicalVolume*
G3VolTable::GetPV(const G4String& name, G4int pcopy){
  if (name == "-") {
    return G3toG4PhysicalMother;
  } else { 
    G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
    G4int npv = pvs->entries();
    _copy = 0;
    for (int i=0; i<npv; i++){
      G4String vn = (*pvs)[i]->GetName();
      if (vn == name) {
	if (pcopy == _copy) {
	  return (*pvs)[i];
	} else {
	  _copy++;
	}
      }
    }
    G4cerr << "G3VolTable::Physical volume name '" << name 
	   << "' copy number " << pcopy << " not found" << endl;
    return 0;
  }
}
