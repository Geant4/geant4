// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.6 1999-05-18 18:43:20 lockman Exp $
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
  : G3toG4LogicalMother(0), G3toG4PhysicalMother(0), _VTE(0) {
    _VTD = new RWTPtrHashDictionary<G4String,VolTableEntry>(RWCString::hash);
};

G3VolTable::~G3VolTable(){
  _VTD->clearAndDestroy();
  delete _VTD;
  G4cout << "Destructing VolTable Hash Dictionary" << endl;
};

void
G3VolTable::SetMother(G4LogicalVolume* LV){
  if (G3toG4LogicalMother == 0) {
    G3toG4LogicalMother = LV;
  }
};

VolTableEntry::VolTableEntry(G4String& Vname, G4String& Shape, 
			     G4double* Rpar,
			     G4int Npar, G4Material* Mat, 
			     G4VSolid* Solid,
			     G4LogicalVolume* LV, G4bool Deferred, 
			     G4bool NegVolPars)
  : _Rpar(0), _Npar(0), _Mat(0), _Solid(0), _LV(0), _Deferred(0), 
    _NegVolPars(0){

    _Vname = Vname;
    _Shape = Shape;
    if (Rpar!=0 && Npar>0) {
      _Rpar = new G4double[Npar];
      for (int i=0; i<Npar; i++) _Rpar[i] = Rpar[i];
    }
    _Npar = Npar;
    _Mat = Mat;
    _Solid = Solid;
    _LV = LV;
    _NegVolPars = NegVolPars;
    _Deferred = Deferred;
};

VolTableEntry::~VolTableEntry(){
  if (_Rpar!=0 && _Npar>0) delete [] _Rpar;
};

VolTableEntry*
G3VolTable::GetVTE(const G4String& Vname){
  _VTE = _VTD->findValue(&Vname);
  return _VTE;
};

void
G3VolTable::PutLV(G4String& Vname, G4String& Shape, G4double* Rpar,
		  G4int Npar, G4Material* Mat, G4VSolid* Solid,
		  G4bool Deferred, G4bool NegVolPars){
  G4LogicalVolume* LV = 0;
  if (GetVTE(Vname) == 0) {
    if (!(Deferred || NegVolPars)) {
      // create logical volume now
      LV = new G4LogicalVolume(Solid, Mat, Vname);
      SetMother(LV);
    }
    // create a VolTableEntry
    _VTE = new VolTableEntry(Vname, Shape, Rpar, Npar, Mat, Solid, LV, 
			     Deferred, NegVolPars);
    G4String* _HashID = new G4String(Vname);
    // insert into dictionary
    _VTD->insertKeyAndValue(_HashID, _VTE);
  }
};

G4LogicalVolume* 
G3VolTable::GetLV(const G4String& name){
  if (name == "-") {
    return G3toG4LogicalMother;
  } else { 
    _VTE = _VTD->findValue(&name);
    if (_VTE != 0) {
      return _VTE->GetLV();
    } else {
      return 0;
    }
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





