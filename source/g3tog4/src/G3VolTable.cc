// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.8 1999-05-22 06:51:00 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3VolTable.hh"
#include "G3Pos.hh"

G3VolTable::G3VolTable() 
  : _VTD(0), G3toG4LogicalMother(0), _FirstKey(0) {
    _VTD = new RWTPtrHashDictionary<G4String,VolTableEntry>(RWCString::hash);
};

G3VolTable::~G3VolTable(){
  _VTD->clearAndDestroy();
  delete _VTD;
  G4cout << "Destructing VolTable Hash Dictionary" << endl;
};

VolTableEntry*
G3VolTable::GetVTE(const G4String& Vname){
  return _VTD->findValue(&Vname);
};

RWTPtrHashDictionary <G4String, VolTableEntry>* 
G3VolTable::GetVTD(){return _VTD;}

void
G3VolTable::PutVTE(const G4String& Vname, const G4String& Shape, 
		   const G4double* Rpar,  const G4int Npar, const G4int Nmed, 
		   const G4Material* Mat, const G4VSolid* Solid,
		   const G4bool Deferred, const G4bool NegVolPars){
  if (GetVTE(Vname) == 0) {

    // create a VolTableEntry
    _VTE = new 
      VolTableEntry(Vname, Shape, Rpar, Npar, Nmed, Mat, Solid, Deferred, 
		    NegVolPars);

    // create a hash key
    G4String* _HashID = new G4String(Vname);

    if (_FirstKey == 0) _FirstKey = _HashID;

    // insert into dictionary
    _VTD->insertKeyAndValue(_HashID, _VTE);
  }
};

VolTableEntry* 
G3VolTable::GetFirstVTE(){
  _VTE = _VTD->findValue(_FirstKey);
  if (_VTE->NPCopies() > 0) {
    _FirstKey = _VTE->GetG3PosCopy(0)->GetMotherVTE()->GetName();
    _VTE = GetFirstVTE();
  }
  return _VTE;
}

G4LogicalVolume*
G3VolTable::GetG3toG4Mother() {
  return GetFirstVTE()->GetLV();
};
