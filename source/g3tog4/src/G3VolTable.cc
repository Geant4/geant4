// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.9 1999-05-26 03:47:33 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3VolTable.hh"
#include "G3Pos.hh"

G3VolTable::G3VolTable() 
  : _VTD(0), G3toG4TopVTE(0), _FirstKey("UnDefined"), _NVTE(0), _NG3Pos(0){
    _VTD = new RWTPtrHashDictionary<G4String,VolTableEntry>(RWCString::hash);
};

G3VolTable::~G3VolTable(){
  _VTD->clearAndDestroy();
  delete _VTD;
  G4cout << "Destructing VolTable Hash Dictionary" << endl;
};

VolTableEntry*
G3VolTable::GetVTE(G4String& Vname) {
  return _VTD->findValue(&Vname);
};

RWTPtrHashDictionary <G4String, VolTableEntry>* 
G3VolTable::GetVTD() {return _VTD;}

void
G3VolTable::PutVTE(G4String& Vname, G4String& Shape, 
		   G4double* Rpar,  G4int Npar, G4int Nmed, 
		   G4Material* Mat, G4VSolid* Solid,
		   G4bool Deferred, G4bool NegVolPars){
  if (GetVTE(Vname) == 0) {

    _NVTE++;

    // create a VolTableEntry
    _VTE = new 
      VolTableEntry(Vname, Shape, Rpar, Npar, Nmed, Mat, Solid, Deferred, 
		    NegVolPars);

    // create a hash key
    G4String* _HashID = new G4String(Vname);

    if (_FirstKey == "UnDefined") _FirstKey = *(_HashID);

    // insert into dictionary
    _VTD->insertKeyAndValue(_HashID, _VTE);
  }
};

void 
G3VolTable::CountG3Pos(){
  _NG3Pos++;
}

void
G3VolTable::SetFirstVTE(){
  G3toG4TopVTE = _VTD->findValue(&_FirstKey);
  if (G3toG4TopVTE->NPCopies() > 0) {
    _FirstKey = G3toG4TopVTE->GetG3PosCopy(0)->GetMotherVTE()->GetName();
    SetFirstVTE();
  }
}

VolTableEntry*
G3VolTable::GetFirstVTE() {
  return G3toG4TopVTE;
};

void
G3VolTable::VTEStat() {
  G4cout << "VTEStat: instantiated " << _NVTE << " logical  and "
	 << _NG3Pos << " positioned volumes." << endl;
}


