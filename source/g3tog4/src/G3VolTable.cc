// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.14 1999-07-29 03:44:38 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <iomanip.h>
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
  // G4cout << "Deleted G3VolTable..." << endl;
};

VolTableEntry*
G3VolTable::GetVTE(const G4String& Vname) {
  return _VTD->findValue(&Vname);
};

void 
G3VolTable::ListVTE(){
  RWTPtrHashDictionaryIterator<G4String, VolTableEntry> iter(*_VTD);
  if (_VTD->entries()>0) {
    for (int i=0;iter();i++){
      _VTE = iter.value();
      G4cout << "G3VolTable element " << setw(3) << i << " name "
	     << _VTE->GetName() << " has " << _VTE->GetNoDaughters() 
	     << " daughters" << endl;
    }
  }
}

RWTPtrHashDictionary <G4String, VolTableEntry>* 
G3VolTable::GetVTD() {return _VTD;}

VolTableEntry*
G3VolTable::PutVTE(VolTableEntry* aVolTableEntry){
  
  if (GetVTE(aVolTableEntry->GetName()) == 0 ){
    
    // create a hash key
    G4String* _HashID = new G4String(aVolTableEntry->GetName());
    
    if (_FirstKey == "UnDefined") _FirstKey = *(_HashID);
    
    // insert into dictionary
    _VTD->insertKeyAndValue(_HashID, aVolTableEntry);
    _NVTE++;
  }
  return GetVTE(aVolTableEntry->GetName());
};

void 
G3VolTable::CountG3Pos(){
  _NG3Pos++;
}

void
G3VolTable::SetFirstVTE(){
  G3toG4TopVTE = _VTD->findValue(&_FirstKey);
  if (G3toG4TopVTE->NPCopies() > 0) {
    _FirstKey = G3toG4TopVTE->GetMother()->GetName();
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


