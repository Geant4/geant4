// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.17 1999-12-05 17:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I.Hrivnacova, 13.10.99

#include "g4std/iomanip"
#include "globals.hh"
#include "G3VolTable.hh"
#include "G3Pos.hh"

typedef G4std::map<G4String, G3VolTableEntry*, less<G4String> >
::iterator VTDiterator;

G3VolTable::G3VolTable() 
  : G3toG4TopVTE(0), _FirstKey("UnDefined"), _NG3Pos(0){
};

G3VolTable::~G3VolTable(){
  if (VTD.size()>0){
    //    G4cout << "Deleting VTD" << endl;
    for (VTDiterator i=VTD.begin(); i != VTD.end(); i++) {
      delete (*i).second;
    }
    VTD.clear();
  }
};

G3VolTableEntry*
G3VolTable::GetVTE(const G4String& Vname) {
  VTDiterator i = VTD.find(Vname);
  return (*i).second;
};

void 
G3VolTable::PrintAll(){
  if (VTD.size()){
    G4int i=0;
    G4cout << "Dump of VTD - " << VTD.size() << " entries:" << endl;
    VTEStat();
    for (VTDiterator v=VTD.begin(); v != VTD.end(); v++){
      G3VolTableEntry* VTE = (*v).second;
      G4cout << "G3VolTable element " << setw(3) << i++ << " name "
	     << VTE->GetName() << " has " << VTE->GetNoDaughters() 
	     << " daughters" << endl;
    }
  }
}

G3VolTableEntry*
G3VolTable::PutVTE(G3VolTableEntry* aG3VolTableEntry){
  
  if (GetVTE(aG3VolTableEntry->GetName()) == 0 ){
    
    // create a hash key
    G4String HashID = aG3VolTableEntry->GetName();
    
    if (_FirstKey == "UnDefined") _FirstKey = HashID;
    
    // insert into dictionary
    VTD[HashID] = aG3VolTableEntry;
  }
  return GetVTE(aG3VolTableEntry->GetName());
};

void 
G3VolTable::CountG3Pos(){
  _NG3Pos++;
}

void
G3VolTable::SetFirstVTE(){
  G3toG4TopVTE = VTD[_FirstKey];
  
  if (G3toG4TopVTE->NPCopies() > 0) {
    _FirstKey = G3toG4TopVTE->GetMother()->GetName();
    SetFirstVTE();
  }
}

G3VolTableEntry*
G3VolTable::GetFirstVTE() {
  return G3toG4TopVTE;
};

void
G3VolTable::VTEStat() {
  G4cout << "Instantiated " << VTD.size() << 
    " volume table entries \n" 
 	 << "                      " << _NG3Pos << " positions." << endl;
}
