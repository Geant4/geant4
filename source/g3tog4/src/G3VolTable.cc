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
// $Id: G3VolTable.cc,v 1.20 2001-07-11 09:58:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I.Hrivnacova, 13.10.99

#include "g4std/iomanip"
#include "globals.hh"
#include "G3VolTable.hh"
#include "G3Pos.hh"

typedef G4std::map<G4String, G3VolTableEntry*, G4std::less<G4String> >
::iterator VTDiterator;

G3VolTable::G3VolTable() 
  : G3toG4TopVTE(0), _FirstKey("UnDefined"), _NG3Pos(0){
}

G3VolTable::~G3VolTable(){
  if (VTD.size()>0){
    //    G4cout << "Deleting VTD" << G4endl;
    for (VTDiterator i=VTD.begin(); i != VTD.end(); i++) {
      delete (*i).second;
    }
    VTD.clear();
  }
}

G3VolTableEntry*
G3VolTable::GetVTE(const G4String& Vname) {
  VTDiterator i = VTD.find(Vname);
  if (i == VTD.end()) return 0;
  else                return (*i).second;
}

void 
G3VolTable::PrintAll(){
  if (VTD.size()){
    G4int i=0;
    G4cout << "Dump of VTD - " << VTD.size() << " entries:" << G4endl;
    VTEStat();
    for (VTDiterator v=VTD.begin(); v != VTD.end(); v++){
      G3VolTableEntry* VTE = (*v).second;
      G4cout << "G3VolTable element " << G4std::setw(3) << i++ << " name "
	     << VTE->GetName() << " has " << VTE->GetNoDaughters() 
	     << " daughters" << G4endl;
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
}

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
}

void
G3VolTable::VTEStat() {
  G4cout << "Instantiated " << VTD.size() << 
    " volume table entries \n" 
 	 << "                      " << _NG3Pos << " positions." << G4endl;
}

void 
G3VolTable::Clear() {
  if (VTD.size()>0){
    for (VTDiterator i=VTD.begin(); i != VTD.end(); i++) {
      delete (*i).second;
    }
    VTD.clear();
  }
  G3toG4TopVTE = 0;
  _FirstKey = "UnDefined";
  _NG3Pos = 0;
}  
