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
// $Id: G3DetTable.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
#include "globals.hh"
#include "G3DetTable.hh"

typedef std::map<G4String, G3DetTableEntry*, std::less<G4String> >
::iterator DTDiterator;

G4String 
G3DetTable::MakeHash(G4String& set, G4String& det){;
  return set+" "+det;
}

G3DetTable::G3DetTable(){
}

G3DetTable::~G3DetTable(){
  if (DTD.size() > 0) {
    //    G4cout << "Deleting DTD" << G4endl;
    for (DTDiterator i=DTD.begin(); i != DTD.end(); i++) {
      delete (*i).second;
    }
    DTD.clear();
  }
}

G4VSensitiveDetector* 
G3DetTable::GetSD(G4String& set, G4String& det){

  // make hash ID
  const G4String ShashID = MakeHash(set, det);

  // search the map
  DTDiterator i = DTD.find(ShashID);
  G3DetTableEntry* DTE = (*i).second;
  if (DTE != 0) {
    return DTE->GetSD();
  } else {
    return 0;
  }  
}

G4int 
G3DetTable::GetID(G4String& set, G4String& det){

  // make hash ID
  G4String ShashID = MakeHash(set, det);

  // search the Hash Dictionary
  DTDiterator i = DTD.find(ShashID);
  G3DetTableEntry* DTE = (*i).second;
  if (DTE != 0) {
    return DTE->GetID();
  } else {
    return 0;
  }
}

void 
G3DetTable::Put(G4String& set, G4String& det, G4int id, 
		G4VSensitiveDetector* D){
  // make hash ID
  G4String ShashID = MakeHash(set, det);
  G3DetTableEntry* DTE = new G3DetTableEntry(set, det, id, D);
  G4cout << "Inserted DTE with id " << ShashID << G4endl;
  DTD[ShashID] = DTE;
}

void
G3DetTable::PrintAll(){
  if (DTD.size()>0){
    G4int count=0;
    G4cout << "Dump of DTD - " << DTD.size() << " entries:" << G4endl;
    for (DTDiterator i=DTD.begin(); i != DTD.end(); i++) {
      count++;
      G3DetTableEntry* DTE = (*i).second;
      G4cout << "DTD entry " << std::setw(3) << count << " sensitive detector name: " 
	     << DTE->GetSD()->GetName() << G4endl;
    }
  }
}





