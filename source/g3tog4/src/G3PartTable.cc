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
// $Id: G3PartTable.cc 67982 2013-03-13 10:36:03Z gcosmo $
//

#include "G4Types.hh"
#include <sstream>
#include <iomanip>
#include "G3PartTable.hh"

typedef std::map<G4String, G4ParticleDefinition*, std::less<G4String> >
::iterator PTDiterator;

G3PartTable::G3PartTable(){
}

G3PartTable::~G3PartTable(){
  if (PTD.size()>0){
    //    G4cout << "Deleting PTD" << G4endl;
    for (PTDiterator i=PTD.begin(); i != PTD.end(); i++) {
      delete (*i).second;
    }
    PTD.clear();
  }
}

G4ParticleDefinition*
G3PartTable::Get(G4int partid){
  G4String ShashID; // static
  HashID(partid, ShashID);
  PTDiterator i = PTD.find(ShashID);
  return (*i).second;
}

void 
G3PartTable::Put(G4int partid, G4ParticleDefinition *partpt){
  G4String ShashID; // static
  HashID(partid, ShashID);
  PTD[ShashID]=partpt;
}

void
G3PartTable::HashID(G4int partid, G4String& theHashID){
  std::ostringstream ostr;
  ostr << "Part" << partid << std::ends;
  theHashID = ostr.str();
}

void 
G3PartTable::HashID(G4int partid, G4String* theHashID){
  HashID(partid, *theHashID);
}

void
G3PartTable::PrintAll(){
  if (PTD.size()>0){
    G4int count=0;
    G4cout << "Dump of PTD - " << PTD.size() << " entries: " << G4endl;
    for (PTDiterator i=PTD.begin(); i != PTD.end(); i++) {
      count++;
      G4ParticleDefinition* aPTD = (*i).second;
      G4cout << "PTD entry " << std::setw(3) << count << " particle name: " 
	     << aPTD->GetParticleName() << G4endl;
    }
  }
}





