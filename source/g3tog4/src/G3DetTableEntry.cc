// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTableEntry.cc,v 1.3 2000-03-02 17:54:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3DetTableEntry.hh"
#include "G4VSensitiveDetector.hh"

G3DetTableEntry::G3DetTableEntry(G4String& set, G4String& det, G4int id, 
				 G4VSensitiveDetector* D){
  _set = set;
  _det = det;
  _id  = id;
  _detpt = D;
}

G3DetTableEntry::~G3DetTableEntry(){;}

G4VSensitiveDetector* 
G3DetTableEntry::GetSD(){
  return _detpt;
}

G4String 
G3DetTableEntry::GetSet(){
  return _set;
}

G4String 
G3DetTableEntry::GetDet(){
  return _det;
}

G4int
G3DetTableEntry::GetID(){
  return _id;
}

