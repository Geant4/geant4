// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ggclos.cc,v 1.5 1999-12-05 17:50:12 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G3toG4.hh"
#include "G3VolTable.hh"

void PG4ggclos(){
  G4ggclos();
}

void G4ggclos(){
  G4cout << "G4ggclos: setting top-level VolTableEntry" << endl;
  G3Vol.SetFirstVTE();
}
