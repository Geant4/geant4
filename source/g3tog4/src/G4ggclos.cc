// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ggclos.cc,v 1.4 1999-05-26 03:49:48 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "G3VolTable.hh"
#include "G3G4Interface.hh"

void PG4ggclos(){
  G4ggclos();
}

void G4ggclos(){
  G4cout << "G4ggclos: setting top-level VolTableEntry" << endl;
  G3Vol.SetFirstVTE();
}
