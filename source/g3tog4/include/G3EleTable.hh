// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3EleTable.hh,v 1.5 2000-11-24 09:50:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// The table of elements.
// In the constructor an array of strings with element name,
// symbol and A are created. 
// The G4Element instances of given Z are created when
// GetEle(G4double Z) is called using the string array.
// For each Z the G4Element is created only once.

// ----------------------

#ifndef G4ELETABLE_HH
#define G4ELETABLE_HH 1

#include "g4rw/ctoken.h"
#include "g4rw/tpordvec.h"
#include "globals.hh"
#include "G4Element.hh"

class G3EleTable
{

public:  // with description

  G3EleTable();
  virtual ~G3EleTable();
  G4Element* GetEle(G4double Z);

private:

  void LoadUp();
  int parse(G4double& Z, char* name, char* sym, G4double& A); 

private:

  char** _EleNames;
  G4Element** _Ele;
  int _MaxEle;

};

extern G3EleTable G3Ele;
#endif
