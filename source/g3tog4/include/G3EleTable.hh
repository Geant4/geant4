// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3EleTable.hh,v 1.4 1999-12-09 01:27:42 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4ELETABLE_HH
#define G4ELETABLE_HH 1

#include "g4rw/ctoken.h"
#include "g4rw/tpordvec.h"
#include "globals.hh"
#include "G4Element.hh"

class G3EleTable{

public:
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
