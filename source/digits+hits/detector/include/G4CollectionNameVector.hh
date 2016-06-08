// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CollectionNameVector.hh,v 1.1 2001/02/08 06:08:53 asaim Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#ifndef G4CollectionNameVector_H 
#define G4CollectionNameVector_H 1

#include "globals.hh"
#include "g4std/vector"

class G4CollectionNameVector : public G4std::vector<G4String>
{
  public:
    G4CollectionNameVector() {;}
    virtual ~G4CollectionNameVector() {;}

    void insert(G4String s) { push_back(s); }
};

#endif

