// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.hh,v 1.7 2001-03-07 14:37:42 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  9th August 1998

#ifndef G4SCENELIST_HH
#define G4SCENELIST_HH

#include "G4Scene.hh"
#include "g4std/vector"

class G4SceneList: public G4std::vector <G4Scene*> {
public:
  void remove(G4Scene*);
};

typedef G4std::vector<G4Scene*>::iterator G4SceneListIterator;
typedef G4std::vector<G4Scene*>::const_iterator G4SceneListConstIterator;

#endif
