#ifndef _G3MEDTABLE_
#define _G3MEDTABLE_ 1

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.hh,v 1.4 1999-07-20 14:16:23 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Maps G3 tracking media indices to G4Material*

#include <rw/tphdict.h>
class G4Material;

class G3MedTable {
private:
  RWTPtrHashDictionary<G4String, G4Material>* _Med;
  void HashID(G4int matid, G4String* _HID);
  void HashID(G4int matid, G4String& _HID);
public:
  G3MedTable();
  ~G3MedTable();
  G4Material* get(G4int medid);
  void put(G4int medid, G4Material* matpt);
};

extern G3MedTable G3Med;
#endif
