#ifndef _G3MATTABLE_
#define _G3MATTABLE_ 1

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.hh,v 1.5 1999-07-21 08:39:43 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Maps G3 material indices to G4Material*

#include <rw/tphdict.h>
class G4Material;

class G3MatTable {
private:
  RWTPtrHashDictionary<G4String, G4Material>* _Mat;
  void HashID(G4int matid, G4String* _HID);
  void HashID(G4int matid, G4String& _HID);
public:
  G3MatTable();
  ~G3MatTable();
  G4Material* get(G4int matid);
  void put(G4int matid, G4Material* matpt);
  void print(G4int matid=-99999);
};

extern G3MatTable G3Mat;
#endif
