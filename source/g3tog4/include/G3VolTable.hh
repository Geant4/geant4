// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.hh,v 1.15 2000-11-24 09:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// G3 volumes table.
// The G3 volumes are represented with the G3VolTableEntry objects
// that are stored in the map (sorted by names).
// In the phase of filling the G3 tables (defining G3 geometry, 
// eg. by parsing the G3 input via clparse.cc)
// a G3 volume can be defined with incomplete parameters
// (negative or none) that have to be retrieved from its mother
// which may be defined later. These parameters are being resolved
// subsequently in the phase of filling G3 tables.
// That's why the G4 object counterparts (solids, logical volumes
// and physical volumes) can be created only after filling the G3 tables
// is finished and all incomplete parameters are resolved. 

// ----------------------
//
// modified by I.Hrivnacova, 13.10.99

#ifndef G3VOLTABLE_HH
#define G3VOLTABLE_HH 1

#include "g4std/map"
#include "G3VolTableEntry.hh"

class G4LogicalVolume;
class G4Material;
class G4VSolid;

class G3VolTable
{

public:  // with description

  G3VolTableEntry* PutVTE(G3VolTableEntry* aVTE);  
  G3VolTableEntry* GetVTE(const G4String& Vname);
  void PrintAll();
  G3VolTable();
  virtual ~G3VolTable();
  G4LogicalVolume* GetG3toG4Mother();
  G3VolTableEntry* GetFirstVTE();
  void SetFirstVTE();
  void VTEStat();
  void CountG3Pos();
  void Clear();

private:

  G3VolTableEntry* G3toG4TopVTE;
  G4String _FirstKey;
  G4std::map<G4String, G3VolTableEntry*, G4std::less<G4String> > VTD;
  G4int _NG3Pos;
};

extern G3VolTable G3Vol;

#endif
