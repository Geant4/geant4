//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G3VolTable.hh 67982 2013-03-13 10:36:03Z gcosmo $
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

#include <map>
#include "G3VolTableEntry.hh"
#include "G3toG4Defs.hh"

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
  std::map<G4String, G3VolTableEntry*, std::less<G4String> > VTD;
  G4int _NG3Pos;
};

extern G3G4DLL_API G3VolTable G3Vol;

#endif
