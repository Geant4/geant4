//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G3PartTable.hh,v 1.8 2001-07-11 09:58:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// G3 particles table.
// Maps G3 particles indices to their G4 particle object counterparts.

// ----------------------

#ifndef G3PARTTABLE_HH
#define G3PARTTABLE_HH 1
#include "g4std/map"
#include "G4ParticleDefinition.hh"

class G3PartTable
{

public:  // with description

  G3PartTable();
  virtual ~G3PartTable();
  G4ParticleDefinition* Get(G4int partid);
  void Put(G4int partid, G4ParticleDefinition* partpt);
  void PrintAll();

private:

  G4std::map<G4String, G4ParticleDefinition*, G4std::less<G4String> > PTD;
  void HashID(G4int partid, G4String* _HID);
  void HashID(G4int partid, G4String& _HID);
};

extern G3PartTable G3Part;

#endif
