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
// $Id: G3PartTable.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// G3 particles table.
// Maps G3 particles indices to their G4 particle object counterparts.

// ----------------------

#ifndef G3PARTTABLE_HH
#define G3PARTTABLE_HH 1
#include <map>
#include "G3toG4Defs.hh"
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

  std::map<G4String, G4ParticleDefinition*, std::less<G4String> > PTD;
  void HashID(G4int partid, G4String* _HID);
  void HashID(G4int partid, G4String& _HID);
};

extern G3G4DLL_API G3PartTable G3Part;

#endif
