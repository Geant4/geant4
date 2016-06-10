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
// $Id: G3DetTable.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// G3DetTable class

#ifndef G3DETTABLE_HH
#define G3DETTABLE_HH 1

#include <map>
#include "G3toG4Defs.hh"
#include "globals.hh"
#include "G3DetTableEntry.hh"

class G3DetTable {
private:
  std::map<G4String, G3DetTableEntry*, std::less<G4String> > DTD;
  G4String MakeHash(G4String& set, G4String& det);

public:
  G3DetTable();
  virtual ~G3DetTable();
  G4int GetID(G4String& set, G4String& det);
  void Put(G4String& set, G4String& det, G4int id, G4VSensitiveDetector* D);
  G4VSensitiveDetector* GetSD(G4String& set, G4String& det); 
  void PrintAll();
};

extern G3G4DLL_API G3DetTable G3Det;
#endif
