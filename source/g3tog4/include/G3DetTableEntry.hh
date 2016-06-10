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
// $Id: G3DetTableEntry.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// G3DetTableEntry class

#ifndef DETTABLEENTRY_HH
#define DETTABLEENTRY_HH 1

#include <map>
#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class G3DetTableEntry {
private:
  G4String _set;
  G4String _det;
  G4int _id;
  G4VSensitiveDetector* _detpt;

public:
  G3DetTableEntry(G4String& set, G4String& det, G4int id, 
		G4VSensitiveDetector* D);
  ~G3DetTableEntry();
  G4VSensitiveDetector* GetSD();
  G4String GetSet();
  G4String GetDet();
  G4int GetID();
};
#endif
