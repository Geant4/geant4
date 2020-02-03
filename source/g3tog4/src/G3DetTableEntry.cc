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
//
#include "globals.hh"
#include "G3DetTableEntry.hh"
#include "G4VSensitiveDetector.hh"

G3DetTableEntry::G3DetTableEntry(G4String& set, G4String& det, G4int id, 
				 G4VSensitiveDetector* D){
  _set = set;
  _det = det;
  _id  = id;
  _detpt = D;
}

G3DetTableEntry::~G3DetTableEntry(){;}

G4VSensitiveDetector* 
G3DetTableEntry::GetSD(){
  return _detpt;
}

G4String 
G3DetTableEntry::GetSet(){
  return _set;
}

G4String 
G3DetTableEntry::GetDet(){
  return _det;
}

G4int
G3DetTableEntry::GetID(){
  return _id;
}

