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
// $Id: G3DetTableEntry.cc,v 1.4 2001-07-11 09:58:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

