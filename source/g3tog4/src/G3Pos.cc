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
// $Id: G3Pos.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// modified by I.Hrivnacova, 13.10.99

#include "globals.hh"
#include "G3Pos.hh"

G3Pos::G3Pos(G4String motherName, G4int C, G4ThreeVector* Position, G4int irot, 
	     G4String Only) 
  : _MotherName(motherName),
    _Copy(C), 
    _Position(Position), 
    _Irot(irot), 
    _Only(Only)
{
  if (_Only == "MANY") {
    // warning when MANY position is created
    G4String text = "G3Pos warning: Not supported MANY option has been encountered.\n";
    text = text +   "               It may cause overlapping volumes.";
    G4cerr << text << G4endl;
  }
}

G3Pos::~G3Pos(){;}

G4bool 
G3Pos::operator == ( const G3Pos& lv) const {
  return (this==&lv) ? true : false;
}

G4String&
G3Pos::GetMotherName() {
  return _MotherName;
}

G4int
G3Pos::GetIrot() {
  return _Irot;
}

G4int
G3Pos::GetCopy() {
  return _Copy;
}

G4ThreeVector* 
G3Pos::GetPos() {
  return _Position;
}

G4String&
G3Pos::GetOnly() {
  return _Only;
}
