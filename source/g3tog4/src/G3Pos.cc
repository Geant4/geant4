// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3Pos.cc,v 1.4 1999-05-28 23:03:31 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3Pos.hh"
#include "VolTableEntry.hh"
#include "G3VolTable.hh"

extern G3VolTable G3Vol;

G3Pos::G3Pos(G4String& LV, G4int C, G4ThreeVector* Position, G4int irot, 
	     G4String& Only) 
  : _VTE(0),_Copy(0), _LV(" "),
    _Position(0),_Irot(0),_Only(" "){
      _VTE = G3Vol.GetVTE(LV);
      _LV = LV;
      _Copy = C;
      _Position = Position;
      _Irot = irot;
      _Only = Only;
};

G3Pos::~G3Pos(){;};

VolTableEntry*
G3Pos::GetVTE() {
  return _VTE;
};

G4bool 
G3Pos::operator == ( const G3Pos& lv) const {
  return (this==&lv) ? true : false;
};

G4String&
G3Pos::GetName() {
  return _LV;
};

G4int
G3Pos::GetIrot() {
  return _Irot;
};

G4int
G3Pos::GetCopy() {
  return _Copy;
};

G4ThreeVector* 
G3Pos::GetPos() {
  return _Position;
};



