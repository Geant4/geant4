// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: VolTableEntry.cc,v 1.5 1999-05-26 05:15:29 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3VolTable.hh"

// This "class" remembers the top-level G3toG4 mother volume. It also returns
// logical and physical volume pointers, given their names and PV copy numbers.
// If no LV/PV name is given, the top-level LV/PV pointer is returned.

#include "globals.hh"
#include "VolTableEntry.hh"
#include "G4LogicalVolume.hh"
#include "G3Pos.hh"

VolTableEntry::VolTableEntry(G4String& Vname, G4String& Shape, 
			     G4double* Rpar,  G4int Npar, 
			     G4int Nmed, 
			     G4Material* Mat, G4VSolid* Solid,
			     G4bool Deferred, G4bool NegVolPars)
  : _Rpar(0), _Npar(0), _Nmed(0), _Mat(0), _Solid(0), _LV(0), _Deferred(0), 
    _NegVolPars(0), _Daughters(0), _G3Pos(0){
      _Nmed = Nmed;
      _Vname = Vname;
      _Shape = Shape;
      if (Rpar!=0 && Npar>0) {
	_Rpar = new G4double[Npar];
	for (int i=0; i<Npar; i++) _Rpar[i] = Rpar[i];
      }
      _Npar = Npar;
      _Mat = Mat;
      _Solid = Solid;
      _NegVolPars = NegVolPars;
      _Deferred = Deferred;
}

void
VolTableEntry::SetLV(G4LogicalVolume* ll){
  _LV = ll;
};

G4String
VolTableEntry::GetName() {
  return _Vname;
};

G4VSolid*
VolTableEntry::GetSolid() {
  return _Solid;
};

G4bool 
VolTableEntry::HasDeferred(){
  return _Deferred;
};

G4String
VolTableEntry::GetShape() {
  return _Shape;
};

G4double* 
VolTableEntry::GetRpar() {
  return _Rpar;
};

G4int 
VolTableEntry::GetNpar() {
  return _Npar;
};

G4Material* 
VolTableEntry::GetMaterial() {
  return _Mat;
}

G4LogicalVolume* 
VolTableEntry::GetLV() {
  return _LV;
};

G4int 
VolTableEntry::GetNmed() {
  return _Nmed;
};

G4bool 
VolTableEntry::HasNegVolPars(){
  return _NegVolPars;
};

G3Pos* 
VolTableEntry::GetG3PosCopy(G4int copy) {
  return _G3Pos[copy];
}

G4int 
VolTableEntry::NPCopies() {
  return _G3Pos.length();
};

void
VolTableEntry::AddG3Pos(G3Pos* aG3Pos){
  G3Vol.CountG3Pos();
  _G3Pos.resize(_G3Pos.length()+1);
  _G3Pos.insert(aG3Pos);
};

void
VolTableEntry::AddDaughter(VolTableEntry* aDaughter){
  if (FindDaughter(aDaughter->GetName()) == 0) {
    _Daughters.resize(_Daughters.length()+1);
    _Daughters.insert(aDaughter);
  }
};

VolTableEntry*
VolTableEntry::FindDaughter(const G4String& Dname){
  for (int idau=0; idau<GetNoDaughters(); idau++){
    if (GetDaughter(idau)->GetName() == Dname) return GetDaughter(idau);
  }
  return 0;
};

G4int
VolTableEntry::GetNoDaughters() {
  return _Daughters.length();
};

VolTableEntry* 
VolTableEntry::GetDaughter(G4int i) {
  return _Daughters[i];
};

inline G4bool 
VolTableEntry::operator == ( const VolTableEntry& lv) const {
  return (this==&lv) ? true : false;
};

VolTableEntry::~VolTableEntry(){
  if (_Rpar!=0 && _Npar>0) delete [] _Rpar;
};








