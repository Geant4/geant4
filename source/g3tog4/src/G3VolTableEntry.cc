// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTableEntry.cc,v 1.6 2000-11-28 12:07:54 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I.Hrivnacova, 13.10.99

#include "globals.hh"
#include "G3VolTableEntry.hh"
#include "G3VolTable.hh"
#include "G4LogicalVolume.hh"
#include "G3Pos.hh"
#include "G3toG4.hh"

G3VolTableEntry::G3VolTableEntry(G4String& vname, G4String& shape, 
			     G4double* rpar, G4int npar, G4int nmed, 
			     G4VSolid* solid, G4bool hasNegPars)
  : fVname(vname), fShape(shape), fRpar(0), fNpar(npar), fNmed(nmed), 
    fSolid(solid), fLV(0), fHasNegPars(hasNegPars), fDivision(0)		  		     
{
  if (npar>0 && rpar!=0) {
    fRpar = new G4double[npar];
    for (G4int i=0; i<npar; i++) fRpar[i] = rpar[i];
  }
  fClones.insert(this);
}

G3VolTableEntry::~G3VolTableEntry(){
  if (fRpar!=0) delete [] fRpar;
  delete fDivision;
}

inline G4bool 
G3VolTableEntry::operator == ( const G3VolTableEntry& lv) const {
  return (this==&lv) ? true : false;
}

void 
G3VolTableEntry::AddG3Pos(G3Pos* aG3Pos){
  G3Vol.CountG3Pos();
  fG3Pos.insert(aG3Pos);
}

void 
G3VolTableEntry::AddDaughter(G3VolTableEntry* aDaughter){
  if (FindDaughter(aDaughter->GetName()) == 0) {
    fDaughters.insert(aDaughter);
  }
}

void 
G3VolTableEntry::AddMother(G3VolTableEntry* itsMother){
  if (FindMother(itsMother->GetName()) == 0) {
    fMothers.insert(itsMother);
  }  
}

void 
G3VolTableEntry::AddClone(G3VolTableEntry* itsClone){
  if (FindClone(itsClone->GetName()) == 0) {
    fClones.insert(itsClone);
  }  
}

void 
G3VolTableEntry::ReplaceDaughter(G3VolTableEntry* vteOld, 
                                 G3VolTableEntry* vteNew)
{
  G4int index = -1;
  for (G4int i=0; i<GetNoDaughters(); i++){
    if (fDaughters[i]->GetName() == vteOld->GetName()) index = i;
  }
  if (index<0) {
    G4Exception(
      "G3VolTableEntry::ReplaceDaughter: old daughter " +
       vteOld->GetName() + " does not exist.");
  }      
  fDaughters[index] = vteNew;
}

void 
G3VolTableEntry::ReplaceMother(G3VolTableEntry* vteOld, 
                               G3VolTableEntry* vteNew)
{
  G4int index = -1;
  for (G4int i=0; i<GetNoMothers(); i++){
    if (fMothers[i]->GetName() == vteOld->GetName()) index = i;
  }
  if (index<0) {
    G4Exception(
      "G3VolTableEntry::ReplaceMother: old mother " +
       vteOld->GetName() + " does not exist.");
  }      
  fMothers[index] = vteNew;
}

G3VolTableEntry*
G3VolTableEntry::FindDaughter(const G4String& Dname){
  for (int idau=0; idau<GetNoDaughters(); idau++){
    if (GetDaughter(idau)->GetName() == Dname) return GetDaughter(idau);
  }
  return 0;
}

G3VolTableEntry*
G3VolTableEntry::FindMother(const G4String& Mname){
  for (G4int i=0; i<GetNoMothers(); i++){
    G3VolTableEntry* mvte = GetMother(i);
    if (mvte->GetName() == Mname) return mvte;
  }
  return 0;
}

G3VolTableEntry*
G3VolTableEntry::FindClone(const G4String& Cname){
  for (G4int i=0; i<GetNoClones(); i++){
    G3VolTableEntry* cvte = GetClone(i);
    if (cvte->GetName() == Cname) return cvte;
  }
  return 0;
}

void G3VolTableEntry::PrintSolidInfo() {
// only parameters related to solid definition
// are printed
  G4cout << "VTE: " << fVname << " " << this << G4endl;
  G4cout << "Solid: " << fSolid << G4endl;
  G4cout << "Parameters (npar = " << fNpar << ") fRpar: ";
  for (G4int i=0; i<fNpar; i++) G4cout << fRpar[i] << " ";
  G4cout << G4endl;
  G4cout << "HasNegPars: " << fHasNegPars << G4endl;
  G4cout << "================================= " << G4endl;
}

void
G3VolTableEntry::SetName(G4String name){
  fVname = name;
}

void
G3VolTableEntry::SetLV(G4LogicalVolume* lv){
  fLV = lv;
}

void 
G3VolTableEntry::SetSolid(G4VSolid* solid){
  fSolid = solid;
}

void G3VolTableEntry::SetNmed(G4int nmed) {
  fNmed = nmed;
}

void G3VolTableEntry::SetNRpar(G4int npar, G4double* Rpar) {
  if (npar != fNpar) {
    fNpar = npar;
    delete [] fRpar;
    fRpar = new G4double[fNpar];
  }      
  for (G4int i=0; i<fNpar; i++) fRpar[i] = Rpar[i];
}  

void G3VolTableEntry::SetHasNegPars(G4bool hasNegPars) {
  fHasNegPars = hasNegPars;
}

void G3VolTableEntry::ClearG3PosCopy(G4int copy) {
  if (fG3Pos.entries()>0 && copy>=0 && copy<fG3Pos.entries()) 
    fG3Pos.removeAt(copy);
}

void G3VolTableEntry::ClearDivision() {
  delete fDivision;
  fDivision = 0;
}

G4String
G3VolTableEntry::GetName() {
  return fVname;
}

G4String
G3VolTableEntry::GetShape() {
  return fShape;
}

G4int
G3VolTableEntry::GetNmed() {
  return fNmed;
}

G4int 
G3VolTableEntry::GetNpar() {
  return fNpar;
}

G4double* 
G3VolTableEntry::GetRpar() {
  return fRpar;
}

G4int 
G3VolTableEntry::NPCopies() {
  return fG3Pos.entries();
}

G3Pos* 
G3VolTableEntry::GetG3PosCopy(G4int copy) {
  if (fG3Pos.entries()>0 && copy>=0)
    return fG3Pos[copy];
  else
    return 0;
}

G4bool 
G3VolTableEntry::HasNegPars(){
  return fHasNegPars;
}

G4VSolid*
G3VolTableEntry::GetSolid() {
  return fSolid;
}

G4LogicalVolume* 
G3VolTableEntry::GetLV() {
  return fLV;
}

G4int
G3VolTableEntry::GetNoDaughters() {
  return fDaughters.entries();
}

G4int
G3VolTableEntry::GetNoMothers() {
  return fMothers.entries();
}

G4int
G3VolTableEntry::GetNoClones() {
  return fClones.entries();
}

G3VolTableEntry* 
G3VolTableEntry::GetDaughter(G4int i) {
  if (i<fDaughters.entries() && i>=0)
    return fDaughters[i];
  else 
    return 0;
}

G3VolTableEntry*
G3VolTableEntry::GetMother(G4int i){
  if (i<fMothers.entries() && i>=0)
    return fMothers[i];
  else
    return 0;
}

// to be removed
G3VolTableEntry*
G3VolTableEntry::GetMother(){
  if (fMothers.entries()>0)
    return fMothers[0];
  else
    return 0;  
}

G3VolTableEntry*
G3VolTableEntry::GetClone(G4int i){
  if (i<fClones.entries() && i>=0)
    return fClones[i];
  else
    return 0;
}

G3VolTableEntry*
G3VolTableEntry::GetMasterClone(){
  G3VolTableEntry* master;
  G4String name = fVname;
  if (name.contains(gSeparator)) {
    name = name(0, name.first(gSeparator));
    master = G3Vol.GetVTE(name); 
  }
  else 
    master = this;

  return master;
}
