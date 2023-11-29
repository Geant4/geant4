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
// modified by I.Hrivnacova, 13.10.99

#include "globals.hh"
#include "G3VolTableEntry.hh"
#include "G3VolTable.hh"
#include "G3RotTable.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G3Pos.hh"
#include "G3toG4.hh"

G3VolTableEntry::G3VolTableEntry(G4String& vname, G4String& shape, 
			         G4double* rpar, G4int npar, G4int nmed, 
			         G4VSolid* solid, G4bool hasNegPars)
  : fVname(vname), fShape(shape), fRpar(0), fNpar(npar), fNmed(nmed), 
    fSolid(solid), fLV(0), fHasNegPars(hasNegPars), fHasMANY(false),
    fDivision(0)		  		     
{
  if (npar>0 && rpar!=0) {
    fRpar = new G4double[npar];
    for (G4int i=0; i<npar; i++) fRpar[i] = rpar[i];
  }
  fClones.push_back(this);
}

G3VolTableEntry::~G3VolTableEntry(){
  if (fRpar!=0) delete [] fRpar;
  delete fDivision;
}

G4bool
G3VolTableEntry::operator == ( const G3VolTableEntry& lv) const {
  return (this==&lv) ? true : false;
}

void 
G3VolTableEntry::AddG3Pos(G3Pos* aG3Pos){

  // insert this position to the vector
  G3Vol.CountG3Pos();
  fG3Pos.push_back(aG3Pos);

  // pass MANY info 
  G4String vonly = aG3Pos->GetOnly();
  if (vonly == "MANY") SetHasMANY(true);
}

void 
G3VolTableEntry::AddDaughter(G3VolTableEntry* aDaughter){
  if (FindDaughter(aDaughter->GetName()) == 0) {
    fDaughters.push_back(aDaughter);
  }
}

void 
G3VolTableEntry::AddMother(G3VolTableEntry* itsMother){
  if (FindMother(itsMother->GetName()) == 0) {
    fMothers.push_back(itsMother);
  }  
}

void 
G3VolTableEntry::AddClone(G3VolTableEntry* itsClone){
  if (FindClone(itsClone->GetName()) == 0) {
    fClones.push_back(itsClone);
  }  
}

void 
G3VolTableEntry::AddOverlap(G3VolTableEntry* overlap){
    fOverlaps.push_back(overlap);
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
    G4String err_message = "Old daughter " + vteOld->GetName()
                         + " does not exist.";
    G4Exception("G3VolTableEntry::ReplaceDaughter()", "G3toG40007",
                FatalException, err_message);
    return;
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
    G4String err_message = "Old mother " + vteOld->GetName()
                         + " does not exist.";
    G4Exception("G3VolTableEntry::ReplaceMother()", "G3toG40008",
                FatalException, err_message);
    return;
  }      
  fMothers[index] = vteNew;
}

G3VolTableEntry*
G3VolTableEntry::FindDaughter(const G4String& Dname){
  for (G4int idau=0; idau<GetNoDaughters(); idau++){
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
  G4cout << "HasMANY: " << fHasMANY << G4endl;
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

void G3VolTableEntry::SetNRpar(G4int npar, G4double* rpar) {
  if (npar != fNpar) {
    fNpar = npar;
    delete [] fRpar;
    fRpar = new G4double[fNpar];
  }      
  for (G4int i=0; i<fNpar; i++) fRpar[i] = rpar[i];
}  

void G3VolTableEntry::SetHasNegPars(G4bool hasNegPars) {
  fHasNegPars = hasNegPars;
}

void G3VolTableEntry::SetHasMANY(G4bool hasMANY) {
  fHasMANY = hasMANY;
}

void G3VolTableEntry::ClearG3PosCopy(G4int copy) {
  if (fG3Pos.size()>0 && copy>=0 && copy<G4int(fG3Pos.size())) {
     std::vector<G3Pos*>::iterator it=fG3Pos.begin();
     for(G4int j=0;j<copy;j++) it++;
     if(it!=fG3Pos.end()) {
         fG3Pos.erase(it);
     }
  }
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
  return (G4int)fG3Pos.size();
}

G3Pos* 
G3VolTableEntry::GetG3PosCopy(G4int copy) {
  if (fG3Pos.size()>0 && copy>=0)
    return fG3Pos[copy];
  else
    return 0;
}

G4bool 
G3VolTableEntry::HasNegPars(){
  return fHasNegPars;
}

G4bool 
G3VolTableEntry::HasMANY(){
  return fHasMANY;
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
  return (G4int)fDaughters.size();
}

G4int
G3VolTableEntry::GetNoMothers() {
  return (G4int)fMothers.size();
}

G4int
G3VolTableEntry::GetNoClones() {
  return (G4int)fClones.size();
}

G4int
G3VolTableEntry::GetNoOverlaps() {
  return (G4int)fOverlaps.size();
}

G3VolTableEntry* 
G3VolTableEntry::GetDaughter(G4int i) {
  if (i<G4int(fDaughters.size()) && i>=0)
    return fDaughters[i];
  else 
    return 0;
}

G3VolTableEntry*
G3VolTableEntry::GetMother(G4int i){
  if (i<G4int(fMothers.size()) && i>=0)
    return fMothers[i];
  else
    return 0;
}

// to be removed
G3VolTableEntry*
G3VolTableEntry::GetMother(){
  if (fMothers.size()>0)
    return fMothers[0];
  else
    return 0;  
}

G3VolTableEntry*
G3VolTableEntry::GetClone(G4int i){
  if (i<G4int(fClones.size()) && i>=0)
    return fClones[i];
  else
    return 0;
}

G3VolTableEntry*
G3VolTableEntry::GetMasterClone(){
  G3VolTableEntry* master;
  G4String name = fVname;
  if (G4StrUtil::contains(name, gSeparator)) {
    name = name.substr(0, name.find(gSeparator));
    master = G3Vol.GetVTE(name); 
  }
  else 
    master = this;

  return master;
}

std::vector<G3VolTableEntry*>*
G3VolTableEntry::GetOverlaps(){
  return &fOverlaps;
}
