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
// $Id: G4IStore.cc,v 1.12 2003/04/03 10:45:51 dressel Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4IStore.cc
//
// ----------------------------------------------------------------------

#include "G4IStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryCell.hh"
#include "G4GeometryCellStepStream.hh"
#include "G4LogicalVolume.hh"

G4IStore::G4IStore(const G4VPhysicalVolume &worldvolume) :
  fWorldVolume(worldvolume)
{}

G4IStore::~G4IStore()
{}

const G4VPhysicalVolume &G4IStore::GetWorldVolume() const 
{
  return fWorldVolume;
}

void G4IStore::SetInternalIterator(const G4GeometryCell &gCell) const
{
  fCurrentIterator = fGeometryCelli.find(gCell);
}

void G4IStore::AddImportanceGeometryCell(G4double importance,
			 const G4GeometryCell &gCell){
  if (importance < 0 ) {
    Error("AddImportanceGeometryCell: invalid importance value given");
  }  
  if (!IsInWorld(gCell.GetPhysicalVolume()) ) {
    Error("AddImportanceGeometryCell: physical volume not in this World");
  }
  SetInternalIterator(gCell);
  if (fCurrentIterator!=fGeometryCelli.end()) {
    Error("AddImportanceGeometryCell: Region allready exists");
  }
  fGeometryCelli[gCell] = importance;
}

void G4IStore::AddImportanceGeometryCell(G4double importance,
				   const G4VPhysicalVolume &aVolume,
				   G4int aRepNum)
{
  AddImportanceGeometryCell(importance,
		      G4GeometryCell(aVolume, aRepNum));
}

void G4IStore::ChangeImportance(G4double importance,
				const G4GeometryCell &gCell){
  if (importance < 0 ) {
    Error("ChangeImportance: Invalid importance value given");
  }
  if (!IsInWorld(gCell.GetPhysicalVolume()) ) {
    Error("ChangeImportance: physical volume not in this World");
  }
  SetInternalIterator(gCell);
  if (fCurrentIterator==fGeometryCelli.end()) {
    Error("ChangeImportance: Region does not exist");
  }
  fGeometryCelli[gCell] = importance;

}
void G4IStore::ChangeImportance(G4double importance,
				const G4VPhysicalVolume &aVolume,
				G4int aRepNum)
{
  ChangeImportance(importance, G4GeometryCell(aVolume, aRepNum));
}

G4double G4IStore::GetImportance(const G4VPhysicalVolume &aVolume,
				 G4int aRepNum) const
{  
  SetInternalIterator(G4GeometryCell(aVolume, aRepNum));
  if (fCurrentIterator==fGeometryCelli.end()) {
    Error("GetImportance: Region does not exist");
  }
  return (*fCurrentIterator).second;
}


G4double G4IStore::GetImportance(const G4GeometryCell &gCell) const
{
  SetInternalIterator(gCell);
  if (fCurrentIterator==fGeometryCelli.end()) {
    G4cout << "PGeometryCell gCell: " << gCell << G4endl;
    G4cout << "Not found in: " << G4endl;
    G4cout << fGeometryCelli << G4endl;
    Error("GetImportance(gCell): Region does not exist");
  }
  return (*fCurrentIterator).second;
}

G4bool G4IStore::IsKnown(const G4GeometryCell &gCell) const {
  G4bool inWorldKnown(IsInWorld(gCell.GetPhysicalVolume()));
		      
  if ( inWorldKnown ) {
    SetInternalIterator(gCell);
    inWorldKnown = (fCurrentIterator!=fGeometryCelli.end());
  }
  return inWorldKnown;
}

G4bool G4IStore::IsInWorld(const G4VPhysicalVolume &aVolume) const
{
  G4bool isIn(true);
  if (!(aVolume == fWorldVolume)) {
    isIn = fWorldVolume.GetLogicalVolume()->IsAncestor(&aVolume);
  }
  return isIn;
}



void G4IStore::Error(const G4String &m) const
{
  G4cout << "ERROR - G4IStore::" << m << G4endl;
  G4Exception("Program aborted.");
}
