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
// $Id: G4IStore.cc,v 1.3 2002-04-09 16:23:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4IStore.cc
//
// ----------------------------------------------------------------------

#include "G4IStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PTouchableKey.hh"
#include "G4PStepStream.hh"


G4IStore::G4IStore(G4VPhysicalVolume &worldvolume) :
  G4VIStore(worldvolume),
  fWorldVolume(worldvolume)
{}

G4IStore::~G4IStore()
{}

G4VPhysicalVolume &G4IStore::GetWorldVolume()
{
  return fWorldVolume;
}

void G4IStore::SetInternalIterator(const G4VPhysicalVolume &aVolume,
				   G4int aRepNum) const
{
  if (!IsInWorld(aVolume)) {
    Error("SetInternalIterator: physical volume not in this World");
  }
  fCurrentIterator = fPtki.find(G4PTouchableKey(aVolume, aRepNum));
}

void G4IStore::AddImportanceRegion(G4double importance,
				   const G4VPhysicalVolume &aVolume,
				   G4int aRepNum)
{
  if (importance <=0 ) {
    Error("AddImportanceRegion: invalid importance value given");
  }
  SetInternalIterator(aVolume, aRepNum);
  if (fCurrentIterator!=fPtki.end()) {
    Error("AddImportanceRegion: Region allready exists");
  }
  fPtki[G4PTouchableKey(aVolume, aRepNum)] = importance;
}

void G4IStore::ChangeImportance(G4double importance,
				const G4VPhysicalVolume &aVolume,
				G4int aRepNum)
{
  if (importance <=0 ) {
    Error("ChangeImportance: Invalid importance value given");
  }
  SetInternalIterator(aVolume, aRepNum);
  if (fCurrentIterator==fPtki.end()) {
    Error("ChangeImportance: Region does not exist");
  }
  fPtki[G4PTouchableKey(aVolume, aRepNum)] = importance;
}

G4double G4IStore::GetImportance(const G4VPhysicalVolume &aVolume,
				 G4int aRepNum) const
{  
  SetInternalIterator(aVolume, aRepNum);
  if (fCurrentIterator==fPtki.end()) {
    Error("GetImportance: Region does not exist");
  }
  return (*fCurrentIterator).second;
}


G4double G4IStore::GetImportance(const G4PTouchableKey &ptk) const
{
  fCurrentIterator = fPtki.find(ptk);
  if (fCurrentIterator==fPtki.end()) {
    G4cout << "PTouchableKey ptk: " << ptk << G4endl;
    G4cout << "Not found in: " << G4endl;
    G4cout << fPtki << G4endl;
    Error("GetImportance(ptk): Region does not exist");
  }
  return (*fCurrentIterator).second;
}

G4bool G4IStore::IsInWorld(const G4VPhysicalVolume &aVolume) const
{
  return true;
  /*
  if (!aVolume) return false;
  if (*aVolume==G4ParallelWorld::GetWorldVolume()) {
    return true;
  }
  return IsInWorld(aVolume->GetMother());
  */
}

void G4IStore::Error(const G4String &m) const
{
  G4cout << "ERROR - G4IStore::" << m << G4endl;
  G4Exception("Program aborted.");
}
