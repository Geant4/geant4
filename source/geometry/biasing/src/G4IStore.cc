#include "G4IStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PTouchableKey.hh"
#include "G4PStepStream.hh"


G4IStore::G4IStore(G4VPhysicalVolume &worldvolume) :
  G4VIStore(worldvolume),
  fWorldVolume(worldvolume){}

G4VPhysicalVolume &G4IStore::GetWorldVolume(){
  return fWorldVolume;
}

void G4IStore::SetInternalIterator(const G4VPhysicalVolume &aVolume,
				   G4int aRepNum) const {
  if (!IsInWorld(aVolume)) {
    Error("SetInternalIterator: physical volume not in this World");
  }
  fCurrentIterator = fPtki.find(G4PTouchableKey(aVolume, aRepNum));
}

void G4IStore::AddImportanceRegion(G4double importance,
				   const G4VPhysicalVolume &aVolume,
				   G4int aRepNum){
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
				G4int aRepNum){
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
				 G4int aRepNum) const {
  
  SetInternalIterator(aVolume, aRepNum);
  if (fCurrentIterator==fPtki.end()) {
    Error("GetImportance: Region does not exist");
  }
  return fCurrentIterator->second;
}


G4double G4IStore::GetImportance(const G4PTouchableKey &ptk) const {
  fCurrentIterator = fPtki.find(ptk);
  if (fCurrentIterator==fPtki.end()) {
    G4cout << "PTouchableKey ptk: " << ptk << G4endl;
    G4cout << "Not found in: " << G4endl;
    G4cout << fPtki << G4endl;
    Error("GetImportance(ptk): Region does not exist");
  }
  return fCurrentIterator->second;
}

G4bool G4IStore::IsInWorld(const G4VPhysicalVolume &aVolume) const {
  return true;
  /*
  if (!aVolume) return false;
  if (*aVolume==G4ParallelWorld::GetWorldVolume()) {
    return true;
  }
  return IsInWorld(aVolume->GetMother());
  */
}



