#ifndef G4IStore_hh
#define G4IStore_hh G4IStore_hh

#include "G4VIStore.hh"
#include "G4PtkImportance.hh"

class G4IStore : public G4VIStore {
public:
  G4IStore(G4VPhysicalVolume &worldvolume);
  ~G4IStore(){};
  void AddImportanceRegion(G4double importance,
			   const G4VPhysicalVolume &,
			   G4int aRepNum = 0);
  void ChangeImportance(G4double importance,
			const G4VPhysicalVolume &,
			G4int aRepNum = 0);
  G4double GetImportance(const G4VPhysicalVolume &,
			 G4int aRepNum = 0) const ;
  G4double GetImportance(const G4PTouchableKey &ptk) const;
  G4VPhysicalVolume &GetWorldVolume();
  
private:
  G4bool IsInWorld(const G4VPhysicalVolume &) const;
  void SetInternalIterator(const G4VPhysicalVolume &,
			   G4int aRepNum) const;
  void Error(const G4String &m) const {
    G4cout << "ERROE: in G4IStore: " << m << G4endl;
    exit(1);
  }
  
  G4VPhysicalVolume &fWorldVolume;
  G4PtkImportance fPtki;


  mutable G4PtkImportance::const_iterator fCurrentIterator;
};

#endif
