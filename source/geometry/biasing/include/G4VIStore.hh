#ifndef G4VIStore_hh
#define G4VIStore_hh  G4VIStore_hh

#include "globals.hh"


class G4PTouchableKey;
class G4VPhysicalVolume;

class  G4VIStore {
public:
  G4VIStore(G4VPhysicalVolume &worldVolume){}
  virtual  ~G4VIStore(){}
  virtual void AddImportanceRegion(G4double importance,
				   const G4VPhysicalVolume &,
				   G4int aRepNum = 0) = 0;
  virtual void ChangeImportance(G4double importance,
				const G4VPhysicalVolume &,
				G4int aRepNum = 0) = 0;
  virtual G4double GetImportance(const G4VPhysicalVolume &,
				 G4int aRepNum = 0) const = 0;
  virtual G4double GetImportance(const G4PTouchableKey &ptk) const = 0;
  virtual G4VPhysicalVolume &GetWorldVolume() = 0;
};

#endif









