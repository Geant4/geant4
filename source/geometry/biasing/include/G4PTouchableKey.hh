#ifndef G4PTouchableKey_hh
#define G4PTouchableKey_hh G4PTouchableKey_hh

#include "globals.hh"
class G4VPhysicalVolume;

class G4PTouchableKey {
public:
  G4PTouchableKey(const G4VPhysicalVolume &aVolume,
		G4int RepNum) :
    fVPhysiclaVolume(&aVolume),
    fRepNum(RepNum){}
  
  const G4VPhysicalVolume *fVPhysiclaVolume;
  G4int fRepNum;
};

class G4PTkComp {
public:
  bool operator() (const G4PTouchableKey &k1, const G4PTouchableKey &k2) 
    const {
    if (k1.fVPhysiclaVolume != k2.fVPhysiclaVolume) {
      return  k1.fVPhysiclaVolume < k2.fVPhysiclaVolume;
    } else {
      return k1.fRepNum < k2.fRepNum;
    }
  }
};

bool operator==(const G4PTouchableKey &k1, const G4PTouchableKey &k2);
bool operator!=(const G4PTouchableKey &k1, const G4PTouchableKey &k2);


#endif
