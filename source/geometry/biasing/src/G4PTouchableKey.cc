#include "G4PTouchableKey.hh"

bool operator==(const G4PTouchableKey &k1, const G4PTouchableKey &k2) {
  if (k1.fVPhysiclaVolume != k2.fVPhysiclaVolume) return false;
  if (k1.fRepNum != k2.fRepNum) return false;
  return true;
}

bool operator!=(const G4PTouchableKey &k1, const G4PTouchableKey &k2) {
  if (k1.fVPhysiclaVolume != k2.fVPhysiclaVolume) return true;
  if (k1.fRepNum != k2.fRepNum) return true;
  return false;  
}
