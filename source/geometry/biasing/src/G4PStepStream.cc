#include "G4PStepStream.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PStep.hh"

ostream& operator<<(ostream &out, const G4PTouchableKey &tk){
  out << "Volume name = " << tk.fVPhysiclaVolume->GetName() << ", ";
  out << "Replica number = " << tk.fRepNum;
  return out;
}
ostream& operator<<(ostream &out, const G4PStep &ps){
  out << "PreTouchableKey : " <<  ps.fPreTouchableKey << " ";
  out << "PostTouchableKey: " <<  ps.fPostTouchableKey << " ";
  out << "CrossBoundary   : " <<  ps.fCrossBoundary << "\n";
  return out;
}
