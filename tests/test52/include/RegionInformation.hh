#ifndef REGIONINFORMATION_HH
#define REGIONINFORMATION_HH

#include "G4VUserRegionInformation.hh"
#include "globals.hh"


class RegionInformation : public G4VUserRegionInformation {
 
 public:
   RegionInformation() : isWorld(false), isTarget(false) {}
   ~RegionInformation() {}

   void Print() const {}

 private:
   G4bool isWorld;   
   G4bool isTarget;

 public:
   inline void FlagRegionAsWorld() { 
     isWorld = true;
     isTarget = false;
   }
   inline void FlagRegionAsTarget() {
     isTarget = true;
     isWorld = false;
   }
   inline G4bool IsWorld() const {
     return isWorld;
   }
   inline G4bool IsTarget() const {
     return isTarget;
   }
};

#endif // REGIONINFORMATION_HH
