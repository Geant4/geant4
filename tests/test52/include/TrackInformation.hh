#ifndef TRACKINFORMATION_HH
#define TRACKINFORMATION_HH

#include "G4VUserTrackInformation.hh"
#include "G4Allocator.hh"
#include "globals.hh"


class TrackInformation : public G4VUserTrackInformation {
 
 public:
   TrackInformation() : hitTarget(false), escapeEnergy(0) {}
   ~TrackInformation() {}
   inline void *operator new(size_t);
   inline void operator delete(void* trackInfo);

   void Print() const {}

 private:
   G4bool hitTarget;
   G4double escapeEnergy;

 public:
   inline void SetHitTarget() { 
     hitTarget = true;
   }
   inline G4bool HitTarget() const {
     return hitTarget;
   }
   inline void SetEscapeEnergy(G4double en) {
     escapeEnergy = en;
   }
   inline G4double EscapeEnergy() {
     return escapeEnergy;
   }
};

extern G4Allocator<TrackInformation> TrackInformationAllocator;

inline void* TrackInformation::operator new(size_t) { 
void* trackInfo;
  trackInfo = (void*) TrackInformationAllocator.MallocSingle();
  return trackInfo;
}

inline void TrackInformation::operator delete(void* trackInfo) { 
  TrackInformationAllocator.FreeSingle((TrackInformation*) trackInfo);
}

#endif // TRACKINFORMATION_HH
