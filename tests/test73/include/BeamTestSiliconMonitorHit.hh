//
#ifndef BEAMTESTSILICONMONITORHIT_HH
#define BEAMTESTSILICONMONITORHIT_HH

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include "G4ParticleDefinition.hh"
#include "G4VHit.hh"
#include "G4StepStatus.hh"

class BeamTestSiliconMonitorHit : public G4VHit {

public:
  
  // Constructors
  BeamTestSiliconMonitorHit();

  // Destructor
  virtual ~BeamTestSiliconMonitorHit();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods
  virtual void Draw();

  virtual void Print();

/*  // Incidence Information, getters and setters 
  inline void SetIncidenceDefinition(G4ParticleDefinition* pd) {fIPD = pd;}
  inline G4ParticleDefinition* GetIncidenceDefinition() const {return fIPD;}

  inline void SetIncidenceKineticEnergy(G4double e) {fIKEnergy = e;}
  inline G4double GetIncidenceKineticEnergy() const {return fIKEnergy;}

  inline void SetIncidencePosition(G4ThreeVector position) {fIPosition = position;}
  inline G4ThreeVector GetIncidencePosition() const {return fIPosition;}

  inline void SetIncidenceMomentumDirection(G4ThreeVector momentum) {fIMomentumD = momentum;}
  inline G4ThreeVector GetIncidenceMomentumDirection() const {return fIMomentumD;
*/
  // Exiting Information, getters and setters
  inline void SetExitDefinition(G4ParticleDefinition* pd) {fEPD = pd;}
  inline G4ParticleDefinition* GetExitDefinition() const {return fEPD;}

  inline void SetExitKineticEnergy(G4double e) {fEKEnergy = e;}
  inline G4double GetExitKineticEnergy() const {return fEKEnergy;}

  inline void SetExitPosition(G4ThreeVector position) {fEPosition = position;}
  inline G4ThreeVector GetExitPosition() const {return fEPosition;}

  inline void SetExitMomentumDirection(G4ThreeVector momentum) {fEMomentumD = momentum;}
  inline void SetExitMomentum(G4ThreeVector momentum) {fEMomentum = momentum;}
  inline G4ThreeVector GetExitMomentumDirection() const {return fEMomentumD;}
  inline G4ThreeVector GetExitMomentum() const {return fEMomentum;}
    inline void SetChamberNumber( G4int val ) { fChamberNum = val; }
    inline G4int GetChamberNumber() const { return fChamberNum; }
    inline void SetTrackId( G4int id ) { fTrackId = id; }
    inline G4int GetTrackId() const { return fTrackId; }
    inline void SetEDep(G4double val) { fEnergy = val; }
    inline G4double GetEDep() const { return fEnergy; }
    inline void SetStepStatus( G4StepStatus val ) { fStatus =val; }
    inline G4StepStatus GetStepStatus() const { return fStatus; }
    inline void SetStepLength( G4double val ) { fStepLength = val; }
    inline G4double GetStepLength() const { return fStepLength; }

   private:
  
  // Data members
/*     G4ParticleDefinition* fIPD;
     G4double fIKEnergy;
     G4ThreeVector fIPosition;
     G4ThreeVector fIMomentumD;
*/
     G4ParticleDefinition* fEPD;
     G4double fEKEnergy;
     G4ThreeVector fEPosition;
     G4ThreeVector fEMomentumD;
     G4ThreeVector fEMomentum;
    G4int fChamberNum;
    G4int fTrackId;
    G4double fEnergy;
    G4double fStepLength;
    G4StepStatus fStatus;
    

};

typedef G4THitsCollection<BeamTestSiliconMonitorHit> BeamTestSiliconMonitorHitsCollection;

extern G4Allocator<BeamTestSiliconMonitorHit> BeamTestSiliconMonitorHitAllocator;

inline void* BeamTestSiliconMonitorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)BeamTestSiliconMonitorHitAllocator.MallocSingle();
  return aHit;
}

inline void BeamTestSiliconMonitorHit::operator delete(void* aHit)
{
   BeamTestSiliconMonitorHitAllocator.FreeSingle((BeamTestSiliconMonitorHit*) aHit);
}

#endif
