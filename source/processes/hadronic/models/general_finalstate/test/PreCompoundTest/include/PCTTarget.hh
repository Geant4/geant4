#ifndef PCTTarget_hh
#define PCTTarget_hh 1

#include <vector>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"

class PCTOptions;

class PCTTarget
{
public:
  PCTTarget(const PCTOptions& opt);
  ~PCTTarget();
  
  inline void SetMomentum(const G4ThreeVector& p);
  
  inline G4LorentzVector GetMomentum() const;
  
  inline G4bool IsValid() const;
  inline G4bool IsNucleus() const;
  inline G4bool IsNatural() const;
    
  G4int GetA();
  G4int GetZ() const;

  inline G4double GetMass() const;

private:
  inline PCTTarget();
  inline PCTTarget(const PCTTarget& target);
  inline const PCTTarget& operator=(const PCTTarget & right);

  void Initialize(const G4int a, const G4int z);
  void Initialize(const std::vector< std::pair<G4int,G4double> >& mat, const G4int z);

  
  void SetParticleMomentum(const G4ThreeVector& p);
  void SetMaterialMomentum(const G4ThreeVector& p);
  
  inline G4LorentzVector GetParticleMomentum() const;
  inline G4LorentzVector GetMaterialMomentum() const;
  
private:
  G4ParticleDefinition * theParticle;
  std::vector< std::pair<G4double,G4ParticleDefinition*> > * theMaterial;
  G4LorentzVector ParticleMomentum;
  std::vector< G4LorentzVector > MaterialMomentum;
  G4int lastIsotope;
};

#include "PCTTarget.icc"

#endif // PCTTarget_hh

